#include "routing.hpp"

// C++ standard libraries
#include <atomic>
#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
#include <omp.h>

//my functions
#include "dependency_graph.hpp"
#include "boundary_exchange.hpp"
#include "partition.hpp"
#include "model_setup.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "models/RHS.hpp"
#include "utils/time.hpp"
#include "utils/level_timing.hpp"

//--------------------------------------------------------------------------------------------------
// Function Definitions
//--------------------------------------------------------------------------------------------------
/**
 * @brief Writes the output results to netCDF files.
 * This function writes the final time step as a snapshot and also writes time series data based on user-defined criteria.
 * @param setup The model setup containing configuration and node information.
 * @param results The results vector containing the integrated values for each link.
 * @param n_steps The number of time steps in the simulation.
 * @param sim_times The vector of time points corresponding to the results.
 * @param q_final The vector to store final results for each link.
 * @param time_string A string representing the current time in the simulation.
 */
void writeOutput(const ModelSetup& setup,
                 const Partition& part,
                 const std::vector<float>& results,
                 size_t n_steps,
                 const std::vector<int>& sim_times,
                 std::vector<float>& q_final,
                 const std::string& time_string,
                 bool is_last_chunk)
{
    // Every array here is in local index space: a rank writes only the links it owns.
    // Local indices are assigned in increasing global order, so a single-rank run walks
    // the same links in the same order as it did before partitioning existed.
    const size_t n_owned = part.n_owned();

    // A distributed run writes one file per rank; a single-rank run keeps the original
    // name, so its output stays directly comparable with earlier runs.
    const std::string suffix = part.distributed()
        ? "_rank" + std::to_string(part.rank) : std::string();

    //------------------------------- SNAPSHOT OUTPUT -----------------------------------------------------------------
    std::cout << "  Writing final time step (snapshot) to netcdf...";
    std::vector<int> stream_ids(n_owned);
    size_t last_step = n_steps - 1;
    for (size_t i_link = 0; i_link < n_owned; ++i_link) {
        q_final[i_link] = results[i_link * n_steps + last_step];
        stream_ids[i_link] = setup.node_map.at(part.global_of[i_link]).stream_id;
    }
    if (!setup.config.snapshot_per_year || is_last_chunk) {
    std::string snapshot_filename = setup.config.snapshot_filepath + "_" + time_string + suffix + ".nc";
    write_snapshot_netcdf(snapshot_filename, q_final.data(), stream_ids.data(), n_owned);
    std::cout << "completed!" << std::endl;

    }
    // --------------------------------- MAXIMUM OUTPUT -----------------------------------------------------------
    if( setup.config.max_output == 1) {
        std::cout << "  Writing maximum values to netcdf...";
        // Find maximum values for each link
        std::vector<float> max_results(n_owned, 0.0f);
        // Parallelize over all links
        #pragma omp parallel for
        for (size_t i_link = 0; i_link < n_owned; ++i_link) {
            float local_max = 0.0f;
            for (size_t t = 0; t < n_steps; ++t) {
                float val = results[i_link * n_steps + t];
                if (val > local_max) {
                    local_max = val;
                }
            }
            max_results[i_link] = local_max;
        }
        std::string max_filename = setup.config.max_output_filepath + "_" + time_string + suffix + ".nc";
        write_snapshot_netcdf(max_filename, max_results.data(), stream_ids.data(), n_owned);
        std::cout << " completed!" << std::endl;
    }


    // --------------------------------- Time series output -----------------------------------------------------------
    if (setup.config.output_flag == 0) {
        std::cout << "No output requested." << std::endl;
        return;
    }
    std::cout << "  Writing time series to netcdf...";
    
    std::vector<size_t> keep_indices;
    std::vector<int> keep_links;
    // keep_indices holds LOCAL indices: a rank can only write series it computed.
    if (setup.config.output_flag == 1) {
        std::cout << "Outputting subset by level >= " << setup.config.min_level << "...";
        for (const auto& [id, node] : setup.node_map) {
            if (node.level >= setup.config.min_level && part.owns(node.index)) {
                keep_indices.push_back(part.local_of[node.index]);
                keep_links.push_back(node.stream_id);
            }
        }
    }
    else if (setup.config.output_flag == 2) {
        std::cout << "Outputting subset by list...";
        for (const auto& [id, node] : setup.node_map) {
            if (node.level > 0 && setup.save_info.stream_ids.count(node.stream_id)
                && part.owns(node.index)) {
                keep_indices.push_back(part.local_of[node.index]);
                keep_links.push_back(node.stream_id);
            }
        }
    }


    // Compact results based on keep_indices
    size_t n_keep_links = keep_indices.size();
    // Number of steps to skip per output
    size_t output_res_steps = static_cast<size_t>(setup.config.output_resolution / setup.config.dt);
    size_t n_saved_steps = (n_steps + output_res_steps - 1) / output_res_steps;
    std::vector<float> compacted_results(n_saved_steps * n_keep_links);
    // Parallel nested loop with fixed indexing
    #pragma omp parallel for
    for (size_t i = 0; i < n_saved_steps; ++i) {
        size_t t = i * output_res_steps;
        if(t >= n_steps) t = n_steps - 1; // clamp to last step
        for (size_t j = 0; j < n_keep_links; ++j) {
            size_t link_index = keep_indices[j];
            compacted_results[i * n_keep_links + j] = results[link_index * n_steps + t];
        }
    }

    std::string series_filename = setup.config.series_filepath + "_" + time_string + suffix + ".nc";
    write_timeseries_netcdf(series_filename,
                           compacted_results.data(),
                           sim_times.data(),
                           keep_links.data(),
                           n_saved_steps,
                           n_keep_links,
                           setup.config.calendar,
                           time_string);

    std::cout << "completed!" << std::endl;
}

/**
 * @brief Solves the ODE for a single link and writes its slice of the results matrix.
 * Every traversal calls SolveLink.
 * Every parent of this link must already have
 * written its slice of @p results before this is entered.
 *
 * Parent inflow is accumulated in routing-table order (node.parents).
 *
 * @param setup The model setup containing configuration and node information.
 * @param runoff The runoff data for the current time chunk.
 * @param results The results vector to store the integrated values for each link.
 * @param node The link to solve.
 * @param n_steps The number of time steps in the simulation.
 * @param total_time_steps The total number of time steps processed so far.
 * @param tc The current time chunk index.
 * @param q_final The vector holding each link's final value from the previous chunk.
 */
static void SolveLink(const ModelSetup& setup,
                      const Partition& part,
                      const BoundaryExchange& ex,
                      const RunoffData& runoff,
                      std::vector<float>& results,
                      const NodeInfo& node,
                      size_t n_steps,
                      size_t total_time_steps,
                      size_t tc,
                      std::vector<float>& q_final)
{
    // results and q_final are in local index space; the link's own slice is here.
    const size_t local = part.local_of[node.index];
    // Prefine items for solver
    auto rk45_dopri_stepper = make_controlled(setup.config.atol,
                            setup.config.rtol,
                            rk45_type());
    double start_time = 0.0;                               // in minutes
    double end_time   = (n_steps-1) * setup.config.dt;         // total time in minutes

    const size_t level = node.level;

    // Initialize the inflow series (y_p_series) for this link. This will be used to store inflow from parent nodes or boundary conditions
    std::vector<float> y_p_series(n_steps, 0.0);
    size_t y_p_resolution = setup.config.dt; // Resolution in minutes for y_p_series, default to dt unless boundary conditions are used

    //Check for BC
    bool has_bc = (setup.config.boundary_conditions_flag == 1) &&
                (setup.boundary_conditions.idToIndex.find(node.stream_id) != setup.boundary_conditions.idToIndex.end());
    if (has_bc) {
        size_t bc_index = setup.boundary_conditions.idToIndex.at(node.stream_id);
        size_t nTime = setup.boundary_conditions.nTime;
        size_t t_start = total_time_steps/setup.config.boundary_conditions_resolution; // Convert total_time_steps to units of boundary conditions
        y_p_resolution = setup.config.boundary_conditions_resolution; // resolution in minutes for y_p_series if boundary conditions are used
        size_t bc_steps = static_cast<size_t>(y_p_resolution / setup.config.dt); // How many ODE steps correspond to one BC time step
        for (size_t step = 0; step < n_steps; ++step) {
            // Map ODE step index to BC index
            size_t bc_time_index = t_start + step / bc_steps;
            if (bc_time_index < nTime) {
                y_p_series[step] = static_cast<double>(
                    setup.boundary_conditions.data[bc_index * nTime + bc_time_index]
                );
            }
        }
    } else if (level > 0) {
        for (size_t parent_index : node.parents) {
            // A parent on another rank contributes the series received in the boundary
            // exchange; a local one contributes its own slice of results. The two are
            // added in the same routing-table order either way, so where a parent lives
            // cannot change the sum.
            const size_t parent_local = part.local_of[parent_index];
            const float* parent_series = nullptr;
            if (parent_local != Partition::NOT_OWNED) {
                parent_series = results.data() + parent_local * n_steps;
            } else {
                auto it = ex.arrived.find(parent_index);
                if (it == ex.arrived.end()) {
                    std::cerr << "Error: link " << node.index << " needs parent "
                              << parent_index << " from rank " << part.rank_of[parent_index]
                              << ", which did not arrive in the boundary exchange."
                              << std::endl;
                    exit(EXIT_FAILURE);
                }
                parent_series = it->second;
            }
            #pragma omp simd //vectorization
            for (size_t t = 0; t < n_steps; ++t) {
                y_p_series[t] += parent_series[t];
            }
        }
    }

    // If reservoir routing is not needed, we can proceed with the ODE integration
    // This is a placeholder for future implementation
    if(setup.config.reservoir_routing_flag == 0){
        //Get initial condition for this link
        double q0;
        if(tc == 0){
            q0 = setup.uini(node.stream_id); // initial condition for this link from user
        }else{
            q0 = q_final[local]; // final step from the previous chunk, in local index space

            if (q0 <= 0.0) {
                std::cerr << "Warning: Initial discharge for link " << node.index << " is non-positive." << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        //  Parameters for the ODE
        const double A_h = node.params[0]; // hillslope area in m^2
        const double lambda_1 = node.params[2];
        const double L_i = node.params[1]; // stream link length in m
        const double v_0 = node.params[3]; // reference channel velocity in m/s
        const double invtau = 60.0 * v_0 / ((1.0 - lambda_1) * L_i);

        //runoff pointer
        const size_t runoff_index = runoff.idToIndex.at(node.stream_id);
        const float* runoff_ptr = &runoff.data[runoff_index * runoff.nTime];

        //solve ODE
        // Callback function to store results
        auto callback = [&](const double& x, const double t) {
            // Convert time t to step index
            size_t step_idx = static_cast<size_t>(t / setup.config.dt);
            // Clamp to last valid index
            if (step_idx >= n_steps) step_idx = n_steps - 1;
            size_t idx = local * n_steps + step_idx;
            results[idx] = std::max(x, 1e-8);
        };
        RHS rhs(runoff_ptr, setup.config.runoff_resolution,
                y_p_series, y_p_resolution,
                A_h,lambda_1,invtau);

        //integrator based on level
        if(level <= setup.config.rk4_level){
            integrate_const(rk4_stepper, rhs, q0, start_time, end_time, setup.config.dt, callback);
        }else{
            integrate_const(rk45_dopri_stepper,
                            rhs,
                            q0,
                            start_time,
                            end_time,
                            setup.config.dt,                       // time step in minutes
                            callback);
        }
    }else{
        // Placeholder for reservoir routing logic
        //exit code with failure
        std::cerr << "Reservoir routing is not implemented yet. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Integrates the ODEs for each link at a given level.
 * This function uses OpenMP to parallelize the integration process for each link.
 * Every link at a level is independent: a link's parents all sit at lower levels and
 * were solved before the barrier that closed the previous level.
 *
 * @param setup The model setup containing configuration and node information.
 * @param runoff The runoff data for the current time chunk.
 * @param results The results vector to store the integrated values for each link.
 * @param level The current level of links being processed.
 * @param nodes_at_level The vector of node indices at the current level.
 * @param n_steps The number of time steps in the simulation.
 * @param total_time_steps The total number of time steps processed so far.
 * @param tc The current time chunk index.
 * @param q_final The vector to store final results for each link.
 * @param profiler Per-level timer. Records the level's wall time and each link's solve
 *                 time when enabled; a no-op otherwise. Does not affect results.
 */
void IntegrateLinksAtLevel(const ModelSetup& setup,
                           const Partition& part,
                           const BoundaryExchange& ex,
                           const RunoffData& runoff,
                           std::vector<float>& results,
                           size_t level,
                           const std::vector<size_t>& nodes_at_level,
                           size_t n_steps,
                           size_t total_time_steps,
                           size_t tc,
                           std::vector<float>& q_final,
                           LevelProfiler& profiler)
{
    profiler.beginLevel(tc, level, nodes_at_level.size());

    // Solve ODE for each link at this level.
    // The schedule is set at runtime from profiling.omp_schedule, so static, dynamic and
    // guided can be compared without a rebuild. It defaults to static, as before.
    #pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < nodes_at_level.size(); ++i) {
        double link_start = profiler.enabled() ? omp_get_wtime() : 0.0;
        size_t link_index = nodes_at_level[i];
        if (!part.owns(link_index)) continue;  // another rank solves this one
        const NodeInfo& node = setup.node_map.at(link_index);  // Safe to access from multiple threads

        SolveLink(setup, part, ex, runoff, results, node, n_steps, total_time_steps, tc, q_final);

        if (profiler.enabled()) profiler.recordLink(omp_get_wtime() - link_start);
    }

    // Outside the loop, so this runs after the implicit barrier and the level's wall
    // time includes the wait.
    profiler.endLevel();
}

//--------------------------------------------------------------------------------------------------
// Dependency-driven traversal
//--------------------------------------------------------------------------------------------------
/**
 * @brief Solves one link, then releases its downstream link if it was the last parent.
 * This is where the dependency is enforced; there is no barrier.
 *
 * Parents may finish in any order. Exactly one of them reads back the old value 1 and
 * spawns the child, so the child is spawned once: no lost wakeup, no double spawn.
 * acq_rel also publishes this link's writes to results and acquires the other parents',
 * so the child reads complete parent series.
 *
 * @param ctx Shared context for the chunk being solved.
 * @param link_index The link to solve.
 */
static void SolveAndRelease(TaskContext ctx, size_t link_index)
{
    const NodeInfo& node = ctx.setup.node_map.at(link_index);
    SolveLink(ctx.setup, ctx.part, ctx.ex, ctx.runoff, ctx.results, node,
              ctx.n_steps, ctx.total_time_steps, ctx.tc, ctx.q_final);

    const size_t child_index = ctx.graph.child[link_index];
    if (child_index == DependencyGraph::NO_CHILD) return;  // an outlet releases nothing
    // A child on another rank is released by the boundary exchange after this chunk, not
    // by a counter here. Its rank does not count this parent in its in-degree.
    if (!ctx.part.owns(child_index)) return;

    if (ctx.pending[child_index].fetch_sub(1, std::memory_order_acq_rel) == 1) {
        // ctx is firstprivate, not shared: the task is deferred, and this frame is gone
        // long before another thread runs it. Sharing it would leave the task reading a
        // dead stack frame -- which segfaults late and at random, not here.
        #pragma omp task firstprivate(ctx, child_index) default(shared)
        SolveAndRelease(ctx, child_index);
    }
}

/**
 * @brief Solves every link in the network by dependency rather than by level.
 * A link becomes eligible as soon as its own upstream links are done, so a finished
 * headwater branch runs downstream while other branches are still upstream.
 *
 * NOTE ON STACK: tasks run nested inside barrier waits, so stack depth grows with
 * concurrency, and Boost's integrate_const is not cheap in stack. Set OMP_STACKSIZE
 * generously. Check here first if this path crashes at width.
 *
 * @param setup The model setup containing configuration and node information.
 * @param runoff The runoff data for the current time chunk.
 * @param results The results vector to store the integrated values for each link.
 * @param graph The dependency graph, built once for the run.
 * @param pending Per-link countdown, reset here for this chunk.
 * @param n_steps The number of time steps in the simulation.
 * @param total_time_steps The total number of time steps processed so far.
 * @param tc The current time chunk index.
 * @param q_final The vector to store final results for each link.
 */
void IntegrateLinksByDependency(const ModelSetup& setup,
                                const Partition& part,
                                const BoundaryExchange& ex,
                                const RunoffData& runoff,
                                std::vector<float>& results,
                                const DependencyGraph& graph,
                                std::vector<std::atomic<int>>& pending,
                                size_t n_steps,
                                size_t total_time_steps,
                                size_t tc,
                                std::vector<float>& q_final)
{
    // Reset the countdown for this chunk. The network never changes, so this is a copy
    // of the in-degrees rather than a rebuild of the graph 
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < pending.size(); ++i) {
        pending[i].store(graph.in_degree[i], std::memory_order_relaxed);
    }

    TaskContext ctx{setup, part, ex, runoff, results, graph, pending.data(), q_final,
                    n_steps, total_time_steps, tc};

    #pragma omp parallel default(shared)
    {
        // One thread seeds the headwaters; every other thread is already at the closing
        // barrier and starts taking tasks immediately.
        #pragma omp single
        {
            for (size_t source_index : graph.sources) {
                #pragma omp task firstprivate(ctx, source_index) default(shared)
                SolveAndRelease(ctx, source_index);
            }
        }
    }
}

/**
 * @brief Processes a single chunk of runoff data.
 * This function reads runoff data from a netCDF file, sets up the time series, initializes the results matrix,
 * and integrates the ODEs for each link at different levels.
 * It also handles the output of results to netCDF files.   
 * @param setup The model setup containing configuration and node information.
 * @param tc The current time chunk index.
 * @param total_time_steps The total number of time steps processed so far.
 * @param q_final The vector to store final results for each link.  
 * @param startIndex The start index for the current chunk.
 */

void ProcessChunk(const ModelSetup& setup,
                  const Partition& part,
                  BoundaryExchange& ex,
                  size_t tc,
                  size_t& total_time_steps,
                  std::vector<float>& q_final,
                  size_t& startIndex,
                  std::vector<float>& results,
                  const DependencyGraph& graph,
                  std::vector<std::atomic<int>>& pending,
                  LevelProfiler& profiler)
{
    std::cout << "Processing chunk/file " << tc + 1 << " of " << setup.runoff_info.nchunks << ":" << std::endl;

    // Compute the start index for this chunk if files change
    if (tc > 0 && setup.runoff_info.filenames[tc] != setup.runoff_info.filenames[tc - 1]) startIndex = 0;
    
    //Compute starttime for this chunk
    std::string time_string = addTimeDelta(setup.config.start_date, setup.config.calendar, total_time_steps); //time string to store the start time for this chunk
    std::cout << "  Start time for this chunk: " << time_string << std::endl;

    // ----------------- RUNOFF DATA --------------------------------------
    std::cout << "  Reading in runoff from netcdf file: " << setup.runoff_info.filenames[tc] << "...";
    auto read_start = std::chrono::high_resolution_clock::now();
    RunoffData runoff = readTotalRunoff(setup.runoff_info.filenames[tc], 
                                        setup.config.runoff_varname, 
                                        setup.config.runoff_id_varname,
                                        startIndex,
                                        setup.config.chunk_size);
    auto read_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_elapsed = read_end - read_start;
    std::cout << "completed!" << std::endl;
    std::cout << "  Total read in time: " << read_elapsed.count() << " seconds" << std::endl;

    //Update start index for the next chunk
    startIndex += setup.config.chunk_size;

    // -------------------- TIME SERIES SETUP --------------------------------------  
    // User defined parameters for simulation time (user input)
    size_t t_final = runoff.nTime * setup.config.runoff_resolution; // minutes in a file chunk from input file
    size_t n_steps = t_final / setup.config.dt; // number of steps based on dt 
    // Times to store results for outputs
    std::vector<int> sim_times; // still in minutes for NetCDF
    size_t output_res_steps = static_cast<size_t>(setup.config.output_resolution / setup.config.dt);
    // output_resolution and dt both in minutes → steps per output
    for (size_t step = 0; step < n_steps; step += output_res_steps) {
        sim_times.push_back(static_cast<int>(step * setup.config.dt)); // store time in minutes
    }

    // Initialize the results matrix
    std::cout << "  Allocating space for results...";
    auto alloc_start = std::chrono::high_resolution_clock::now();
    size_t required_size = n_steps * part.n_owned();
    // Only resize if vector is smaller than needed
    if (results.size() < required_size) {
        results.resize(required_size);
    }
    auto alloc_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> alloc_elapsed = alloc_end - alloc_start;
    std::cout << "completed!" << std::endl;
    std::cout << "  Total allocation time: " << alloc_elapsed.count() << " seconds" << std::endl;

    //  -------------- SOLVING ODEs ----------------------------------
    std::cout << "  Starting integration for each link..."  << std::flush;
    auto solve_start = std::chrono::high_resolution_clock::now();
    // Every series this rank needs from upstream ranks, before any link is solved. The
    // call blocks, so a waiting rank sleeps rather than spinning.
    ReceiveBoundaries(ex, n_steps);
    if (setup.config.traversal == "counter") {
        // Dependency-driven: a link runs as soon as its own upstream links are done.
        IntegrateLinksByDependency(setup, part, ex, runoff, results, graph, pending, n_steps, total_time_steps, tc, q_final);
    } else {
        // Level-synchronous: loop through each level and process nodes.
        for (const auto& [level, nodes_at_level] : setup.level_groups) {
            IntegrateLinksAtLevel(setup, part, ex, runoff, results, level, nodes_at_level, n_steps, total_time_steps, tc, q_final, profiler);
        }
    }
    // Hand this rank's cut-edge series to the ranks downstream of it.
    SendBoundaries(ex, part, results, n_steps);
    std::cout << "completed!" << std::endl;
    auto solve_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> solve_elapsed = solve_end - solve_start ;
    std::cout << "  Total integration time: " << solve_elapsed.count() << " seconds" << std::endl;

    // Append this chunk's per-level rows now, so a job that is killed part way through a
    // multi-chunk run still leaves the levels it did finish on disk.
    profiler.writeChunk();

    // Update total time steps for the next chunk
    total_time_steps += t_final; //time in minutes for this chunk

    // -----------OUTPUT --------------------------------------------
    auto write_start = std::chrono::high_resolution_clock::now();
    bool is_last_chunk = (tc + 1 == static_cast<size_t>(setup.runoff_info.nchunks));
    writeOutput(setup, part, results, n_steps, sim_times, q_final, time_string, is_last_chunk);
    auto write_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> write_elapsed = write_end - write_start;
    std::cout << "  Total write time: " << write_elapsed.count() << " seconds" << std::endl;
    
}

//-------------------------------------------------------------------------------------------------
// MAIN FUNCTION
//-------------------------------------------------------------------------------------------------
/**
 * @brief Runs the routing process for the given model setup.
 * This function processes runoff data in time chunks, integrates the ODEs for each link,
 * and outputs the results to netCDF files.
 * @param setup The model setup containing configuration, node information, and other parameters.
 */
void runRouting(const ModelSetup& setup, int rank, int n_ranks){
    std::cout << "_________________STARTING ROUTING_________________ \n" << std::endl;
    // Which links this rank owns. Without mpi.partition_file this is one rank owning
    // everything, with local index == global index, i.e. the behaviour before
    // partitioning existed.
    const Partition part = LoadPartition(setup.n_links, setup.config.mpi_partition_file, rank, n_ranks);

    std::vector<float> q_final(part.n_owned()); //final value per owned link, local indices

    std::cout << "  Traversal: " << setup.config.traversal
              << (setup.config.traversal == "counter" ? " (dependency-driven)"
                                                      : " (level-synchronous)") << std::endl;

    // Per-level timing, off unless profiling.level_timing is 1
    LevelProfiler profiler;
    if (setup.config.profile_level_timing == 1) {
        if (setup.config.traversal == "counter") {
            // The counter traversal has no levels to time
            std::cout << "  Note: profiling.level_timing does not apply to the counter "
                      << "traversal. Per-level timing is off for this run." << std::endl;
        } else {
            profiler.enable(setup.config.profile_filepath);
        }
    }

    // Built once: the network is fixed for the whole run, so only the countdown is reset
    // per chunk. Costs nothing on the level path beyond the setup pass.
    // The graph is needed by the counter traversal and by the exchange (for child[]), so
    // it is built whenever either applies.
    DependencyGraph graph;
    std::vector<std::atomic<int>> pending;
    if (setup.config.traversal == "counter" || part.distributed()) {
        graph = BuildDependencyGraph(setup, part);
    }
    if (setup.config.traversal == "counter") {
        pending = std::vector<std::atomic<int>>(setup.n_links);
    }
    BoundaryExchange ex = BuildBoundaryExchange(part, graph);
    ex.lookahead = setup.config.mpi_lookahead_chunks;
    if (part.distributed()) {
        std::cout << "  Pipelining: lookahead_chunks = " << ex.lookahead
                  << (ex.lookahead ? " (a rank runs a chunk ahead)" : " (lockstep)")
                  << std::endl;
    }

    //reserving max size up front
    size_t max_size = static_cast<size_t>(setup.config.chunk_size * setup.config.runoff_resolution * part.n_owned() / setup.config.dt);
    std::vector<float> results;          // declare the vector
    results.reserve(max_size);           // reserve memory upfront

    // process chunks
    size_t total_time_steps = 0; // keep tract of total simulation time
    size_t startIndex = 0; // start index for the first chunk
    for(int tc = 0; tc < setup.runoff_info.nchunks; ++tc){
        ProcessChunk(setup, part, ex, tc, total_time_steps, q_final, startIndex, results, graph, pending, profiler);
    }
    FinishBoundaries(ex);   // no message may still be in flight at MPI_Finalize
    std::cout << "__________________________________________________ \n" << std::endl;

    profiler.printSummary();
}
// End of file: Tiger_HLM_Routing/src/routing.cpp