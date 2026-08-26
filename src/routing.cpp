#include "routing.hpp"

// C++ standard libraries
#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;
#include <omp.h>

//my functions
#include "model_setup.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "models/RHS.hpp"
#include "utils/time.hpp"

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
                 const std::vector<float>& results,
                 size_t n_steps,
                 const std::vector<int>& sim_times,
                 std::vector<float>& q_final,
                 const std::string& time_string,
                 bool is_last_chunk) 
{   
        //------------------------------- SNAPSHOT OUTPUT -----------------------------------------------------------------
    // q_final MUST be computed every chunk — it's the state handoff to the next chunk.
    size_t last_step = n_steps - 1;
    #pragma omp parallel for
    for (size_t i_link = 0; i_link < setup.n_links; ++i_link) {
        q_final[i_link] = results[i_link * n_steps + last_step];
    }

    // stream_ids is shared by both the snapshot and max_output blocks below,
    // so it's computed once here regardless of which flags are set.
    std::vector<int> stream_ids(setup.n_links);
    for (size_t i_link = 0; i_link < setup.n_links; ++i_link) {
        stream_ids[i_link] = setup.node_map.at(i_link).stream_id;
    }

    //------------------------------- SNAPSHOT OUTPUT -----------------------------------------------------------------
    if (!setup.config.snapshot_per_year || is_last_chunk) {
        std::cout << "  Writing final time step (snapshot) to netcdf ("
                  << (setup.config.snapshot_per_year ? "per-year" : "per-chunk")
                  << ")...";
        std::string snapshot_filename = setup.config.snapshot_filepath + "_" + time_string + ".nc";
        write_snapshot_netcdf(snapshot_filename, q_final.data(), stream_ids.data(), setup.n_links);
        std::cout << "completed!" << std::endl;
    }

    // --------------------------------- MAXIMUM OUTPUT -----------------------------------------------------------
    // Deliberately NOT gated by snapshot_per_year — writes every chunk regardless,
    // since it's a running max, independent of the snapshot checkpoint cadence.
    if (setup.config.max_output == 1) {
        std::cout << "  Writing maximum values to netcdf...";
        std::vector<float> max_results(setup.n_links, 0.0f);
        #pragma omp parallel for
        for (size_t i_link = 0; i_link < setup.n_links; ++i_link) {
            float local_max = 0.0f;
            for (size_t t = 0; t < n_steps; ++t) {
                float val = results[i_link * n_steps + t];
                if (val > local_max) {
                    local_max = val;
                }
            }
            max_results[i_link] = local_max;
        }
        std::string max_filename = setup.config.max_output_filepath + "_" + time_string + ".nc";
        write_snapshot_netcdf(max_filename, max_results.data(), stream_ids.data(), setup.n_links);
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
    if (setup.config.output_flag == 1) {
        std::cout << "Outputting subset by level >= " << setup.config.min_level << "...";
        for (const auto& [id, node] : setup.node_map) {
            if (node.level >= setup.config.min_level) {
                keep_indices.push_back(node.index);
                keep_links.push_back(node.stream_id);
            }
        }
    }
    else if (setup.config.output_flag == 2) {
        std::cout << "Outputting subset by list...";
        for (const auto& [id, node] : setup.node_map) {
            if (node.level > 0 && setup.save_info.stream_ids.count(node.stream_id)) {
                keep_indices.push_back(node.index);
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

    std::string series_filename = setup.config.series_filepath + "_" + time_string + ".nc";
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
 * @brief Integrates the ODEs for each link at a given level.
 * This function uses OpenMP to parallelize the integration process for each link.
 * It initializes the inflow series for each link, checks for boundary conditions, and integrates the ODEs.
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
 */
void IntegrateLinksAtLevel(const ModelSetup& setup,
                           const RunoffData& runoff,
                           std::vector<float>& results,
                           size_t level,
                           const std::vector<size_t>& nodes_at_level,
                           size_t n_steps,
                           size_t total_time_steps,
                           size_t tc,
                           std::vector<float>& q_final)
{   
    // Prefine items for solver
    auto rk45_dopri_stepper = make_controlled(setup.config.atol, 
                            setup.config.rtol, 
                            rk45_type());
    double start_time = 0.0;                               // in minutes
    double end_time   = (n_steps-1) * setup.config.dt;         // total time in minutes
    
    // Solve ODE for each link at this level
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < nodes_at_level.size(); ++i) {
        size_t link_index = nodes_at_level[i];
        const NodeInfo& node = setup.node_map.at(link_index);  // Safe to access from multiple threads

        // // Initialize the inflow series (y_p_series) for this link. This will be used to store inflow from parent nodes or boundary conditions
        // std::vector<float> y_p_series(n_steps, 0.0);
        // Reused per-thread buffer, persists across calls on the same OS
        // thread, avoiding a fresh heap allocation per node/chunk.
        thread_local std::vector<float> y_p_series;
        y_p_series.assign(n_steps, 0.0f);
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
                #pragma omp simd //vectorization
                for (size_t t = 0; t < n_steps; ++t) {
                    y_p_series[t] += results[parent_index * n_steps + t];
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
                q0 = q_final[node.index]; // final step from results which uses node.index

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
                size_t idx = node.index * n_steps + step_idx;
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
}

/**
 * @brief Per-node integration body — identical to the loop body inside
 * IntegrateLinksAtLevel, extracted so IntegrateWithTasks can call it
 * from each case of its parent-count switch.
 */
template <typename Rk45Stepper>
inline void integrate_one_node_body(const ModelSetup& setup,
                                     const RunoffData& runoff,
                                     std::vector<float>& results,
                                     size_t level,
                                     size_t link_index,
                                     size_t n_steps,
                                     size_t total_time_steps,
                                     size_t tc,
                                     std::vector<float>& q_final,
                                     Rk45Stepper& rk45_dopri_stepper,
                                     double start_time,
                                     double end_time)
{
    const NodeInfo& node_local = setup.node_map.at(link_index);

    thread_local std::vector<float> y_p_series;
    y_p_series.assign(n_steps, 0.0f);
    size_t y_p_resolution = setup.config.dt;

    bool has_bc = (setup.config.boundary_conditions_flag == 1) &&
                  (setup.boundary_conditions.idToIndex.find(node_local.stream_id)
                   != setup.boundary_conditions.idToIndex.end());

    if (has_bc) {
        size_t bc_index = setup.boundary_conditions.idToIndex.at(node_local.stream_id);
        size_t nTime = setup.boundary_conditions.nTime;
        size_t t_start = total_time_steps / setup.config.boundary_conditions_resolution;
        y_p_resolution = setup.config.boundary_conditions_resolution;
        size_t bc_steps = static_cast<size_t>(y_p_resolution / setup.config.dt);
        for (size_t step = 0; step < n_steps; ++step) {
            size_t bc_time_index = t_start + step / bc_steps;
            if (bc_time_index < nTime) {
                y_p_series[step] = static_cast<double>(
                    setup.boundary_conditions.data[bc_index * nTime + bc_time_index]);
            }
        }
    } else if (level > 0) {
        for (size_t parent_index : node_local.parents) {
            #pragma omp simd
            for (size_t t = 0; t < n_steps; ++t) {
                y_p_series[t] += results[parent_index * n_steps + t];
            }
        }
    }

    if (setup.config.reservoir_routing_flag == 0) {
        double q0;
        if (tc == 0) {
            q0 = setup.uini(node_local.stream_id);
        } else {
            q0 = q_final[node_local.index];
            if (q0 <= 0.0) {
                std::cerr << "Warning: Initial discharge for link "
                          << node_local.index << " is non-positive." << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        const double A_h = node_local.params[0];
        const double lambda_1 = node_local.params[2];
        const double L_i = node_local.params[1];
        const double v_0 = node_local.params[3];
        const double invtau = 60.0 * v_0 / ((1.0 - lambda_1) * L_i);

        const size_t runoff_index = runoff.idToIndex.at(node_local.stream_id);
        const float* runoff_ptr = &runoff.data[runoff_index * runoff.nTime];

        auto callback = [&](const double& x, const double t) {
            size_t step_idx = static_cast<size_t>(t / setup.config.dt);
            if (step_idx >= n_steps) step_idx = n_steps - 1;
            size_t idx = node_local.index * n_steps + step_idx;
            results[idx] = std::max(x, 1e-8);
        };
        RHS rhs(runoff_ptr, setup.config.runoff_resolution,
                y_p_series, y_p_resolution,
                A_h, lambda_1, invtau);

        if (level <= setup.config.rk4_level) {
            integrate_const(rk4_stepper, rhs, q0, start_time, end_time,
                            setup.config.dt, callback);
        } else {
            integrate_const(rk45_dopri_stepper, rhs, q0, start_time, end_time,
                            setup.config.dt, callback);
        }
    } else {
        std::cerr << "Reservoir routing is not implemented yet. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Testing task-scheduling approach
/**
 * @brief Dataflow-scheduled integration using OpenMP task depend().
 * Alternative to IntegrateLinksAtLevel: runs the entire chunk in one
 * parallel region, emits one task per node with dependencies declared
 * on parents' output slices. Nodes start the moment their parents
 * finish, regardless of "level" per-level barriers. The old functionality
 * of per level barriers should be retained and maybe have the user specify which
 * function to use. Added use_task_scheduling flag to the config.
 * Uses explicit switch over parent count (max 5 hardcoded for testing).
 * Level is extracted to a plain local to avoid the icx limitation on
 * capturing structured bindings inside OpenMP tasks.
 */
void IntegrateWithTasks(const ModelSetup& setup,
                         const RunoffData& runoff,
                         std::vector<float>& results,
                         size_t n_steps,
                         size_t total_time_steps,
                         size_t tc,
                         std::vector<float>& q_final)
{
    auto rk45_dopri_stepper = make_controlled(setup.config.atol,
                                                setup.config.rtol,
                                                rk45_type());
    double start_time = 0.0;
    double end_time = (n_steps - 1) * setup.config.dt;

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (const auto& kv : setup.level_groups) {
                size_t level_val = kv.first;                       // plain local, safe to capture in tasks
                const std::vector<size_t>& nodes_at_level = kv.second;

                for (size_t link_index : nodes_at_level) {
                    const NodeInfo& node = setup.node_map.at(link_index);
                    float* my_out = &results[node.index * n_steps];
                    size_t n_parents = node.parents.size();

                    // Pointers to parent output slices, one per possible parent.
                    float* p0 = (n_parents > 0) ? &results[node.parents[0] * n_steps] : my_out;
                    float* p1 = (n_parents > 1) ? &results[node.parents[1] * n_steps] : my_out;
                    float* p2 = (n_parents > 2) ? &results[node.parents[2] * n_steps] : my_out;
                    float* p3 = (n_parents > 3) ? &results[node.parents[3] * n_steps] : my_out;
                    float* p4 = (n_parents > 4) ? &results[node.parents[4] * n_steps] : my_out;

                    // Switch on parent count max is 5 per CONUS East diagnostic.
                    switch (n_parents) {
                    case 0:
                        #pragma omp task depend(out: my_out[0]) \
                            firstprivate(link_index, level_val)
                        {
                            integrate_one_node_body(setup, runoff, results, level_val,
                                                     link_index, n_steps, total_time_steps,
                                                     tc, q_final, rk45_dopri_stepper,
                                                     start_time, end_time);
                        }
                        break;
                    case 1:
                        #pragma omp task depend(in: p0[0]) depend(out: my_out[0]) \
                            firstprivate(link_index, level_val)
                        {
                            integrate_one_node_body(setup, runoff, results, level_val,
                                                     link_index, n_steps, total_time_steps,
                                                     tc, q_final, rk45_dopri_stepper,
                                                     start_time, end_time);
                        }
                        break;
                    case 2:
                        #pragma omp task depend(in: p0[0], p1[0]) depend(out: my_out[0]) \
                            firstprivate(link_index, level_val)
                        {
                            integrate_one_node_body(setup, runoff, results, level_val,
                                                     link_index, n_steps, total_time_steps,
                                                     tc, q_final, rk45_dopri_stepper,
                                                     start_time, end_time);
                        }
                        break;
                    case 3:
                        #pragma omp task depend(in: p0[0], p1[0], p2[0]) depend(out: my_out[0]) \
                            firstprivate(link_index, level_val)
                        {
                            integrate_one_node_body(setup, runoff, results, level_val,
                                                     link_index, n_steps, total_time_steps,
                                                     tc, q_final, rk45_dopri_stepper,
                                                     start_time, end_time);
                        }
                        break;
                    case 4:
                        #pragma omp task depend(in: p0[0], p1[0], p2[0], p3[0]) depend(out: my_out[0]) \
                            firstprivate(link_index, level_val)
                        {
                            integrate_one_node_body(setup, runoff, results, level_val,
                                                     link_index, n_steps, total_time_steps,
                                                     tc, q_final, rk45_dopri_stepper,
                                                     start_time, end_time);
                        }
                        break;
                    case 5:
                        #pragma omp task depend(in: p0[0], p1[0], p2[0], p3[0], p4[0]) depend(out: my_out[0]) \
                            firstprivate(link_index, level_val)
                        {
                            integrate_one_node_body(setup, runoff, results, level_val,
                                                     link_index, n_steps, total_time_steps,
                                                     tc, q_final, rk45_dopri_stepper,
                                                     start_time, end_time);
                        }
                        break;
                    default:
                        std::cerr << "ERROR: node " << node.index << " has " << n_parents
                                  << " parents, exceeding hardcoded max of 5. "
                                  << "Increase switch cases in IntegrateWithTasks." << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        } // end omp single
    } // end omp parallel (implicit barrier waits for all tasks)
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
                  size_t tc, 
                  size_t& total_time_steps, 
                  std::vector<float>& q_final,
                  size_t& startIndex,
                  std::vector<float>& results)
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
    size_t required_size = n_steps * setup.n_links;
    // Only resize if vector is smaller than needed
    if (results.size() < required_size) {
        results.resize(setup.n_links * n_steps);
    }
    auto alloc_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> alloc_elapsed = alloc_end - alloc_start;
    std::cout << "completed!" << std::endl;
    std::cout << "  Total allocation time: " << alloc_elapsed.count() << " seconds" << std::endl;

    // //  -------------- SOLVING ODEs ----------------------------------
    // std::cout << "  Starting integration for each link..."  << std::flush;
    // auto solve_start = std::chrono::high_resolution_clock::now();
    // // Loop through each level and process nodes
    // // for (const auto& [level, nodes_at_level] : setup.level_groups) {
    // //     IntegrateLinksAtLevel(setup, runoff, results, level, nodes_at_level, n_steps, total_time_steps, tc, q_final);
    // // }
    // for (const auto& [level, nodes_at_level] : setup.level_groups) {
    //     double level_start = omp_get_wtime();
    //     IntegrateLinksAtLevel(setup, runoff, results, level, nodes_at_level, n_steps, total_time_steps, tc, q_final);
    //     double level_end = omp_get_wtime();
    //     std::cout << "    [LEVEL] chunk=" << tc
    //                << " level=" << level
    //                << " nodes=" << nodes_at_level.size()
    //                << " time=" << (level_end - level_start) << "s"
    //                << std::endl;
    // }
    // std::cout << "completed!" << std::endl;

        //  -------------- SOLVING ODEs ----------------------------------
    std::cout << "  Starting integration for each link..."  << std::flush;
    auto solve_start = std::chrono::high_resolution_clock::now();

    if (setup.config.use_task_scheduling) {
        std::cout << " [task-scheduled dataflow mode]" << std::flush;
        double dfs_start = omp_get_wtime();
        IntegrateWithTasks(setup, runoff, results, n_steps, total_time_steps, tc, q_final);
        double dfs_end = omp_get_wtime();
        std::cout << "\n    [TASKS] chunk=" << tc
                  << " total_time=" << (dfs_end - dfs_start) << "s"
                  << std::endl;
    } else {
        std::cout << " [level-barrier mode]" << std::flush;
        for (const auto& [level, nodes_at_level] : setup.level_groups) {
            double level_start = omp_get_wtime();
            IntegrateLinksAtLevel(setup, runoff, results, level, nodes_at_level, n_steps, total_time_steps, tc, q_final);
            double level_end = omp_get_wtime();
            std::cout << "    [LEVEL] chunk=" << tc
                       << " level=" << level
                       << " nodes=" << nodes_at_level.size()
                       << " time=" << (level_end - level_start) << "s"
                       << std::endl;
        }
    }
    std::cout << "completed!" << std::endl;

    auto solve_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> solve_elapsed = solve_end - solve_start ;
    std::cout << "  Total integration time: " << solve_elapsed.count() << " seconds" << std::endl;

    // Update total time steps for the next chunk
    total_time_steps += t_final; //time in minutes for this chunk

    // -----------OUTPUT --------------------------------------------
    bool is_last_chunk = (tc + 1 == static_cast<size_t>(setup.runoff_info.nchunks));
    auto write_start = std::chrono::high_resolution_clock::now();
    writeOutput(setup, results, n_steps, sim_times, q_final, time_string, is_last_chunk);
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
void runRouting(const ModelSetup& setup){
    std::cout << "_________________STARTING ROUTING_________________ \n" << std::endl;
    std::vector<float> q_final(setup.n_links); //define q_final to store final results for each link

    //reserving max size up front
    size_t max_size = static_cast<size_t>(setup.config.chunk_size * setup.config.runoff_resolution * setup.n_links / setup.config.dt);
    std::vector<float> results;          // declare the vector
    results.reserve(max_size);           // reserve memory upfront

    // process chunks
    size_t total_time_steps = 0; // keep tract of total simulation time
    size_t startIndex = 0; // start index for the first chunk
    for(int tc = 0; tc < setup.runoff_info.nchunks; ++tc){
        ProcessChunk(setup, tc, total_time_steps, q_final, startIndex, results);
    }
    std::cout << "__________________________________________________ \n" << std::endl;
}
// End of file: Tiger_HLM_Routing/src/routing.cpp