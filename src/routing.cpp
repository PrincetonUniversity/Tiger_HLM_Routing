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
                 const std::string& time_string) 
{   
    //------------------------------- SNAPSHOT OUTPUT -----------------------------------------------------------------
    std::cout << "  Writing final time step (snapshot) to netcdf...";
    std::vector<int> stream_ids(setup.n_links);
    size_t last_step = n_steps - 1;
    for (size_t i_link = 0; i_link < setup.n_links; ++i_link) {
        q_final[i_link] = results[i_link * n_steps + last_step];
        stream_ids[i_link] = setup.node_map.at(i_link).stream_id;
    }
    std::string snapshot_filename = setup.config.snapshot_filepath + "_" + time_string + ".nc";
    write_snapshot_netcdf(snapshot_filename, q_final.data(), stream_ids.data(), setup.n_links);
    std::cout << "completed!" << std::endl;

    // --------------------------------- MAXIMUM OUTPUT -----------------------------------------------------------
    if( setup.config.max_output == 1) {
        std::cout << "  Writing maximum values to netcdf...";
        // Find maximum values for each link
        std::vector<float> max_results(setup.n_links, 0.0f);
        // Parallelize over all links
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
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < nodes_at_level.size(); ++i) {
        size_t link_index = nodes_at_level[i];
        const NodeInfo& node = setup.node_map.at(link_index);  // Safe to access from multiple threads

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

    //  -------------- SOLVING ODEs ----------------------------------
    std::cout << "  Starting integration for each link..."  << std::flush;
    auto solve_start = std::chrono::high_resolution_clock::now();
    // Loop through each level and process nodes
    for (const auto& [level, nodes_at_level] : setup.level_groups) {
        IntegrateLinksAtLevel(setup, runoff, results, level, nodes_at_level, n_steps, total_time_steps, tc, q_final);
    }
    std::cout << "completed!" << std::endl;
    auto solve_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> solve_elapsed = solve_end - solve_start ;
    std::cout << "  Total integration time: " << solve_elapsed.count() << " seconds" << std::endl;

    // Update total time steps for the next chunk
    total_time_steps += t_final; //time in minutes for this chunk

    // -----------OUTPUT --------------------------------------------
    auto write_start = std::chrono::high_resolution_clock::now();
    writeOutput(setup, results, n_steps, sim_times, q_final, time_string);
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