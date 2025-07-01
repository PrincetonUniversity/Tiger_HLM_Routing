#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath> // Required for pow
#include <boost/numeric/odeint.hpp>
#include <omp.h>

//my functions
#include "I_O/node_info.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "models/RHS.hpp"


// C++ standard libraries
using namespace std;
using namespace boost::numeric::odeint;


// Tasks to be completed:
// 1. Move RHS to a separate file.
// 2. Add index for resolution to the RHS function to avoid recalculating the index for each time step.
// 3. Add time calculation for dates from start time
// 4. Add a res function to handle the reservoir routing (place holder).
// 5. Add ability to output from list of LinkIDs (stream IDs) to a netcdf file.
// 6. Add reading user inputs from yaml. 
// 7. Add checks for:
//   - If the initial conditions are valid for the given links.
//   - If the boundary conditions are valid for the entire simulation time.
// 8. Move all parts of main to functions to improve readability and maintainability.


// ROUTING ODES--------------------------------------------------
// struct RHS {
//     const float* runoff_series; // runoff, hourly
//     const std::vector<double>& y_p_series; // inflow from parent nodes, hourly (for now)
//     const double A_h; // hillslope area in m^2
//     const double lambda_1; // exponent
//     const double invtau; // constant term

//     RHS(const float* runoff_series,
//         const std::vector<double>& y_p_series,
//         const double A_h, 
//         const double lambda_1, 
//         const double invtau
//         )
//         : runoff_series(runoff_series),
//         y_p_series(y_p_series),
//         A_h(A_h),
//         lambda_1(lambda_1),
//         invtau(invtau){}

//     inline void operator()(const double q, double& dQdt, const double t) const {
//         // t is in minutes, inputs are in hours 
//         size_t hour_idx = t/60;  // assumes t in minutes (!!!!need to pass index for runoff and y_p series)

//         // runoff_series is hourly 
//         double runoff = runoff_series[hour_idx]*(0.001/60); // convert from mm/h to m/min (assuming runoff is in mm/h)

//         // y_p is hourly — 60 minutes per step
//         double y_p = y_p_series[hour_idx];

//         // Nonlinear routing ODE
//         double q_safe = std::max(q, 1e-8); //set a safe minimum to avoid division by zero
//         dQdt = invtau * pow(q_safe, lambda_1) * (-q_safe + (runoff * A_h / 60.0) + y_p);
//     }
// };


// // RK45 Solver and methodology
// // state_type = value_type = deriv_type = time_type = double
// typedef runge_kutta_dopri5< double , double , double , double , vector_space_algebra , default_operations , never_resizer > stepper_type;

int main()
{   
    // ----------------- HEADER --------------------------------------
    std::cout << "___________________BUILD INFO_____________________ \n" << std::endl;
    std::cout << "Tiger HLM routing version 0.1" << std::endl;
    std::cout << "Tiger HLM routing build date: " << __DATE__ << " " << __TIME__ << std::endl;
    std::cout << "Tiger HLM routing build by: am2192" << std::endl;
    std::cout << "__________________________________________________ \n" << std::endl;

    // ----------------- SETUP OMP --------------------------------------
    std::cout << "___________________OpenMP INFO____________________ \n" << std::endl;
    const char* env_threads = std::getenv("OMP_NUM_THREADS");
    if (env_threads) {
        // Environment variable is set — let OpenMP use it
        std::cout << "OMP_NUM_THREADS is set to " << env_threads << "! \n" << std::endl;
    } else {
        // Not set — manually set to 1
        omp_set_num_threads(1);
        std::cout << "OMP_NUM_THREADS not set — defaulting to 1 thread." << std::endl;    
        std::cout << "If you are running this on a SLURM cluster, please set the number of threads using the SLURM_CPUS_PER_TASK variable." << std::endl;
        std::cout << "For example, you can set it in your SLURM script as follows: export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK" << std::endl;
    }
    std::cout << "__________________________________________________ \n" << std::endl;

    // ----------------- STARTING ROUTING --------------------------------------
    std::cout << "_________________STARTING ROUTING_________________ \n" << std::endl;

    // ----------------- INPUTS --------------------------------------

    // Read node levels from CSV file
    std::cout << "Reading node levels and parameters from CSV...";
    std::unordered_map<size_t, NodeInfo> node_map;
    std::map<size_t, std::vector<size_t>> level_groups;
    std::string node_levels_filename = "../data/node_levels_params.csv"; //user input 
    read_node_levels(node_levels_filename, node_map, level_groups);
    size_t n_links = node_map.size(); //number of links used for allocating results
    std::cout << "completed!" << std::endl;


    // Initial conditions (optional, can be set to a constant value)
    std::cout << "Loading initial conditions...";
    int initial_conditions_flag = 0; // 0 for constant value, 1 for reading from file
    float initial_value = 1.0; // constant value for q0
    // If reading from file, set the filename, variable name, and ID variable name
    std::string initial_conditions_filename = "/scratch/gpfs/GVILLARI/am2192/snapshot.nc"; //user input 
    std::string initial_conditions_varname = "snapshot"; //user input 
    std::string initial_conditions_id_varname = "LinkID"; //user input 
    auto uini = loadInitialConditions(initial_conditions_flag,
                                     initial_value, initial_conditions_filename, 
                                     initial_conditions_varname, 
                                     initial_conditions_id_varname);
    std::cout << "completed!" << std::endl;

    //read boundary conditions from file if they exist
    std::cout << "Reading boundary conditions from file...";
    int boundary_conditions_flag = 0; // 0 for no boundary conditions, 1 for reading from file for whole simulation
    int boundary_conditions_resolution = 60; // resolution in minutes (user input)
    std::string boundary_conditions_filename = "/scratch/gpfs/GVILLARI/am2192/routing/BC.nc"; //user input
    std::string boundary_conditions_varname = "BC"; //user input
    std::string boundary_conditions_id_varname = "LinkID"; //user input
    BoundaryConditions  boundary_conditions; // to store boundary conditions if they exist
    if(boundary_conditions_flag == 1){
        boundary_conditions = readBoundaryConditions(boundary_conditions_filename, 
                                                    boundary_conditions_varname, 
                                                    boundary_conditions_id_varname); 
    }
    std::cout << "completed!" << std::endl;

    // TIME CHUNKING STARTS HERE ------------------------------------------
    int input_flag = 1; // 0 for single file with time chunks, 1 for multiple files without time chunks (user input)
    int runoff_resolution = 60; // resolution in minutes (user input)
    size_t chunk_size = 2000; // size of each time chunk in hours (user input)
    std::string path = "/scratch/gpfs/GVILLARI/am2192/routing/data/"; ///scratch/gpfs/GVILLARI/am2192/routing/total_runoff_test.nc
    std::string runoff_varname = "ro"; //user input 
    std::string runoff_id_varname = "LinkID"; //user input 

    // Get runoff chunk info
    RunoffChunkInfo runoff_info = getRunoffChunkInfo(path, input_flag, chunk_size);

    //define q_final to store final results for each link
    std::vector<float> q_final(n_links); //updated each loop

    // Users parameters for simulation time
    int simulation_resolution = 60; // resolution in minutes (user input)
    double dt = 1.0; //user input for solving ODE time step in minutes

    // keep tract of total simulation time
    size_t total_time_steps = 0; // total time steps across all chunks

    //need to define q_final to store final results for each link
    for(int tc = 0; tc < runoff_info.nchunks; ++tc){ // Loop over time chunks or multiple files
        std::cout << "Processing chunk/file " << tc + 1 << " of " << runoff_info.nchunks << ":" << std::endl;
        size_t startIndex = tc * chunk_size; // start index for this chunk (need to set default chunk size for when no time chunking)


        // ----------------- RUNOFF DATA --------------------------------------

        // Runoff data
        std::cout << "  Read in runoff from netcdf file..." << runoff_info.filenames[tc] << "...";
        RunoffData runoff = readTotalRunoff(input_flag,
                                            runoff_info.filenames[tc], 
                                            runoff_varname, 
                                            runoff_id_varname,
                                            startIndex,
                                            chunk_size);
        std::cout << "completed!" << std::endl;

        // -------------------- TIME SERIES SETUP --------------------------------------
        
        // User defined parameters for simulation time (user input)
        double tf = runoff.nTime * runoff_resolution; // minutes in a file chunk from input file 
        
        //times to store results
        size_t n_steps = static_cast<size_t>(tf / simulation_resolution);
        std::vector<int> times(n_steps);
        for(size_t i = 0; i < n_steps; ++i){
            times[i] = i * simulation_resolution; // time in minutes
        }     

        // Initialize the results matrix
        std::cout << "  Allocating space for results...";
        std::vector<float> results(n_steps * n_links, 0.0);
        std::cout << "completed!" << std::endl;

        //  -------------- SOLVING ODEs ----------------------------------
        std::cout << "  Starting integration for each link..."  << std::flush;
        auto start = std::chrono::high_resolution_clock::now();

        // Loop through each level and process nodes
        for (const auto& [level, nodes_at_level] : level_groups) {
            // Solve ODE for each link at this level
            #pragma omp parallel for
            for (size_t i = 0; i < nodes_at_level.size(); ++i) {
                size_t link_index = nodes_at_level[i];
                const NodeInfo& node = node_map.at(link_index);  // Safe to access from multiple threads
                

                // Initialize the inflow series (y_p_series) for this link
                // This will be used to store inflow from parent nodes or boundary conditions
                std::vector<double> y_p_series(n_steps, 0.0);

                bool has_bc = (boundary_conditions_flag == 1) &&
                            (boundary_conditions.idToIndex.find(node.stream_id) != boundary_conditions.idToIndex.end());

                if (has_bc) {
                    size_t bc_index = boundary_conditions.idToIndex.at(node.stream_id);
                    size_t nTime = boundary_conditions.nTime;
                    size_t t_start = total_time_steps/boundary_conditions_resolution; // Convert total_time_steps to hours

                    for (size_t t = t_start; t < n_steps; ++t) {
                        y_p_series[t] = static_cast<double>(
                            boundary_conditions.data[bc_index * nTime + t]
                        );
                    }
                } else if (level > 0) {
                    for (size_t parent_index : node.parents) {
                        for (size_t t = 0; t < n_steps; ++t) {
                            y_p_series[t] += results[t * n_links + parent_index];
                        }
                    }
                }

                //Get initial condition for this link
                double q0;
                if(tc == 0){
                    q0 = uini(node.stream_id); // initial condition for this link from user
                }else{
                    q0 = q_final[node.index]; // final step from results which uses node.index
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
                    size_t idx = static_cast<size_t>(t / simulation_resolution);
                    if (idx < n_steps)
                        results[idx * n_links + node.index] = x;
                };

                RHS rhs(runoff_ptr, y_p_series,A_h,lambda_1,invtau);
                auto stepper = make_controlled(1E-9, 1E-6, stepper_type());
                integrate_times(stepper, rhs, q0, times.begin(), times.end(), dt, callback);

                // //Solving for res logic
                // if(has_res == 0){
                //     RHS rhs(runoff_template, y_p_series,A_h, lambda_1, L_i, v_0, invtau);

                //     auto stepper = make_controlled(1E-9, 1E-6, stepper_type());
                //     integrate_times(stepper, rhs, q0, times.begin(), times.end(), dt, callback);
                // }
                // else{
                //     res_func()
                // }

            }
        }
        std::cout << "completed!" << std::endl;

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "  Total integration time: " << elapsed.count() << " seconds" << std::endl;


        // -----------OUTPUT --------------------------------------------

        // Process results for output
        std::cout << "  Processing results for output...";

        // Get final time step for all links
        std::vector<int> stream_ids(n_links);
        size_t last_step = n_steps - 1;
        for(size_t i_link = 0; i_link < n_links; ++i_link) {
            q_final[i_link] = results[last_step * n_links + i_link];
            stream_ids[i_link] = node_map.at(i_link).stream_id;
        }

        // subset for level 2 and above:
        int min_level = 4; // user input for minimum level to keep links
        if (min_level < 2) {
            int min_level = 2; // Ensure minimum level is at least 2
        }
        // Step 1: Collect indices of links to keep (those with level >= min_level)
        std::vector<size_t> keep_indices;
        std::vector<int> keep_links;
        for (const auto& [id, node] : node_map) {
            if (node.level >= min_level) { //user input 
                keep_indices.push_back(node.index);
                keep_links.push_back(node.stream_id);
            }
        }
        size_t n_keep_links = keep_indices.size();

        // Step 2: Compact results in-place
        size_t write_pos = 0;
        for (size_t t = 0; t < n_steps; ++t) {
            for (size_t link_index : keep_indices) {
                if (link_index < n_links) {
                    size_t read_pos = t * n_links + link_index;
                    results[write_pos++] = results[read_pos];
                }
            }
        }
        results.resize(n_steps * keep_indices.size()); //Resize results to new size
        std::cout << "  completed!" << std::endl;


        // Save to netcdf
        std::cout << "  Writing results to NetCDF...";

        std::string series_filename = "/scratch/gpfs/GVILLARI/am2192/routing/output/series"; // user input
        series_filename += std::to_string(tc + 1) + ".nc"; // append chunk number to filename

        std::string snapshot_filename = "/scratch/gpfs/GVILLARI/am2192/routing/output/snapshot"; // user input
        snapshot_filename += std::to_string(tc + 1) + ".nc"; // append chunk number to filename
        
        int compression_level=0;
        // Dump time series to netcdf41321
        write_timeseries_netcdf(series_filename,
                            results.data(),
                            times.data(),
                            keep_links.data(),
                            n_steps,
                            n_keep_links,
                            compression_level);    
        // Snapshot output
        write_snapshot_netcdf(snapshot_filename, 
                              q_final.data(), 
                              stream_ids.data(),
                              n_links,
                              compression_level);

        std::cout << "completed!" << std::endl;

        // Update total time steps
        total_time_steps += tf; //time in minutes for this chunk


    }

    // ----------------- END OF ROUTING --------------------------------------
    std::cout << "_________________END OF ROUTING__________________ \n" << std::endl;
    std::cout << "Routing completed successfully!" << std::endl;
    std::cout << "Thank you for using Tiger HLM routing!" << std::endl;

    
    return 0;
}