
#include "model_setup.hpp"  // Your struct and setupModel() declaration

// C++ standard libraries
#include <iostream>         // For std::cout, std::endl

// External libraries
#include "I_O/node_info.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "I_O/config_loader.hpp"


/**
 * @brief Sets up the model by loading configuration, node levels, initial conditions,
 * boundary conditions, runoff data, and output options.    
 */
ModelSetup setupModel(const char* config_path) {
    ModelSetup setup;

    // ----------------- SETUP --------------------------------------
    std::cout << "_________________MODEL SET UP_____________________ \n" << std::endl;

    // INPUTS --------------------------------------
    // Read user inputs from YAML file
    std::cout << "Loading user inputs from YAML file...";
    setup.config = ConfigLoader::loadConfig(config_path);
    std::cout << "completed!" << std::endl;

    // Read node levels from CSV file
    std::cout << "Loading network parameters...";
    read_node_levels(setup.config.parameters_file, setup.node_map, setup.level_groups);
    setup.n_links = setup.node_map.size(); //number of links used for allocating results
    std::cout << "completed!" << std::endl;

    // Initial conditions (optional, can be set to a constant value)
    std::cout << "Loading initial conditions...";
    setup.uini = loadInitialConditions(setup.config.initial_conditions_flag,
                                     setup.config.initial_value,
                                     setup.config.initial_conditions_filename,
                                     setup.config.initial_conditions_varname,
                                     setup.config.initial_conditions_id_varname);
    std::cout << "completed!" << std::endl;

    //read boundary conditions from file if they exist
    std::cout << "Loading boundary condition...";
    if(setup.config.boundary_conditions_flag == 1){
        setup.boundary_conditions = readBoundaryConditions(setup.config.boundary_conditions_filename,
                                                          setup.config.boundary_conditions_varname,
                                                          setup.config.boundary_conditions_id_varname);
    }
    std::cout << "completed!" << std::endl;

    // Check if reservoir routing is needed
    // This is a placeholder for future implementation
    std::cout << "Checking if reservoir routing is needed (placeholder)...";
    std::cout << "completed!" << std::endl;


    //Get runoff chunk info
    std::cout << "Loading runoff data chunk info...";
    setup.runoff_info = getRunoffChunkInfo(setup.config.runoff_path,
                                           setup.config.runoff_varname,
                                           setup.config.chunk_size);
    std::cout << "completed!" << std::endl;

    // OUTPUT OPTIONS------------------------
    std::cout << "Setting up output options...";
    if(setup.config.output_flag == 2) setup.save_info = readSaveList(setup.config.link_list_filename); //read the save list from file
    std::cout << "completed!" << std::endl;
    std::cout << "__________________________________________________ \n" << std::endl;

    
    // Print out memory information ----------------------
    // Kills system above 1TB
    size_t input_memory = setup.n_links * setup.config.chunk_size * sizeof(float); // Memory for input runoff data
    size_t results_memory = setup.n_links * setup.config.chunk_size * setup.config.runoff_resolution * sizeof(float); // Memory for results vector
    size_t output_memory = 0.5 * (setup.n_links * setup.config.chunk_size * setup.config.runoff_resolution * sizeof(float)) / setup.config.output_resolution; // Memory for output results
    size_t total_memory = (input_memory + results_memory + output_memory) / (1024.0 * 1024.0 * 1024.0); // Assuming 1.1x for intermediate results and overhead
    if(total_memory > 500.0){
        std::cout << "Warning: Estimated memory usage (saving level 1 and above) is high (" << total_memory << " GB). Consider reducing chunk size or ensure saving out levels higher than 1 only. Please monitor your system's memory usage." << std::endl;
    }
    return setup;
}