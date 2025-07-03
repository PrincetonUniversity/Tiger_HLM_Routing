
#include "model_setup.hpp"  // Your struct and setupModel() declaration

// C++ standard libraries
#include <iostream>         // For std::cout, std::endl
#include <unordered_map>   // For std::unordered_map
#include <map>             // For std::map
#include <vector>          // For std::vector
#include <string>          // For std::string

// External libraries
#include "I_O/node_info.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "I_O/config_loader.hpp"



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
                                                     setup.config.input_flag,
                                                     setup.config.chunk_size);
    std::cout << "completed!" << std::endl;

    // OUTPUT OPTIONS------------------------
    std::cout << "Setting up output options...";
    if(setup.config.output_flag == 2) setup.save_info = readSaveList(setup.config.link_list_filename); //read the save list from file
    std::cout << "completed!" << std::endl;
    std::cout << "__________________________________________________ \n" << std::endl;


    return setup;
}