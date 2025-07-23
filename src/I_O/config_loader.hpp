#pragma once

#include <string>
#include <vector>
#include <map>


/**
 * @brief Configuration structure for the model.
 * This structure holds all the necessary parameters
 * and settings required to run the model simulation.
 */
struct ModelConfig {
    
    // Time period
    std::string start_date;
    std::string calendar;

    // Solver
    double dt;
    double rtol;
    double atol;
       
    // Parameters
    std::string parameters_file;

    // Initial conditions
    int initial_conditions_flag;
    float initial_value; // constant value for initial conditions
    std::string initial_conditions_filename;
    std::string initial_conditions_varname;
    std::string initial_conditions_id_varname;

    // Boundary conditions
    int boundary_conditions_flag;
    int boundary_conditions_resolution; // resolution in minutes
    std::string boundary_conditions_filename;
    std::string boundary_conditions_varname;
    std::string boundary_conditions_id_varname;

    // Reservoirs
    int reservoir_routing_flag; // 0 for no reservoir routing, 1 for reservoir routing
    std::string reservoir_file; // file containing reservoir data if reservoir routing is used

    // Runoff
    int input_flag; // 0 for single file with time chunks, 1 for multiple files without time chunks (user input)
    int runoff_resolution; // resolution in minutes (user input)
    size_t chunk_size; // size of each time chunk in hours (user input)
    std::string runoff_path; // path to runoff data files
    std::string runoff_varname; // variable name for runoff data
    std::string runoff_id_varname; // ID variable name for runoff data

    // Output
    int output_flag;
    int min_level;
    int output_resolution;
    std::string link_list_filename;
    std::string series_filepath; 
    std::string snapshot_filepath; 

};

/**
 * @brief YAML parser for configuration files.
 * This class provides methods to parse a YAML configuration file
 * and retrieve configuration parameters in a structured way.   
 */
// Simple YAML parser class
class SimpleYamlParser {
public:
    // Public interface
    void parseFile(const std::string& filename);
    
    // Getter methods
    std::string getString(const std::string& key, const std::string& defaultValue = "");
    int getInt(const std::string& key, int defaultValue = 0);
    double getDouble(const std::string& key, double defaultValue = 0.0);

private:
    // Member variables
    std::map<std::string, std::string> keyValueMap;
    std::map<std::string, std::vector<std::map<std::string, std::string>>> arrayMap;
    std::map<std::string, std::vector<std::string>> simpleArrayMap;
    
    // Private helper methods
    void parseLines(const std::vector<std::string>& lines);
    std::string getSectionKey(const std::vector<std::string>& path);
    
    // Static utility functions
    static std::string trim(const std::string& str);
    static std::string removeQuotes(const std::string& str);
    static bool isInlineArray(const std::string& str);
    static std::vector<std::string> parseInlineArray(const std::string& str);
    static int getIndentLevel(const std::string& line);
    static bool isArrayItem(const std::string& line);
    static bool isComment(const std::string& line);
};

/**
 * @brief Configuration loader class.
 * This class provides a method to load the model configuration
 * from a YAML file and return a ModelConfig object.
 */
class ConfigLoader {
public:
    static ModelConfig loadConfig(const std::string& filename);
};