#pragma once
//
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>



/**
 * @brief Structure to hold boundary conditions for the routing model.
 *  * This structure contains the boundary condition data, the number of links, 
 * the number of time steps, the IDs of the links, and a mapping from link IDs to their corresponding indices in the data vector.
 * This could be combined with RunoffData, but kept separate for clarity.
 */
struct BoundaryConditions {
    std::vector<float> data; // Boundary condition data
    size_t nLink;            // Number of links
    size_t nTime;            // Number of time steps
    std::vector<int> ids;    // Link IDs
    std::unordered_map<int, size_t> idToIndex; // Map from link ID to index in data vector
};
BoundaryConditions readBoundaryConditions(const std::string& filename, 
                            const std::string& varname,
                            const std::string& id_varname);


/**
 * @brief Structure to hold information about runoff chunks.
 * This structure contains the number of chunks and their filenames.
 * It is used to manage runoff data that may be split into multiple files or time chunks.
 */
struct RunoffChunkInfo{
    int nchunks = 1; // Number of chunks
    std::vector<std::string> filenames;
};

RunoffChunkInfo getRunoffChunkInfo(const std::string& path, 
                            const int flag = 0, // 0 for single file with time chunks, 1 for multiple files without time chunks
                            const int chunk_size = 0 // Size of each chunk in hours, 0 for no chunking
                        );
size_t GetNCTimeSize(const std::string& filename);



/**
 * @brief Structure to hold runoff data read from a NetCDF file.
 * This structure contains the runoff data, the number of links,
 * the number of time steps, the IDs of the links, and a mapping from
 * link IDs to their corresponding indices in the data vector.
 */
struct RunoffData {
    std::vector<float> data;
    size_t nLink;
    size_t nTime;
    std::vector<int> ids;
    std::unordered_map<int, size_t> idToIndex;
};

RunoffData readTotalRunoff(const int flag,
                           const std::string& filename, 
                           const std::string& varname, \
                           const std::string& id_varname,
                           const size_t startIndex,
                           const size_t chunk_size);


/**
 * @brief Loads initial conditions for the routing model.
 * This function can either return a constant value for q0 or read initial conditions from a file.
 * 
 * @param flag Determines the behavior of the function:
 *             - 0: Returns a constant function with the given initial_value.
 *             - 1: Reads initial conditions from a file (not implemented in this snippet).
 * @param initial_value The constant value to return if flag is 0.
 * @param filename The name of the file to read initial conditions from (if flag is 1).
 * @param varname The variable name in the file containing the initial conditions.
 * @param id_varname The variable name in the file containing link IDs.
 * @return A function that takes an integer (link ID) and returns the corresponding initial condition.
 */
std::function<float(int)> loadInitialConditions(const int flag = 0, 
                                                const float initial_value = 1.0, 
                                                const std::string& filename = "", 
                                                const std::string& varname = "snapshot", 
                                                const std::string& id_varname="LinkID");