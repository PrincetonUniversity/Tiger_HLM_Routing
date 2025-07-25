#include "inputs.hpp"
//
#include <netcdf.h>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <set>

#define ERR(e) { std::cerr << "NetCDF error: " << nc_strerror(e) << " at " << __FILE__ << ":" << __LINE__ << std::endl; exit(EXIT_FAILURE); }

/**
 * @brief Reads boundary conditions from a NetCDF file. 
 * * This function reads a 2D variable (time, link) from the specified NetCDF file,
 * @param filename The path to the NetCDF file.
 * @param varname The name of the variable containing boundary condition data.
 * @param id_varname The name of the variable containing link IDs.
 * 
 */
BoundaryConditions readBoundaryConditions(const std::string& filename, 
                           const std::string& varname,
                           const std::string& id_varname){
    int ncid, varid, retval, idVarId;

    // Open the NetCDF file
    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
        ERR(retval);

    // Inquire the variable ID
    if ((retval = nc_inq_varid(ncid, varname.c_str(), &varid)))
        ERR(retval);

    // Inquire the variable dimensions
    // Note: The variable is expected to be 2D (Link, Time)
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    if ((retval = nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr)))
        ERR(retval);

    size_t dim_sizes[2];
    for (int i = 0; i < 2; ++i)
        if ((retval = nc_inq_dimlen(ncid, dimids[i], &dim_sizes[i])))
            ERR(retval);

    // Read the variable data
    size_t nLink = dim_sizes[0];
    size_t nTime = dim_sizes[1];
    std::vector<float> data(nLink * nTime);

    //Read in entire dataset
    if ((retval = nc_get_var_float(ncid, varid, data.data())))
        ERR(retval);


    // Now get the link ID variable: coordinate variable with same name as link dimension
    if ((retval = nc_inq_varid(ncid, id_varname.c_str(), &idVarId)))
        ERR(retval);

    // Read the ID values (assuming they're integers)
    std::vector<int> ids(nLink);
    if ((retval = nc_get_var_int(ncid, idVarId, ids.data())))
        ERR(retval);

    // Create mapping from ID value to index
    std::unordered_map<int, size_t> idToIndex;
    for (size_t i = 0; i < nLink; ++i)
        idToIndex[ids[i]] = i;

    // Close the NetCDF file
    if ((retval = nc_close(ncid)))
        ERR(retval);

    return {data, nLink, nTime, ids, idToIndex};
}



/**
 * @brief Gets information about runoff chunks based on the provided path.
 * @param path The path to the runoff data files or directory.
 * @param varname The name of the variable to read from the files.
 * @param chunk_size Size of each chunk in temporal resolution. If 0, no chunking is applied.
 * @return A RunoffChunkInfo structure containing the number of chunks and their filenames.
 */

RunoffChunkInfo getRunoffChunkInfo(const std::string& path, 
                                   const std::string& varname,
                                   const int chunk_size){
    RunoffChunkInfo info;
    
    // Get files from the specified path
    std::set<std::filesystem::path> sorted_by_name;
    for (auto &entry : std::filesystem::directory_iterator(path))
        sorted_by_name.insert(entry.path());
    
    // Ensure the path is valid and contains files
    if(sorted_by_name.empty()){
        std::cerr << "Error: No files found in the specified path: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    // push all filenames into the info struct
    for (auto &filename : sorted_by_name){
        // If chunk_size is 0, treat all files as a single chunk
        if(chunk_size == 0){
            info.filenames.push_back(filename.c_str()); 
        } 
        // Chunk files based on user chunk size
        else {
            // Get number of time steps in the file
            size_t nTimeSteps = GetNCTimeSize(filename, varname);
            if(nTimeSteps <= chunk_size){
                // If chunk size is larger than number of time steps or zero, treat as single file
                info.filenames.push_back(filename.c_str()); // Add the single file path
            }else{
                //chunk size plus one for the last chunk
                size_t nchunks = nTimeSteps / chunk_size + 1;
                for(int i=0; i < nchunks; ++i){
                    info.filenames.push_back(filename.c_str());
                }
            }
        }

    }

    //Chunks are number of files
    info.nchunks = info.filenames.size(); // number of files

    // Return empty info
    return info;
};

/**
 * @brief Gets the number of time steps in a NetCDF file.
 * @param filename The path to the NetCDF file.
 * @return The number of time steps in the file.
 */
size_t GetNCTimeSize(const std::string& filename,
                     const std::string& varname){
    
    int ncid, varid, retval;

    // Open the NetCDF file
    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
        ERR(retval);

    // Inquire the variable ID
    if ((retval = nc_inq_varid(ncid, varname.c_str(), &varid)))
        ERR(retval);

    // Inquire the variable dimensions
    // Note: The variable is expected to be 2D (link,time)
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    if ((retval = nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr)))
        ERR(retval);

    if (ndims < 2) {
        std::cerr << "Error: Variable '" << varname << "' is not 2D as expected.\n";
        exit(EXIT_FAILURE);
    }
        
    size_t dim_sizes[2];
    for (int i = 0; i < 2; ++i)
        if ((retval = nc_inq_dimlen(ncid, dimids[i], &dim_sizes[i])))
            ERR(retval);
    
    // Close the NetCDF file
    if ((retval = nc_close(ncid)))
        ERR(retval);

    return {dim_sizes[1]};
}



/**
 * @brief Reads total runoff data from a NetCDF file.
 * @param filename The path to the NetCDF file.
 * @param varname The name of the variable containing runoff data.
 * @param id_varname The name of the variable containing link IDs.
 * @return A RunoffData structure containing the runoff data, number of links, number of time steps, link IDs, and a mapping from ID to index.
 */
RunoffData readTotalRunoff(const std::string& filename, 
                           const std::string& varname, 
                           const std::string& id_varname,
                           const size_t startIndex,
                           const size_t chunk_size){
    int ncid, varid, retval, idVarId;

    // Open the NetCDF file
    if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
        ERR(retval);

    // Inquire the variable ID
    if ((retval = nc_inq_varid(ncid, varname.c_str(), &varid)))
        ERR(retval);

    // Inquire the variable dimensions
    // Note: The variable is expected to be 2D (link,time)
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    if ((retval = nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr)))
        ERR(retval);

    size_t dim_sizes[2];
    for (int i = 0; i < 2; ++i)
        if ((retval = nc_inq_dimlen(ncid, dimids[i], &dim_sizes[i])))
            ERR(retval);

    // Read the variable data
    size_t nLink = dim_sizes[0];
    size_t nTime = dim_sizes[1];
    std::vector<float> data(nLink * nTime);

    if(chunk_size == 0)
    {
        //Read in entire dataset
        if ((retval = nc_get_var_float(ncid, varid, data.data())))
            ERR(retval);

    } else {
        // Read a subset of the data based on start and size
        size_t edge_case = nTime - startIndex;
        nTime = std::min(chunk_size, edge_case); // Ensure size does not exceed available time steps
        data.resize(nLink * nTime);              // Resize the outer data vector

        // Define start and count arrays for subsetting
        size_t start[2] = {0, startIndex};
        size_t count[2] = {nLink, nTime};

        // get the total runoff variable
        if ((retval = nc_get_vara_float(ncid, varid, start, count, data.data())))
            ERR(retval);
    } 

    // Now get the link ID variable: coordinate variable with same name as link dimension
    if ((retval = nc_inq_varid(ncid, id_varname.c_str(), &idVarId)))
        ERR(retval);

    // Read the ID values (assuming they're integers)
    std::vector<int> ids(nLink);
    if ((retval = nc_get_var_int(ncid, idVarId, ids.data())))
        ERR(retval);

    // Create mapping from ID value to index
    std::unordered_map<int, size_t> idToIndex;
    for (size_t i = 0; i < nLink; ++i)
        idToIndex[ids[i]] = i;

    // Close the NetCDF file
    if ((retval = nc_close(ncid)))
        ERR(retval);

    return {data, nLink, nTime, ids, idToIndex};
}


/**
 * @brief Reads initial conditions for the routing model.
 * @param flag Indicates how to read initial conditions:
 *             0 - Constant value for q0 (initial_value must be provided).
 *             1 - Read from a file (filename, varname, id_varname must be provided).
 * @param initial_value The constant value for q0 if flag is 0.
 * @param filename The path to the file containing initial conditions if flag is 1.
 * @param varname The name of the variable in the file containing initial conditions.
 * @param id_varname The name of the variable in the file containing link IDs.
 * @return A structure containing initial conditions for each link.
 */

std::function<float(int)> loadInitialConditions(const int flag, 
                                        const float initial_value, 
                                        const std::string& filename, 
                                        const std::string& varname, 
                                        const std::string& id_varname){
    

    // If flag is 0: return constant function
    if (flag == 0) {
        return [initial_value](int) { return initial_value; };
    }

    // If flag is 1: read from netcdf file
    if (flag == 1){

        int ncid, varid, retval, idVarId;

        // Open the NetCDF file read-only
        if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid))) {
            ERR(retval);
        }

        // Get variable ID for snapshot variable (float array)
        if ((retval = nc_inq_varid(ncid, varname.c_str(), &varid))) {
            nc_close(ncid);
            ERR(retval);
        }

        // Query variable dimensions (expect 1D)
        int ndims;
        int dimids[NC_MAX_VAR_DIMS];
        if ((retval = nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dimids, nullptr))) {
            nc_close(ncid);
            ERR(retval);
        }

        // Get dimension length
        size_t dim_size;
        if ((retval = nc_inq_dimlen(ncid, dimids[0], &dim_size))) {
            nc_close(ncid);
            ERR(retval);
        }

        // Read snapshot data (float array)
        std::vector<float> data(dim_size);
        if ((retval = nc_get_var_float(ncid, varid, data.data()))) {
            nc_close(ncid);
            ERR(retval);
        }

        // Get variable ID for ID variable (assumed int)
        if ((retval = nc_inq_varid(ncid, id_varname.c_str(), &idVarId))) {
            nc_close(ncid);
            ERR(retval);
        }

        // Read the ID variable data
        std::vector<int> ids(dim_size);
        if ((retval = nc_get_var_int(ncid, idVarId, ids.data()))) {
            nc_close(ncid);
            ERR(retval);
        }

        // Close NetCDF file now that data is loaded
        if ((retval = nc_close(ncid))) {
            std::cerr << "Warning: NetCDF file close error: " << nc_strerror(retval) << std::endl;
        }

        // Build map from ID to snapshot value
        std::unordered_map<int, float> map;
        for (size_t i = 0; i < dim_size; ++i) {
            map[ids[i]] = data[i];
        }

        // Return lookup lambda
        return [map, initial_value](int key) {
            auto it = map.find(key);
            return (it != map.end()) ? it->second : initial_value;
        };
    }

    // Fallback if file not found or empty
    return [initial_value](int) { return initial_value; };
};



