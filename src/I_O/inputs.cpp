#include "inputs.hpp"
//
#include <netcdf.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <filesystem>
#include <set>
#include <cctype>

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
        size_t nTimeSteps = GetNCTimeSize(filename, varname);
        // If chunk_size is 0, treat all files as a single chunk
        if(chunk_size == 0){
            info.filenames.push_back(filename.c_str());
            info.ntime.push_back(nTimeSteps);
        }
        // Chunk files based on user chunk size
        else {
            if(nTimeSteps <= chunk_size){
                // If chunk size is larger than number of time steps or zero, treat as single file
                info.filenames.push_back(filename.c_str()); // Add the single file path
                info.ntime.push_back(nTimeSteps);
            }else{
                //chunk size plus one for the last chunk if required
                size_t nchunks = nTimeSteps / chunk_size;
                if (nTimeSteps % chunk_size != 0) {
                    nchunks += 1;
                }
                for(int i=0; i < nchunks; ++i){
                    info.filenames.push_back(filename.c_str());
                    info.ntime.push_back(std::min<size_t>(chunk_size, nTimeSteps - i * chunk_size));
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

    return dim_sizes[1];
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
 * @brief Adds the (LinkID, value) pairs of one snapshot file to map.
 * @param reject_duplicates Exit if a LinkID is already in map; used when combining
 *                          per-rank files, where each link must come from exactly one file.
 */
static void readSnapshotInto(const std::string& filename,
                             const std::string& varname,
                             const std::string& id_varname,
                             std::unordered_map<int, float>& map,
                             bool reject_duplicates){

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

    for (size_t i = 0; i < dim_size; ++i) {
        if (!reject_duplicates) {
            map[ids[i]] = data[i];
        } else if (!map.emplace(ids[i], data[i]).second) {
            std::cerr << "Error: LinkID " << ids[i] << " in " << filename
                      << " also appears in another per-rank snapshot file. Remove stale"
                      << " rank files left by an earlier run." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}


/**
 * @brief Lists the per-rank snapshot files <stem>_rank<N>.nc for a filename <stem>.nc, sorted.
 */
static std::vector<std::string> findRankSnapshots(const std::string& filename){
    namespace fs = std::filesystem;
    const fs::path target(filename);
    const fs::path dir = target.has_parent_path() ? target.parent_path() : fs::path(".");
    const std::string prefix = target.stem().string() + "_rank";
    const std::string ext = target.extension().string();

    std::vector<std::string> files;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(dir, ec)) {
        const std::string name = entry.path().filename().string();
        if (name.size() <= prefix.size() + ext.size()
            || name.compare(0, prefix.size(), prefix) != 0
            || name.compare(name.size() - ext.size(), ext.size(), ext) != 0) {
            continue;
        }
        const std::string rank = name.substr(prefix.size(), name.size() - prefix.size() - ext.size());
        if (std::all_of(rank.begin(), rank.end(), [](unsigned char c) { return std::isdigit(c); })) {
            files.push_back(entry.path().string());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}


/**
 * @brief Reads initial conditions for the routing model.
 * @param flag Indicates how to read initial conditions:
 *             0 - Constant value for q0 (initial_value must be provided).
 *             1 - Read from a file (filename, varname, id_varname must be provided).
 *                 If filename does not exist, the per-rank files <stem>_rank<N>.nc written
 *                 by a distributed run are read instead.
 * @param initial_value The constant value for q0 if flag is 0, and for links missing from the file.
 * @param filename The path to the file containing initial conditions if flag is 1.
 * @param varname The name of the variable in the file containing initial conditions.
 * @param id_varname The name of the variable in the file containing link IDs.
 * @param n_links Number of links in the network; a warning is printed if the file covers fewer.
 * @return A function that takes a link ID and returns its initial condition.
 */

std::function<float(int)> loadInitialConditions(const int flag,
                                        const float initial_value,
                                        const std::string& filename,
                                        const std::string& varname,
                                        const std::string& id_varname,
                                        const size_t n_links){


    // If flag is 0: return constant function
    if (flag == 0) {
        return [initial_value](int) { return initial_value; };
    }

    // If flag is 1: read from netcdf file
    if (flag == 1){

        std::unordered_map<int, float> map;
        if (std::filesystem::exists(filename)) {
            readSnapshotInto(filename, varname, id_varname, map, false);
        } else {
            // Every rank reads all per-rank files, so the previous run's rank count
            // does not need to match this run's.
            const std::vector<std::string> parts = findRankSnapshots(filename);
            if (parts.empty()) {
                std::cerr << "Error: initial condition file " << filename
                          << " not found, and no per-rank files (<name>_rank<N>.nc) next to it."
                          << std::endl;
                exit(EXIT_FAILURE);
            }
            for (const auto& part : parts) {
                readSnapshotInto(part, varname, id_varname, map, true);
            }
            std::cout << "read " << parts.size() << " per-rank files...";
        }

        if (map.size() < n_links) {
            std::cerr << "Warning: initial conditions cover " << map.size() << " links but the"
                      << " network has " << n_links << "; the rest start at " << initial_value
                      << "." << std::endl;
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



