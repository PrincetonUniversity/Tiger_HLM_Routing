#include "output_series.hpp"
#include <netcdf.h>
#include <string>
#include <iostream>
#include <vector>


#define NC_CHECK(call) \
do { \
    int status = (call); \
    if (status != NC_NOERR) { \
        std::cerr << "NetCDF error: " << nc_strerror(status) << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return; \
    } \
} while (0)


/**
 * @brief Write a streamflow array to a NetCDF file with optional compression.
 */
void write_timeseries_netcdf(const std::string& filename,
                        const float* results,
                        const int* time_vals,
                        const int* linkid_vals,
                        int n_steps,
                        int n_links,
                        int compression_level) {

    int ncid, sys_dimid, time_dimid;
    int sys_varid, time_varid, results_varid;

    // Create file
    NC_CHECK(nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid));

    // Define dimensions
    NC_CHECK(nc_def_dim(ncid, "LinkID", n_links, &sys_dimid));
    NC_CHECK(nc_def_dim(ncid, "time", n_steps, &time_dimid));

    // Define coordinate variables
    NC_CHECK(nc_def_var(ncid, "LinkID", NC_INT, 1, &sys_dimid, &sys_varid));
    NC_CHECK(nc_def_var(ncid, "time", NC_INT, 1, &time_dimid, &time_varid));

    // Define main data variable
    int dims[2] = {time_dimid,sys_dimid};
    NC_CHECK(nc_def_var(ncid, "outputs", NC_FLOAT, 2, dims, &results_varid));

    // Set compression if requested
    if (compression_level > 0) {
        NC_CHECK(nc_def_var_deflate(ncid, results_varid, 1, 1, compression_level));
    }

    // Add attributes
    NC_CHECK(nc_put_att_text(ncid, sys_varid, "long_name", 36, "ID associated with each stream link"));
    NC_CHECK(nc_put_att_text(ncid, time_varid, "long_name", 5, "Time"));
    NC_CHECK(nc_put_att_text(ncid, time_varid, "units", 34, "minutes since start of simulation"));
    NC_CHECK(nc_put_att_text(ncid, results_varid, "long_name", 10, "Discharge"));
    NC_CHECK(nc_put_att_text(ncid, results_varid, "units", 6, "m^3/s"));


    // End define mode
    NC_CHECK(nc_enddef(ncid));

    // Write coordinate variables
    NC_CHECK(nc_put_var_int(ncid, sys_varid, linkid_vals));
    NC_CHECK(nc_put_var_int(ncid, time_varid, time_vals));

    // Write main data
    NC_CHECK(nc_put_var_float(ncid, results_varid, results));

    // Close file
    NC_CHECK(nc_close(ncid));
}





/**
 * @brief Write only the final time step of a results 2D array to a NetCDF file (no time dimension). 
 * 
 */
void write_snapshot_netcdf(const std::string& filename,
                        const float* q_final,
                        const int* linkid_vals,
                        int n_links,
                        int compression_level) {

    // NetCDF identifiers
    int ncid, sys_dimid;
    int sys_varid, final_varid;

    // Create file
    NC_CHECK(nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid));

    // Define dimensions

    NC_CHECK(nc_def_dim(ncid, "LinkID", n_links, &sys_dimid));

    // Define coordinate variables
    NC_CHECK(nc_def_var(ncid, "LinkID", NC_INT, 1, &sys_dimid, &sys_varid));
    NC_CHECK(nc_put_att_text(ncid, sys_varid, "long_name", 36, "ID associated with each stream link"));

    // Define main data variable
    int dims[1] = {sys_dimid};
    NC_CHECK(nc_def_var(ncid, "snapshot", NC_FLOAT, 1, dims, &final_varid));
    NC_CHECK(nc_put_att_text(ncid, final_varid, "long_name", 22, "Final discharge state"));
    NC_CHECK(nc_put_att_text(ncid, final_varid, "units", 6, "m^3/s"));

    // Set compression if requested
    if (compression_level > 0) {
        NC_CHECK(nc_def_var_deflate(ncid, final_varid, 1, 1, compression_level));
    }

    // End define mode
    NC_CHECK(nc_enddef(ncid));

    // Write coordinate variables
    NC_CHECK(nc_put_var_int(ncid, sys_varid, linkid_vals));

    // Write main data
    NC_CHECK(nc_put_var_float(ncid, final_varid, q_final));

    // Close file
    NC_CHECK(nc_close(ncid));
}
