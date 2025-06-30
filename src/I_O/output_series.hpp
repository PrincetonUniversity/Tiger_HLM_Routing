#pragma once

#include <netcdf.h>
#include <string>
#include <iostream>
#include <vector>


/**
 * @brief Write a streamflow array to a NetCDF file with optional compression.
 * 
 * @param filename The name of the output NetCDF file.
 * @param results Pointer to the 2D array of results (time steps x links).
 * @param time_vals Pointer to the 1D array of time values.
 * @param linkid_vals Pointer to the 1D array of link IDs.
 * @param n_steps Number of time steps in the results array.
 * @param n_links Number of links in the results array.
 * @param compression_level Compression level for the NetCDF file (0 for no compression).
 */
void write_timeseries_netcdf(const std::string& filename,
                        const float* results,
                        const int* time_vals,
                        const int* linkid_vals,
                        int n_steps,
                        int n_links,
                        int compression_level = 0);

/**
 * @brief Write only the final time step of a results 2D array to a NetCDF file (no time dimension).
 * 
 * @param filename The name of the output NetCDF file.
 * @param q_final Pointer to the 1D array of final discharge values.
 * @param linkid_vals Pointer to the 1D array of link IDs.
 * @param n_links Number of links in the results array.
 * @param compression_level Compression level for the NetCDF file (0 for no compression).       
 * */
void write_snapshot_netcdf(const std::string& filename,
                        const float* q_final,
                        const int* linkid_vals,         
                        int n_links,
                        int compression_level = 0);