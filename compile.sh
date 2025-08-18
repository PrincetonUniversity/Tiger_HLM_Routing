#!/bin/bash

# Load environment modules 
module purge
module load intel/2024.2
module load intel-oneapi/2024.2    # Optional, but some HDF5/NetCDF modules require it
module load hdf5/oneapi-2024.2/1.14.4   
module load netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# Build the project
make clean
make