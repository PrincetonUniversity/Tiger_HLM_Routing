#!/bin/bash

module purge
module load intel/2024.2
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load hdf5/oneapi-2024.2/1.14.4
module load netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# Uncomment for GPU build:
module load cudatoolkit/12.9
make clean
make USE_GPU=1

# CPU-only build (default):
# make clean
# make