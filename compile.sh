#!/bin/bash

set -e

module purge
module load intel/2024.2
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load hdf5/oneapi-2024.2/1.14.4
module load netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2

# CPU build
make clean
make
mv bin/routing routing_cpu.tmp  # so make clean doesn't remove

# GPU build
module load cudatoolkit/12.9
make clean
make USE_GPU=1
mv bin/routing bin/routing_gpu
mv routing_cpu.tmp bin/routing_cpu