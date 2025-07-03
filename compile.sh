#!/bin/bash

# Load environment modules 
module load hdf5/gcc/1.14.4
module load netcdf/gcc/hdf5-1.14.4/4.9.2

# Build the project
make clean
make