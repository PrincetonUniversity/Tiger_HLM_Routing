#!/bin/bash
#SBATCH --job-name=profile_routing
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/gpfs/GVILLARI/mb6477/Tiger_HLM_Routing/profile_%j.log

module purge
module load intel/2024.2
module load intel-oneapi/2024.2
module load hdf5/oneapi-2024.2/1.14.4
module load netcdf/oneapi-2024.2/hdf5-1.14.4/4.9.2
module load intel-advisor/oneapi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd /scratch/gpfs/GVILLARI/mb6477/Tiger_HLM_Routing

advisor --collect=survey --project-dir=./profiler_outputs -- \
    ./bin/routing /scratch/gpfs/GVILLARI/mb6477/AORC_1km_CONUS/Tiger_HLM_tutorial/config_routing_tutorial.yaml
