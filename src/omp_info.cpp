#include "omp_info.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <omp.h>  // Required here because omp_set_num_threads is used

/**
 * @brief Sets up OpenMP by checking the environment variable OMP_NUM_THREADS.
 * If the variable is set, it uses that value; otherwise, it defaults to 1 thread.
 * It also provides instructions for setting the variable in a SLURM cluster environment.
 */
void setupOpenMP() {
    std::cout << "___________________OpenMP INFO____________________ \n" << std::endl;

    const char* env_threads = std::getenv("OMP_NUM_THREADS");
    if (env_threads) {
        std::cout << "OMP_NUM_THREADS is set to " << env_threads << "! \n" << std::endl;
    } else {
        omp_set_num_threads(1);
        std::cout << "OMP_NUM_THREADS not set — defaulting to 1 thread." << std::endl;
        std::cout << "If you are running this on a SLURM cluster, please set the number of threads using the SLURM_CPUS_PER_TASK variable." << std::endl;
        std::cout << "For example, you can set it in your SLURM script as follows: export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK" << std::endl;
    }

    std::cout << "__________________________________________________ \n" << std::endl;
}

/**
 * @brief Sets the schedule used by the link loop, which is declared schedule(runtime).
 * Lets one binary compare static, dynamic and guided without a rebuild, so that
 * intra-level load imbalance can be separated from barrier idle.
 */
void applyOmpSchedule(const std::string& kind, int chunk)
{
    omp_sched_t sched = omp_sched_static;
    if (kind == "dynamic") {
        sched = omp_sched_dynamic;
    } else if (kind == "guided") {
        sched = omp_sched_guided;
    } else if (kind != "static") {
        std::cerr << "Warning: unknown profiling.omp_schedule '" << kind
                  << "'. Falling back to static." << std::endl;
    }

    omp_set_schedule(sched, chunk);
}
