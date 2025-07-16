#pragma once

/**
 * @brief Sets up OpenMP by checking the environment variable OMP_NUM_THREADS.
 * If the variable is set, it uses that value; otherwise, it defaults to 1 thread.
 * It also provides instructions for setting the variable in a SLURM cluster environment.
 */
void setupOpenMP();