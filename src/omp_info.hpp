#pragma once

#include <string>

/**
 * @brief Sets up OpenMP by checking the environment variable OMP_NUM_THREADS.
 * If the variable is set, it uses that value; otherwise, it defaults to 1 thread.
 * It also provides instructions for setting the variable in a SLURM cluster environment.
 */
void setupOpenMP();

/**
 * @brief Sets the schedule used by the link loop, which is declared schedule(runtime).
 * Lets one binary compare static, dynamic and guided without a rebuild, so that
 * intra-level load imbalance can be separated from barrier idle.
 * @param kind "static", "dynamic" or "guided". Anything else falls back to static.
 * @param chunk Chunk size for that schedule; 0 selects the OpenMP default.
 */
void applyOmpSchedule(const std::string& kind, int chunk);