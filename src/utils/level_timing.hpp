#pragma once

#include <cstddef>
#include <string>
#include <vector>

/**
 * @brief Timing summary for one topological level within one time chunk.
 *
 * The level-synchronous path solves every link at a level concurrently and then waits on
 * the implicit barrier at the end of the parallel loop, so wall covers both the work and
 * the wait. The time the barrier costs is wall * threads - busy_sum.
 */
struct LevelStat {
    size_t chunk;      // index of the runoff time chunk
    size_t level;      // topological level, 0 at the headwaters
    size_t n_links;    // links solved at this level
    int    threads;    // threads that solved at least one link
    double wall;       // wall time for the level, including the closing barrier
    double busy_sum;   // summed time inside link solves, over all threads
    double busy_max;   // time inside link solves for the busiest thread
    double busy_min;   // time inside link solves for the least busy thread that took a link
    double link_min;   // fastest single link solve at this level
    double link_max;   // slowest single link solve at this level
};

/**
 * @brief Per-level timing for the level-synchronous routing path.
 *
 * Separates the three ways the makespan can be spent: useful work (busy_sum), imbalance
 * within a level (busy_max - busy_min), and idle waiting on the barrier (the remainder).
 * Disabled by default, in which case every call is a single predicted branch.
 *
 * Accumulators are one per thread and padded to a cache line, so recording a link never
 * makes two threads write the same line. Without that padding the profiler would
 * manufacture the very contention it exists to measure.
 */
class LevelProfiler {
public:
    // Turns profiling on and sets the CSV path. Rows are appended per chunk so a job
    // that is killed part way still leaves usable data behind.
    void enable(const std::string& filepath);
    bool enabled() const { return enabled_; }

    // Called outside the parallel region, around each level.
    void beginLevel(size_t chunk, size_t level, size_t n_links);
    void endLevel();

    // Called once per link from inside the parallel region.
    void recordLink(double seconds);

    // Appends the levels recorded since the last call, then clears them.
    void writeChunk();

    // Whole-run totals: parallel efficiency, and how much of the wall time went to
    // levels too narrow to fill the machine.
    void printSummary() const;

private:
    struct alignas(64) ThreadAccum {
        double busy;
        double link_min;
        double link_max;
        size_t n_links;
    };

    bool                     enabled_        = false;
    bool                     header_written_ = false;
    std::string              filepath_;
    std::vector<ThreadAccum> accum_;
    std::vector<LevelStat>   stats_;

    size_t chunk_       = 0;
    size_t level_       = 0;
    size_t n_links_     = 0;
    double level_start_ = 0.0;

    // Running totals, kept separately because stats_ is cleared as it is written out.
    double total_wall_    = 0.0;
    double total_busy_    = 0.0;
    double narrow_wall_   = 0.0;  // wall spent in levels holding fewer links than threads
    size_t narrow_levels_ = 0;
    size_t total_levels_  = 0;
    size_t total_links_   = 0;
    int    max_threads_   = 0;
};
