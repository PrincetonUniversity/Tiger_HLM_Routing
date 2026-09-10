#include "level_timing.hpp"

// C++ standard libraries
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>

//--------------------------------------------------------------------------------------------------
// Function Definitions
//--------------------------------------------------------------------------------------------------
/**
 * @brief Turns profiling on and sets the CSV path.
 * Sizes one padded accumulator per thread up front, so recordLink never allocates.
 * @param filepath Path to write the per-level CSV to.
 */
void LevelProfiler::enable(const std::string& filepath)
{
    enabled_     = true;
    filepath_    = filepath;
    max_threads_ = omp_get_max_threads();
    accum_.assign(static_cast<size_t>(max_threads_), ThreadAccum{});
}

/**
 * @brief Starts timing a level. Called outside the parallel region.
 * @param chunk The current time chunk index.
 * @param level The topological level about to be solved.
 * @param n_links The number of links at this level.
 */
void LevelProfiler::beginLevel(size_t chunk, size_t level, size_t n_links)
{
    if (!enabled_) return;

    chunk_   = chunk;
    level_   = level;
    n_links_ = n_links;

    for (auto& a : accum_) {
        a.busy     = 0.0;
        a.link_min = std::numeric_limits<double>::max();
        a.link_max = 0.0;
        a.n_links  = 0;
    }

    level_start_ = omp_get_wtime();
}

/**
 * @brief Records one link solve. Called from inside the parallel region.
 * Writes only to this thread's own cache line.
 * @param seconds Wall time the link's integration took.
 */
void LevelProfiler::recordLink(double seconds)
{
    if (!enabled_) return;

    ThreadAccum& a = accum_[static_cast<size_t>(omp_get_thread_num())];
    a.busy += seconds;
    a.link_min = std::min(a.link_min, seconds);
    a.link_max = std::max(a.link_max, seconds);
    ++a.n_links;
}

/**
 * @brief Closes out a level and folds its numbers into the running totals.
 * Called outside the parallel region, after the implicit barrier, so wall includes
 * the barrier wait.
 */
void LevelProfiler::endLevel()
{
    if (!enabled_) return;

    LevelStat s{};
    s.chunk    = chunk_;
    s.level    = level_;
    s.n_links  = n_links_;
    s.wall     = omp_get_wtime() - level_start_;
    s.busy_sum = 0.0;
    s.busy_max = 0.0;
    s.busy_min = std::numeric_limits<double>::max();
    s.link_min = std::numeric_limits<double>::max();
    s.link_max = 0.0;
    s.threads  = 0;

    for (const auto& a : accum_) {
        if (a.n_links == 0) continue;  // thread took no link at this level
        ++s.threads;
        s.busy_sum += a.busy;
        s.busy_max = std::max(s.busy_max, a.busy);
        s.busy_min = std::min(s.busy_min, a.busy);
        s.link_min = std::min(s.link_min, a.link_min);
        s.link_max = std::max(s.link_max, a.link_max);
    }
    if (s.threads == 0) {  // empty level, leave the mins as zero rather than max double
        s.busy_min = 0.0;
        s.link_min = 0.0;
    }

    total_wall_ += s.wall;
    total_busy_ += s.busy_sum;
    total_links_ += s.n_links;
    ++total_levels_;
    if (s.n_links < static_cast<size_t>(max_threads_)) {
        narrow_wall_ += s.wall;
        ++narrow_levels_;
    }

    stats_.push_back(s);
}

/**
 * @brief Appends the levels recorded since the last call, then clears them.
 * Appending per chunk keeps memory flat over long runs and leaves usable data behind
 * if the job is killed.
 */
void LevelProfiler::writeChunk()
{
    if (!enabled_ || stats_.empty()) return;

    std::ofstream out(filepath_, header_written_ ? std::ios::app : std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "Warning: could not open profiling file " << filepath_
                  << ", per-level timing not written." << std::endl;
        stats_.clear();
        return;
    }

    if (!header_written_) {
        out << "chunk,level,n_links,threads,wall,busy_sum,busy_max,busy_min,"
               "link_min,link_max\n";
        header_written_ = true;
    }

    out << std::setprecision(9);
    for (const auto& s : stats_) {
        out << s.chunk << ',' << s.level << ',' << s.n_links << ',' << s.threads << ','
            << s.wall << ',' << s.busy_sum << ',' << s.busy_max << ',' << s.busy_min
            << ',' << s.link_min << ',' << s.link_max << '\n';
    }

    stats_.clear();
}

/**
 * @brief Prints whole-run totals: where the makespan went.
 * The three numbers that matter are the parallel efficiency, the share of wall time
 * spent in levels too narrow to fill the machine, and the imbalance within levels.
 */
void LevelProfiler::printSummary() const
{
    if (!enabled_) return;

    const double capacity = total_wall_ * static_cast<double>(max_threads_);
    const double idle     = capacity - total_busy_;

    std::cout << "\n________________LEVEL TIMING SUMMARY______________ \n" << std::endl;
    std::cout << "  Threads:                  " << max_threads_ << std::endl;
    std::cout << "  Levels solved:            " << total_levels_ << std::endl;
    std::cout << "  Link solves:              " << total_links_ << std::endl;
    std::cout << "  Wall in integration:      " << total_wall_ << " s" << std::endl;
    std::cout << "  Useful work:              " << total_busy_ << " CPU-s" << std::endl;

    if (capacity > 0.0) {
        std::cout << "  Idle:                     " << idle << " CPU-s ("
                  << 100.0 * idle / capacity << "%)" << std::endl;
        std::cout << "  Parallel efficiency:      " << 100.0 * total_busy_ / capacity
                  << "%" << std::endl;
    }
    if (total_wall_ > 0.0) {
        std::cout << "  Levels below " << max_threads_ << " links:     " << narrow_levels_
                  << " of " << total_levels_ << ", holding "
                  << 100.0 * narrow_wall_ / total_wall_ << "% of the wall time"
                  << std::endl;
    }
    std::cout << "\n  Per-level detail written to: " << filepath_ << std::endl;
    std::cout << "__________________________________________________ \n" << std::endl;
}
