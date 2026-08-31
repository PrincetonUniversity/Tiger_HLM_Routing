#pragma once

#include <cstddef>
#include <string>
#include <vector>

/**
 * @brief Which links this rank owns, and where they sit in its local arrays.
 *
 * A rank owns a scattered set of global link indices, so `results` cannot be indexed
 * globally: at CONUS scale that would mean allocating the whole 99.7 GiB array on every
 * rank to use a fraction of it. Every rank therefore keeps a dense local index space
 * 0..n_owned-1, and this maps between the two.
 *
 * Local indices are assigned in increasing global order, so a single-rank run is the
 * identity map and behaves exactly as before.
 */
struct Partition {
    static constexpr size_t NOT_OWNED = static_cast<size_t>(-1);

    int                 rank = 0;
    int                 n_ranks = 1;
    std::vector<int>    rank_of;    // global link index -> owning rank
    std::vector<size_t> local_of;   // global link index -> local index, or NOT_OWNED
    std::vector<size_t> global_of;  // local index -> global link index

    size_t n_owned() const { return global_of.size(); }
    bool   owns(size_t global_index) const { return local_of[global_index] != NOT_OWNED; }
    bool   distributed() const { return n_ranks > 1; }
};

/**
 * @brief Loads a partition file, or builds the single-rank identity partition.
 *
 * An empty path gives one rank owning every link, with local index == global index, which
 * is the behaviour of every run before partitioning existed.
 *
 * The file is produced by tools/partition.py; its format is documented there.
 *
 * @param n_links Number of links in the routing table.
 * @param path Partition file, or "" for the single-rank identity partition.
 * @param rank This process's rank.
 * @param n_ranks Total ranks in the run.
 * @return The partition, or exits on a malformed file or a link-count mismatch.
 */
Partition LoadPartition(size_t n_links,
                        const std::string& path,
                        int rank,
                        int n_ranks);
