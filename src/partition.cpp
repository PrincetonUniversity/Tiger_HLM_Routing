#include "partition.hpp"

// C++ standard libraries
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

namespace {

// Matches the header written by tools/partition.py.
const char PART_MAGIC[8] = {'H', 'L', 'M', 'P', 'A', 'R', 'T', '1'};

/**
 * @brief Builds the local index space from an already-filled rank_of.
 *
 * Local indices are handed out in increasing global order. That is what makes a
 * single-rank run the identity map, so it stays comparable to a pre-partition run.
 */
void BuildLocalIndices(Partition& part, size_t n_links)
{
    part.local_of.assign(n_links, Partition::NOT_OWNED);
    part.global_of.clear();
    for (size_t i = 0; i < n_links; ++i) {
        if (part.rank_of[i] == part.rank) {
            part.local_of[i] = part.global_of.size();
            part.global_of.push_back(i);
        }
    }
}

} // namespace

Partition LoadPartition(size_t n_links,
                        const std::string& path,
                        int rank,
                        int n_ranks)
{
    Partition part;
    part.rank = rank;
    part.n_ranks = n_ranks;

    if (path.empty()) {
        // No partition file: one rank owns everything and local index == global index.
        if (n_ranks != 1) {
            std::cerr << "Error: " << n_ranks << " ranks but no mpi.partition_file. "
                      << "A distributed run needs a partition; build one with "
                      << "tools/partition.py." << std::endl;
            exit(EXIT_FAILURE);
        }
        part.rank_of.assign(n_links, 0);
        BuildLocalIndices(part, n_links);
        std::cout << "  Partition: single rank, all " << n_links << " links."
                  << std::endl;
        return part;
    }

    std::ifstream file(path, std::ios::binary);
    if (!file) {
        std::cerr << "Error: failed to open partition file: " << path << std::endl;
        exit(EXIT_FAILURE);
    }

    char magic[8];
    file.read(magic, 8);
    if (!file || std::memcmp(magic, PART_MAGIC, 8) != 0) {
        std::cerr << "Error: " << path << " is not a partition file (bad magic)."
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    uint64_t n_links_file = 0, n_ranks_file = 0;
    file.read(reinterpret_cast<char*>(&n_links_file), sizeof(n_links_file));
    file.read(reinterpret_cast<char*>(&n_ranks_file), sizeof(n_ranks_file));
    if (!file) {
        std::cerr << "Error: " << path << " is truncated in its header." << std::endl;
        exit(EXIT_FAILURE);
    }

    // The partition is indexed by link index, so a table of a different size is a
    // different network and every index in the file would mean something else.
    if (n_links_file != n_links) {
        std::cerr << "Error: " << path << " was built for " << n_links_file
                  << " links but the routing table has " << n_links
                  << ". Wrong partition for this network." << std::endl;
        exit(EXIT_FAILURE);
    }
    if (static_cast<int>(n_ranks_file) != n_ranks) {
        std::cerr << "Error: " << path << " was built for " << n_ranks_file
                  << " ranks but the run has " << n_ranks << "." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<int32_t> raw(n_links);
    file.read(reinterpret_cast<char*>(raw.data()),
              static_cast<std::streamsize>(raw.size() * sizeof(int32_t)));
    if (!file) {
        std::cerr << "Error: " << path << " is truncated: expected " << n_links
                  << " rank entries." << std::endl;
        exit(EXIT_FAILURE);
    }

    part.rank_of.assign(n_links, 0);
    for (size_t i = 0; i < n_links; ++i) {
        if (raw[i] < 0 || raw[i] >= n_ranks) {
            std::cerr << "Error: " << path << " assigns link " << i << " to rank "
                      << raw[i] << ", outside 0.." << n_ranks - 1 << "." << std::endl;
            exit(EXIT_FAILURE);
        }
        part.rank_of[i] = raw[i];
    }

    BuildLocalIndices(part, n_links);
    std::cout << "  Partition: rank " << rank << " of " << n_ranks << " owns "
              << part.n_owned() << " of " << n_links << " links ("
              << (100.0 * part.n_owned() / n_links) << "%)." << std::endl;
    return part;
}
// End of file: Tiger_HLM_Routing/src/partition.cpp
