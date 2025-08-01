#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <map>

/** 
 * @struct NodeInfo
 * @brief Structure to hold information about a node in the network.
 * 
 * This structure contains:
 * - index: Unique identifier for the node.
 * - stream_id: Identifier for the stream associated with the node.
 * - level: The level of the node in the hierarchy.
 * - parents: A vector of indices representing parent nodes.
 * - params: A vector of parameters associated with the node.
 */
struct NodeInfo {
    size_t index;
    int stream_id;
    size_t level;
    std::vector<size_t> parents;
    std::vector<double> params;
};

// Function declaration
void read_node_levels(
    const std::string& filename,
    std::unordered_map<size_t, NodeInfo>& node_map,
    std::map<size_t, std::vector<size_t>>& level_groups
);
