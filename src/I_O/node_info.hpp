#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <map>

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
