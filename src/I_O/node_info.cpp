#include "node_info.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <stdexcept>

void read_node_levels(
    const std::string& filename,
    std::unordered_map<size_t,NodeInfo>& node_map,
    std::map<size_t,std::vector<size_t>>& level_groups
) {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Failed to open: " + filename);

    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string token;

        NodeInfo node;

        // index
        std::getline(ss, token, ',');
        node.index = std::stoul(token);

        // stream_id
        std::getline(ss, token, ',');
        node.stream_id = std::stoi(token);

        // level
        std::getline(ss, token, ',');
        node.level = std::stoul(token);

        // parents (semicolon separated)
        std::getline(ss, token, ',');
        if (!token.empty()) {
            std::istringstream ps(token);
            while (std::getline(ps, token, ';')) {
                if (!token.empty())
                    node.parents.push_back(std::stoul(token));
            }
        }

        // params (semicolon separated)
        std::getline(ss, token);
        if (!token.empty()) {
            std::istringstream ps(token);
            while (std::getline(ps, token, ';')) {
                if (!token.empty())
                    node.params.push_back(std::stod(token));
            }
        }

        node_map[node.index] = node;
        level_groups[node.level].push_back(node.index);
    }
}
