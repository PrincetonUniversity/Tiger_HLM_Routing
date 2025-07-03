#pragma once

// C++ standard libraries
#include <unordered_map>
#include <map>
#include <vector>
#include <functional>
#include <string>

// External libraries
#include "I_O/node_info.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "I_O/config_loader.hpp"


struct ModelSetup {
    ModelConfig config;
    std::unordered_map<size_t, NodeInfo> node_map;
    std::map<size_t, std::vector<size_t>> level_groups;
    size_t n_links = 0;
    std::function<float(int)> uini;
    BoundaryConditions boundary_conditions;
    RunoffChunkInfo runoff_info;
    SaveInfo save_info;
};

ModelSetup setupModel(const char* config_path);