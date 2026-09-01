#pragma once
#include <vector>
#include "model_setup.hpp"
#include "I_O/inputs.hpp"

void IntegrateLevel0GPU(const ModelSetup& setup,
                        const RunoffData& runoff,
                        std::vector<float>& results,
                        const std::vector<size_t>& nodes_at_level,
                        size_t n_steps,
                        size_t tc,
                        std::vector<float>& q_final);