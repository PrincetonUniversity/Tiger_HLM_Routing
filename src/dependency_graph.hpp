#pragma once

#include <atomic>
#include <cstddef>
#include <vector>

#include "model_setup.hpp"
#include "I_O/inputs.hpp"

/**
 * @brief The routing network expressed as what each link is waiting on.
 * Built once. The network is fixed for the run, so only the countdown is reset per chunk.
 */
struct DependencyGraph {
    static constexpr size_t NO_CHILD = static_cast<size_t>(-1);

    std::vector<size_t> child;      // the one downstream link, NO_CHILD at an outlet
    std::vector<int>    in_degree;  // upstream links that must finish first
    std::vector<size_t> sources;    // in_degree 0, i.e. the headwaters: the seed set
    size_t              outlets = 0;
};

/**
 * @brief Everything a link task needs, bundled so the task clause stays readable.
 * Lives in IntegrateLinksByDependency's frame, so it outlives every task referring to it.
 * Owns nothing; pending is an interior pointer, as in RHS::runoff_series.
 */
struct TaskContext {
    const ModelSetup&       setup;
    const RunoffData&       runoff;
    std::vector<float>&     results;
    const DependencyGraph&  graph;
    std::atomic<int>*       pending;
    std::vector<float>&     q_final;
    const size_t            n_steps;
    const size_t            total_time_steps;
    const size_t            tc;
};

/**
 * @brief Inverts the parent lists into the child/in-degree form the countdown needs.
 * Derived from node_map; no new input file.
 *
 * @param setup The model setup containing the node map.
 * @return The dependency graph for the whole network.
 */
DependencyGraph BuildDependencyGraph(const ModelSetup& setup);
