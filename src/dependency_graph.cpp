#include "dependency_graph.hpp"

// C++ standard libraries
#include <cstdlib>
#include <iostream>

/**
 * @brief Inverts the parent lists into the child/in-degree form the countdown needs.
 * Derived from node_map; no new input file.
 *
 * child[] is one index, not a list, because every link drains to exactly one downstream
 * link. Checked here rather than assumed: a link with two children would lose a release
 * and deadlock the traversal.
 *
 * @param setup The model setup containing the node map.
 * @return The dependency graph for the whole network.
 */
DependencyGraph BuildDependencyGraph(const ModelSetup& setup)
{
    DependencyGraph graph;
    graph.child.assign(setup.n_links, DependencyGraph::NO_CHILD);
    graph.in_degree.assign(setup.n_links, 0);

    for (const auto& [id, node] : setup.node_map) {
        if (node.index >= setup.n_links) {
            std::cerr << "Error: link index " << node.index << " is outside 0.."
                      << setup.n_links - 1 << ". The routing table is not densely indexed."
                      << std::endl;
            exit(EXIT_FAILURE);
        }
        graph.in_degree[node.index] = static_cast<int>(node.parents.size());

        for (size_t parent_index : node.parents) {
            if (parent_index >= setup.n_links) {
                std::cerr << "Error: link " << node.index << " lists parent " << parent_index
                          << ", which is outside 0.." << setup.n_links - 1 << "." << std::endl;
                exit(EXIT_FAILURE);
            }
            if (graph.child[parent_index] != DependencyGraph::NO_CHILD) {
                std::cerr << "Error: link " << parent_index << " drains to both "
                          << graph.child[parent_index] << " and " << node.index
                          << ". The dependency traversal requires one downstream link per link."
                          << std::endl;
                exit(EXIT_FAILURE);
            }
            graph.child[parent_index] = node.index;
        }
    }

    for (size_t i = 0; i < setup.n_links; ++i) {
        if (graph.in_degree[i] == 0) graph.sources.push_back(i);
        if (graph.child[i] == DependencyGraph::NO_CHILD) ++graph.outlets;
    }

    std::cout << "  Dependency graph: " << graph.sources.size() << " headwaters, "
              << graph.outlets << " outlets." << std::endl;
    return graph;
}
// End of file: Tiger_HLM_Routing/src/dependency_graph.cpp
