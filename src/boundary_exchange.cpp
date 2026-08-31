#include "boundary_exchange.hpp"

// C++ standard libraries
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>

#include <mpi.h>

namespace {

// One message per rank pair, so a single tag suffices.
constexpr int BOUNDARY_TAG = 7001;

} // namespace

BoundaryExchange BuildBoundaryExchange(const Partition& part, const DependencyGraph& graph)
{
    BoundaryExchange ex;
    if (!part.distributed()) return ex;   // nothing crosses a boundary in a 1-rank run

    // std::map keeps peers in ascending rank order, which is what makes the
    // receive-from-lower-then-send-to-higher ordering below safe.
    std::map<int, std::vector<size_t>> incoming, outgoing;

    for (size_t i = 0; i < part.rank_of.size(); ++i) {
        const size_t child = graph.child[i];
        if (child == DependencyGraph::NO_CHILD) continue;
        const int owner_i = part.rank_of[i];
        const int owner_c = part.rank_of[child];
        if (owner_i == owner_c) continue;   // wholly inside one rank

        if (owner_i == part.rank) {
            outgoing[owner_c].push_back(i);         // I compute it, they need it
        } else if (owner_c == part.rank) {
            incoming[owner_i].push_back(i);         // they compute it, I need it
        }
        // edges between two other ranks are not this rank's concern
    }

    // The partitioner orders ranks topologically, so a cut edge always runs from a lower
    // rank to a higher one. If that ever fails the exchange order below would deadlock,
    // so check rather than assume.
    for (const auto& [peer, links] : incoming) {
        if (peer >= part.rank) {
            std::cerr << "Error: rank " << part.rank << " expects input from rank " << peer
                      << ", which is not upstream of it. The partition's rank graph is not "
                      << "topologically ordered and the exchange would deadlock." << std::endl;
            exit(EXIT_FAILURE);
        }
        BoundaryExchange::Peer p;
        p.rank = peer;
        p.links = links;
        std::sort(p.links.begin(), p.links.end());
        ex.recv_from.push_back(std::move(p));
    }
    for (const auto& [peer, links] : outgoing) {
        if (peer <= part.rank) {
            std::cerr << "Error: rank " << part.rank << " would send to rank " << peer
                      << ", which is not downstream of it. The partition's rank graph is not "
                      << "topologically ordered and the exchange would deadlock." << std::endl;
            exit(EXIT_FAILURE);
        }
        BoundaryExchange::Peer p;
        p.rank = peer;
        p.links = links;
        std::sort(p.links.begin(), p.links.end());
        ex.send_to.push_back(std::move(p));
    }

    size_t n_in = 0, n_out = 0;
    for (const auto& p : ex.recv_from) n_in += p.links.size();
    for (const auto& p : ex.send_to)   n_out += p.links.size();
    std::cout << "  Boundary exchange: rank " << part.rank << " receives " << n_in
              << " series from " << ex.recv_from.size() << " rank(s), sends " << n_out
              << " to " << ex.send_to.size() << "." << std::endl;
    return ex;
}

void ReceiveBoundaries(BoundaryExchange& ex, size_t n_steps)
{
    ex.arrived.clear();
    for (auto& peer : ex.recv_from) {
        peer.buffer.resize(peer.links.size() * n_steps);
        MPI_Recv(peer.buffer.data(), static_cast<int>(peer.buffer.size()), MPI_FLOAT,
                 peer.rank, BOUNDARY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Both sides pack in ascending global link index, so position determines identity.
        for (size_t k = 0; k < peer.links.size(); ++k) {
            ex.arrived[peer.links[k]] = peer.buffer.data() + k * n_steps;
        }
    }
}

void SendBoundaries(BoundaryExchange& ex,
                    const Partition& part,
                    const std::vector<float>& results,
                    size_t n_steps)
{
    // The previous chunk's messages must have left before their buffers are refilled.
    // Waiting here rather than at the send is the whole point: the rank spent the
    // intervening chunk solving, so by now the transfer has almost certainly finished.
    if (ex.sends_in_flight) {
        for (auto& peer : ex.send_to) {
            MPI_Wait(&peer.request, MPI_STATUS_IGNORE);
        }
        ex.sends_in_flight = false;
    }

    for (auto& peer : ex.send_to) {
        peer.buffer.resize(peer.links.size() * n_steps);
        for (size_t k = 0; k < peer.links.size(); ++k) {
            const size_t local = part.local_of[peer.links[k]];
            std::copy(results.begin() + static_cast<std::ptrdiff_t>(local * n_steps),
                      results.begin() + static_cast<std::ptrdiff_t>((local + 1) * n_steps),
                      peer.buffer.begin() + static_cast<std::ptrdiff_t>(k * n_steps));
        }
        if (ex.lookahead > 0) {
            MPI_Isend(peer.buffer.data(), static_cast<int>(peer.buffer.size()), MPI_FLOAT,
                      peer.rank, BOUNDARY_TAG, MPI_COMM_WORLD, &peer.request);
        } else {
            MPI_Send(peer.buffer.data(), static_cast<int>(peer.buffer.size()), MPI_FLOAT,
                     peer.rank, BOUNDARY_TAG, MPI_COMM_WORLD);
        }
    }
    if (ex.lookahead > 0) ex.sends_in_flight = true;
}

void FinishBoundaries(BoundaryExchange& ex)
{
    if (!ex.sends_in_flight) return;
    for (auto& peer : ex.send_to) {
        MPI_Wait(&peer.request, MPI_STATUS_IGNORE);
    }
    ex.sends_in_flight = false;
}
// End of file: Tiger_HLM_Routing/src/boundary_exchange.cpp
