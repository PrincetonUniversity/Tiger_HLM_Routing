#include "boundary_exchange.hpp"

// C++ standard libraries
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <map>
#include <climits>

#include <mpi.h>

namespace {

// One message per rank pair, so a single tag suffices until the message has to be split.
constexpr int BOUNDARY_TAG = 7001;

// Messages longer than an MPI int count are sent in pieces, one tag each
constexpr size_t MAX_MPI_COUNT = static_cast<size_t>(INT_MAX);

void SendChunked(const float* data, size_t count, int peer_rank, int tag_base)
{
    size_t offset = 0;
    int tag = BOUNDARY_TAG + tag_base;
    while (offset < count) {
        const size_t piece = std::min(count - offset, MAX_MPI_COUNT);
        MPI_Send(data + offset, static_cast<int>(piece), MPI_FLOAT,
                 peer_rank, tag, MPI_COMM_WORLD);
        offset += piece;
        ++tag;
    }
}

void RecvChunked(float* data, size_t count, int peer_rank, int tag_base)
{
    size_t offset = 0;
    int tag = BOUNDARY_TAG + tag_base;
    while (offset < count) {
        const size_t piece = std::min(count - offset, MAX_MPI_COUNT);
        MPI_Recv(data + offset, static_cast<int>(piece), MPI_FLOAT,
                 peer_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        offset += piece;
        ++tag;
    }
}

// Non-blocking forms of the two above: every piece is posted at once and its request
// appended to `requests`, which the caller waits on before the data may be touched.
void IsendChunked(const float* data, size_t count, int peer_rank, int tag_base,
                  std::vector<MPI_Request>& requests)
{
    size_t offset = 0;
    int tag = BOUNDARY_TAG + tag_base;
    while (offset < count) {
        const size_t piece = std::min(count - offset, MAX_MPI_COUNT);
        requests.emplace_back();
        MPI_Isend(data + offset, static_cast<int>(piece), MPI_FLOAT,
                  peer_rank, tag, MPI_COMM_WORLD, &requests.back());
        offset += piece;
        ++tag;
    }
}

void IrecvChunked(float* data, size_t count, int peer_rank, int tag_base,
                  std::vector<MPI_Request>& requests)
{
    size_t offset = 0;
    int tag = BOUNDARY_TAG + tag_base;
    while (offset < count) {
        const size_t piece = std::min(count - offset, MAX_MPI_COUNT);
        requests.emplace_back();
        MPI_Irecv(data + offset, static_cast<int>(piece), MPI_FLOAT,
                  peer_rank, tag, MPI_COMM_WORLD, &requests.back());
        offset += piece;
        ++tag;
    }
}

void WaitPieces(std::vector<MPI_Request>& requests)
{
    if (requests.empty()) return;
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    requests.clear();
}

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
        // blocking send reads this group from `results`,
        // GroupLocalIndicesByDestination ensures contiguity and index order
        p.first_local = part.local_of[p.links[0]];
        for (size_t k = 0; k < p.links.size(); ++k) {
            if (part.local_of[p.links[k]] != p.first_local + k) {
                std::cerr << "Error: rank " << part.rank << "'s links for rank " << peer
                          << " are not contiguous in local index order. Local indices must be "
                          << "grouped by destination before the exchange is built." << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        ex.send_to.push_back(std::move(p));
    }

    size_t n_in = 0, n_out = 0;
    for (const auto& p : ex.recv_from) n_in += p.links.size();
    for (const auto& p : ex.send_to)   n_out += p.links.size();
    // if sending most of a rank's links, skip copy and double-buffer `results`
    ex.direct_send = n_out * 2 > part.n_owned();
    std::cout << "  Boundary exchange: rank " << part.rank << " receives " << n_in
              << " series from " << ex.recv_from.size() << " rank(s), sends " << n_out
              << " to " << ex.send_to.size() << "."
              << (ex.direct_send ? " Sends straight from a double-buffered results." : "")
              << std::endl;
    return ex;
}

int RankGraphDepth(const Partition& part, const DependencyGraph& graph)
{
    if (!part.distributed()) return 1;

    // Collect distinct rank pairs. Repeats do not deepen the chain and would grow the
    // next pass quadratically.
    std::vector<std::vector<int>> succ(part.n_ranks);
    {
        std::vector<std::vector<bool>> seen(part.n_ranks, std::vector<bool>(part.n_ranks, false));
        for (size_t i = 0; i < part.rank_of.size(); ++i) {
            const size_t child = graph.child[i];
            if (child == DependencyGraph::NO_CHILD) continue;
            const int u = part.rank_of[i];
            const int v = part.rank_of[child];
            if (u == v || seen[u][v]) continue;
            seen[u][v] = true;
            succ[u].push_back(v);
        }
    }

    // Ranks are topological (partition.py guarantees it, BuildBoundaryExchange checks),
    // so one forward sweep suffices.
    std::vector<int> longest(part.n_ranks, 1);
    int best = 1;
    for (int u = 0; u < part.n_ranks; ++u) {
        for (int v : succ[u]) {
            if (longest[u] + 1 > longest[v]) {
                longest[v] = longest[u] + 1;
                if (longest[v] > best) best = longest[v];
            }
        }
    }
    return best;
}

namespace {

// With lookahead two chunks can be in flight per rank pair, so each parity gets its own
// tag range. Every send and receive path derives its tags here.
constexpr int PIECE_STRIDE = 4096;   // pieces per chunk; enough for a year of CONUS level 0

int TagBase(const BoundaryExchange& ex, size_t chunk)
{
    return ex.lookahead > 0 ? static_cast<int>(chunk % 2) * PIECE_STRIDE : 0;
}

// Posts the receives for `chunk` into that parity's buffer without waiting.
void PostReceives(BoundaryExchange& ex, size_t n_steps, size_t chunk)
{
    const int parity = static_cast<int>(chunk % 2);
    for (auto& peer : ex.recv_from) {
        if (peer.buffers.size() < 2) peer.buffers.resize(2);
        std::vector<float>& buffer = peer.buffers[parity];
        buffer.resize(peer.links.size() * n_steps);
        IrecvChunked(buffer.data(), buffer.size(), peer.rank, TagBase(ex, chunk),
                     peer.piece_requests[parity]);
    }
}

} // namespace

// Returns the seconds spent blocked
double ReceiveBoundaries(BoundaryExchange& ex, size_t n_steps, size_t chunk, size_t next_n_steps)
{
    const auto wait_start = std::chrono::high_resolution_clock::now();
    ex.arrived.clear();
    const int parity = static_cast<int>(chunk % 2);

    if (ex.lookahead > 0) {
        // Wait for this chunk's receives, then post the next chunk's so it arrives during the solve
        if (chunk == 0) PostReceives(ex, n_steps, chunk);
        for (auto& peer : ex.recv_from) WaitPieces(peer.piece_requests[parity]);
        if (next_n_steps > 0) PostReceives(ex, next_n_steps, chunk + 1);
    } else {
        for (auto& peer : ex.recv_from) {
            if (peer.buffers.empty()) peer.buffers.resize(1);
            std::vector<float>& buffer = peer.buffers[0];
            buffer.resize(peer.links.size() * n_steps);
            RecvChunked(buffer.data(), buffer.size(), peer.rank, 0);
        }
    }

    for (auto& peer : ex.recv_from) {
        const std::vector<float>& buffer = peer.buffers[ex.lookahead > 0 ? parity : 0];
        // Both sides pack in ascending global link index, so position determines identity.
        for (size_t k = 0; k < peer.links.size(); ++k) {
            ex.arrived[peer.links[k]] = buffer.data() + k * n_steps;
        }
    }
    const std::chrono::duration<double> waited =
        std::chrono::high_resolution_clock::now() - wait_start;
    return waited.count();
}

void WaitSends(BoundaryExchange& ex, int parity)
{
    for (auto& peer : ex.send_to) WaitPieces(peer.piece_requests[parity]);
}

void SendBoundaries(BoundaryExchange& ex,
                    const std::vector<float>& results,
                    size_t n_steps,
                    size_t chunk)
{
    const size_t slots = static_cast<size_t>(std::max(1, ex.lookahead));
    const size_t slot = chunk % slots;
    const int tag_base = TagBase(ex, chunk);

    for (auto& peer : ex.send_to) {
        const size_t count = peer.links.size() * n_steps;
        const float* group = results.data() + peer.first_local * n_steps;

        if (ex.double_buffered()) {
            // Doesn't block, no copy
            IsendChunked(group, count, peer.rank, tag_base, peer.piece_requests[chunk % 2]);
        } else if (ex.lookahead > 0 && count <= MAX_MPI_COUNT) {
            // Doesn't block: copy, since the next chunk overwrites `results` before this send completes
            if (peer.buffers.empty()) {
                peer.buffers.resize(slots);
                peer.requests.assign(slots, MPI_REQUEST_NULL);
            }
            MPI_Wait(&peer.requests[slot], MPI_STATUS_IGNORE);   // reclaim this slot
            std::vector<float>& buffer = peer.buffers[slot];
            buffer.assign(group, group + count);
            MPI_Isend(buffer.data(), static_cast<int>(count), MPI_FLOAT,
                      peer.rank, BOUNDARY_TAG + tag_base, MPI_COMM_WORLD, &peer.requests[slot]);
        } else {
            // blocks: send directly from `results`, which nothing overwrites until this returns
            if (ex.lookahead > 0 && !peer.oversize_warned) {
                std::cerr << "Warning: boundary message to rank " << peer.rank << " is "
                          << count << " floats, past the " << MAX_MPI_COUNT
                          << " MPI can count. Sending it blocking, so lookahead is off"
                          << " for this peer." << std::endl;
                peer.oversize_warned = true;
            }
            SendChunked(group, count, peer.rank, tag_base);
        }
    }
}

void FinishBoundaries(BoundaryExchange& ex)
{
    for (auto& peer : ex.send_to) {
        for (auto& request : peer.requests) {
            MPI_Wait(&request, MPI_STATUS_IGNORE);
        }
        WaitPieces(peer.piece_requests[0]);
        WaitPieces(peer.piece_requests[1]);
    }
    // No receives should be outstanding here
    for (auto& peer : ex.recv_from) {
        WaitPieces(peer.piece_requests[0]);
        WaitPieces(peer.piece_requests[1]);
    }
}
// End of file: Tiger_HLM_Routing/src/boundary_exchange.cpp
