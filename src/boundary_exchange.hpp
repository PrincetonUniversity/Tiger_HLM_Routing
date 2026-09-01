#pragma once

#include <cstddef>
#include <mpi.h>
#include <unordered_map>
#include <vector>

#include "dependency_graph.hpp"
#include "partition.hpp"

/**
 * @brief The links whose discharge has to cross between ranks, and the buffers for them.
 *
 * A cut edge is a link whose downstream link is owned by a different rank. The upstream
 * rank computes the series; the downstream rank needs it to form that child's inflow.
 *
 * partition.py assigns ranks in topological order, so every cut edge runs low rank ->
 * high rank and its validator rejects backward edges. Receiving from all lower ranks
 * then sending to all higher ones therefore cannot deadlock.
 *
 * All series crossing a given rank pair travel in ONE message, packed in increasing global
 * link index. Both sides derive that order from the partition, so no keys or tags per link
 * are needed and the two sides cannot disagree about the layout.
 */
struct BoundaryExchange {
    struct Peer {
        int                 rank = 0;   // the other end
        std::vector<size_t> links;      // global link indices, ascending
        std::vector<float>  buffer;     // links.size() * n_steps
        MPI_Request         request = MPI_REQUEST_NULL;  // in-flight send, lookahead only
    };

    std::vector<Peer> recv_from;  // lower ranks; their links are parents of links I own
    std::vector<Peer> send_to;    // higher ranks; my links are parents of links they own

    // Global link index -> its received series. Only remote parents appear here, so it is
    // small: one entry per cut edge arriving at this rank.
    std::unordered_map<size_t, const float*> arrived;

    // 0: send blocking, so a rank cannot start the next chunk until its downstream
    //    neighbour has taken this one. 1: send non-blocking and wait only when the buffer
    //    is next needed, so a rank runs a chunk ahead of the rank below it.
    int lookahead = 1;
    bool sends_in_flight = false;

    bool empty() const { return recv_from.empty() && send_to.empty(); }
};

/**
 * @brief Works out which links cross which rank boundary.
 *
 * Derived from the partition and the graph, both already built. Nothing is read from disk
 * and no communication happens here.
 *
 * @param part This rank's partition.
 * @param graph The dependency graph, for child[].
 * @return The exchange plan for this rank.
 */
BoundaryExchange BuildBoundaryExchange(const Partition& part, const DependencyGraph& graph);

/**
 * @brief Receives every incoming boundary series for this chunk, before solving starts.
 *
 * Blocking, so a rank waiting on an upstream neighbour sleeps rather than spinning.
 *
 * After this returns, every remote parent's series is in `arrived`, so a link with only
 * remote parents is genuinely ready and can be seeded like any headwater.
 *
 * @param ex The exchange plan; its buffers are filled.
 * @param n_steps Steps in this chunk.
 */
double ReceiveBoundaries(BoundaryExchange& ex, size_t n_steps);

/**
 * @brief Sends this rank's outgoing boundary series once the chunk is solved.
 *
 * With lookahead the send is an MPI_Isend and this returns immediately, so the rank moves
 * on to the next chunk while the message is still in flight. The wait happens here on the
 * NEXT call, before the buffer is refilled -- one chunk of lookahead, which is what fills
 * the pipeline along the rank chain.
 *
 * Without lookahead it is a plain MPI_Send. That is correct but lockstep, and above Intel
 * MPI's eager threshold (~256 KB, reached by a 7-day chunk) MPI_Send blocks until the
 * receiver posts its receive, so the ranks advance one at a time.
 *
 * @param ex The exchange plan.
 * @param part This rank's partition, to find each link's local slice.
 * @param results This rank's solved series, in local index space.
 * @param n_steps Steps in this chunk.
 */
void SendBoundaries(BoundaryExchange& ex,
                    const Partition& part,
                    const std::vector<float>& results,
                    size_t n_steps);

/**
 * @brief Waits for any outstanding non-blocking sends. Call once after the last chunk.
 *
 * MPI_Finalize must not be reached with requests still in flight, and the send buffers
 * must outlive their messages.
 */
void FinishBoundaries(BoundaryExchange& ex);
