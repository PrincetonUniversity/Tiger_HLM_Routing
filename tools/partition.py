#!/usr/bin/env python3
"""
Partition a Tiger-HLM streamflow routing network across MPI ranks.

The partitions follow parent-link subtrees, allowing MPI ranks to run independently.
Ranks own connected pieces: a piece is a subtree minus any subtrees already cut
away upstream of it -- pieces tile the network and each has exactly one outlet.

ALGORITHM
    decompose   Accumulate subtree size in routing parameter file order.
                When a link's accumulated size reaches the target,
                cut the edge below it: this becomes the piece's outlet.
                Small outlet trees that never reach the target stay whole.
    order       Ensure topological order by sorting pieces by outlet index.
    assign      Split the sorted list of pieces into R contiguous runs of roughly equal link count.

    Known limitation:
    Balancing link count does not necessarily balance wall-clock time among ranks,
    so a rank holding mostly headwaters finishes early and idles.
    This is because when rk4_level: 1 in the config,
    levels 0 and 1 are solved with RK4 rather than RK45.

USAGE
    partition.py build     <routing_table.csv> <n_ranks> <out.part> [granularity=16]
    partition.py buildmany <routing_table.csv> <r1,r2,..> <out_prefix> [granularity=16]
    partition.py check     <routing_table.csv> <out.part>

    routing_table.csv   The routing parameter table. Only the parents column is read.
    n_ranks             How many MPI ranks to divide the network across.
    out.part            Where to write the partition. A .txt summary is written beside it.
    out_prefix          As out.part, but buildmany appends _r<N>.part for each rank count.
    granularity         Pieces per rank. At granularity 1 there is one piece per rank.
                        Additional pieces allow more balanced packing at the cost of more cut edges.

`buildmany` reads the table once and emits <out_prefix>_r<N>.part for each rank count,
and is more efficient than calling `build` repeatedly.

`check` verifies a partition file already on disk: every link is assigned,
every rank is in range, and the rank graph has no backward edges. Exits non-zero on failure.


OUTPUT FORMAT (.part)
    Little-endian binary, designed to be read in one fread from C++:
        magic   char[8]  "HLMPART1"
        n_links uint64
        n_ranks uint64
        rank    int32[n_links]   rank owning each link, indexed by link index
    A sidecar <out.part>.txt carries the human-readable summary.
"""

import sys
import struct
from array import array

MAGIC = b"HLMPART1"
NO_CHILD = -1


def read_children(csv_path):
    """
    Return (n_links, child) where child is a flat array of length n_links,
    and child[i] is the single downstream link of i, or NO_CHILD.

    Only column 4 (parents) is parsed. The index column is the row number, and params are
    irrelevant here. Each parent p of row i contributes the edge p -> i, so one pass over
    the file builds the whole child array without ever holding a parent list.
    """
    child = array("q")
    n = 0
    with open(csv_path, "r") as fh:
        fh.readline()  # header
        for line in fh:
            # index,node,level,parents,params -- parents is the 4th field
            f = line.split(",", 4)
            if len(f) < 4:
                continue
            child.append(NO_CHILD)
            n += 1
    # Second pass fills the edges. Two passes avoid growing and indexing at once, and the
    # file is read from page cache the second time.
    with open(csv_path, "r") as fh:
        fh.readline()
        i = 0
        for line in fh:
            f = line.split(",", 4)
            if len(f) < 4:
                continue
            par = f[3]
            if par:
                for p in par.split(";"):
                    if p:
                        pi = int(p)
                        if not (0 <= pi < n):
                            sys.exit(f"ERROR: link {i} lists parent {pi}, outside 0..{n-1}")
                        if child[pi] != NO_CHILD:
                            sys.exit(f"ERROR: link {pi} drains to both {child[pi]} and {i}. "
                                     "The traversal requires one downstream link per link.")
                        child[pi] = i
            i += 1
    return n, child


def decompose(n, child, target):
    """
    Cut the forest into connected pieces of at least `target` links where possible.

    One forward sweep. acc[i] is the number of links in the piece being built at i,
    counting i and every upstream link not already cut away. Because child index > parent
    index always, a link's accumulated size is final by the time the sweep reaches it.

    Returns piece[i], the piece id owning each link, and a list of each piece's outlet.
    """
    acc = array("q", [1]) * n
    cut = bytearray(n)

    for i in range(n):
        c = child[i]
        if acc[i] >= target or c == NO_CHILD:
            # This link ends its piece: either it is big enough, or it is a true outlet
            # with nothing downstream to merge into.
            cut[i] = 1
        elif c != NO_CHILD:
            acc[c] += acc[i]

    # Assign piece ids by sweeping backwards: a link belongs to its child's piece unless
    # it was cut. Reverse file order visits every child before its parents.
    piece = array("q", [-1]) * n
    outlets = []
    for i in range(n - 1, -1, -1):
        if cut[i]:
            piece[i] = len(outlets)
            outlets.append(i)
        else:
            piece[i] = piece[child[i]]
    return piece, outlets


def assign_ranks(n, child, piece, outlets, n_ranks):
    """
    Pack pieces into n_ranks, keeping the rank graph acyclic.

    Pieces are ordered by outlet index, which is a topological order,
    then split into contiguous runs of roughly equal link count. Contiguity means every
    cross-rank edge points from a lower rank to a higher one, so no cycle can exist.
    """
    n_pieces = len(outlets)
    size = [0] * n_pieces
    for i in range(n):
        size[piece[i]] += 1

    order = sorted(range(n_pieces), key=lambda p: outlets[p])
    sizes = [size[p] for p in order]

    # Split the ordered list into at most n_ranks contiguous runs, minimising the largest
    # run. Binary search the load limit: the smallest limit that still fits is optimal, and
    # a greedy left-to-right fill is the feasibility test for a given limit.
    #
    # A previous version advanced ranks by comparing a running total against rank
    # boundaries, which skipped ranks entirely whenever one piece exceeded a rank's share
    # and left them with no links at all (Raritan R=8 produced an empty rank). Searching
    # for the limit has no such failure mode: every rank is filled before the next opens.
    def runs_needed(limit):
        used, load = 1, 0
        for s in sizes:
            if load + s <= limit:
                load += s
            else:
                used += 1
                load = s
        return used

    lo, hi = max(sizes), n          # a run can never be smaller than the biggest piece
    while lo < hi:
        mid = (lo + hi) // 2
        if runs_needed(mid) <= n_ranks:
            hi = mid
        else:
            lo = mid + 1

    # Fill left to right under that limit. The second clause is what stops ranks at the
    # end being left empty: greedy under an optimal limit can finish in fewer than n_ranks
    # runs, so once only as many pieces remain as ranks, every one of them opens a new
    # rank. Without it, Raritan R=16 handed the last ranks nothing at all.
    piece_rank = [0] * n_pieces
    r, load = 0, 0
    for idx, p in enumerate(order):
        pieces_left = n_pieces - idx
        ranks_left = n_ranks - r
        if r < n_ranks - 1 and load > 0 and (load + sizes[idx] > lo or pieces_left <= ranks_left):
            r += 1
            load = 0
        piece_rank[p] = r
        load += sizes[idx]

    rank = array("i", [0]) * n
    for i in range(n):
        rank[i] = piece_rank[piece[i]]
    return rank, size, piece_rank


def summarize(n, child, rank, n_ranks, size, piece_rank):
    """Cut edges, per-rank balance, and a check that the rank graph really is acyclic."""
    counts = [0] * n_ranks
    for i in range(n):
        counts[rank[i]] += 1

    cut_edges = 0
    backward = 0
    edges = set()
    for i in range(n):
        c = child[i]
        if c != NO_CHILD and rank[c] != rank[i]:
            cut_edges += 1
            edges.add((rank[i], rank[c]))
            if rank[c] < rank[i]:
                backward += 1

    lines = []
    lines.append(f"links          {n}")
    lines.append(f"ranks          {n_ranks}")
    lines.append(f"pieces         {len(size)}")
    lines.append(f"cut edges      {cut_edges}   (links whose child is on another rank)")
    lines.append(f"rank pairs     {len(edges)}  (distinct sender->receiver pairs)")
    lines.append(f"backward edges {backward}   (MUST be 0: nonzero means a possible cycle)")
    lo, hi = min(counts), max(counts)
    lines.append(f"balance        min {lo}  max {hi}  imbalance {hi/float(max(lo,1)):.3f}x")
    lines.append("")
    lines.append("rank    links   share")
    for r in range(n_ranks):
        lines.append(f"{r:4d} {counts[r]:10d}  {100.0*counts[r]/n:5.2f}%")
    return "\n".join(lines), backward


def build_one(n, child, n_ranks, out_path, granularity):
    """Partition an already-parsed network. Returns 1 if the result looks unsafe."""
    # Pieces are sized well below one rank's share. Piece size is what limits how evenly
    # the packer can fill ranks, and a piece cannot be split once formed, so sizing them at
    # N/n_ranks leaves barely more pieces than ranks and balance collapses (Raritan R=4 hit
    # 6.4x). Each extra piece costs at most one cut edge, and a cut edge carries ~5.7 KB
    # per 1-day chunk, so finer pieces are close to free.
    target = max(1, n // (n_ranks * granularity))
    piece, outlets = decompose(n, child, target)
    print(f"  decomposed into {len(outlets)} connected pieces (target {target} links)")
    rank, size, piece_rank = assign_ranks(n, child, piece, outlets, n_ranks)
    text, backward = summarize(n, child, rank, n_ranks, size, piece_rank)
    print(text)

    with open(out_path, "wb") as fh:
        fh.write(MAGIC)
        fh.write(struct.pack("<QQ", n, n_ranks))
        rank.tofile(fh)
    with open(out_path + ".txt", "w") as fh:
        fh.write(text + "\n")
    print(f"\n  wrote {out_path} and {out_path}.txt")
    return 1 if backward else 0


def cmd_build(csv_path, n_ranks, out_path, granularity=16):
    n, child = read_children(csv_path)
    print(f"  read {n} links")
    return build_one(n, child, n_ranks, out_path, granularity)


def cmd_buildmany(csv_path, rank_list, out_prefix, granularity=16):
    n, child = read_children(csv_path)
    print(f"  read {n} links (once, for {len(rank_list)} rank counts)")
    status = 0
    for r in rank_list:
        print(f"\n----- {r} ranks -----")
        status |= build_one(n, child, r, f"{out_prefix}_r{r}.part", granularity)
    return status


def cmd_check(csv_path, part_path):
    """Re-derive everything from the table and confirm the stored partition is sound."""
    with open(part_path, "rb") as fh:
        if fh.read(8) != MAGIC:
            sys.exit(f"ERROR: {part_path} is not an HLMPART1 file")
        n_stored, n_ranks = struct.unpack("<QQ", fh.read(16))
        rank = array("i")
        rank.fromfile(fh, n_stored)

    n, child = read_children(csv_path)
    if n != n_stored:
        sys.exit(f"ERROR: table has {n} links, partition file has {n_stored}")

    bad = sum(1 for i in range(n) if not (0 <= rank[i] < n_ranks))
    if bad:
        sys.exit(f"ERROR: {bad} links carry a rank outside 0..{n_ranks-1}")

    text, backward = summarize(n, child, rank, n_ranks, [0], None)
    print(text)
    if backward:
        print("\nFAIL - backward edges present, the rank graph may contain a cycle")
        return 1
    print("\nPASS - every link assigned, rank graph is acyclic")
    return 0


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    if argv[1] == "build":
        if len(argv) not in (5, 6):
            sys.exit("usage: partition.py build <table.csv> <n_ranks> <out.part> [granularity]")
        return cmd_build(argv[2], int(argv[3]), argv[4],
                         int(argv[5]) if len(argv) > 5 else 16)
    if argv[1] == "buildmany":
        if len(argv) not in (5, 6):
            sys.exit("usage: partition.py buildmany <table.csv> <r1,r2,..> <prefix> [gran]")
        ranks = [int(x) for x in argv[3].split(",") if x]
        return cmd_buildmany(argv[2], ranks, argv[4],
                             int(argv[5]) if len(argv) > 5 else 16)
    if argv[1] == "check":
        if len(argv) != 4:
            sys.exit("usage: partition.py check <table.csv> <out.part>")
        return cmd_check(argv[2], argv[3])
    sys.exit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
