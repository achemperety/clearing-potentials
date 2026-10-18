"""Constructors for the digraph families used in the paper."""
import itertools
import random

from .dvg import Digraph

__all__ = ["directed_cycle", "blowup", "cycle_with_pendants", "cycle_with_legs",
           "book_with_pendants", "blowup_with_pendants", "chain_into_cycle",
           "random_tournament", "transitive_tournament", "random_digraph", "with_universal_vertex",
           "random_scc_symmetric"]


def directed_cycle(m):
    return Digraph(m, [(i, (i + 1) % m) for i in range(m)])


def blowup(sizes):
    """C_m(s_0,...,s_{m-1}): class S_i of s_i twins, every arc from S_i to S_{i+1}.  Returns (D, classes)."""
    classes, nxt = [], 0
    for s in sizes:
        classes.append(list(range(nxt, nxt + s)))
        nxt += s
    m = len(sizes)
    arcs = [(a, b) for i in range(m) for a in classes[i] for b in classes[(i + 1) % m]]
    return Digraph(nxt, arcs), classes


def _add_pendants(n, arcs, at, k):
    arcs = list(arcs)
    for _ in range(k):
        arcs += [(at, n), (n, at)]
        n += 1
    return n, arcs


def _add_leg(n, arcs, at, length):
    arcs, prev = list(arcs), at
    for _ in range(length):
        arcs += [(prev, n), (n, prev)]
        prev, n = n, n + 1
    return n, arcs


def cycle_with_pendants(m, counts):
    """Directed m-cycle c_0..c_{m-1} (vertices 0..m-1) with counts[i] pendant vertices at c_i."""
    n, arcs = m, [(i, (i + 1) % m) for i in range(m)]
    for i, k in enumerate(counts):
        n, arcs = _add_pendants(n, arcs, i, k)
    return Digraph(n, arcs)


def cycle_with_legs(m, lengths):
    """Directed m-cycle with a symmetric path ("leg") of lengths[i] vertices hanging at c_i."""
    n, arcs = m, [(i, (i + 1) % m) for i in range(m)]
    for i, L in enumerate(lengths):
        n, arcs = _add_leg(n, arcs, i, L)
    return Digraph(n, arcs)


def book_with_pendants(r, k):
    """Book B_r = C_3(r,1,1): spine 0 -> 1, pages 2..r+1 with 1 -> p -> 0; k pendant vertices at 0."""
    arcs = [(0, 1)] + [(1, p) for p in range(2, r + 2)] + [(p, 0) for p in range(2, r + 2)]
    n, arcs = _add_pendants(r + 2, arcs, 0, k)
    return Digraph(n, arcs)


def blowup_with_pendants(sizes, pendants):
    """Blow-up C_m(sizes) with pendants[i] pendant vertices at the first vertex of class i."""
    D0, classes = blowup(sizes)
    n, arcs = D0.n, D0.arcs
    for i, k in enumerate(pendants):
        n, arcs = _add_pendants(n, arcs, classes[i][0], k)
    return Digraph(n, arcs)


def chain_into_cycle(m, k):
    """Directed m-cycle with a directed path f_1 -> ... -> f_k -> c_0 attached (f_1 a source)."""
    arcs = [(i, (i + 1) % m) for i in range(m)]
    for j in range(k):
        arcs.append((m + j, m + j + 1) if j < k - 1 else (m + j, 0))
    return Digraph(m + k, arcs)


def transitive_tournament(n):
    return Digraph(n, [(i, j) for i in range(n) for j in range(i + 1, n)])


def random_tournament(n, rng):
    return Digraph(n, [(i, j) if rng.random() < 0.5 else (j, i) for i, j in itertools.combinations(range(n), 2)])


def random_digraph(n, p, rng):
    return Digraph(n, [(i, j) for i in range(n) for j in range(n) if i != j and rng.random() < p])


def with_universal_vertex(n, p, rng):
    """A random digraph on n vertices in which vertex 0 is adjacent (in some direction) to all others."""
    arcs = [(i, j) for i in range(1, n) for j in range(1, n) if i != j and rng.random() < p]
    for v in range(1, n):
        r = rng.random()
        if r < 1 / 3:
            arcs.append((0, v))
        elif r < 2 / 3:
            arcs.append((v, 0))
        else:
            arcs += [(0, v), (v, 0)]
    return Digraph(n, arcs)


def random_scc_symmetric(n, rng, blocks=3, p=0.5, q=0.4):
    """Random SCC-symmetric digraph: symmetric random graphs on consecutive blocks, joined by
    one-way arcs from earlier to later blocks."""
    cuts = sorted(rng.sample(range(1, n), min(blocks - 1, n - 1)))
    bounds = [0] + cuts + [n]
    block = [None] * n
    for b in range(len(bounds) - 1):
        for v in range(bounds[b], bounds[b + 1]):
            block[v] = b
    arcs = []
    for i, j in itertools.combinations(range(n), 2):
        if block[i] == block[j]:
            if rng.random() < p:
                arcs += [(i, j), (j, i)]
        elif rng.random() < q:
            arcs.append((i, j))
    return Digraph(n, arcs)
