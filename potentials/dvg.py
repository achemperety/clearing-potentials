"""Directed vertex geography (DVG) on small digraphs, with exact resolvents.

Vertices are 0..n-1 and vertex sets are bitmasks.  A position is a pair (W, v) with v in W,
W the set of unvisited vertices including the token; its options are (W - v, w) for the
out-neighbours w of v in W - v.  The root resolvent satisfies

        1 / R(W, v) = z - sum_{w} R(W - v, w)            (Section 2.1 of the paper).

A position is reachable when some directed path ends at v and has vertex set (V - W) + v.
"""
from fractions import Fraction
import itertools

from .qz import RatFunc, Z

__all__ = ["Digraph", "bits", "popcount"]


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def popcount(mask):
    return bin(mask).count("1")


class Digraph:
    def __init__(self, n, arcs):
        self.n = n
        self.arcs = sorted({(int(a), int(b)) for a, b in arcs})
        for a, b in self.arcs:
            if a == b or not (0 <= a < n and 0 <= b < n):
                raise ValueError("bad arc %r" % ((a, b),))
        self.full = (1 << n) - 1
        self.out = [0] * n
        self.inn = [0] * n
        for a, b in self.arcs:
            self.out[a] |= 1 << b
            self.inn[b] |= 1 << a
        self._ends = None
        self._R = {}
        self._Rq = {}

    # ------------------------------------------------------------ structure
    def reach(self, v, adj=None):
        adj = self.out if adj is None else adj
        seen = todo = 1 << v
        while todo:
            x = (todo & -todo).bit_length() - 1
            todo &= todo - 1
            new = adj[x] & ~seen
            seen |= new
            todo |= new
        return seen

    def is_strongly_connected(self):
        return self.n == 0 or (self.reach(0) == self.full and self.reach(0, self.inn) == self.full)

    def is_symmetric(self):
        return self.out == self.inn

    def is_unilateral(self):
        r = [self.reach(v) for v in range(self.n)]
        return all((r[u] >> v) & 1 or (r[v] >> u) & 1 for u in range(self.n) for v in range(self.n))

    def is_scc_symmetric(self):
        """Every arc that lies on a directed cycle has its reverse."""
        r = [self.reach(v) for v in range(self.n)]
        return all(not ((r[b] >> a) & 1) or ((self.out[b] >> a) & 1) for a, b in self.arcs)

    def comparable(self, u, v):
        return bool((self.reach(u) >> v) & 1 or (self.reach(v) >> u) & 1)

    def strong_component(self, v):
        return self.reach(v) & self.reach(v, self.inn)

    # ------------------------------------------------------------ paths
    def ends(self):
        """ends[S] = bitmask of the vertices at which some directed path with vertex set S ends."""
        if self._ends is None:
            e = [0] * (1 << self.n)
            for v in range(self.n):
                e[1 << v] = 1 << v
            for S in range(1, 1 << self.n):
                for v in bits(e[S]):
                    for w in bits(self.out[v] & ~S):
                        e[S | (1 << w)] |= 1 << w
            self._ends = e
        return self._ends

    def is_path_set(self, S):
        return S != 0 and self.ends()[S] != 0

    def path_with(self, S, end=None):
        """An explicit directed path with vertex set S (ending at `end` if given), or None."""
        e = self.ends()
        if not S or not e[S] or (end is not None and not (e[S] >> end) & 1):
            return None
        cur = end if end is not None else (e[S] & -e[S]).bit_length() - 1
        rev, T = [cur], S
        while T != (1 << cur):
            T2 = T & ~(1 << cur)
            prev = next(p for p in bits(e[T2] & self.inn[cur]))
            rev.append(prev)
            cur, T = prev, T2
        return rev[::-1]

    def is_directed_path(self, P):
        return len(set(P)) == len(P) and all((self.out[P[i]] >> P[i + 1]) & 1 for i in range(len(P) - 1))

    def directed_paths(self):
        """All directed paths (as lists), each listed once."""
        stack = [[v] for v in range(self.n)]
        while stack:
            P = stack.pop()
            yield P
            used = sum(1 << v for v in P)
            for w in bits(self.out[P[-1]] & ~used):
                stack.append(P + [w])

    def back(self, P):
        """Number of arcs t_j -> t_i with i < j along the path P."""
        return sum(1 for j in range(len(P)) for i in range(j) if (self.out[P[j]] >> P[i]) & 1)

    def reachable_positions(self):
        seen, stack = set(), [(self.full, s) for s in range(self.n)]
        while stack:
            W, v = stack.pop()
            if (W, v) in seen:
                continue
            seen.add((W, v))
            W2 = W & ~(1 << v)
            for w in bits(self.out[v] & W2):
                stack.append((W2, w))
        return seen

    # ------------------------------------------------------------ resolvents
    def R(self, W, v):
        """Root resolvent of the position (W, v), exactly in Q(z)."""
        key = (W, v)
        r = self._R.get(key)
        if r is None:
            W2 = W & ~(1 << v)
            s = Z
            for w in bits(self.out[v] & W2):
                s = s - self.R(W2, w)
            r = s.inv()
            self._R[key] = r
        return r

    def Rq(self, W, v, z0):
        """The same resolvent evaluated exactly at a rational point z0 (see Lemma 2.2)."""
        key = (W, v, z0)
        r = self._Rq.get(key)
        if r is None:
            W2 = W & ~(1 << v)
            s = Fraction(z0)
            for w in bits(self.out[v] & W2):
                s -= self.Rq(W2, w, z0)
            if s == 0:
                raise ZeroDivisionError("pole at %s" % z0)
            r = 1 / s
            self._Rq[key] = r
        return r

    def gamma(self, P, z0=None):
        """Product of the resolvents along the play P started at (V, P[0])."""
        if not self.is_directed_path(P):
            raise ValueError("not a directed path: %r" % (P,))
        W = self.full
        g = RatFunc.const(1) if z0 is None else Fraction(1)
        for v in P:
            g = g * (self.R(W, v) if z0 is None else self.Rq(W, v, z0))
            W &= ~(1 << v)
        return g

    # ------------------------------------------------------------ obstructions
    def back_arc_conflict(self):
        """Two directed paths with equal vertex sets and different back counts, or None."""
        seen = {}
        for P in self.directed_paths():
            S = sum(1 << v for v in P)
            b = self.back(P)
            if S in seen and seen[S][0] != b:
                return seen[S][1], P
            seen.setdefault(S, (b, P))
        return None

    def gamma_conflict(self, z0=Fraction(7, 3)):
        """Two directed paths with equal vertex sets whose products Gamma differ at z0, or None."""
        seen = {}
        for P in self.directed_paths():
            S = sum(1 << v for v in P)
            g = self.gamma(P, z0)
            if S in seen and seen[S][0] != g:
                return seen[S][1], P, seen[S][0], g
            seen.setdefault(S, (g, P))
        return None

    def set_potential(self):
        """M(W)/M(V) on all unvisited sets of reachable positions if Gamma depends only on vertex
        sets (exactly, in Q(z)); otherwise None."""
        order = sorted(self.reachable_positions(), key=lambda p: -popcount(p[0]))
        M = {self.full: RatFunc.const(1)}
        for (W, v) in order:
            val = self.R(W, v) * M[W]
            W2 = W & ~(1 << v)
            if W2 in M:
                if M[W2] != val:
                    return None
            else:
                M[W2] = val
        return M

    # ------------------------------------------------------------ misc
    def relabel(self, perm):
        return Digraph(self.n, [(perm[a], perm[b]) for a, b in self.arcs])

    def is_isomorphic(self, other):
        if self.n != other.n or len(self.arcs) != len(other.arcs):
            return False
        target = set(other.arcs)
        outdeg = sorted(popcount(x) for x in self.out)
        if outdeg != sorted(popcount(x) for x in other.out):
            return False
        for p in itertools.permutations(range(self.n)):
            if all((p[a], p[b]) in target for a, b in self.arcs):
                return True
        return False

    @staticmethod
    def from_directg(line):
        t = list(map(int, line.split()))
        n, m = t[0], t[1]
        return Digraph(n, [(t[2 + 2 * i], t[3 + 2 * i]) for i in range(m)])

    def to_directg(self):
        return " ".join(map(str, [self.n, len(self.arcs)] + [x for a in self.arcs for x in a]))
