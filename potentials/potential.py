"""Weighted matching polynomials and local weighted matching potentials (Section 2 of the paper).

For a graph U on V (a set of frozenset edges) and weights a : V -> ring,
    mu_a(U[S]) = sum over matchings M of U[S] of (-1)^{|M|} prod_{x in S not covered by M} a_x,
computed by the vertex recurrence  mu_a(H) = a_v mu_a(H - v) - sum_{w ~ v} mu_a(H - v - w).
A potential is a pair (U, a) with M(W) != 0 and M(W - v) = R(W, v) M(W) at every reachable position,
where M(S) = mu_a(U[S]).
"""
from .dvg import bits
from .qz import RatFunc

__all__ = ["MatchingPoly", "check_potential", "directed_collapse_potential", "symmetric_part"]


class MatchingPoly:
    def __init__(self, n, U, a, one):
        self.adj = [0] * n
        for e in U:
            i, j = tuple(e)
            self.adj[i] |= 1 << j
            self.adj[j] |= 1 << i
        self.a, self.one = a, one
        self.memo = {0: one}

    def __call__(self, S):
        r = self.memo.get(S)
        if r is None:
            v = (S & -S).bit_length() - 1
            rest = S & ~(1 << v)
            r = self.a[v] * self(rest)
            for w in bits(self.adj[v] & rest):
                r = r - self(rest & ~(1 << w))
            self.memo[S] = r
        return r


def check_potential(D, U, a):
    """Exact check of the potential identities in Q(z) at every reachable position.
    Returns (ok, first failing position or None)."""
    M = MatchingPoly(D.n, U, a, RatFunc.const(1))
    for (W, v) in sorted(D.reachable_positions()):
        MW = M(W)
        if MW.is_zero():
            return False, (W, v)
        if M(W & ~(1 << v)) != D.R(W, v) * MW:
            return False, (W, v)
    return True, None


def symmetric_part(D):
    return {frozenset((a, b)) for a, b in D.arcs if (D.out[b] >> a) & 1}


def directed_collapse_potential(D):
    """The potential of Theorem 2.9 on an SCC-symmetric digraph: U = symmetric part,
    a_x = z - eta_x, eta_x = sum of R(V, y) over out-neighbours y of x outside the strong component of x."""
    assert D.is_scc_symmetric()
    z = RatFunc.var()
    a = []
    for x in range(D.n):
        comp = D.strong_component(x)
        eta = RatFunc.const(0)
        for y in bits(D.out[x] & ~comp):
            eta = eta + D.R(D.full, y)
        a.append(z - eta)
    return symmetric_part(D), a
