"""Clearing positions and normality propagation (Section 3 of the paper).

Definition 3.4.  For N a set of vertices, an N-clearing set for x is a path set S containing x,
such that some directed path with vertex set S ends at x, and every vertex w outside S with
x -/-> w satisfies
   (a) S + w is a path set, or
   (b) V - (S + w) is contained in N, and some path set X has V - X contained in N.
N_0 = {}, N_{k+1} = N_k + {x : x has an N_k-clearing set}, N_inf = union.
D is fully clearable when N_inf = V.

Lemma 3.5 shows a_x = z + O(1/z) for x in N_inf in every potential over the algebraic closure
of Q(z).  This module computes N_inf with a certificate and re-checks certificates path by path.
"""
from .dvg import bits

__all__ = ["normal_closure", "check_certificate", "fully_clearable"]


def _has_cover(ends, full, N):
    """Some path set X with V - X contained in N (the global clause of (b))."""
    for X in range(1, full + 1):
        if ends[X] and ((full & ~X) & ~N) == 0:
            return X
    return None


def normal_closure(D):
    """Return (N_inf, certificate).

    The certificate is a list of rounds; round k lists triples (x, S, reasons) with S an N_k-clearing
    set for x and reasons[w] in {'a', 'b'} for each w outside S with x -/-> w.  Each round also records
    the path set X used by clause (b), if any.
    """
    ends, full = D.ends(), D.full
    N, rounds = 0, []
    while True:
        X = _has_cover(ends, full, N)
        found = []
        for x in range(D.n):
            if (N >> x) & 1:
                continue
            for S in range(1, full + 1):
                if not (ends[S] >> x) & 1:
                    continue
                rest = full & ~S
                reasons, ok = {}, True
                for w in bits(rest & ~D.out[x]):
                    if ends[S | (1 << w)]:
                        reasons[w] = "a"
                    elif X is not None and ((rest & ~(1 << w)) & ~N) == 0:
                        reasons[w] = "b"
                    else:
                        ok = False
                        break
                if ok:
                    found.append((x, S, reasons))
                    break
        if not found:
            return N, rounds
        rounds.append({"cover": X, "steps": found})
        for x, _, _ in found:
            N |= 1 << x


def fully_clearable(D):
    return normal_closure(D)[0] == D.full


def check_certificate(D, rounds):
    """Re-check a certificate with explicit paths.  Returns the certified set."""
    full, N = D.full, 0
    for rnd in rounds:
        X = rnd["cover"]
        if X is not None:
            P = D.path_with(X)
            assert P is not None and D.is_directed_path(P) and sum(1 << v for v in P) == X
            assert ((full & ~X) & ~N) == 0
        for x, S, reasons in rnd["steps"]:
            P = D.path_with(S, end=x)
            assert P is not None and D.is_directed_path(P) and P[-1] == x and sum(1 << v for v in P) == S
            rest = full & ~S
            for w in bits(rest & ~D.out[x]):
                r = reasons.get(w)
                if r == "a":
                    Q = D.path_with(S | (1 << w))
                    assert Q is not None and D.is_directed_path(Q) and sum(1 << v for v in Q) == S | (1 << w)
                elif r == "b":
                    assert X is not None and ((rest & ~(1 << w)) & ~N) == 0
                else:
                    raise AssertionError("vertex %d is not justified for %d" % (w, x))
        for x, _, _ in rnd["steps"]:
            N |= 1 << x
    return N
