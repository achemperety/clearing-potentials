#!/usr/bin/env python3
"""Report for "Clearing positions and the rigidity of complex matching potentials on digraphs".

Every check is exact (Fractions, Q[z], Q(z)).  Sections follow the paper.

    python3 scripts/verify_all.py            # about one minute
    python3 scripts/verify_all.py --quick    # a few seconds

The exhaustive census (Theorem 5.1) is run by csrc/census.c (see README); this script re-reads its outputs
in data/ and logs/ and certifies the survivors.  The Groebner computations of Section 6 are in groebner/.
"""
import argparse
import itertools
import json
import os
import random
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, ROOT)

from potentials import (Digraph, RatFunc, Z, Poly, bits, popcount, normal_closure, check_certificate,  # noqa: E402
                        check_potential, directed_collapse_potential, MatchingPoly)
from potentials.families import (directed_cycle, blowup, cycle_with_pendants, cycle_with_legs,  # noqa: E402
                                 book_with_pendants, blowup_with_pendants, chain_into_cycle,
                                 random_tournament, transitive_tournament, random_digraph,
                                 with_universal_vertex, random_scc_symmetric)

OK = []


def check(name, cond, detail=""):
    OK.append(bool(cond))
    print(("  [ok]   " if cond else "  [FAIL] ") + name + (("  " + detail) if detail else ""))


def hdr(s):
    print("\n" + s + "\n" + "-" * len(s))


# ----------------------------------------------------------------------------------------------- Section 2
def section2(rng, trials):
    hdr("Section 2: resolvents and potentials")
    ok_asym = ok_int = True
    for _ in range(trials):
        D = random_digraph(rng.randint(1, 6), rng.random(), rng)
        for (W, v) in D.reachable_positions():
            R = D.R(W, v)
            deg, co = R.laurent_at_infinity(3)
            d = popcount(D.out[v] & W & ~(1 << v))
            ok_asym &= (deg == -1 and co[0] == 1 and co[1] == 0 and co[2] == d)
            ok_int &= all(c.denominator == 1 for c in R.num.c + R.den.c) and R.num.lead() == 1 and R.den.lead() == 1
            ok_int &= R.den.deg() == R.num.deg() + 1
    check("Lemma 2.1: R(W,v) = 1/z + 0/z^2 + d^+_W(v)/z^3 + ... at every reachable position", ok_asym)
    check("Lemma 2.2: numerator and denominator monic in Z[z], deg den = deg num + 1", ok_int)

    ok_coll = ok_cav = ok_inv = True
    for _ in range(trials):
        D = random_scc_symmetric(rng.randint(2, 7), rng)
        U, a = directed_collapse_potential(D)
        good, _ = check_potential(D, U, a)
        ok_coll &= good
        M = MatchingPoly(D.n, U, a, RatFunc.const(1))
        NU = [0] * D.n
        for e in U:
            i, j = tuple(e)
            NU[i] |= 1 << j
            NU[j] |= 1 << i
        for (W, v) in D.reachable_positions():
            W2 = W & ~(1 << v)
            rhs = RatFunc.const(0)
            for w in bits(NU[v] & ~D.out[v] & W2):
                rhs = rhs + M(W2 & ~(1 << w)) / M(W2)
            for t in bits(D.out[v] & ~NU[v] & W2):
                rhs = rhs - D.R(W2, t)
            ok_cav &= (a[v] - Z == rhs)
        for P in itertools.islice(D.directed_paths(), 200):
            X = sum(1 << x for x in P)
            ok_inv &= (D.gamma(P) == M(D.full & ~X) / M(D.full))
    check("Theorem 2.7: the directed-collapse potential satisfies every identity (random SCC-symmetric digraphs)", ok_coll)
    check("Lemma 2.5: cavity identity at every reachable position of these potentials", ok_cav)
    check("Lemma 2.3: Gamma(P) = M(V - V(P)) / M(V) along plays", ok_inv)


# ----------------------------------------------------------------------------------------------- Section 3
def section3(rng, trials):
    hdr("Section 3: clearing positions")
    # Lemma 3.1 (terminal rigidity): a vertex that ends a directed Hamiltonian path has weight z in every potential
    ok_term, nterm = True, 0
    for _ in range(trials):
        D = random_scc_symmetric(rng.randint(2, 7), rng)
        U, a = directed_collapse_potential(D)
        reach = set(D.reachable_positions())
        for x in bits(D.ends()[D.full]):
            nterm += 1
            ok_term &= a[x] == Z and ((1 << x), x) in reach and D.R(1 << x, x) == RatFunc.const(1) / Z
    check("Lemma 3.1: every end x of a directed Hamiltonian path has a_x = z (directed-collapse potentials); "
          "({x},x) is reachable with R = 1/z", ok_term and nterm > 0, "%d vertices" % nterm)
    ok_inc = ok_cover = ok_cert = True
    for _ in range(trials):
        D = random_digraph(rng.randint(2, 7), rng.random(), rng)
        N, cert = normal_closure(D)
        ok_cert &= (check_certificate(D, cert) == N)
        for u in range(D.n):
            for x in range(u + 1, D.n):
                if not D.comparable(u, x):
                    ok_inc &= not ((N >> u) & 1) and not ((N >> x) & 1)
        if N != D.full:
            ends = D.ends()
            ok_cover &= not any(ends[X] and ((D.full & ~X) & ~N) == 0 for X in range(1, D.full + 1))
    check("certificates of N_inf re-checked path by path (random digraphs)", ok_cert)
    check("Lemma 3.5: incomparable vertices are never cleared", ok_inc)
    check("Proposition 3.8: if N_inf != V, no path set contains V - N_inf", ok_cover)
    # Example 3.7: two isolated vertices, U = {xy}, a = alpha with alpha^2 = z alpha + 1
    # identities: a_y = M(V)/z and a_x = M(V)/z with M(V) = alpha^2 - 1;  alpha^2 - 1 = z alpha.
    # Work in Q(z)[alpha]/(alpha^2 - z alpha - 1): elements c0 + c1 alpha with c0, c1 in Q(z).
    def mul(p, q):
        c0 = p[0] * q[0] + p[1] * q[1]               # alpha^2 = z alpha + 1
        c1 = p[0] * q[1] + p[1] * q[0] + p[1] * q[1] * Z
        return (c0, c1)
    one, zero = RatFunc.const(1), RatFunc.const(0)
    alpha = (zero, one)
    MV = mul(alpha, alpha)
    MV = (MV[0] - one, MV[1])
    lhs = (MV[0] / Z, MV[1] / Z)
    check("Example 3.7: alpha = (alpha^2 - 1)/z holds in Q(z)[alpha]/(alpha^2 - z alpha - 1)", lhs == alpha)


# ----------------------------------------------------------------------------------------------- Section 4
def section4(rng, trials, quick):
    hdr("Section 4: families")
    cnt = 0
    ok_bl = True
    mmax = 6 if quick else 8
    for m in range(3, mmax + 1):
        for sizes in itertools.product(range(1, 4), repeat=m):
            if sum(sizes) > 12 or tuple(sizes) != min(tuple(sizes[i:] + sizes[:i]) for i in range(m)):
                continue
            D, _ = blowup(list(sizes))
            N, cert = normal_closure(D)
            ok_bl &= (N == D.full and check_certificate(D, cert) == N and D.is_strongly_connected()
                      and not D.is_scc_symmetric())
            cnt += 1
    check("Theorem 4.3: blow-ups C_m(s), m <= %d, s_i <= 3, n <= 12, up to rotation, are fully clearable" % mmax,
          ok_bl, "(%d blow-ups)" % cnt)
    ok_nh = True
    for m in (4, 6):
        for s0, s1 in [(2, 1), (1, 2), (3, 1), (3, 2)]:
            D, _ = blowup([s0, s1] * (m // 2))
            if D.n <= 12:
                ok_nh &= (D.ends()[D.full] == 0 and D.set_potential() is not None)
    check("Remark 4.4: two-periodic blow-ups with s0 != s1 have a set potential and no Hamiltonian path", ok_nh)
    ok_t = True
    for _ in range(trials):
        n = rng.randint(2, 8)
        T = random_tournament(n, rng)
        ok_t &= normal_closure(T)[0] == T.full
    for n in range(2, 8):
        T = transitive_tournament(n)
        U, a = directed_collapse_potential(T)
        ok_t &= check_potential(T, U, a)[0] and len(U) == 0
    check("Theorem 4.5 / Corollary 1.4(a): random tournaments are fully clearable; transitive ones carry the "
          "collapse potential with U empty", ok_t)
    ok_h = True
    for _ in range(trials):
        n = rng.randint(2, 7)
        perm = list(range(n))
        rng.shuffle(perm)
        arcs = [(perm[i], perm[i + 1]) for i in range(n - 1)]
        arcs += [(i, j) for i in range(n) for j in range(n) if i != j and rng.random() < 0.3]
        D = Digraph(n, arcs)
        ok_h &= normal_closure(D)[0] == D.full
    check("Theorem 4.5: digraphs with a planted directed Hamiltonian path are fully clearable", ok_h)
    ok_u = True
    for _ in range(trials):
        D = with_universal_vertex(rng.randint(3, 7), rng.random() * 0.5, rng)
        if D.is_strongly_connected():
            ok_u &= normal_closure(D)[0] == D.full
    check("Corollary 4.2(b): strongly connected digraphs with a universal vertex are fully clearable", ok_u)
    ok_cls = True
    for r in range(2, 6):
        D = book_with_pendants(r, r - 1)
        ok_cls &= normal_closure(D)[0] == D.full and D.set_potential() is not None
    D = blowup_with_pendants([2, 1, 1, 1], [0, 1, 1, 0])
    ok_cls &= normal_closure(D)[0] == D.full and D.set_potential() is not None and D.ends()[D.full] == 0
    for m in (3, 4, 5):
        for k in (1, 2, 3):
            E = chain_into_cycle(m, k)
            ok_cls &= normal_closure(E)[0] == E.full and E.set_potential() is not None
    check("Example 4.6: books with r-1 pendants (r <= 5), C4(2,1,1,1) with pendants, chains into cycles: "
          "fully clearable, with set potentials", ok_cls)
    # Table 1
    arcs = [(0, 2), (1, 2), (2, 3), (3, 4), (4, 0), (4, 1), (2, 5), (5, 2), (3, 6), (6, 3)]
    D = Digraph(7, arcs)
    table = [(0, {0, 2, 3, 4}, {1: "a", 5: "a", 6: "a"}), (1, {1, 2, 3, 4}, {0: "a", 5: "a", 6: "a"}),
             (2, {0, 2, 3, 4}, {1: "a", 6: "a"}), (3, {2, 3}, {0: "a", 1: "a", 5: "a"}),
             (4, {0, 2, 3, 4}, {5: "a", 6: "a"}), (5, {0, 2, 3, 4, 5, 6}, {1: "b"}),
             (6, {6}, {0: "b", 1: "b", 2: "b", 4: "b", 5: "b"})]
    N, ok1 = 0, True
    Xm = sum(1 << v for v in (2, 3, 5, 6))
    for x, S, reasons in table:
        Sm = sum(1 << v for v in S)
        ok1 &= D.path_with(Sm, end=x) is not None
        need = sorted(bits(D.full & ~Sm & ~D.out[x]))
        ok1 &= need == sorted(reasons)
        for w, r in reasons.items():
            if r == "a":
                ok1 &= D.is_path_set(Sm | (1 << w))
            else:
                ok1 &= ((D.full & ~Sm & ~(1 << w)) & ~N) == 0 and D.is_directed_path([5, 2, 3, 6]) \
                    and ((D.full & ~Xm) & ~N) == 0
        N |= 1 << x
    check("Table 1: the clearing certificate of C4(2,1,1,1) with pendant vertices", ok1 and N == D.full)


# ----------------------------------------------------------------------------------------------- Section 5
def section5():
    hdr("Section 5: small digraphs")
    logs = os.path.join(ROOT, "logs")
    data = os.path.join(ROOT, "data")
    tot = {}
    for fn in ("census_n7_part0.err", "census_n7_part1.err"):
        p = os.path.join(logs, fn)
        if os.path.exists(p):
            txt = open(p).read()
            for key in ("digraphs", "strongly connected", "non-symmetric", "with Hamiltonian path",
                        "not fully clearable", "of these passing the back-arc test"):
                v = int(txt.split(key + " ")[1].split()[0])
                tot[key] = tot.get(key, 0) + v
    if tot:
        check("Theorem 5.1(2): census totals on seven vertices",
              tot["strongly connected"] == 705422362 and tot["non-symmetric"] == 705421509
              and tot["not fully clearable"] == 1041 and tot["of these passing the back-arc test"] == 72,
              str(tot))
    for n, (sc, ns) in {3: (5, 3), 4: (83, 77), 5: (5048, 5027), 6: (1047008, 1046896)}.items():
        p = os.path.join(logs, "census_n%d.err" % n)
        if os.path.exists(p):
            txt = open(p).read()
            check("Theorem 5.1(1): n = %d, all %d non-symmetric strongly connected classes fully clearable" % (n, ns),
                  ("strongly connected %d " % sc) in txt and ("non-symmetric %d " % ns) in txt
                  and "not fully clearable 0 " in txt)
    p = os.path.join(data, "n7_survivors.txt")
    rows = [l.strip() for l in open(p) if l.strip()]
    ok = len(rows) == 72
    for line in rows:
        D = Digraph.from_directg(line)
        ok &= D.is_strongly_connected() and not D.is_symmetric() and normal_closure(D)[0] != D.full
        ok &= D.back_arc_conflict() is None and D.gamma_conflict(Fraction(7, 3)) is not None
    check("Theorem 5.1(2): each of the 72 survivors has a Gamma-conflict at z = 7/3", ok)
    p = os.path.join(data, "unilateral6_survivors.txt")
    rows = [l.strip() for l in open(p) if l.strip()]
    nconf = nset = 0
    for line in rows:
        D = Digraph.from_directg(line)
        if D.gamma_conflict(Fraction(7, 3)) is not None:
            nconf += 1
        elif D.set_potential() is not None:
            nset += 1
    check("Proposition 5.2: 145 unilateral survivors = 132 Gamma-conflicts + 13 set potentials",
          len(rows) == 145 and nconf == 132 and nset == 13)
    try:
        import networkx as nx
        from networkx.generators.atlas import graph_atlas_g
        stats, bad = {}, []
        for G in graph_atlas_g():
            n = G.number_of_nodes()
            if 1 <= n <= 6 and nx.is_connected(G):
                D = Digraph(n, [(a, b) for a, b in G.edges()] + [(b, a) for a, b in G.edges()])
                full = normal_closure(D)[0] == D.full
                stats.setdefault(n, [0, 0])
                stats[n][0] += 1
                stats[n][1] += full
                if not full:
                    bad.append(G)
        check("Remark 5.4: connected graphs on <= 5 vertices fully clearable; 108 of 112 on six; the four "
              "exceptions are trees that are neither paths nor stars",
              all(stats[n][0] == stats[n][1] for n in range(1, 6)) and stats[6] == [112, 108]
              and all(nx.is_tree(G) and max(d for _, d in G.degree()) not in (2, 5) for G in bad), str(stats))
    except ImportError:
        print("  (networkx not installed: Remark 5.4 skipped)")


# ----------------------------------------------------------------------------------------------- Section 6
def section6():
    hdr("Section 6: where clearing stops (Table 2)")
    rows = [("C4 pendants (1,1,1,1)", cycle_with_pendants(4, [1, 1, 1, 1]), 8, [4, 5, 6, 7]),
            ("C4 pendants (2,0,2,0)", cycle_with_pendants(4, [2, 0, 2, 0]), 8, [4, 5, 6, 7]),
            ("C4 legs of length 2 at c0,c2", cycle_with_legs(4, [2, 0, 2, 0]), 8, list(range(8))),
            ("C6 pendants (1,0,1,0,1,0)", cycle_with_pendants(6, [1, 0, 1, 0, 1, 0]), 9, [6, 7, 8]),
            ("C3 pendants (2,2,2)", cycle_with_pendants(3, [2, 2, 2]), 9, [3, 4, 5, 6, 7, 8]),
            ("C3 legs of length 2", cycle_with_legs(3, [2, 2, 2]), 9, list(range(9))),
            ("C5(2,1,1,1,1) pendants at x1,x2,x3", blowup_with_pendants([2, 1, 1, 1, 1], [0, 1, 1, 1, 0]), 9, [6, 7, 8]),
            ("C5 pendants (1,1,1,1,1)", cycle_with_pendants(5, [1, 1, 1, 1, 1]), 10, [5, 6, 7, 8, 9]),
            ("C4(3,1,1,1) two pendants at x1,x2", blowup_with_pendants([3, 1, 1, 1], [0, 2, 2, 0]), 10, [6, 7, 8, 9])]
    for name, D, n, outside in rows:
        N, _ = normal_closure(D)
        out = [v for v in range(D.n) if not (N >> v) & 1]
        check("%-38s n=%2d, strongly connected, set potential, no Ham. path, outside N_inf as listed" % (name, n),
              D.n == n and D.is_strongly_connected() and D.set_potential() is not None and D.ends()[D.full] == 0
              and out == outside)
    for name, D in [("C3 with one pendant per vertex", cycle_with_pendants(3, [1, 1, 1])),
                    ("C4 with pendants at two opposite vertices", cycle_with_pendants(4, [1, 0, 1, 0]))]:
        check(name + ": fully clearable (Section 5, after Theorem 5.1)", normal_closure(D)[0] == D.full)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    trials = 30 if args.quick else 150
    t0 = time.time()
    section2(rng, trials)
    section3(rng, trials)
    section4(rng, trials, args.quick)
    section5()
    section6()
    print("\n%d checks, %d failed, %.1f s" % (len(OK), OK.count(False), time.time() - t0))
    sys.exit(0 if all(OK) else 1)


if __name__ == "__main__":
    main()
