#!/usr/bin/env python3
"""Certify the census survivors of Theorem 5.1(2) and Proposition 5.2.

Input: digraphs in nauty's `directg -T` text format, one per line (the stdout of csrc/census.c).
For each digraph the script re-checks, with code independent of census.c, that it is strongly connected
(or unilateral, with --unilateral) and not (SCC-)symmetric, that it is not fully clearable and that it
passes the back-arc test; it then searches for two directed paths with the same vertex set whose products
of resolvents differ at z = 7/3 (exact rational arithmetic; Lemma 2.2 makes one evaluation a proof).
Digraphs without such a pair are tested for a set potential exactly in Q(z).

    python3 scripts/certify_survivors.py data/n7_survivors.txt --out data/n7_gamma_certificates.json
    python3 scripts/certify_survivors.py data/unilateral6_survivors.txt --unilateral
"""
import argparse
import json
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from potentials import Digraph, normal_closure  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--unilateral", action="store_true", help="survivors of csrc/census_unilateral.c")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rows = [l.strip() for l in open(args.infile) if l.strip()]
    certs, setpot = [], []
    for line in rows:
        D = Digraph.from_directg(line)
        if args.unilateral:
            assert D.is_unilateral() and not D.is_scc_symmetric(), line
        else:
            assert D.is_strongly_connected() and not D.is_symmetric(), line
        N, _ = normal_closure(D)
        assert N != D.full, "fully clearable: " + line
        assert D.back_arc_conflict() is None, "fails the back-arc test: " + line
        c = D.gamma_conflict(Fraction(7, 3))
        if c is not None:
            P1, P2, g1, g2 = c
            assert sorted(P1) == sorted(P2) and g1 != g2
            certs.append({"digraph": line, "path1": P1, "path2": P2,
                          "gamma1_at_7/3": str(g1), "gamma2_at_7/3": str(g2)})
        else:
            assert D.set_potential() is not None, "neither a conflict nor a set potential: " + line
            setpot.append(line)
    print("%d survivors: %d certified by Gamma-conflicts at z = 7/3, %d with an exact set potential"
          % (len(rows), len(certs), len(setpot)))
    for line in setpot:
        print("  set potential:", line)
    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"conflicts": certs, "set_potentials": setpot}, fh, indent=1)


if __name__ == "__main__":
    main()
