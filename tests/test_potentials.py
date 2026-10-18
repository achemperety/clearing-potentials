"""Fast tests for the `potentials` package (python3 -m pytest -q, a few seconds)."""
import os
import random
import shutil
import subprocess
import sys
from fractions import Fraction

import pytest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
from potentials import (Digraph, RatFunc, Z, MatchingPoly, check_potential, directed_collapse_potential,  # noqa
                        normal_closure, check_certificate, bits)
from potentials.families import (blowup, cycle_with_pendants, directed_cycle, random_scc_symmetric,  # noqa: E402
                                 random_tournament, transitive_tournament)


def test_ratfunc_arithmetic():
    one = RatFunc.const(1)
    f = (Z * Z - one) / (Z - one)
    assert f == Z + one
    assert (one / Z).laurent_at_infinity(2)[0] == -1
    assert (f * f.inv()) == one


def test_resolvent_asymptotics():
    rng = random.Random(1)
    for _ in range(20):
        n = rng.randint(1, 5)
        D = Digraph(n, [(i, j) for i in range(n) for j in range(n) if i != j and rng.random() < 0.5])
        for (W, v) in D.reachable_positions():
            deg, co = D.R(W, v).laurent_at_infinity(3)
            assert deg == -1 and co[0] == 1 and co[1] == 0


def test_directed_collapse_potential():
    rng = random.Random(2)
    for _ in range(10):
        D = random_scc_symmetric(rng.randint(2, 6), rng)
        U, a = directed_collapse_potential(D)
        assert check_potential(D, U, a)[0]


def test_blowups_fully_clearable():
    for sizes in ([2, 1, 1], [2, 1, 2, 1], [3, 1, 1, 1], [1, 2, 1, 2]):
        D, _ = blowup(sizes)
        N, cert = normal_closure(D)
        assert N == D.full and check_certificate(D, cert) == N


def test_tournaments():
    rng = random.Random(3)
    for n in range(3, 8):
        T = random_tournament(n, rng)
        assert normal_closure(T)[0] == T.full
        assert normal_closure(transitive_tournament(n))[0] == transitive_tournament(n).full


def test_pendant_cycle_not_clearable():
    D = cycle_with_pendants(4, [1, 1, 1, 1])
    N, _ = normal_closure(D)
    assert N != D.full
    assert sorted(v for v in range(D.n) if not (N >> v) & 1) == [4, 5, 6, 7]
    assert D.set_potential() is not None


def test_incomparable_never_cleared():
    D = Digraph(3, [(0, 1), (0, 2)])          # 1 and 2 are incomparable
    N, _ = normal_closure(D)
    assert not (N >> 1) & 1 and not (N >> 2) & 1


def test_gamma_certificate_first_survivor():
    line = open(os.path.join(ROOT, "data", "n7_survivors.txt")).readline().strip()
    D = Digraph.from_directg(line)
    assert D.is_strongly_connected() and not D.is_symmetric()
    assert normal_closure(D)[0] != D.full
    assert D.back_arc_conflict() is None
    c = D.gamma_conflict(Fraction(7, 3))
    assert c is not None and sorted(c[0]) == sorted(c[1]) and c[2] != c[3]


def test_terminal_rigidity_on_cycles():
    for m in range(3, 7):
        C = directed_cycle(m)
        assert C.ends()[C.full] == C.full
        assert all(C.R(1 << x, x) == RatFunc.const(1) / Z for x in range(m))


@pytest.mark.skipif(shutil.which("cc") is None and shutil.which("gcc") is None, reason="no C compiler")
def test_census_c_compiles(tmp_path):
    cc = shutil.which("cc") or shutil.which("gcc")
    out = tmp_path / "census"
    r = subprocess.run([cc, "-O2", "-o", str(out), os.path.join(ROOT, "csrc", "census.c")], capture_output=True)
    assert r.returncode == 0, r.stderr.decode()
