#!/usr/bin/env python3
"""Modular cross-check of Corollary 6.6 for the 4-cycle with one pendant vertex at every cycle vertex.

Vertices: cycle c_i = i (arcs c_i -> c_{i+1}), pendant p_i = 4 + i (arcs c_i <-> p_i).
By Proposition 6.5 every potential has U[P] equal, up to rotation, to {}, {p0p2} or {p0p2,p1p3}, with pendant
weights z or roots of t^2 - z t - 1; by Proposition 6.4 the cycle weights are z + (|N_U(c_i) & P| - 1)/z.
U[P] = {} gives a real potential, impossible by [ALC, Thm 6.10].  For the other three weight patterns this
script runs over all 2^22 graphs on the pairs meeting the cycle (up to the rotations fixing the pattern) and
evaluates every identity den(R) Phi(W - v) - num(R) Phi(W) at a point (z0, alpha0) of F_p, p = 2^61 - 1, with
alpha0^2 = z0 alpha0 + 1.  Reduction mod p followed by evaluation is a ring homomorphism from
Z_(p)[z, 1/z, alpha]/(alpha^2 - z alpha - 1), so a nonzero value proves that the identity fails.

A positive control runs the same evaluator on the symmetric 4-cycle with pendants, whose standard potential
(U = underlying graph, a = z) must pass, while a perturbed U must fail.
"""
import json, random, sys, time
from scf_core import Digraph, cycle_with_pendants, popcount

p = (1 << 61) - 1


def sqrt_mod(a):
    r = pow(a % p, (p + 1) // 4, p)   # p = 3 mod 4
    return r if r * r % p == a % p else None


def poly_at(poly, x):
    s = 0
    for c in reversed(poly):
        s = (s * x + c.numerator * pow(c.denominator, p - 2, p)) % p
    return s


def identities_at(D, z0):
    pos = sorted(D.reachable_positions(), key=lambda q: popcount(q[0]))
    return [(W, v, poly_at(D.R(W, v).num, z0), poly_at(D.R(W, v).den, z0)) for (W, v) in pos]


def satisfies_all(n, eqs, U, a):
    adj = [0] * n
    for (i, j) in U:
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    memo = {0: 1}

    def mu(S):
        r = memo.get(S)
        if r is not None:
            return r
        v = (S & -S).bit_length() - 1
        rest = S & ~(1 << v)
        val = a[v] * mu(rest)
        nb = adj[v] & rest
        while nb:
            w = (nb & -nb).bit_length() - 1
            nb &= nb - 1
            val -= mu(rest & ~(1 << w))
        val %= p
        memo[S] = val
        return val
    for (W, v, num, den) in eqs:
        if (den * mu(W & ~(1 << v)) - num * mu(W)) % p:
            return False
    return True


def choose_point(seed):
    rng = random.Random(seed)
    while True:
        z0 = rng.randrange(2, p - 1)
        s = sqrt_mod(z0 * z0 + 4)
        if s is not None:
            al = (z0 + s) * pow(2, p - 2, p) % p
            assert (al * al - z0 * al - 1) % p == 0
            return z0, al


def positive_control(z0):
    arcs = []
    for i in range(4):
        arcs += [(i, (i + 1) % 4), ((i + 1) % 4, i), (i, 4 + i), (4 + i, i)]
    D = Digraph(8, arcs)
    eqs = identities_at(D, z0)
    G = [(i, (i + 1) % 4) for i in range(4)] + [(i, 4 + i) for i in range(4)]
    good = satisfies_all(8, eqs, G, [z0] * 8)
    bad = satisfies_all(8, eqs, G[1:], [z0] * 8)
    return good, bad


def main(seed=1):
    z0, al = choose_point(seed)
    good, bad = positive_control(z0)
    print('positive control: standard potential of the symmetric digraph passes: %s; perturbed U passes: %s' % (good, bad))
    D = cycle_with_pendants(4, [1, 1, 1, 1])
    eqs = identities_at(D, z0)
    zinv = pow(z0, p - 2, p)
    CC = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    CP = [(i, 4 + j) for i in range(4) for j in range(4)]
    pairs = CC + CP
    idx = {e: k for k, e in enumerate(pairs)}

    def rot(e, k):
        f = lambda x: (x + k) % 4 if x < 4 else 4 + (x - 4 + k) % 4
        a, b = f(e[0]), f(e[1])
        return (min(a, b), max(a, b))
    cases = [('U[P] = {p0p2}, a = (alpha, z, alpha, z) on P', [(4, 6)], {4: al, 5: z0, 6: al, 7: z0}, [0, 2]),
             ('U[P] = {p0p2, p1p3}, a = (alpha, z - alpha, alpha, z - alpha)', [(4, 6), (5, 7)], {4: al, 5: (z0 - al) % p, 6: al, 7: (z0 - al) % p}, [0, 2]),
             ('U[P] = {p0p2, p1p3}, a = alpha on P', [(4, 6), (5, 7)], {4: al, 5: al, 6: al, 7: al}, [0, 1, 2, 3])]
    out = {}
    for name, H, pw, stab in cases:
        t = time.time()
        seen = bytearray(1 << len(pairs))
        maps = [[idx[rot(e, k)] for e in pairs] for k in stab]
        tested = 0
        surv = []
        for mask in range(1 << len(pairs)):
            if seen[mask]:
                continue
            bl = [k for k in range(len(pairs)) if mask >> k & 1]
            for mp in maps:
                m2 = 0
                for k in bl:
                    m2 |= 1 << mp[k]
                seen[m2] = 1
            U = H + [pairs[k] for k in bl]
            a = [0] * 8
            for i in range(4):
                d = sum(1 for (x, y) in U if (x == i and y >= 4) or (y == i and x >= 4))
                a[i] = (z0 + (d - 1) * zinv) % p
            for v, w in pw.items():
                a[v] = w
            tested += 1
            if satisfies_all(8, eqs, U, a):
                surv.append(U)
        out[name] = {'graphs': tested, 'survivors': surv, 'seconds': round(time.time() - t, 1)}
        print('%-62s graphs %8d  survivors %d  (%.0fs)' % (name, tested, len(surv), time.time() - t), flush=True)
    return out


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
