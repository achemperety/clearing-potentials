"""Core engine for Sections 25-29 (strongly connected frontier of complex potentials).

Everything here is exact and uses only the Python standard library:
  * polynomials and rational functions over Q in one variable z (class RF);
  * directed vertex geography (DVG) on small digraphs encoded by bitmasks;
  * resolvents R(W, v) in Q(z) and at rational points;
  * the set-potential test of [ALC, Thm 10.7 / 13.1] and explicit obstructions;
  * the clearing criterion of Section 3 with normality propagation, and an
    independent checker for its certificates.

Notation follows [ALC, Section 10].  A position is a pair (W, v) with v in W, W the set
of unvisited vertices; its options are (W - v, w) for out-neighbours w of v in W - v.
A set X is a *path set* if it is the vertex set of a directed path.
"""
from fractions import Fraction
import itertools

# ----------------------------------------------------------------------------------------
# Exact univariate polynomials and rational functions over Q
# ----------------------------------------------------------------------------------------

def _trim(c):
    c = list(c)
    while c and c[-1] == 0:
        c.pop()
    return tuple(c)


def p_add(a, b):
    n = max(len(a), len(b))
    return _trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def p_neg(a):
    return tuple(-x for x in a)


def p_mul(a, b):
    if not a or not b:
        return ()
    r = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            r[i + j] += x * y
    return _trim(r)


def p_divmod(a, b):
    if not b:
        raise ZeroDivisionError('polynomial division by zero')
    a = list(a)
    q = [Fraction(0)] * max(len(a) - len(b) + 1, 1)
    lb = b[-1]
    while len(a) >= len(b) and a:
        k = len(a) - len(b)
        c = a[-1] / lb
        q[k] = c
        for i, y in enumerate(b):
            a[i + k] -= c * y
        a = list(_trim(a))
    return _trim(q), tuple(a)


def p_monic(a):
    if not a:
        return a
    l = a[-1]
    return tuple(x / l for x in a)


def p_gcd(a, b):
    while b:
        _, r = p_divmod(a, b)
        a, b = b, r
    return p_monic(a)


def p_eval(a, x):
    s = Fraction(0)
    for c in reversed(a):
        s = s * x + c
    return s


ONE_P = (Fraction(1),)
Z_P = (Fraction(0), Fraction(1))


class RF:
    """Element of Q(z) in lowest terms with monic denominator."""
    __slots__ = ('num', 'den')

    def __init__(self, num, den=ONE_P, reduce=True):
        num = _trim(Fraction(x) for x in num)
        den = _trim(Fraction(x) for x in den)
        if not den:
            raise ZeroDivisionError('zero denominator')
        if reduce:
            if not num:
                den = ONE_P
            else:
                g = p_gcd(num, den)
                if len(g) > 1:
                    num, _ = p_divmod(num, g)
                    den, _ = p_divmod(den, g)
                l = den[-1]
                num = tuple(x / l for x in num)
                den = tuple(x / l for x in den)
        self.num, self.den = num, den

    @staticmethod
    def const(c):
        return RF((Fraction(c),))

    @staticmethod
    def z():
        return RF(Z_P)

    def __add__(self, o):
        if not isinstance(o, RF):
            o = RF.const(o)
        return RF(p_add(p_mul(self.num, o.den), p_mul(o.num, self.den)), p_mul(self.den, o.den))

    __radd__ = __add__

    def __neg__(self):
        return RF(p_neg(self.num), self.den, reduce=False)

    def __sub__(self, o):
        if not isinstance(o, RF):
            o = RF.const(o)
        return self + (-o)

    def __rsub__(self, o):
        return RF.const(o) - self

    def __mul__(self, o):
        if not isinstance(o, RF):
            o = RF.const(o)
        return RF(p_mul(self.num, o.num), p_mul(self.den, o.den))

    __rmul__ = __mul__

    def inv(self):
        if not self.num:
            raise ZeroDivisionError('inverse of zero')
        return RF(self.den, self.num)

    def __truediv__(self, o):
        if not isinstance(o, RF):
            o = RF.const(o)
        return self * o.inv()

    def __eq__(self, o):
        if not isinstance(o, RF):
            o = RF.const(o)
        return self.num == o.num and self.den == o.den

    def __hash__(self):
        return hash((self.num, self.den))

    def is_zero(self):
        return not self.num

    def degree(self):
        """Degree at z = infinity: deg num - deg den (None for 0)."""
        if not self.num:
            return None
        return (len(self.num) - 1) - (len(self.den) - 1)

    def leading(self):
        return self.num[-1] / self.den[-1]

    def __call__(self, x):
        d = p_eval(self.den, x)
        if d == 0:
            raise ZeroDivisionError('pole')
        return p_eval(self.num, x) / d

    def __repr__(self):
        def ps(p):
            if not p:
                return '0'
            t = []
            for i, c in enumerate(p):
                if c == 0:
                    continue
                t.append(('%s' % c) + ('' if i == 0 else ('*z' if i == 1 else '*z^%d' % i)))
            return ' + '.join(t)
        if self.den == ONE_P:
            return ps(self.num)
        return '(%s)/(%s)' % (ps(self.num), ps(self.den))


def path_poly(k):
    """Matching polynomial of the path P_k with vertex weights z: mu(P_k)(z)."""
    a, b = ONE_P, Z_P  # P_0, P_1
    if k == 0:
        return a
    for _ in range(k - 1):
        a, b = b, p_add(p_mul(Z_P, b), p_neg(a))
    return b


# ----------------------------------------------------------------------------------------
# Digraphs and directed vertex geography
# ----------------------------------------------------------------------------------------

def bits(m):
    while m:
        b = m & -m
        yield b.bit_length() - 1
        m ^= b


def popcount(m):
    return bin(m).count('1')


class Digraph:
    def __init__(self, n, arcs):
        self.n = n
        self.arcs = sorted(set((int(a), int(b)) for a, b in arcs))
        assert all(a != b and 0 <= a < n and 0 <= b < n for a, b in self.arcs)
        self.full = (1 << n) - 1
        self.out = [0] * n
        self.inn = [0] * n
        for a, b in self.arcs:
            self.out[a] |= 1 << b
            self.inn[b] |= 1 << a
        self._R = {}
        self._Rz = {}
        self._can = None

    # ---------------- structure ----------------
    def reach(self, v, adj=None):
        adj = self.out if adj is None else adj
        seen = 1 << v
        frontier = seen
        while frontier:
            x = (frontier & -frontier).bit_length() - 1
            frontier &= frontier - 1
            new = adj[x] & ~seen
            seen |= new
            frontier |= new
        return seen

    def is_strong(self):
        return self.n == 0 or (self.reach(0) == self.full and self.reach(0, self.inn) == self.full)

    def is_symmetric(self):
        return self.out == self.inn

    def is_unilateral(self):
        r = [self.reach(v) for v in range(self.n)]
        return all((r[u] >> v & 1) or (r[v] >> u & 1) for u in range(self.n) for v in range(self.n))

    def is_scc_symmetric(self):
        r = [self.reach(v) for v in range(self.n)]
        return all(not (r[b] >> a & 1) or (self.out[b] >> a & 1) for a, b in self.arcs)

    def adjacent(self, u, v):
        return bool((self.out[u] | self.inn[u]) >> v & 1)

    # ---------------- paths ----------------
    def path_table(self):
        """can[S] = bitmask of vertices v such that some directed path has vertex set S and ends at v."""
        if self._can is None:
            n = self.n
            can = [0] * (1 << n)
            for v in range(n):
                can[1 << v] |= 1 << v
            for S in range(1, 1 << n):
                e = can[S]
                while e:
                    v = (e & -e).bit_length() - 1
                    e &= e - 1
                    o = self.out[v] & ~S
                    while o:
                        w = (o & -o).bit_length() - 1
                        o &= o - 1
                        can[S | (1 << w)] |= 1 << w
            self._can = can
        return self._can

    def is_path_set(self, S):
        return self.path_table()[S] != 0

    def a_path(self, S, end=None):
        """An explicit directed path with vertex set S (ending at `end` if given), or None."""
        can = self.path_table()
        if not can[S] or (end is not None and not (can[S] >> end & 1)):
            return None
        if end is None:
            end = (can[S] & -can[S]).bit_length() - 1
        path = [end]
        cur, T = end, S
        while popcount(T) > 1:
            T2 = T & ~(1 << cur)
            cand = can[T2] & self.inn[cur]
            assert cand, 'inconsistent path table'
            prev = (cand & -cand).bit_length() - 1
            path.append(prev)
            cur, T = prev, T2
        return path[::-1]

    def hamiltonian_path(self):
        return self.a_path(self.full)

    def is_path(self, P):
        return len(set(P)) == len(P) and all(self.out[P[i]] >> P[i + 1] & 1 for i in range(len(P) - 1))

    def reachable_positions(self):
        """All reachable positions (W, v) as bitmask pairs."""
        seen = set()
        stack = [(self.full, s) for s in range(self.n)]
        while stack:
            W, v = stack.pop()
            if (W, v) in seen:
                continue
            seen.add((W, v))
            W2 = W & ~(1 << v)
            o = self.out[v] & W2
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                stack.append((W2, w))
        return seen

    # ---------------- resolvents ----------------
    def R(self, W, v):
        """Root resolvent of the position (W, v) in Q(z): 1/R = z - sum R(W - v, w)."""
        key = (W, v)
        r = self._R.get(key)
        if r is None:
            W2 = W & ~(1 << v)
            s = RF.z()
            o = self.out[v] & W2
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                s = s - self.R(W2, w)
            r = s.inv()
            self._R[key] = r
        return r

    def Rz(self, W, v, z0):
        """The same resolvent evaluated exactly at the rational point z0 (ZeroDivisionError at a pole)."""
        key = (W, v, z0)
        r = self._Rz.get(key)
        if r is None:
            W2 = W & ~(1 << v)
            s = Fraction(z0)
            o = self.out[v] & W2
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                s -= self.Rz(W2, w, z0)
            if s == 0:
                raise ZeroDivisionError('pole at z0')
            r = 1 / s
            self._Rz[key] = r
        return r

    def gamma(self, P, z0=None):
        """Product of resolvents along the play P (a directed path), in Q(z) or at z0."""
        assert self.is_path(P)
        W = self.full
        g = RF.const(1) if z0 is None else Fraction(1)
        for v in P:
            g = g * (self.R(W, v) if z0 is None else self.Rz(W, v, z0))
            W &= ~(1 << v)
        return g

    # ---------------- set potentials ----------------
    def set_potential(self):
        """Phi normalised by Phi(V) = 1 on the family of unvisited sets, or None if none exists.
        Exact over Q(z): Phi(W - v) = R(W, v) Phi(W) at every reachable position [ALC, Section 13.1]."""
        pos = sorted(self.reachable_positions(), key=lambda p: -popcount(p[0]))
        Phi = {self.full: RF.const(1)}
        for (W, v) in pos:
            val = self.R(W, v) * Phi[W]
            W2 = W & ~(1 << v)
            if W2 in Phi:
                if Phi[W2] != val:
                    return None
            else:
                Phi[W2] = val
        return Phi

    def gamma_conflict(self, z0=Fraction(7, 3)):
        """Two directed paths with the same vertex set whose products Gamma differ at z0, or None.
        A difference at one rational point proves that the rational functions differ, so a returned
        pair certifies that D has no set potential."""
        best = {}
        n = self.n
        result = [None]

        def dfs(path, mask, g, W):
            v = path[-1]
            g = g * self.Rz(W, v, z0)
            if mask in best:
                if best[mask][0] != g:
                    result[0] = (best[mask][1], list(path), best[mask][0], g)
                    return True
            else:
                best[mask] = (g, list(path))
            W2 = W & ~(1 << v)
            o = self.out[v] & W2
            while o:
                w = (o & -o).bit_length() - 1
                o &= o - 1
                path.append(w)
                if dfs(path, mask | (1 << w), g, W2):
                    return True
                path.pop()
            return False

        for s in range(n):
            if dfs([s], 1 << s, Fraction(1), self.full):
                break
        return result[0]


# ----------------------------------------------------------------------------------------
# Clearing (Section 3)
# ----------------------------------------------------------------------------------------

def clearing(D, refine=True):
    """Normality propagation of Section 3.

    Returns a dict with
      'normal': bitmask N of vertices proved to satisfy a_x = z + O(1/z) in every potential;
      'steps' : list of (x, S, kind) in the order found, where S is the vertex set of a directed path
                ending at x (the play reaching the clearing position (W, x), W = (V - S) + x), and
                kind[w] in {'path', 'normal'} records why each non-out-neighbour w in V - S is harmless;
      'delta' : None, or the path set X used to certify deg Phi(V) = n (V - X inside N at that time).
    """
    n, full = D.n, D.full
    can = D.path_table()
    normal = 0
    steps = []
    delta = full if can[full] else None  # Hamiltonian path: X = V
    changed = True
    while changed:
        changed = False
        if refine and delta is None:
            for S in range(1, full + 1):
                if can[S] and ((full & ~S) & ~normal) == 0:
                    delta = S
                    break
        for x in range(n):
            if normal >> x & 1:
                continue
            for S in range(1, full + 1):
                if not (can[S] >> x & 1):
                    continue
                rest = full & ~S
                kinds = {}
                ok = True
                for w in bits(rest & ~D.out[x]):
                    if can[S | (1 << w)]:
                        kinds[w] = 'path'
                    elif refine and delta is not None and ((rest & ~(1 << w)) & ~normal) == 0:
                        kinds[w] = 'normal'
                    else:
                        ok = False
                        break
                if ok:
                    normal |= 1 << x
                    steps.append((x, S, kinds))
                    changed = True
                    break
    return {'normal': normal, 'steps': steps, 'delta': delta, 'full': normal == full}


def check_clearing_certificate(D, steps, delta):
    """Independent check of a clearing certificate.  Returns the certified normal set.

    Each step (x, S, kinds) must satisfy: S is the vertex set of a directed path ending at x (an explicit
    path is reconstructed and checked arc by arc); every w outside S with x -/-> w is either covered by an
    explicit directed path with vertex set S + w ('path'), or ('normal') the certificate has already
    certified deg Phi(V) = n and every vertex of (V - S) - w is already certified.
    """
    n, full = D.n, D.full
    normal = 0
    delta_ok = False
    for (x, S, kinds) in steps:
        if not delta_ok and delta is not None:
            P = D.a_path(delta)
            if P is not None and D.is_path(P) and ((full & ~delta) & ~normal) == 0:
                delta_ok = True
        P = D.a_path(S, end=x)
        assert P is not None and D.is_path(P) and P[-1] == x and sum(1 << v for v in P) == S
        rest = full & ~S
        for w in bits(rest & ~D.out[x]):
            k = kinds.get(w)
            if k == 'path':
                Q = D.a_path(S | (1 << w))
                assert Q is not None and D.is_path(Q) and sum(1 << v for v in Q) == S | (1 << w)
            elif k == 'normal':
                assert delta_ok and ((rest & ~(1 << w)) & ~normal) == 0
            else:
                raise AssertionError('vertex %d not justified at step for %d' % (w, x))
        normal |= 1 << x
    return normal


# ----------------------------------------------------------------------------------------
# Families
# ----------------------------------------------------------------------------------------

def directed_cycle(m):
    return Digraph(m, [(i, (i + 1) % m) for i in range(m)])


def blowup(sizes):
    m = len(sizes)
    cls, nxt = [], 0
    for s in sizes:
        cls.append(list(range(nxt, nxt + s)))
        nxt += s
    arcs = [(a, b) for i in range(m) for a in cls[i] for b in cls[(i + 1) % m]]
    return Digraph(nxt, arcs), cls


def with_pendants(n, arcs, at, k):
    arcs = list(arcs)
    for _ in range(k):
        arcs += [(at, n), (n, at)]
        n += 1
    return n, arcs


def with_leg(n, arcs, at, length):
    arcs = list(arcs)
    prev = at
    for _ in range(length):
        arcs += [(prev, n), (n, prev)]
        prev = n
        n += 1
    return n, arcs


def cycle_with_pendants(m, counts):
    n, arcs = m, [(i, (i + 1) % m) for i in range(m)]
    for i, c in enumerate(counts):
        n, arcs = with_pendants(n, arcs, i, c)
    return Digraph(n, arcs)


def cycle_with_legs(m, lengths):
    n, arcs = m, [(i, (i + 1) % m) for i in range(m)]
    for i, L in enumerate(lengths):
        n, arcs = with_leg(n, arcs, i, L)
    return Digraph(n, arcs)


def book_with_pendants(r, k):
    """Book B_r (spine 0 -> 1, pages 2..r+1 with 1 -> p -> 0) with k pendant vertices at the spine tail 0."""
    arcs = [(0, 1)] + [(1, p) for p in range(2, r + 2)] + [(p, 0) for p in range(2, r + 2)]
    n, arcs = with_pendants(r + 2, arcs, 0, k)
    return Digraph(n, arcs)


def one_blown_with_pendants(m, r):
    """C_m(r,1,...,1) with r-1 pendant vertices at each of x_1..x_{m-2} [ALC, Thm 13.12 with e = 0]."""
    D0, cls = blowup([r] + [1] * (m - 1))
    n, arcs = D0.n, D0.arcs
    for i in range(1, m - 1):
        n, arcs = with_pendants(n, arcs, cls[i][0], r - 1)
    return Digraph(n, arcs)


def chain_into_cycle(m, k):
    arcs = [(i, (i + 1) % m) for i in range(m)]
    prev = None
    for j in range(k):
        v = m + j
        if prev is not None:
            arcs.append((prev, v))
        prev = v
    if k:
        arcs.append((prev, 0))
    return Digraph(m + k, arcs)


def parse_directg_line(line):
    t = list(map(int, line.split()))
    n, ne = t[0], t[1]
    return n, [(t[2 + 2 * i], t[3 + 2 * i]) for i in range(ne)]


def isomorphic(D1, D2):
    """Brute-force isomorphism test for small digraphs."""
    if D1.n != D2.n or len(D1.arcs) != len(D2.arcs):
        return False
    s2 = set(D2.arcs)
    for p in itertools.permutations(range(D1.n)):
        if all((p[a], p[b]) in s2 for a, b in D1.arcs):
            return True
    return False


def sturm_count(coeffs):
    """Number of distinct real roots of a polynomial with rational coefficients (low degree first)."""
    def trim(c):
        c = list(c)
        while c and c[-1] == 0:
            c.pop()
        return c

    def rem(a, b):
        a = list(a)
        while len(a) >= len(b) and a:
            k = len(a) - len(b)
            q = a[-1] / b[-1]
            for i, y in enumerate(b):
                a[i + k] -= q * y
            a = trim(a)
        return a
    p0 = trim([Fraction(c) for c in coeffs])
    p1 = trim([i * p0[i] for i in range(1, len(p0))])
    seq = [p0, p1]
    while seq[-1] and len(seq[-1]) > 1:
        r = rem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def changes(signs):
        s = [x for x in signs if x != 0]
        return sum(1 for i in range(len(s) - 1) if s[i] * s[i + 1] < 0)
    at_pinf = [1 if q[-1] > 0 else -1 for q in seq if q]
    at_minf = [(1 if q[-1] > 0 else -1) * (-1 if (len(q) - 1) % 2 else 1) for q in seq if q]
    return changes(at_minf) - changes(at_pinf)
