"""Exact arithmetic in Q[z] and Q(z).

Polynomials are tuples of Fractions, lowest degree first, with no trailing zeros
(the zero polynomial is the empty tuple).  Rational functions are kept in lowest
terms with a monic denominator, so equality of RatFunc objects is equality in Q(z).
Only the Python standard library is used.
"""
from fractions import Fraction

__all__ = ["Poly", "RatFunc", "Z", "ONE", "path_matching_poly"]


def _trim(c):
    c = [Fraction(x) for x in c]
    while c and c[-1] == 0:
        c.pop()
    return tuple(c)


class Poly:
    """Univariate polynomial over Q."""

    __slots__ = ("c",)

    def __init__(self, coeffs=()):
        self.c = _trim(coeffs)

    # -- basic structure -------------------------------------------------
    def deg(self):
        return len(self.c) - 1  # -1 for the zero polynomial

    def is_zero(self):
        return not self.c

    def lead(self):
        return self.c[-1]

    def __eq__(self, other):
        return isinstance(other, Poly) and self.c == other.c

    def __hash__(self):
        return hash(self.c)

    # -- ring operations -------------------------------------------------
    def __add__(self, o):
        a, b = self.c, o.c
        n = max(len(a), len(b))
        return Poly([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])

    def __neg__(self):
        return Poly([-x for x in self.c])

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        a, b = self.c, o.c
        if not a or not b:
            return Poly()
        r = [Fraction(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] += x * y
        return Poly(r)

    def scale(self, s):
        return Poly([s * x for x in self.c])

    def divmod(self, o):
        if o.is_zero():
            raise ZeroDivisionError("polynomial division by zero")
        r = list(self.c)
        q = [Fraction(0)] * max(len(r) - len(o.c) + 1, 1)
        lo = o.c[-1]
        while len(r) >= len(o.c) and r:
            k = len(r) - len(o.c)
            f = r[-1] / lo
            q[k] = f
            for i, y in enumerate(o.c):
                r[i + k] -= f * y
            r = list(_trim(r))
        return Poly(q), Poly(r)

    def monic(self):
        if self.is_zero():
            return self
        return self.scale(1 / self.c[-1])

    def gcd(self, o):
        a, b = self, o
        while not b.is_zero():
            a, b = b, a.divmod(b)[1]
        return a.monic()

    def __call__(self, x):
        s = Fraction(0)
        for coef in reversed(self.c):
            s = s * x + coef
        return s

    def __repr__(self):
        if not self.c:
            return "0"
        parts = []
        for i, x in enumerate(self.c):
            if x:
                parts.append(str(x) + ("" if i == 0 else ("*z" if i == 1 else "*z^%d" % i)))
        return " + ".join(parts)


class RatFunc:
    """Element of Q(z), stored in lowest terms with monic denominator."""

    __slots__ = ("num", "den")

    def __init__(self, num, den=None, _reduced=False):
        num = num if isinstance(num, Poly) else Poly(num)
        den = Poly((1,)) if den is None else (den if isinstance(den, Poly) else Poly(den))
        if den.is_zero():
            raise ZeroDivisionError("zero denominator")
        if not _reduced:
            if num.is_zero():
                den = Poly((1,))
            else:
                g = num.gcd(den)
                if g.deg() > 0:
                    num = num.divmod(g)[0]
                    den = den.divmod(g)[0]
                lc = den.lead()
                num, den = num.scale(1 / lc), den.scale(1 / lc)
        self.num, self.den = num, den

    @staticmethod
    def const(c):
        return RatFunc(Poly((Fraction(c),)))

    @staticmethod
    def var():
        return RatFunc(Poly((0, 1)))

    def _coerce(self, o):
        return o if isinstance(o, RatFunc) else RatFunc.const(o)

    def __add__(self, o):
        o = self._coerce(o)
        return RatFunc(self.num * o.den + o.num * self.den, self.den * o.den)

    __radd__ = __add__

    def __neg__(self):
        return RatFunc(-self.num, self.den, _reduced=True)

    def __sub__(self, o):
        return self + (-self._coerce(o))

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __mul__(self, o):
        o = self._coerce(o)
        return RatFunc(self.num * o.num, self.den * o.den)

    __rmul__ = __mul__

    def inv(self):
        if self.num.is_zero():
            raise ZeroDivisionError("inverse of zero")
        return RatFunc(self.den, self.num)

    def __truediv__(self, o):
        return self * self._coerce(o).inv()

    def __rtruediv__(self, o):
        return self._coerce(o) * self.inv()

    def __eq__(self, o):
        o = self._coerce(o)
        return self.num == o.num and self.den == o.den

    def __hash__(self):
        return hash((self.num, self.den))

    def is_zero(self):
        return self.num.is_zero()

    def degree(self):
        """Exponent of the leading term at z = infinity (None for 0)."""
        if self.num.is_zero():
            return None
        return self.num.deg() - self.den.deg()

    def laurent_at_infinity(self, terms):
        """First `terms` coefficients of the expansion in powers of 1/z, starting at z^degree."""
        # num/den = z^(dn-dd) * (sum a_i z^-i)/(sum b_i z^-i); invert the series
        a = list(reversed(self.num.c))
        b = list(reversed(self.den.c))
        out = []
        rem = a + [Fraction(0)] * terms
        for k in range(terms):
            c = rem[k] / b[0]
            out.append(c)
            for i, y in enumerate(b):
                if k + i < len(rem):
                    rem[k + i] -= c * y
        return self.degree(), out

    def __call__(self, x):
        d = self.den(x)
        if d == 0:
            raise ZeroDivisionError("pole")
        return self.num(x) / d

    def __repr__(self):
        if self.den.c == (1,):
            return "(%r)" % self.num
        return "(%r)/(%r)" % (self.num, self.den)


Z = RatFunc.var()
ONE = RatFunc.const(1)


def path_matching_poly(k):
    """mu(P_k)(z), the matching polynomial of the path with k vertices (mu(P_0) = 1)."""
    a, b = Poly((1,)), Poly((0, 1))
    if k == 0:
        return a
    for _ in range(k - 1):
        a, b = b, Poly((0, 1)) * b - a
    return b
