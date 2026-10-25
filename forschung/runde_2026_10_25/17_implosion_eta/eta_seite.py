"""Fatou-Koordinaten von exp_eta an e (beide Seiten), Graph-Normierung.

Z = 2e/(e-y), Phi(Z) = Z - (1/3) log|Z| + sum c_k Z^{-k}   (gleiche c_k auf beiden Seiten)
untere Seite (y<e):  slog_minus(y) = Phi(Z(y)) - C_T        (T_eta(0)=1)
obere Seite (y>e):   slog_plus(y)  = Phi(Z(y)) - C_ch       (cheta(0)=e^2)
"""
import sys
from fractions import Fraction
import mpmath as mp
import sympy as sp

K = 14


def _mul(a, b, M):
    c = [Fraction(0)] * M
    for i, x in enumerate(a):
        if x:
            for j in range(M - i):
                c[i + j] += x * b[j]
    return c


def koeffizienten(K=K):
    """c_k aus q' = (1-exp(-2q))/2 = q*r(q), q = 1/Z, exakt mit Bruechen:
    1/q' - 1/q - 1 + (1/3) log r + sum c_k (q'^k - q^k) = 0, Ordnung fuer Ordnung."""
    M = K + 3
    from math import factorial
    # r = q'/q = sum_{n>=0} (-1)^n 2^n q^n/(n+1)!
    r = [Fraction((-1) ** n * 2 ** n, factorial(n + 1)) for n in range(M)]
    # 1/r
    inv = [Fraction(0)] * M
    inv[0] = Fraction(1)
    for n in range(1, M):
        inv[n] = -sum(r[j] * inv[n - j] for j in range(1, n + 1))
    # 1/q' - 1/q = (1/r - 1)/q  -> Koeffizienten inv[n+1] bei q^n
    A = [inv[n + 1] for n in range(M - 1)] + [Fraction(0)]
    A[0] -= 1
    # log r = sum, via log r = int r'/r
    rp = [(n + 1) * r[n + 1] for n in range(M - 1)] + [Fraction(0)]
    d = _mul(rp, inv, M)
    logr = [Fraction(0)] + [d[n - 1] / n for n in range(1, M)]
    for n in range(M):
        A[n] += logr[n] / 3
    qs = [Fraction(0)] + r[:M - 1]                      # q' als Reihe
    pot = [Fraction(1)] + [Fraction(0)] * (M - 1)
    cs = []
    pots = []
    for k in range(1, K + 1):
        pot = _mul(pot, qs, M)
        pots.append(pot)
    # Gleichung: A + sum_k c_k (pots[k-1] - q^k) = 0; pots[k-1]-q^k beginnt bei q^{k+1} mit Koeffizient k*r1 = -k
    c = []
    for k in range(1, K + 1):
        n = k + 1
        s_ = A[n] + sum(c[j - 1] * (pots[j - 1][n] - (1 if n == j else 0)) for j in range(1, k))
        lead = pots[k - 1][n]
        c.append(-s_ / lead)
    # Probe: Koeffizient q^0, q^1 muessen 0 sein
    assert A[0] == 0 and A[1] == 0, (A[0], A[1])
    return c


_C = None


def Phi(Z):
    global _C
    if _C is None:
        _C = [mp.mpf(c.numerator) / c.denominator for c in koeffizienten()]
    s = Z - mp.log(abs(Z)) / 3
    for k, c in enumerate(_C, 1):
        s += c * Z**(-k)
    return s


def Zv(y):
    e = mp.e
    return 2 * e / (e - y)


def phi_unten(y, n=4000):
    """Phi(Z(y)) fuer y<e: vorwaerts iterieren."""
    a = 1 / mp.e
    for _ in range(n):
        y = mp.exp(a * y)
    return Phi(Zv(y)) - n


def phi_oben(y, n=4000):
    """Phi(Z(y)) fuer y>e: rueckwaerts iterieren (y -> e ln y)."""
    for _ in range(n):
        y = mp.e * mp.log(y)
    return Phi(Zv(y)) + n


if __name__ == '__main__':
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    cs = koeffizienten()
    print('c_k =', cs[:6])
    CT = phi_unten(mp.mpf(1), n)
    CT2 = phi_unten(mp.mpf(1), 2 * n)
    Cch = phi_oben(mp.e**2, n)
    Cch2 = phi_oben(mp.e**2, 2 * n)
    print('C_T  =', mp.nstr(CT, mp.mp.dps - 5), ' Gegenprobe 2n:', mp.nstr(CT - CT2, 3))
    print('C_ch =', mp.nstr(Cch, mp.mp.dps - 5), ' Gegenprobe 2n:', mp.nstr(Cch - Cch2, 3))
    print('C_ch - C_T =', mp.nstr(Cch - CT, mp.mp.dps - 5))
