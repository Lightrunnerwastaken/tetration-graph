"""Universelle Koenigs-Reihen fuer Tuerme (Runde 25.10.2026, Gebiet 30).

Kern: in der Koordinate y = L (x - alpha) wird x -> u*exp(L x) zu f(y) = lam*(exp(y) - 1), lam = L*alpha.
sigma_lam(y) = sum_k s_k(lam) y^k  (Koenigs-Funktion von f), s_1 = 1,
   s_k (1 - lam^(k-1)) = sum_{j<k} s_j lam^(j-1) j! S(k,j) / k!      (S = Stirling 2. Art)
Alle Koeffizienten exakt als Bruchzahlen (fractions.Fraction).
"""
from fractions import Fraction as Fr
from math import factorial
from functools import lru_cache


@lru_cache(maxsize=None)
def stirling2(n, k):
    if n == k:
        return 1
    if k == 0 or k > n:
        return 0
    return k * stirling2(n - 1, k) + stirling2(n - 1, k - 1)


def smul(a, b, deg):
    out = [Fr(0)] * (deg + 1)
    for i, ai in enumerate(a[:deg + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[:deg + 1 - i]):
            if bj:
                out[i + j] += ai * bj
    return out


def sk_series(K, J):
    """s_k(lam) fuer k = 1..K als Potenzreihen in lam bis Grad J (Liste von Listen)."""
    s = {1: [Fr(1)] + [Fr(0)] * J}
    for k in range(2, K + 1):
        rhs = [Fr(0)] * (J + 1)
        kf = factorial(k)
        for j in range(1, k):
            c = Fr(factorial(j) * stirling2(k, j), kf)
            sj = s[j]
            for d in range(0, J + 1 - (j - 1)):
                if sj[d]:
                    rhs[d + j - 1] += c * sj[d]
        # / (1 - lam^(k-1)) = * sum_i lam^((k-1) i)
        out = [Fr(0)] * (J + 1)
        for d in range(J + 1):
            if rhs[d]:
                e = d
                while e <= J:
                    out[e] += rhs[d]
                    e += k - 1
        s[k] = out
    return s


def phi_series(M):
    """Phi(lam) = sigma_lam(-lam) als Potenzreihe in lam bis Grad M."""
    s = sk_series(M, M)
    out = [Fr(0)] * (M + 1)
    for k in range(1, M + 1):
        sg = (-1) ** k
        for d in range(0, M + 1 - k):
            if s[k][d]:
                out[d + k] += sg * s[k][d]
    return out


def lambert_series(M):
    """Y(L) = sum (n+1)^(n-1) L^n / n!  (Y = exp(L Y))."""
    return [Fr((n + 1) ** (n - 1), factorial(n)) if n >= 1 else Fr(1) for n in range(M + 1)]


def compose(outer, inner, M):
    """outer(inner(L)), inner[0] = 0."""
    assert inner[0] == 0
    res = [Fr(0)] * (M + 1)
    pw = [Fr(1)] + [Fr(0)] * M
    for m in range(M + 1):
        if outer[m]:
            for d in range(M + 1):
                res[d] += outer[m] * pw[d]
        pw = smul(pw, inner, M)
    return res


def K_series(M):
    """K(L) = Y(L) * Phi(L Y(L)): kappa = K(log a) fuer a ≡ 1 (Turmstart a^^(-1) = 0)."""
    Y = lambert_series(M)
    lam = [Fr(0)] + Y[:M]
    Ph = phi_series(M)
    return smul(Y, compose(Ph, lam, M), M)


# ---------------- p-adische Auswertung rationaler Reihen -----------------

def vp(x, p):
    if x == 0:
        return 10 ** 9
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def vfr(c, p):
    return vp(c.numerator, p) - vp(c.denominator, p)


def eval_padic(coeffs, z, p, N):
    """sum c_m z^m p-adisch; z ganze Zahl (Vertreter mod hoher p-Potenz).
    Rueckgabe (R, E) mit Wert = R / p^E, R mod p^(N+E) (absolute Genauigkeit p^N)."""
    vz = vp(z, p)
    E = 0
    for m, c in enumerate(coeffs):
        if c:
            E = max(E, -(vfr(c, p) + m * vz))
    q = p ** (N + E)
    R = 0
    zm = 1
    W = p ** (N + E + 60)
    for m, c in enumerate(coeffs):
        if c:
            n, d = c.numerator, c.denominator
            e = vp(d, p)
            dd = d // p ** e
            t = n * zm * p ** E
            assert t % p ** e == 0
            R = (R + (t // p ** e) * pow(dd, -1, q)) % q
        zm = zm * z % W
    return R % q, E
