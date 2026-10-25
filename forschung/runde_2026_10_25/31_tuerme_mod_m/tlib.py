"""Bibliothek Gebiet 31: Tuerme a^^n mod m, Stabilisierungshoehe, Carmichael-Kette.

Konvention (wie Runde 18.10., Bericht 04): a^^(-1) = 0, a^^0 = 1, a^^(n+1) = a^(a^^n).
t(a,m) = kleinstes t >= -1 mit a^^n ≡ a^^t (mod m) fuer alle n >= t.
l(m)   = Laenge der Carmichael-Kette m -> lambda(m) -> ... -> 1 (l(1)=0, l(2)=1).
"""
from functools import lru_cache
from math import gcd

LIM = 10 ** 6
_spf = list(range(LIM + 1))
for i in range(2, int(LIM ** 0.5) + 1):
    if _spf[i] == i:
        for j in range(i * i, LIM + 1, i):
            if _spf[j] == j:
                _spf[j] = i


def fac(n):
    out = []
    while n > 1:
        p = _spf[n] if n <= LIM else _slowp(n)
        k = 0
        while n % p == 0:
            n //= p; k += 1
        out.append((p, k))
    return out


def _slowp(n):
    d = 2
    while d * d <= n:
        if n % d == 0:
            return d
        d += 1
    return n


@lru_cache(maxsize=None)
def lam(m):
    if m <= 2:
        return 1
    l = 1
    for p, k in fac(m):
        if p == 2:
            v = 1 if k == 1 else (2 if k == 2 else 1 << (k - 2))
        else:
            v = (p - 1) * p ** (k - 1)
        l = l * v // gcd(l, v)
    return l


@lru_cache(maxsize=None)
def ell(m):
    return 0 if m == 1 else 1 + ell(lam(m))


def nu(m):
    return max((k for _, k in fac(m)), default=0)


def Lam(m):
    """kgV der Kette m, lambda(m), lambda^2(m), ..."""
    L = 1
    while m > 1:
        L = L * m // gcd(L, m); m = lam(m)
    return L


# ---------------- direkte Rechnung (verallgemeinerter Euler) ----------------
CAP = 200


def capped(a, n):
    """min(a^^n, CAP)"""
    if n == -1:
        return 0
    if n == 0:
        return 1
    e = capped(a, n - 1)
    if e >= CAP:
        return CAP if a >= 2 else 1
    if a == 1:
        return 1
    if a == 0:
        return 1 if e == 0 else 0
    v = 1
    for _ in range(e):
        v *= a
        if v >= CAP:
            return CAP
    return v


def tower_mod(a, n, m):
    if m == 1:
        return 0
    if n == -1:
        return 0
    if n == 0:
        return 1 % m
    E = capped(a, n - 1)
    if E < CAP:
        return pow(a, E, m)
    L = lam(m)
    e = tower_mod(a, n - 1, L)
    e2 = e + L * (CAP // L + 1)  # e2 ≡ a^^(n-1) mod L und e2 >= CAP >= nu(m)
    return pow(a, e2, m)


def t_direct(a, m, extra=8):
    H = ell(m) + extra
    u = [tower_mod(a, n, m) for n in range(-1, H + 1)]
    assert len(set(u[-4:])) == 1, (a, m, u)
    t = H
    while t - 1 >= -1 and u[t - 1 + 1] == u[-1]:
        t -= 1
    return t


# ---------------- Rekursion ueber Ordnungen (Satz R) ----------------
def ordmod(a, q):
    """Ordnung von a mod q (a Einheit)."""
    L = lam(q)
    o = L
    for p, _ in fac(L):
        while o % p == 0 and pow(a, o // p, q) == 1:
            o //= p
    return o


def term_value(a, p, k, large):
    """t(a,p^k) fuer p | a."""
    q = p ** k
    if a % q == 0:
        return 1
    v = 0; x = a
    while x % p == 0:
        x //= p; v += 1
    if large:
        return 2
    n = 2
    while v * capped(a, n - 1) < k:
        n += 1
    return n


def t_rec(a, m, large=False):
    """t(a,m) ueber die Ordnungs-Rekursion. large=True: Endknoten mit p|a wie fuer grosse a (Wert 1 oder 2)."""
    if m == 1:
        return -1
    best = -10
    for p, k in fac(m):
        if a % p == 0:
            v = term_value(a, p, k, large)
        else:
            v = 1 + t_rec(a, ordmod(a % p ** k, p ** k), large)
        best = max(best, v)
    return best
