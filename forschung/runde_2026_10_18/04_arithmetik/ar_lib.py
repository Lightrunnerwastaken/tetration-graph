"""Exakte Ganzzahl-Werkzeuge fuer Tuerme modulo m und p-adisch (Runde 18.10.2026, 04_arithmetik).

Nur Python-Ganzzahlen; keine Gleitkommazahlen. Alle Funktionen sind deterministisch.
"""
from functools import lru_cache
from math import gcd
import sys
sys.setrecursionlimit(100000)
from sympy import factorint

BIG = 1 << 64  # Schwelle: groessere Tuerme gelten als "riesig" (> jeder Exponentenschwelle)


@lru_cache(maxsize=None)
def fac(m):
    return tuple(sorted(factorint(m).items()))


@lru_cache(maxsize=None)
def carmichael(m):
    if m == 1:
        return 1
    l = 1
    for p, k in fac(m):
        if p == 2:
            v = 1 if k == 1 else (2 if k == 2 else 1 << (k - 2))
        else:
            v = (p - 1) * p ** (k - 1)
        l = l * v // gcd(l, v)
    return l


def nu(m):
    """groesster Primzahlexponent von m (Schwelle fuer a^e mod m)."""
    return max((k for _, k in fac(m)), default=0)


@lru_cache(maxsize=None)
def tower_exact(a, n):
    """a^^n exakt, falls < BIG, sonst None. a^^0 = 1, a^^(-1) = 0."""
    if n == -1:
        return 0
    if n == 0:
        return 1
    if a >= 2 and n >= 6:
        return None  # 2^^5 = 2^65536 > BIG
    e = tower_exact(a, n - 1)
    if e is None:
        return None if a >= 2 else 1
    if a <= 1:
        return a ** e
    if e * a.bit_length() > 70:
        return None
    v = a ** e
    return v if v < BIG else None


def pow_big(a, e_mod, m):
    """a^e mod m, wobei e riesig (>= nu(m)) und e ≡ e_mod (mod carmichael(m))."""
    if m == 1:
        return 0
    lam = carmichael(m)
    e = e_mod % lam + lam * (nu(m) + 1)
    return pow(a, e, m)


@lru_cache(maxsize=200000)
def tower_mod(a, n, m):
    """a^^n mod m exakt (n >= -1). Iterativ ueber die Carmichael-Kette m, lambda(m), lambda^2(m), ..."""
    if m == 1:
        return 0
    ex = tower_exact(a, n)
    if ex is not None:
        return ex % m
    chain = [m]
    # Stufe i: a^^(n-i) mod chain[i]; abbrechen, wenn Modul 1 oder Turm exakt bekannt
    while True:
        i = len(chain) - 1
        if chain[i] == 1 or tower_exact(a, n - i - 1) is not None:
            break
        chain.append(carmichael(chain[i]))
    i = len(chain) - 1
    if chain[i] == 1:
        r = 0  # a^^(n-i) mod 1
        start = i - 1
        val = r
        for k in range(start, -1, -1):
            # a^^(n-k) mod chain[k], Exponent a^^(n-k-1) riesig, Rest val mod chain[k+1]
            val = pow_big(a, val, chain[k])
        return val
    e = tower_exact(a, n - i - 1)
    val = pow(a, e, chain[i])
    for k in range(i - 1, -1, -1):
        val = pow_big(a, val, chain[k])
    return val


def tower_mod_height_big(a, H, m):
    """a^^H mod m fuer eine (beliebig grosse, auch nur symbolisch bekannte) Hoehe H >= stab_height+1.
    Nutzt: ab der Stabilisierungshoehe ist der Rest konstant."""
    t = stab_height(a, m)
    return tower_mod(a, max(t, 0), m) if H >= t else tower_mod(a, H, m)


@lru_cache(maxsize=None)
def stab_height(a, m):
    """kleinstes t >= -1 mit a^^j ≡ a^^t (mod m) fuer alle j >= t.
    Beweisbar korrekt: konstant ab Hoehe h, sobald Exponent mod carmichael(m) konstant ist
    und die Exponenten riesig sind; wir pruefen bis zur rekursiven Schranke."""
    if m == 1:
        return -1
    # obere Schranke: Stabilisierungshoehe mod lambda(m) plus 1, plus Schwelle fuer riesige Exponenten
    tl = stab_height(a, carmichael(m))
    hi = max(tl + 1, 0)
    # ab hi: Exponent a^^(j-1) mod lambda konstant; ist er auch >= nu(m)? sicherstellen
    while True:
        e = tower_exact(a, hi - 1) if hi >= 1 else None
        if hi >= 1 and (e is None or e >= nu(m)):
            break
        hi += 1
    # ab hi konstant; minimale Hoehe rueckwaerts suchen
    v = tower_mod(a, hi, m)
    t = hi
    while t - 1 >= -1 and tower_mod(a, t - 1, m) == v:
        t -= 1
    return t


def vp(x, p):
    if x == 0:
        return 10 ** 9
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def order_mod(a, p):
    o, x = 1, a % p
    while x != 1:
        x = x * a % p
        o += 1
    return o


# ---------------- p-adische Hilfen (Arithmetik modulo p^N mit Bewertungsbuchhaltung) -----------

def teichmueller(a, p, N):
    q = p ** N
    if p == 2:
        return 1 if a % 4 == 1 else q - 1
    return pow(a, p ** N, q)


def padic_log_unit1(u, p, N):
    """log_p(u) mod p^N fuer u ≡ 1 mod p (p ungerade) bzw. mod 4 (p=2). Rueckgabe als ganze Zahl mod p^N."""
    # sum (-1)^{n+1} (u-1)^n / n ; Glied n hat Bewertung >= n*v - v_p(n)
    v = vp(u - 1, p)
    extra = 2 * (N.bit_length() + 5)
    M = N + extra
    q = p ** M
    x = (u - 1) % q
    s = 0
    xn = 1
    n = 0
    while True:
        n += 1
        xn = xn * x % q
        if n * v - vp(n, p) >= N + 2 and n * v - (n.bit_length() + 2) > N + 2:
            break
        s_n = vp(n, p)
        num = xn
        # num ist durch p^{n v} teilbar; Division durch n = p^s * w
        w = n // p ** s_n
        term = (num // p ** s_n) * pow(w, -1, q) % q
        s = s + term if n % 2 == 1 else s - term
    return s % p ** N


def padic_log(a, p, N):
    """Iwasawa-log: log_p(a) = log(<a>) mit a = omega(a)<a>."""
    if p == 2:
        u = a % (1 << (N + 8))
        if u % 4 == 3:
            u = (-a) % (1 << (N + 8))
        return padic_log_unit1(u, 2, N)
    q = p ** (N + 4)
    u = pow(a, p - 1, q)
    L = padic_log_unit1(u, p, N + 4)
    return L * pow(p - 1, -1, q) % p ** N


def lambert_tower(y, p, N):
    """sum_{n>=0} (n+1)^{n-1} y^n / n! mod p^N, fuer v_p(y) >= 1 (p ungerade) bzw. >= 2 (p=2).
    Loest Y = exp(y Y)."""
    v = vp(y, p) if y % p ** N else 10 ** 6
    extra = 4 * N // max(1, (p - 1)) + 20
    M = N + extra
    q = p ** M
    s = 1
    yn = 1
    nf_unit = 1  # Einheitsteil von n!
    nf_v = 0
    n = 0
    while True:
        n += 1
        yn = yn * y % q
        k = n
        while k % p == 0:
            k //= p
            nf_v += 1
        nf_unit = nf_unit * k % q
        if n * v - nf_v > N + 5 and n * min(v, 50) - n // (p - 1 if p > 2 else 1) > N + 5:
            break
        num = pow(n + 1, n - 1, q) * yn % q
        # num/p^{nf_v}
        assert num % p ** nf_v == 0 or num == 0
        term = (num // p ** nf_v) * pow(nf_unit, -1, q) % q
        s = (s + term) % q
    return s % p ** N


def alpha_series(a, p, N):
    """p-adischer Turmgrenzwert ueber die Lambert-Reihe (p teilt a nicht)."""
    q = p ** N
    L = padic_log(a, p, N + 5)
    if p == 2:
        u = 1 if a % 4 == 1 else -1
    else:
        c = tower_mod(a, 60, p - 1)  # stabiler Rest mod p-1 (Hoehe 60 >> Stabilisierung)
        u = pow(teichmueller(a, p, N + 5), c, p ** (N + 5))
    y = u * L % p ** (N + 5)
    Y = lambert_tower(y, p, N + 3)
    return u * Y % q


def tower_run(a, p, N, H0, H1):
    """Liste [a^^H mod p^N fuer H = H0..H1], schnell: a^x mod p^N haengt nur von x mod (p-1) und x mod p^(N-1) ab.
    Exakt, sofern H0-1 >= stab_height(a, p-1) und die Exponenten riesig sind (H0 >= 6)."""
    q = p ** N
    assert H0 >= 6 and H0 - 1 >= stab_height(a, max(p - 1, 1))
    x = tower_mod(a, H0, q)
    out = [x]
    if p == 2:
        mod_e = 1 << max(N - 2, 1)
        for H in range(H0 + 1, H1 + 1):
            e = x % mod_e + mod_e * (N + 1)
            x = pow(a, e, q); out.append(x)
        return out
    c = tower_mod(a, H0, p - 1)  # stabiler Rest der Exponenten mod p-1
    pe = p ** (N - 1)
    lam = (p - 1) * pe
    inv = pow(pe, -1, p - 1) if p > 2 else 0
    for H in range(H0 + 1, H1 + 1):
        r = x % pe
        # CRT: e ≡ r (mod p^(N-1)), e ≡ c (mod p-1)
        k = ((c - r) * inv) % (p - 1)
        e = r + pe * k + lam * (N + 1)
        x = pow(a, e, q); out.append(x)
    return out


def alpha_fast(a, p, N):
    H0 = max(6, stab_height(a, max(p - 1, 1)) + 1)
    seq = tower_run(a, p, N, H0, H0 + N + 5)  # Bewertung waechst mindestens um 1 pro Stufe
    return seq[-1]
