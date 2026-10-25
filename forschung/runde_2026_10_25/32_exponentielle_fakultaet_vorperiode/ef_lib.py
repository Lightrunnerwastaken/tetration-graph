"""Bibliothek: exponentielle Fakultaet EF(n) = n^EF(n-1), EF(1) = 1, modulo m.

Zwei unabhaengige Wege:
  (A) EF_mod(n, m): direkter Turm-Algorithmus ueber die Carmichael-Kette (oder Euler-phi),
      ohne jede Theorie ueber die Vorperiode.
  (B) delta(p, n): exakte Rekursion fuer delta_p(n) = v_p(EF(n) - E^_p(n)) (Satz D im Bericht),
      daraus N_formel(m) = 1 + max{n : es gibt p^k || m mit delta_p(n) < k}.
Nur Ganzzahlen. Python: C:/Users/janis/miniconda3/python.exe
"""
from math import gcd
from functools import lru_cache

EF_EXACT = {1: 1, 2: 2, 3: 9, 4: 262144}   # EF(5) = 5^262144 hat 183231 Stellen
INF = 10 ** 40      # echte Gleichheit in Z_p (delta = unendlich)
HUGE = 10 ** 15     # Platzhalter fuer Werte >= EF(5) (> 10^183000); jedes praktische k ist kleiner


def lcm(a, b):
    return a * b // gcd(a, b)


def factor(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


@lru_cache(None)
def carmichael(m):
    if m == 1:
        return 1
    L = 1
    for p, k in factor(m).items():
        if p == 2:
            l = 1 if k == 1 else (2 if k == 2 else 2 ** (k - 2))
        else:
            l = (p - 1) * p ** (k - 1)
        L = lcm(L, l)
    return L


@lru_cache(None)
def phi(m):
    r = m
    for p in factor(m):
        r = r // p * (p - 1)
    return r


def EF_mod(n, m, red=carmichael):
    """EF(n) mod m, direkt. red = carmichael oder phi (Gegenprobe)."""
    if m == 1:
        return 0
    if n <= 4:
        return EF_EXACT[n] % m
    L = red(m)
    e = EF_mod(n - 1, L, red)
    # EF(n-1) >= EF(4) = 262144 > jeder Exponent v_p(m) (m < 2^262144): Exponent kongruent mod L und gross genug
    e2 = e + L * (1 + 64 // L)
    return pow(n, e2, m)


def tilde(x):
    while x % 4 == 0:
        x //= 2
    return x


@lru_cache(None)
def Lam_t(m):
    """Kleinste Periode (Bericht 04 der Runde 18.10.): kgV der Carmichael-Kette mit gekapptem 2-Anteil."""
    L = 1
    while m > 1:
        L = lcm(L, tilde(m))
        m = carmichael(m)
    return L


@lru_cache(None)
def Lam_full(m):
    L = 1
    while m > 1:
        L = lcm(L, m)
        m = carmichael(m)
    return L


def N_brute(m, nmax=40, P=None, red=carmichael):
    """Vorperiode direkt: 1 + groesstes n <= nmax mit EF(n) != EF(n+P) mod m."""
    if P is None:
        P = Lam_t(m)
    last = 0
    for n in range(1, nmax + 1):
        if EF_mod(n, m, red) != EF_mod(n + P, m, red):
            last = n
    return last + 1


# ---------- (B) delta-Rekursion ----------

def vp(x, p):
    if x == 0:
        return INF
    c = 0
    while x % p == 0:
        x //= p
        c += 1
    return c


def ordmod(a, p, fac_pm1):
    """Ordnung von a modulo Primzahl p (p teilt a nicht); fac_pm1 = Faktorisierung von p-1."""
    o = p - 1
    for q in fac_pm1:
        while o % q == 0 and pow(a, o // q, p) == 1:
            o //= q
    return o


def EF_big(n):
    """EF(n) als exakte Zahl fuer n <= 4, sonst HUGE."""
    return EF_EXACT[n] if n <= 4 else HUGE


def mulcap(a, b):
    r = a * b
    return min(r, HUGE) if r < INF else INF


class Delta:
    """delta_p(n) fuer alle Primzahlen p <= pmax und n <= nmax, zeilenweise in n."""

    def __init__(self, pmax, nmax):
        self.pmax, self.nmax = pmax, nmax
        sieve = bytearray([1]) * (pmax + 1)
        sieve[0] = sieve[1] = 0
        for i in range(2, int(pmax ** 0.5) + 1):
            if sieve[i]:
                sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        self.primes = [i for i in range(2, pmax + 1) if sieve[i]]
        # kleinster Primteiler fuer Faktorisierung von p-1
        spf = list(range(pmax + 1))
        for i in range(2, int(pmax ** 0.5) + 1):
            if spf[i] == i:
                for j in range(i * i, pmax + 1, i):
                    if spf[j] == j:
                        spf[j] = i
        self.spf = spf
        self.d = {}   # (p, n) -> delta
        for p in self.primes:
            self.d[(p, 1)] = INF
        for n in range(2, nmax + 1):
            for p in self.primes:
                self.d[(p, n)] = self._delta(p, n)

    def fac(self, x):
        f = {}
        while x > 1:
            q = self.spf[x]
            f[q] = f.get(q, 0) + 1
            x //= q
        return f

    def _delta(self, p, n):
        if n % p == 0:
            return mulcap(vp(n, p), EF_big(n - 1))
        prev = self.d[(p, n - 1)]
        if p == 2:   # n ungerade >= 3: v2(n^2-1) - 1 + delta_2(n-1), sofern delta_2(n-1) >= 1 (immer)
            if prev < 1:
                return 0
            return min(vp(n * n - 1, 2) - 1 + prev, INF)
        f = self.fac(p - 1)
        o = ordmod(n % p, p, f)
        for r, a in self.fac(o).items():
            if self.d[(r, n - 1)] < a:
                return 0
        M = p ** 12
        t = (pow(n, p - 1, M) - 1) % M          # v_p(n^(p-1)-1), exakt solange < 12
        V = vp(t, p) if t else vp(pow(n, p - 1) - 1, p)
        return min(V + prev, INF) if prev < INF else INF

    def bad_max(self, p, k):
        """groesstes n <= nmax mit delta_p(n) < k (0, falls keines)."""
        last = 0
        for n in range(1, self.nmax + 1):
            if self.d[(p, n)] < k:
                last = n
        return last


def N_formel(m, D):
    best = 0
    for p, k in factor(m).items():
        best = max(best, D.bad_max(p, k))
    return best + 1
