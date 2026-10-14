"""Etappe 1, Teil A1: Kandidatenkonstanten auf hoher Genauigkeit, auf zwei Motorstufen (dps 120 und 160) zum Abgleich.

Aufruf:  python a_konstanten.py 120   bzw.   python a_konstanten.py 160
Ausgabe: daten/konstanten_<dps>.json  (Zahlen als Zeichenketten mit dps+10 Stellen)

Motor: fatou.gp-Fork; komplexe Argumente nur für |Im z| ≤ 0.2 (darüber bricht der Fork ab dps 60 ab, eigener Befund),
deshalb Cauchy-Kreis mit Radius 0.2 und N = 192 Knoten (Aliasfehler (0.2/2)^192 = 1e−192).
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

from motor import gpzahl, werte

DPS = int(sys.argv[1])
mp.mp.dps = DPS + 30
HIER = Path(__file__).resolve().parent
OUT = {}
t0 = time.time()


def log(*a):
    print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)


def sexp(basis, zs):
    return werte(basis, [f"sexp({gpzahl(z, DPS + 10)})" for z in zs], DPS)


def slog(basis, ys):
    return werte(basis, [f"slog({gpzahl(y, DPS + 10)})" for y in ys], DPS)


def taylor0(basis, r=mp.mpf("0.2"), N=192, K=40):
    zs = [r * mp.expjpi((2 * mp.mpf(j) + 1) / N) for j in range(N)]       # halbe Schrittweite: kein Knoten auf der reellen Achse
    w = sexp(basis, zs)
    a = []
    for k in range(K + 1):
        s = mp.fsum(w[j] * mp.expjpi(-k * (2 * mp.mpf(j) + 1) / N) for j in range(N)) / N
        a.append(mp.re(s) / r ** k)
    return a


def gl_knoten(n):
    x, wts = [], []
    for i in range(1, n + 1):
        z = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(200):
            p0, p1 = mp.mpf(1), z
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * z * p1 - (k - 1) * p0) / k
            dp = n * (z * p1 - p0) / (z * z - 1)
            dz = p1 / dp
            z -= dz
            if abs(dz) < mp.mpf(10) ** (-mp.mp.dps + 3):
                break
        x.append(z)
        wts.append(2 / ((1 - z * z) * dp * dp))
    return x, wts


def integral(basis, a, b, n):
    x, wts = gl_knoten(n)
    h, m = (b - a) / 2, (a + b) / 2
    v = sexp(basis, [m + h * xi for xi in x])
    return h * mp.fsum(wi * mp.re(vi) for wi, vi in zip(wts, v))


def S(name, v):
    OUT[name] = mp.nstr(v, DPS + 10)
    log(name, mp.nstr(v, 30))


# ------------------------------------------------------------------ Basis e
B = "exp(1)"
a = taylor0(B)
for k in range(16):
    S(f"e.taylor0_a{k}", a[k])
S("e.T(1/2)", mp.re(sexp(B, [mp.mpf(1) / 2])[0]))
S("e.T(-1/2)", mp.re(sexp(B, [-mp.mpf(1) / 2])[0]))
S("e.T(1/3)", mp.re(sexp(B, [mp.mpf(1) / 3])[0]))
S("e.slog(1/2)", mp.re(slog(B, [mp.mpf(1) / 2])[0]))
S("e.slog(2)", mp.re(slog(B, [mp.mpf(2)])[0]))
S("e.slog(e/2)", mp.re(slog(B, [mp.e / 2])[0]))
# Ableitungen der slog bei 0 aus den Taylor-Koeffizienten: T'(−1) = a1, T''(−1) = 2a2 − a1²
a1, a2, a3 = a[1], a[2], a[3]
S("e.slog'(0)", 1 / a1)
T2m1 = 2 * a2 - a1 ** 2
S("e.slog''(0)", -T2m1 / a1 ** 3)
S("e.T''(0)", 2 * a2)
# Integrale (Gauss–Legendre, zwei Knotenzahlen)
for (lo, hi, nm) in ((0, 1, "int_0^1"), (-1, 0, "int_-1^0")):
    for n in (110, 140):
        S(f"e.{nm}_GL{n}", integral(B, mp.mpf(lo), mp.mpf(hi), n))

# D und exponentielle Fakultät
Y = 262144 * mp.log(5) + mp.log(mp.log(6))
k = 0
while Y > 1:
    Y, k = mp.log(Y), k + 1
D = 2 + k + mp.re(slog(B, [Y])[0]) - 6
S("e.D_expfak", D)


def Tser(s):
    return mp.polyval(a[::-1], s)


def sexp_gross(t):
    """sexp_e(t) für t > −1 über T(s) mit s ∈ [−½, ½) aus der Taylorreihe (|s| klein genug)."""
    m = int(mp.floor(t + mp.mpf(1) / 2))
    s = t - m
    assert abs(s) < mp.mpf("0.15"), s
    v = Tser(s)
    for _ in range(m):
        v = mp.exp(v)
    return v


def ef_fort(x, kk):
    L = sexp_gross(x + kk + D - 2) - mp.log(mp.log(x + kk))
    for j in range(kk - 1, 0, -1):
        v = L / mp.log(x + j)
        if j == 1:
            return v
        L = mp.log(v)


# Taylorreihe nur auf |s| < 0.15 genau, aber volle Länge K = 40 reicht dort nicht für 160 Stellen → längere Reihe
a = taylor0(B, K=180)
OUT["_hinweis_taylor_K"] = 180
for kk in (5, 6):                     # innerster Index x + kk = 6 bzw. 7; sexp(x+kk+D−2) bleibt darstellbar
    S(f"e.expfak_F(1)_k{kk}", ef_fort(mp.mpf(1), kk))
    S(f"e.expfak_F(2)_k{kk - 1}", ef_fort(mp.mpf(2), kk - 1))
    S(f"e.expfak_0!2=F'(1)_k{kk}", mp.diff(lambda x: ef_fort(x, kk), mp.mpf(1)))

# ------------------------------------------------------------------ Basis 2
B2 = "2"
a2s = taylor0(B2, K=12)
for k in range(8):
    S(f"2.taylor0_a{k}", a2s[k])
S("2.T(1/2)", mp.re(sexp(B2, [mp.mpf(1) / 2])[0]))
S("2.T(-1/2)", mp.re(sexp(B2, [-mp.mpf(1) / 2])[0]))
for n in (110, 140):
    S(f"2.int_0^1_GL{n}", integral(B2, mp.mpf(0), mp.mpf(1), n))
    S(f"2.int_-1^0_GL{n}", integral(B2, mp.mpf(-1), mp.mpf(0), n))

(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / f"konstanten_{DPS}.json").write_text(json.dumps(OUT, indent=1), encoding="utf-8")
log("fertig")
