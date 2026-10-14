"""Modell der Singularität von θ bei q = 1 und Vergleich mit den gemessenen a_n.

Herleitung (Bericht, Abschnitt 2): nahe z = 0 ist  θ(z) = 2 − z + h(log log sexp z)  (Hauptzweige),
und  log log sexp z = c + log(1 − q) + O(1 − q),  c = ln(A/2π) + iπ/2,  A = sexp'(0),  q = e^{2πiz}.
Modellfunktion  Φ(q) = h(c + log(1 − q))  ist im Einheitskreis holomorph (h ist auf der oberen Halbebene holomorph).
b_n = [q^n] Φ(q) per Cauchy-Integral auf |q| = e^{−2πy} (kein Kneser-Motor nötig, nur h und A).
Zusätzlich Transfer-Entwicklung  a_n ≈ (1/n) Σ_k g_k (−1)^k h^{(k)}(c − ln n),  1/Γ(s) = Σ g_k s^k.
Aufruf: python p3_modell.py <basis e|2> <dps> <y> <M> <datei mit a_n> [A]
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
from kern import Fix

B, DPS, Y, M, DATEI = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
mp.mp.dps = DPS + 20
fx = Fix(mp.e if B == "e" else mp.mpf(B), K=40)
t0 = time.time()

dat = json.loads(Path(DATEI).read_text(encoding="utf-8"))
a = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in dat["a"]]
NMAX = min(len(a) - 1, M // 2)

# A = sexp'(0) aus dem Motor (Ableitung per Cauchy auf kleinem Kreis), oder als Argument
if len(sys.argv) > 6:
    A = mp.mpf(sys.argv[6])
else:
    from kern import sexp
    basis = "exp(1)" if B == "e" else B
    K = 32
    rho = mp.mpf("0.1")
    pts = [rho * mp.expjpi(2 * mp.mpf(j) / K) for j in range(K)]
    vals = sexp(basis, pts, 60)
    A = mp.re(mp.fsum(v * mp.expjpi(-2 * mp.mpf(j) / K) for j, v in enumerate(vals)) / K / rho)
print("A = sexp'(0) =", mp.nstr(A, 30), flush=True)
lb = fx.lb
# Basis b: sexp(z) ≈ 1 + A z, log_b(sexp z) ≈ A z/ln b, log_b log_b sexp z = log(A z/ln b)/ln b
# allgemein: ℓ2 = log_b(log_b(sexp z)),  log_b(sexp z) = A z/ln b·(1+O(z)),  z = i(1−q)/(2π)(1+O)
c = mp.log(1j * A / (2 * mp.pi * lb))          # ℓ2·ln b = c + log(1−q)


def Phi(q):
    return fx.h((c + mp.log(1 - q)) / lb)


y = mp.mpf(Y)
r = mp.exp(-2 * mp.pi * y)


def fft(v):
    n = len(v)
    if n == 1:
        return v[:]
    ev, od = fft(v[0::2]), fft(v[1::2])
    out = [None] * n
    for k in range(n // 2):
        tw = mp.expjpi(-2 * mp.mpf(k) / n) * od[k]
        out[k] = ev[k] + tw
        out[k + n // 2] = ev[k] - tw
    return out


vals = [Phi(r * mp.expjpi(2 * mp.mpf(j) / M)) for j in range(M)]
per = 2j * mp.pi / mp.log(fx.lam)
for j in range(1, M):
    m = min(range(-5, 6), key=lambda m: abs(vals[j] + m * per - vals[j - 1]))
    vals[j] += m * per
print(f"[{time.time()-t0:.0f}s] Φ fertig, Randsprung {mp.nstr(min(abs(vals[0]+m*per-vals[-1]) for m in range(-5,6)),3)}", flush=True)
Fq = fft(vals)
b = [Fq[n] / M / r ** n for n in range(NMAX + 1)]

# Transfer-Entwicklung mit Ableitungen von h (Cauchy auf Kreis Radius 0.5 um v0)
g = mp.taylor(lambda s: mp.rgamma(s), 0, 8)    # g[k] Koeffizienten von 1/Γ


def hders(v0, kmax=8, rad=mp.mpf("0.5"), P=64):
    pts = [v0 + rad * mp.expjpi(2 * mp.mpf(j) / P) for j in range(P)]
    hv = [fx.h(p) for p in pts]
    for j in range(1, P):
        m = min(range(-5, 6), key=lambda m: abs(hv[j] + m * per - hv[j - 1]))
        hv[j] += m * per
    return [mp.fsum(hv[j] * mp.expjpi(-2 * mp.mpf(j * k) / P) for j in range(P)) / P / rad ** k * mp.factorial(k)
            for k in range(kmax + 1)]


zeilen = []
print(" n | |a_n| | lokaler Exp. | |a_n-b_n|/|a_n| | n·(a_n-b_n)/a_n | Transfer(1 Term)/a_n | Transfer(8 Terme)/a_n | arg a_n")
for n in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
    if n > NMAX:
        continue
    loc = -(mp.log(abs(a[min(n + 1, NMAX)])) - mp.log(abs(a[n - 1]))) / (mp.log(min(n + 1, NMAX)) - mp.log(n - 1)) if n > 1 else mp.nan
    rel = abs(a[n] - b[n]) / abs(a[n])
    # Transfer: Φ(λ) = h((c − λ)/ln b), λ = −log(1−q);  Φ^{(k)}(λ) = (−1/ln b)^k h^{(k)}((c−λ)/ln b)
    v0 = (c - mp.log(n)) / lb
    hd = hders(v0)
    t1 = (-1 / lb) * hd[1] / n
    t8 = mp.fsum(g[k] * (-1 / lb) ** k * hd[k] for k in range(1, 9)) / n
    zeilen.append({"n": n, "abs_a": mp.nstr(abs(a[n]), 15), "lok_exp": mp.nstr(loc, 8), "rel_ab": mp.nstr(rel, 6),
                   "n_rel": mp.nstr(n * (a[n] - b[n]) / a[n], 8), "T1": mp.nstr(t1 / a[n], 10), "T8": mp.nstr(t8 / a[n], 12),
                   "arg": mp.nstr(mp.arg(a[n]), 8), "b": mp.nstr(b[n], 20), "a": mp.nstr(a[n], 20)})
    z = zeilen[-1]
    print(f"{n:5d} {z['abs_a']:>14} {z['lok_exp']:>10} {z['rel_ab']:>10} {z['n_rel']:>30} {z['T1']:>30} {z['T8']:>32} {z['arg']}", flush=True)
out = {"A": mp.nstr(A, 40), "c": mp.nstr(c, 40), "y": Y, "M": M, "zeilen": zeilen,
       "b": [[mp.nstr(mp.re(v), 30), mp.nstr(mp.im(v), 30)] for v in b]}
fn = Path("daten") / f"modell_{B}_{DPS}_{Y}_{M}.json"
fn.write_text(json.dumps(out, indent=0), encoding="utf-8")
print(f"[{time.time()-t0:.0f}s] geschrieben {fn}")
