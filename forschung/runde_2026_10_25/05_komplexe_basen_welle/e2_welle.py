"""E2: Fourier-Koeffizienten der oberen und unteren Welle für eine komplexe Basis, und Vergleich mit dem Koeffizientengesetz.

θ_U(z) = h_+(sexp z) − z = Σ a_n q^n,   q = e^{2πiz},  abgetastet auf Im z = +y
θ_L(z) = h_−(sexp z) − z = Σ d_n p^n,   p = e^{−2πiz}, abgetastet auf Im z = −y
Modell (Singularität bei z = −1 bzw. 0, q = 1):
  b_n = [q^n] h_+((c_U + log(1−q))/ln b),  c_U = log(+iA/(2π ln b))
  e_n = [p^n] h_−((c_L + log(1−p))/ln b),  c_L = log(−iA/(2π ln b))
Aufruf: python e2_welle.py "<basis>" <dps> <y> <N> [m]   (m: Verschiebung, h(sexp z) = h(sexp(z−m)) + m; Fixpunkt-Startwerte aus sexp(±4i))
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
from kern import Fix, sexp, schliessen, _tag, HIER, fft, entwirren

basis, DPS, Y, N = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
MSH = int(sys.argv[5]) if len(sys.argv) > 5 else 0   # Verschiebung: h(sexp z) = h(sexp(z-m)) + m
WELLEN = sys.argv[6] if len(sys.argv) > 6 else "UL"   # welche Wellen
mp.mp.dps = DPS + 20
t0 = time.time()
y = mp.mpf(Y)
v = sexp(basis, [mp.mpc(1), mp.mpc(0, 4), mp.mpc(0, -4)], DPS)
b = v[0]
lb = mp.log(b)
sp, sm = v[1], v[2]
fp, fm = Fix(b, sp), Fix(b, sm)
print(f"b = {mp.nstr(b, 15)}  L+ = {mp.nstr(fp.L, 15)} |λ+| = {mp.nstr(abs(fp.lam), 8)} k+={fp.k}   L- = {mp.nstr(fm.L, 15)} |λ-| = {mp.nstr(abs(fm.lam), 8)} k-={fm.k}")
print(f"   |sexp(4i)-L+| = {mp.nstr(abs(v[1]-fp.L), 3)}   |sexp(-4i)-L-| = {mp.nstr(abs(v[2]-fm.L), 3)}")
# A = sexp'(0)
K = 32
rho = mp.mpf("0.1")
cv = sexp(basis, [rho * mp.expjpi(2 * mp.mpf(j) / K) for j in range(K)], DPS)
A = mp.fsum(cv[j] * mp.expjpi(-2 * mp.mpf(j) / K) for j in range(K)) / K / rho
print("A = sexp'(0) =", mp.nstr(A, 20), flush=True)

ts = [mp.mpf(j) / N - mp.mpf(1) / 2 for j in range(N)]
res = {}
for nm, f, sg in (("U", fp, 1), ("L", fm, -1)):
    if nm not in WELLEN:
        continue
    zs = [t + sg * 1j * y for t in ts]
    ws = sexp(basis, [z - MSH for z in zs], DPS)
    th = entwirren([f.h(w) + MSH - z for w, z in zip(ws, zs)], f.per)
    schritt = max(abs(th[j] - th[j - 1]) for j in range(1, N))
    rand = min(abs(th[0] + m * f.per - th[-1]) for m in range(-6, 7))
    F = fft(th)
    # θ_U: e^{-2πi n t_j} -> a_n = (−1)^n F[n]/N e^{2πny};  θ_L in p = e^{−2πiz}: Koeffizient bei e^{+2πi n t} -> F[N−n]
    if sg == 1:
        a = [(-1) ** n * F[n] / N * mp.exp(2 * mp.pi * n * y) for n in range(N // 2)]
    else:
        a = [(-1) ** n * F[(N - n) % N] / N * mp.exp(2 * mp.pi * n * y) for n in range(N // 2)]
    # Modell
    c = mp.log(sg * 1j * A / (2 * mp.pi * lb))
    r = mp.exp(-2 * mp.pi * y)
    ph = entwirren([f.h((c + mp.log(1 - r * mp.expjpi(2 * mp.mpf(j) / N))) / lb) for j in range(N)], f.per)
    G = fft(ph)
    bm = [G[n] / N / r ** n for n in range(N // 2)]
    print(f"[{time.time()-t0:.0f}s] Welle {nm}: max Schritt {mp.nstr(schritt, 3)}, Randsprung {mp.nstr(rand, 3)}, a_0 = {mp.nstr(a[0], 15)}", flush=True)
    print("    n |  |a_n|          | lok.Exp | |a_n-b_n|/|a_n| | n ln n |a_n| | arg a_n")
    nmax_ok = int(DPS * mp.log(10) / (2 * mp.pi * y)) - 2
    for n in [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 160, 200, 250]:
        if n >= N // 2 or n > nmax_ok:
            continue
        loc = -(mp.log(abs(a[n + 1])) - mp.log(abs(a[n - 1]))) / (mp.log(n + 1) - mp.log(n - 1)) if n > 1 else mp.nan
        print(f"  {n:4d} {mp.nstr(abs(a[n]), 10):>14} {mp.nstr(loc, 5):>8} {mp.nstr(abs(a[n]-bm[n])/abs(a[n]), 4):>12} "
              f"{mp.nstr(n*mp.log(n)*abs(a[n]), 6):>10} {mp.nstr(mp.arg(a[n]), 5)}")
    res[nm] = {"L": mp.nstr(f.L, 30), "lam": mp.nstr(f.lam, 20), "c": mp.nstr(c, 20), "schritt": mp.nstr(schritt, 5), "rand": mp.nstr(rand, 5),
               "a": [[mp.nstr(mp.re(x), 25), mp.nstr(mp.im(x), 25)] for x in a],
               "b": [[mp.nstr(mp.re(x), 25), mp.nstr(mp.im(x), 25)] for x in bm]}
out = {"basis": basis, "dps": DPS, "y": Y, "N": N, "b": mp.nstr(b, 30), "A": mp.nstr(A, 30), "wellen": res}
(HIER / "daten").mkdir(exist_ok=True)
fn = HIER / "daten" / f"welle_{_tag(basis)}_{DPS}_{Y}_{N}_m{MSH}_{WELLEN}.json"
fn.write_text(json.dumps(out, indent=0), encoding="utf-8")
print(f"[{time.time()-t0:.0f}s] geschrieben {fn.name}")
schliessen()
