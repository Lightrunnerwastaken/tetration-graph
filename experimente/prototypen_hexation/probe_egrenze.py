"""e-Grenze des Fork-Verfahrens (knhexinit): Stellen, Bogen, Reihengrad, Reichweite der Bildpunkte und nächste Singularität der Abel-Reihe.
 Aufruf: probe_egrenze.py basis dps   (basis als GP-Ausdruck, z. B. exp(1) oder 2)
 Singularität: a_k reell, konjugiertes Paar w* = c + R e^{±iφ} → b_{k+1} ≈ s b_k + t b_{k−1}, Fit über Fenster der Koeffizienten.
"""
import sys, time
sys.path.insert(0, r"D:/Users/janis/Documents/Tetration/gp-tetration/src")
import mpmath as mp
from fatou_backend import FatouGP, find_default_gp_exe

basis, dps = sys.argv[1], int(sys.argv[2])
mp.mp.dps = dps + 10
t0 = time.time()
gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=dps)
INIT = "hexinit(sexp(1),0,1)"
try:
    namen = [INIT, "hexreg", "hexerr", "kp", "klam", "knc", "knr", "knrho", "knres", "poldegree(knA)"] + [f"knpar[{i}]" for i in range(1, 11)]
    w = dict(zip(namen, gp.eval_batch(basis, namen)))
    print(f"Basis {basis}, dps {dps}: {time.time() - t0:.1f} s, hexreg {w['hexreg'].real}", flush=True)
    for n in namen[2:]:
        print(f"  {n:16s} = {mp.nstr(w[n], 15)}", flush=True)
    rad = w["knpar[3]"].real * w["knr"].real
    plan = gp.eval_batch(basis, [f"knplan({mp.nstr(rad, 20)}, {dps - 3})[{i}]" for i in range(1, 8)])
    print("  knplan(rad):", [mp.nstr(v.real, 8) for v in plan], " (Kosten, K, N, rad, hmin, htop, Reichweite)", flush=True)
    grad = int(w["poldegree(knA)"].real)
    a = [v.real for v in gp.eval_batch(basis, [f"polcoeff(knA,{k})" for k in range(grad + 1)])]
finally:
    gp.close()
c, rho, knr = w["knc"].real, w["knrho"].real, w["knr"].real
b = [a[k] / rho ** k for k in range(grad + 1)]
for k in range(0, grad + 1, max(1, grad // 10)):
    print(f"  k = {k:4d}: |b_k| = {mp.nstr(abs(b[k]), 6)},  |b_k|^(-1/k) = {mp.nstr(abs(b[k]) ** (-mp.mpf(1) / k), 6) if k and b[k] else '—'}")
for lo, hi in ((grad // 6, grad // 3), (grad // 4, grad // 2), (grad // 3, 2 * grad // 3)):
    if hi - lo < 6:
        continue
    A = mp.matrix([[b[k], b[k - 1]] for k in range(lo, hi)])
    y = mp.matrix([b[k + 1] for k in range(lo, hi)])
    s, t = mp.lu_solve(A.T * A, A.T * y)
    r1 = (s + mp.sqrt(s * s + 4 * t)) / 2
    print(f"  Fenster k = {lo}..{hi}: Singularität bei w* ≈ {mp.nstr(c + 1 / r1, 8)}, |w* − c| = {mp.nstr(abs(1 / r1), 6)} = {mp.nstr(abs(1 / r1) / knr, 4)} Im p")
