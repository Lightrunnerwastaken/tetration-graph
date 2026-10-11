"""X07 — Rechteckige Abschnitte für KNESER selbst: T_{N+d} J_∞(f_K)|_{P_N}, m = f_K'∘f_K = e^x/f_K'(x) (exakt).

Koeffizienten aus gemeinsam/kneser_halbexp_koeff_4000.json (≈47 Stellen) werden als exakte Daten eines Nachbar-f behandelt
(kohärente Störung, wie C c27 / D d19); Arbeitspräzision 2.6N + 500 Bit. Norm ℓ²_ρ.
Aufruf: python x07_rechteck_kneser.py N d1,d2,…
"""
import json
import os
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
DL = sorted(int(x) for x in sys.argv[2].split(","))
ctx.prec = int(2.6 * N) + 500
t0 = time.time()
HIER = os.path.dirname(os.path.abspath(__file__))
a = json.load(open(os.path.join(HIER, "..", "gemeinsam", "kneser_halbexp_koeff_4000.json")))["a"]
M = N + max(DL)
f = arb_poly([arb(s) for s in a[: M + 1]])
fp = f.derivative()
ex = arb_poly([1 / arb(1).__class__(1) if k == 0 else arb(1) / arb.fac_ui(k) for k in range(M + 1)])
m = (ex * fp.inv_series(M + 1)).truncate(M + 1) if hasattr(fp, "inv_series") else None
if m is None:
    # Reihen-Division per Hand
    inv = [arb(0)] * (M + 1)
    inv[0] = 1 / fp[0]
    for k in range(1, M + 1):
        s = arb(0)
        for i in range(1, k + 1):
            if i <= fp.degree():
                s += fp[i] * inv[k - i]
        inv[k] = -s / fp[0]
    m = (ex * arb_poly(inv)).truncate(M + 1)
Lr = arb("0.31813150520476413531265425158766451720351761387139986692237")
Li = arb("1.33723570143068940890116214319371061253950213846051241687")
rho = (Lr * Lr + Li * Li).sqrt()
rk = [arb(1)]
for k in range(1, M + 1):
    rk.append(rk[-1] * rho)
A = arb_mat(M + 1, N + 1)
pot = arb_poly([1])
for j in range(N + 1):
    col = (pot + m * arb_poly([0] * j + [1])).truncate(M + 1)
    inv_rj = 1 / rk[j]
    for k in range(min(col.degree(), M) + 1):
        A[k, j] = col[k] * rk[k] * inv_rj
    pot = (pot * f).truncate(M + 1)
print(f"Kneser N={N}: Matrix ({time.time() - t0:.0f}s), prec {ctx.prec}, m(0) = {float(m[0].mid()):.6f}", flush=True)
JN = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        JN[i, j] = A[i, j]
G = JN.transpose() * JN
for d in range(0, max(DL) + 1):
    if d > 0:
        row = arb_mat(1, N + 1)
        for j in range(N + 1):
            row[0, j] = A[N + d, j]
        G = G + row.transpose() * row
    if d not in DL:
        continue
    Gi = G.inv()
    Gf = np.array([[float(Gi[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
    ew = np.sort(np.linalg.eigvalsh((Gf + Gf.T) / 2))[::-1]
    s = 1 / np.sqrt(ew[:4])
    print(f"   d={d}: σ_min = " + ", ".join(f"{x:.5f}" for x in s) + f"   [{time.time() - t0:.0f}s]", flush=True)

# ---- Projektionen: quadratisch Π_N f_K = J_N⁻¹T_N g und rechteckig Π̂ = R⁺ T_{N+d} g, g = J_∞f_K = e^x + m·f_K (skaliert)
import mpmath as mp
mp.mp.dps = 40
g = (ex + m * f).truncate(M + 1)
Lc = mp.mpc(str(Lr.mid().str(40, radius=False)), str(Li.mid().str(40, radius=False)))
fK = [mp.mpf(s) for s in a]
def fK_at(z):
    s_ = mp.mpc(0)
    for c in reversed(fK):
        s_ = s_ * z + c
    return s_
kreis = [0.98 * abs(Lc) * mp.expjpi(2 * mp.mpf(t) / 360) for t in range(360)]
fwerte = [fK_at(z) for z in kreis]
def bewerte(coef, name):
    c = [mp.mpf(x.mid().str(45, radius=False)) for x in coef]
    def P(z):
        s_ = mp.mpc(0)
        for cc in reversed(c):
            s_ = s_ * z + cc
        return s_
    e = max(abs(P(z) - v) for z, v in zip(kreis, fwerte))
    eL = abs(P(Lc) - Lc)
    print(f"   {name}: max_(|z|=0.98ρ)|P − f_K| = {mp.nstr(e, 4)},  |P(L) − L| = {mp.nstr(eL, 4)}", flush=True)
rhs = arb_mat(M + 1, 1)
for k in range(M + 1):
    rhs[k, 0] = g[k] * rk[k]
rq = arb_mat(N + 1, 1)
for k in range(N + 1):
    rq[k, 0] = rhs[k, 0]
xq = JN.solve(rq)
bewerte([xq[j, 0] / rk[j] for j in range(N + 1)], "quadratisch Π_N f_K")
for dd in [x for x in DL if x > 0]:
    R = arb_mat(N + 1 + dd, N + 1)
    rr = arb_mat(N + 1 + dd, 1)
    for i in range(N + 1 + dd):
        rr[i, 0] = rhs[i, 0]
        for j in range(N + 1):
            R[i, j] = A[i, j]
    Rt = R.transpose()
    xr = (Rt * R).solve(Rt * rr)
    bewerte([xr[j, 0] / rk[j] for j in range(N + 1)], f"rechteckig d={dd}")
print(f"   [{time.time() - t0:.0f}s]", flush=True)
