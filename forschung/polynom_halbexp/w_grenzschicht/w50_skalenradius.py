"""W50 — Wachstum der Eckmoden in Abhängigkeit vom Skalierungsradius: σ₁(J_N⁻¹) in ℓ²_{ρ'}, ρ' = |L|(1 − s).
Idee für die Wachstumsschranke (I): verschwindet das exponentielle Wachstum für ein ρ' < |L| (Fixpunkt ausserhalb des Kreises), so gilt in ℓ²_ρ
σ₁ ≤ (ρ/ρ')^N · σ₁(ρ') (Normäquivalenz auf P_N), also c' = −ln(1 − s), sofern σ₁(ρ') beschränkt bleibt.
Aufruf: python w50_skalenradius.py N1,N2,...
"""
import sys
import time

import numpy as np
from flint import arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
SS = [float(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 else [-0.02, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.1, 0.15]
for N in [int(v) for v in sys.argv[1].split(",")]:
    t0 = time.time()
    n = N + 1
    ctx.prec = int(1.2 * N) + 300
    a = W.kubik_a(N, coef)
    J = W.J_mat(a, N)
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    lg = np.zeros((n, n))
    sg = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            m = Ji[i, j].mid()
            if m == 0:
                lg[i, j], sg[i, j] = -np.inf, 0.0
            else:
                lg[i, j] = float(abs(m).log())
                sg[i, j] = 1.0 if m > 0 else -1.0
    k = np.arange(n)
    out = []
    for s in SS:
        lr = np.log(W.RHO * (1 - s))
        L2 = lg + lr * (k[:, None] - k[None, :])
        mx = L2[np.isfinite(L2)].max()
        G = sg * np.exp(L2 - mx)
        sv = np.linalg.svd(G, compute_uv=False)
        out.append((s, np.log(sv[0]) + mx, np.log(sv[1]) + mx))
    print(f"N={N} ({time.time()-t0:.0f}s): " + " | ".join(f"s={s:+.2f}: lnσ₁={a1:7.2f} lnσ₂={a2:7.2f}" for s, a1, a2 in out), flush=True)
