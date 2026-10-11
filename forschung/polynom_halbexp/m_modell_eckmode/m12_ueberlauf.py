"""M12 — Überlauf der Eckmoden: wohin bildet J_∞ die Eck-Ausgänge ab?
x := J_N⁻¹ v_i (exakt in arb, v_i aus m07, float) ist ein Polynom vom Grad N mit T_N J_∞ x = v_i und ‖x‖ ≈ σ_i.
Gemessen: Koeffizienten von J_∞ x in den Zeilen N+1 … N+D (exakt), ihr Gewicht relativ zu ‖x‖, und das „Überlauf-Polynom“ q
(J_∞ x ≈ z^{N+1} q). Normiertes radiales Modell, Exp-Lage.
Aufruf: python m12_ueberlauf.py "N1,N2" [D]
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
D = int(sys.argv[2]) if len(sys.argv) > 2 else 24
eps, thg = 0.31382, 76.618
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)


def bild_koeff(x, N, M):
    """Koeffizienten 0..M von J_∞ x = x∘f + m x (x Liste arb, Grad N)."""
    f = arb_poly([arb(eps), arb(c0), arb(eps)])
    mp_ = arb_poly([arb(c0)]) + arb(2 * eps) * f
    # Horner: x∘f
    acc = arb_poly([0])
    for c in reversed(x):
        acc = (acc * f + arb_poly([c])).truncate(M + 1)
    tot = acc + (mp_ * arb_poly(x)).truncate(M + 1)
    return [tot[j] for j in range(M + 1)]


import m_kern as K

for N in NS:
    t0 = time.time()
    d = np.load(f"m07_N{N}_e{eps}_t{thg}.npz")
    S, U, Vt = d["S"], d["U"], d["Vt"]
    ctx.prec = int(1.3 * N) + 300
    J = K.jmat([eps, c0, eps], N, reell=True)
    nmod = int(np.sum(S[:8] > 1.6))
    B = arb_mat(N + 1, nmod)
    for i in range(nmod):
        for k in range(N + 1):
            B[k, i] = arb(float(Vt[i][k]))
    X = J.solve(B, algorithm="approx")
    print(f"N={N}: {nmod} Eckmoden, σ = {S[:nmod]} ({time.time()-t0:.0f}s)", flush=True)
    for i in range(nmod):
        x = [arb(X[k, i].mid()) for k in range(N + 1)]
        xn = np.sqrt(sum(float(c.mid()) ** 2 for c in x))
        img = bild_koeff(x, N, N + D)
        low = np.array([float(c.mid()) for c in img[: N + 1]])
        hi = np.array([float(c.mid()) for c in img[N + 1:]])
        cum = np.sqrt(np.cumsum(hi ** 2))
        print(f"  Mode {i+1}: ‖x‖ = {xn:.4g} (σ = {S[i]:.4g}), ‖T_N J x‖ = {np.linalg.norm(low):.3e}; "
              f"Überlauf ‖[Jx]_(N+1..N+j)‖/‖x‖ für j = 1,2,3,4,6,8,12,{D}: "
              + ", ".join(f"{cum[j-1]/xn:.4f}" for j in (1, 2, 3, 4, 6, 8, 12, D)), flush=True)
        print("     |[Jx]_{N+j}|/‖x‖, j = 1…12: " + " ".join(f"{abs(h)/xn:.2e}" for h in hi[:12]), flush=True)
