"""W03 — Randprofil der adjungierten Eckmoden auf |x| = ρ nahe L (kubisch), direkt aus v_i:
Λ_N(φ) := (Σ_{i≤2} σ_i²|V_i(ρe^{iφ})|²)^{1/2}; Raten (log Λ_{N2} − log Λ_{N1})/(N2 − N1).
Aufruf: python w03_profil.py 300,400,500,600
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
argL = np.degrees(np.angle(W.L))
taus = np.array([-8, -4, -2, -1, -0.5, 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64, 80, 103.38])
res = {}
for N in NS:
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, Vt = d["S"], d["Vt"]
    z = W.RHO * np.exp(1j * np.deg2rad(argL + taus))
    lam = np.sqrt(sum(S[i] ** 2 * np.abs(W.V_auswerten(Vt[i], z)) ** 2 for i in range(2)))
    res[N] = lam
    print(f"N={N} σ1={S[0]:.3f}: " + " ".join(f"{t}:{np.log10(l):.2f}" for t, l in zip(taus, lam)))
print("\nRaten zwischen aufeinanderfolgenden N (je Grad), minus Rate in L:")
for N1, N2 in zip(NS[:-1], NS[1:]):
    r = (np.log(res[N2]) - np.log(res[N1])) / (N2 - N1)
    r0 = r[list(taus).index(0)]
    print(f"{N1}→{N2}: Rate L = {r0:.5f}; λ(τ): " + " ".join(f"{t}:{x - r0:+.5f}" for t, x in zip(taus, r)))
