"""W08 — Asymptotische Trennung des Randprofils: log Λ_N(τ) = a(τ) + κ(τ) log N + λ(τ) N  (Fit über N),
Λ_N(φ) := (Σ_{i≤2} σ_i²|V_i(ρe^{iφ})|²)^{1/2} (kubisch, aus w01). Danach lokale Exponenten von λ(τ) − λ(0).
Aufruf: python w08_profil_fit.py 600,700,800,900,1000
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
argL = np.degrees(np.angle(W.L))
taus = np.array([-4, -2, -1, -0.5, -0.25, 0, 0.05, 0.1, 0.2, 0.25, 0.35, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64, 103.38])
logs = []
for N in NS:
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, Vt = d["S"], d["Vt"]
    z = W.RHO * np.exp(1j * np.deg2rad(argL + taus))
    rr_ = max(1, sum(1 for s in S[:Vt.shape[0]] if s > 1.40)); lam = np.sqrt(sum(S[i] ** 2 * np.abs(W.V_auswerten(Vt[i], z)) ** 2 for i in range(rr_)))
    logs.append(np.log(lam))
logs = np.array(logs)            # [N, τ]
Ns = np.array(NS, dtype=float)
A = np.vstack([np.ones_like(Ns), np.log(Ns), Ns]).T
coef, res, *_ = np.linalg.lstsq(A, logs, rcond=None)
a, kap, lam = coef
pred = A @ coef
rms = np.sqrt(np.mean((pred - logs) ** 2, axis=0))
i0 = list(taus).index(0)
print("τ(°)    κ(τ)     λ(τ)       λ−λ(0)     rms")
for j, t in enumerate(taus):
    print(f"{t:7.2f} {kap[j]:8.3f} {lam[j]:10.6f} {lam[j]-lam[i0]:+10.6f} {rms[j]:.1e}")
print("\nlokale Exponenten d ln(λ−λ0)/d ln τ (τ > 0):")
pos = [j for j, t in enumerate(taus) if t > 0 and lam[j] - lam[i0] > 0]
for j1, j2 in zip(pos[:-1], pos[1:]):
    e = np.log((lam[j2] - lam[i0]) / (lam[j1] - lam[i0])) / np.log(taus[j2] / taus[j1])
    print(f"  {taus[j1]:.2f}→{taus[j2]:.2f}: {e:.3f}")
