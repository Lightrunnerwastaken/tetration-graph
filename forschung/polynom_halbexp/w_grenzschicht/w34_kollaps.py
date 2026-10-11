"""W34 — Skalenkollaps des Eckprofils nahe L: f_N(τ) := ln(Λ_N(τ)/Λ_N(0)), Λ_N(φ) = (Σ_{σ_i>1.40} σ_i²|V_i(ρe^{iφ})|²)^{1/2} (kubisch, w01).
Hypothese: f_N(τ) = F(τ·N^β) (+ κ-Terme) mit β = 1/p (p = Profil-Exponent): p = 3/2 ⇔ β = 2/3 (Airy-artige Innenschicht).
Test: Streuung von f_N über N bei festem ξ = τN^β für β ∈ [0.3, 1.2]; zusätzlich die Variante mit Subtraktion des Eckterms ln N.
Aufruf: python w34_kollaps.py 400,500,...,1000
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
argL = np.angle(W.L)
taus = np.geomspace(0.01, 8.0, 120)            # Grad, Γin-Seite
tneg = -np.geomspace(0.01, 8.0, 120)           # Γout-Seite
F = {}
Fneg = {}
for N in NS:
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, Vt = d["S"], d["Vt"]
    r = max(1, sum(1 for s in S[:Vt.shape[0]] if s > 1.40))
    def lam(tt):
        z = W.RHO * np.exp(1j * (argL + np.deg2rad(tt)))
        return np.sqrt(sum(S[i] ** 2 * np.abs(W.V_auswerten(Vt[i], z)) ** 2 for i in range(r)))
    l0 = lam(np.array([0.0]))[0]
    F[N] = np.log(lam(taus) / l0)
    Fneg[N] = np.log(lam(tneg) / l0)
    print(f"N={N}: r={r}, Λ(0) = {l0:.4g}, f(0.1°)={np.interp(0.1, taus, F[N]):+.3f}, f(1°)={np.interp(1, taus, F[N]):+.3f}, "
          f"f(4°)={np.interp(4, taus, F[N]):+.3f}, f(−1°)={np.interp(1, -tneg, Fneg[N]):+.3f}", flush=True)
print("\nStreuung (rms über N) von f_N(ξ/N^β) auf ξ-Gitter, β-Scan:")
best = None
for beta in np.arange(0.3, 1.21, 0.05):
    xis = np.geomspace(0.3, 200, 40)
    vals = []
    for xi in xis:
        row = []
        for N in NS:
            t = xi / N ** beta
            if taus[0] <= t <= taus[-1]:
                row.append(np.interp(np.log(t), np.log(taus), F[N]))
        if len(row) == len(NS):
            vals.append(np.std(row))
    sc = np.sqrt(np.mean(np.square(vals))) if vals else np.nan
    if best is None or sc < best[1]:
        best = (beta, sc)
    print(f"  β = {beta:.2f}: rms {sc:.4f} ({len(vals)} ξ-Punkte)")
print(f"bestes β = {best[0]:.2f}")
beta = best[0]
print(f"\nKollabierte Kurve F(ξ), ξ = τ·N^{beta:.2f} (τ in Grad):")
for xi in (0.5, 1, 2, 5, 10, 20, 50, 100):
    print(f"  ξ={xi}: " + " ".join(f"{np.interp(np.log(xi/N**beta), np.log(taus), F[N]):+.3f}" for N in NS))
