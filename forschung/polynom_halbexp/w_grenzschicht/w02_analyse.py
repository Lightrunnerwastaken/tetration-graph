"""W02 — Gestalt der Eckmoden (kubisch): Indexprofil, Winkelprofil, Phasenraum (Gauss-Fenster).
Aufruf: python w02_analyse.py N
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
d = np.load(W.HIER + f"/w01_N{N}.npz")
S, U, Vt = d["S"], d["U"], d["Vt"]
n = N + 1
k = np.arange(n)
argL = np.degrees(np.angle(W.L))
print(f"N={N}: σ1..4 = {S[:4]}")
for i in range(2):
    v = Vt[i]
    u = U[:, i]
    print(f"\n=== Mode {i+1}: σ = {S[i]:.4f}")
    # Indexprofil
    lv = np.log10(np.abs(v) + 1e-300)
    lu = np.log10(np.abs(u) + 1e-300)
    imax = np.argmax(np.abs(v))
    print(f" v: max |v_k| bei k/N = {imax/N:.4f}; Quantile Masse 5/50/95%: "
          + " ".join(f"{k[np.searchsorted(np.cumsum(v**2), q)]/N:.3f}" for q in (0.05, 0.5, 0.95)))
    print(" v: log10|v_k| an k/N =", " ".join(f"{x:.2f}:{lv[int(x*N)]:.1f}" for x in (0.5, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 1.0)))
    # Glättung: Hüllkurve als max über Fenster 5
    print(" u: Quantile Masse 5/50/95%: " + " ".join(f"{k[np.searchsorted(np.cumsum(u**2), q)]/N:.3f}" for q in (0.05, 0.5, 0.95)))
    print(" u: log10|u_k| an k/N =", " ".join(f"{x:.2f}:{lu[int(x*N)]:.1f}" for x in (0.1, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1.0)))
    # Winkelprofile
    phis = np.array([0, 30, 60, 70, 74, 76, argL, 77.5, 79, 82, 86, 90, 100, 120, 150, 180])
    Vc = W.V_auswerten(v, W.RHO * np.exp(1j * np.deg2rad(phis)))
    Uc = W.U_auswerten(u, W.RHO * np.exp(1j * np.deg2rad(phis)))
    print(" log10|V(ρe^{iφ})|:", " ".join(f"{p:.1f}:{np.log10(abs(x)+1e-300):.2f}" for p, x in zip(phis, Vc)))
    print(" log10|U(ρe^{iφ})|:", " ".join(f"{p:.1f}:{np.log10(abs(x)+1e-300):.2f}" for p, x in zip(phis, Uc)))
    # Phasenraum: Gauss-Fenster in k, Breite w; für Fensterzentren k0 das Spektrum über φ
    wbr = max(8, int(0.03 * N))
    print(f" Phasenraum (Fenster σ_k = {wbr}): Zentrum k0/N → Winkel des Maximums, log10 Max, Breite")
    phig = np.deg2rad(np.linspace(0, 360, 3601))
    E = np.exp(-1j * np.outer(k, phig))
    for x0 in (0.6, 0.7, 0.8, 0.85, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 0.995):
        g = np.exp(-0.5 * ((k - x0 * N) / wbr) ** 2)
        spec = np.abs((v * g) @ E)
        j = np.argmax(spec[:1801])   # obere Hälfte
        hm = spec[:1801] > 0.5 * spec[j]
        lo, hi = np.degrees(phig[:1801][hm][[0, -1]])
        print(f"   k0/N={x0:.3f}: φ_max = {np.degrees(phig[j]):6.1f}°, log10 = {np.log10(spec[j]):6.2f}, Halbwertsbereich {lo:.0f}…{hi:.0f}°")
