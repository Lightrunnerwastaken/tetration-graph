"""W04 — Nullstellen und Potential (1/N)log|·| der Eckmoden u₁ (Polynom) und V₁ (Laurent-Polynom), kubisches Modell.
Aufruf: python w04_nullstellen.py N
"""
import sys

import numpy as np
import mpmath as mp

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
d = np.load(W.HIER + f"/w01_N{N}.npz")
S, U, Vt = d["S"], d["U"], d["Vt"]
u, v = U[:, 0], Vt[0]
# komplexe Kombination der beiden Eckmoden: wähle Kombination, die bei L̄ minimal ist
zL = W.L / W.RHO
for name, c in (("u1", u), ("v1", v)):
    # Polynom in t: u: Σ u_k t^k (t = z/ρ); v: t^N V = Σ v_j t^{N-j}
    coeffs = c[::-1] if name == "u1" else c       # numpy roots: höchster Grad zuerst
    mp.mp.dps = 50
    r = np.roots(coeffs)
    ab = np.abs(r)
    ang = np.degrees(np.angle(r))
    print(f"\n{name}: {len(r)} Nullstellen (t = z/ρ); |t|-Quantile 10/50/90%: {np.quantile(ab, [0.1, 0.5, 0.9])}")
    # Histogramm über Winkel und Radius
    H, xe, ye = np.histogram2d(ang, ab, bins=[np.arange(-180, 181, 20), [0, 0.5, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.5, 2, 5, 1e9]])
    print("   Winkel\\|t|: " + " ".join(f"{y:>5.2f}" for y in ye[:-1]))
    for i in range(len(xe) - 1):
        print(f"   {xe[i]:5.0f}…{xe[i+1]:4.0f}: " + " ".join(f"{int(h):5d}" for h in H[i]))
# Potential (1/N) log|V(z)| auf Radien
print("\n(1/N) log|V1(ρ r e^{iφ})| (Zeilen r, Spalten φ):")
phis = np.array([0, 40, 60, 70, 76.6, 80, 90, 110, 130, 150, 180])
print("  r\\φ   " + " ".join(f"{p:7.1f}" for p in phis))
for r in (0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4):
    z = W.RHO * r * np.exp(1j * np.deg2rad(phis))
    val = np.log(np.abs(W.V_auswerten(v, z))) / N
    print(f"  {r:4.2f}  " + " ".join(f"{x:7.4f}" for x in val))
print("\n(1/N) log|U1(ρ r e^{iφ})|:")
for r in (0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4):
    z = W.RHO * r * np.exp(1j * np.deg2rad(phis))
    val = np.log(np.abs(W.U_auswerten(u, z))) / N
    print(f"  {r:4.2f}  " + " ".join(f"{x:7.4f}" for x in val))
