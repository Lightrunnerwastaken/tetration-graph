"""M04 — Phasenraum-Bild (gefensterte Fourier-Transformation) der Eckvektoren u1, v1 aus m02 (npz).
|Σ_k w((k−k0)/W) c_k e^{ikφ}| für k0 = κN, Gauss-Fenster der Breite W, auf einem φ-Raster.
Aufruf: python m04_phasenraum.py datei.npz [W]
"""
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
d = np.load(sys.argv[1])
W = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
U, Vt = d["U"], d["Vt"]
N = U.shape[0] - 1
k = np.arange(N + 1)
phis = np.array([0, 40, 60, 70, 74, 76, 76.6, 77, 78, 80, 85, 90, 100, 120, 150, 180])
for name, vec in (("u1", U[:, 0]), ("v1", Vt[0, :])):
    print(f"\n{name} (N={N}, Fensterbreite {W}): log10|Fenster-FT|, Zeilen κ = k0/N, Spalten φ in Grad")
    print("   κ    " + " ".join(f"{p:6.1f}" for p in phis))
    for kap in np.arange(0.5, 1.0001, 0.025):
        k0 = kap * N
        w = np.exp(-0.5 * ((k - k0) / W) ** 2)
        row = []
        for p in phis:
            z = np.exp(1j * np.deg2rad(p) * k)
            row.append(np.log10(abs(np.sum(w * vec * z)) + 1e-300))
        print(f"  {kap:.3f} " + " ".join(f"{x:6.1f}" for x in row))
