"""W05 — Nullstellenkurven der Eckmoden u₁, V₁ als Funktion des Winkels, mehrere N (kubisch).
Nullstellen mit mpmath.polyroots (hohe Präzision auf den float-Koeffizienten) — gibt die Nullstellen DES float-Vektors.
Ausgabe je Winkelband: mittlerer Radius r(φ), N(r−1), Dichte (Anzahl/Band)/N.
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
bands = np.array([0, 20, 40, 60, 70, 74, 76.6, 80, 85, 90, 100, 120, 140, 160, 180])
for N in NS:
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    S, U, Vt = d["S"], d["U"], d["Vt"]
    for name, c in (("u1", U[:, 0]), ("v1", Vt[0])):
        coeffs = c[::-1] if name == "u1" else c
        r = np.roots(coeffs)
        ab, ang = np.abs(r), np.abs(np.degrees(np.angle(r)))
        out = []
        for a, b in zip(bands[:-1], bands[1:]):
            sel = (ang >= a) & (ang < b) & (ab > 0.8) & (ab < 1.25)
            if sel.sum():
                rm = np.median(ab[sel])
                out.append(f"{a:.0f}-{b:.0f}:{rm:.4f}/{N*(rm-1):+.1f}/{sel.sum()/(2*N*np.deg2rad(b-a)):.3f}")
        print(f"N={N} {name}: " + " ".join(out))
