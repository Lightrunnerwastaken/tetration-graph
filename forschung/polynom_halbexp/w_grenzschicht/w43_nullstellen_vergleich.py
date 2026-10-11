"""W43 — Vergleich der Nullstellen von y₁ (w42) mit der Modellvorhersage (w42_modell_kurven.npz aus dem Referenzmodell):
Lage (Abstand zu B bzw. zur Nulllinie der Linse), Anteil innen (Modell 0.4728) und Winkelverteilung (Dichte |∇H|/2π längs der Kurven).
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
M = np.load("w42_modell_kurven.npz")
B, Z, gB, gZ = M["B"], M["Z"], M["gB"], M["gZ"]
B = np.concatenate([np.conj(B[::-1]), B[1:]])
gB = np.concatenate([gB[::-1], gB[1:]])
Z = np.concatenate([Z, np.conj(Z[::-1])[1:]])
gZ = np.concatenate([gZ, gZ[::-1][1:]])


def cdf_modell(curve, g):
    ds = np.abs(np.diff(curve))
    m = np.concatenate([[0], np.cumsum(0.5 * (g[1:] + g[:-1]) * ds)])
    ang = np.unwrap(np.angle(curve))
    return ang, m / m[-1]


aB, FB = cdf_modell(B, gB)
aZ, FZ = cdf_modell(Z, gZ)
R = json.load(open("w42_nullstellen.json"))
for N in sorted(R, key=int):
    x = np.array(R[N]["re"]) + 1j * np.array(R[N]["im"])
    inn = np.abs(x) < 1
    dB = np.array([np.min(np.abs(z - B)) for z in x[inn]])
    dZ = np.array([np.min(np.abs(z - Z)) for z in x[~inn]])
    # Winkel-KS-Abstand
    ai = np.sort(np.angle(x[inn]))
    ks_i = np.max(np.abs(np.interp(ai, aB, FB) - (np.arange(1, len(ai) + 1) - 0.5) / len(ai)))
    ao = np.sort(np.mod(np.angle(x[~inn]), 2 * np.pi))
    ks_o = np.max(np.abs(np.interp(ao, np.mod(aZ, 2 * np.pi), FZ) - (np.arange(1, len(ao) + 1) - 0.5) / len(ao)))
    print(f"N={N}: innen {inn.mean():.4f} (Modell 0.4728); Abstand zu B: Median {np.median(dB):.4f}, max {dB.max():.4f}; "
          f"zur Linsen-Nulllinie: Median {np.median(dZ):.4f}, max {dZ.max():.4f}; Winkel-KS innen {ks_i:.3f}, aussen {ks_o:.3f}")
