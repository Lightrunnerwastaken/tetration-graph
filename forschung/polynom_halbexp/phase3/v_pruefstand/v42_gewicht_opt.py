"""Strang V — VE4, Optimierung des Diagonalgewichts auf gespeicherten UL-Faktoren (v38 …npz): min_x log‖U e^{−X}‖ + log‖e^{X} L‖ (konvex in x),
Subgradient g_j = u_B,j² − v_A,j² mit führenden Singulärvektoren aus warm gestarteter Potenz-Iteration (nur numpy).
Ausgabe: min-Produkt, Faktoren, Schur-Summen beim Optimum, Gestalt von log W_k (Steigungen, Krümmung), Datei …_opt.json.
Aufruf: python v42_gewicht_opt.py DATEI.npz [Iterationen]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
f = Path(sys.argv[1])
IT = int(sys.argv[2]) if len(sys.argv) > 2 else 400
D = np.load(f)
U, L = D["U"], D["L"]
n = U.shape[0]
X = U @ L
nX = np.linalg.norm(X, 2)


def top(M, v0, it=60):
    v = v0 / np.linalg.norm(v0)
    for _ in range(it):
        u = M @ v
        s = np.linalg.norm(u)
        u /= s
        v = M.T @ u
        v /= np.linalg.norm(v)
    u = M @ v
    s = np.linalg.norm(u)
    return s, u / s, v


x = np.zeros(n)
va = np.ones(n)
vb = np.ones(n)
best = None
eta = 0.6
hist = []
for it in range(IT):
    w = np.exp(x)
    A = U / w[None, :]
    B = L * w[:, None]
    a, ua, va = top(A, va, 30 if it else 200)
    b, ub, vb = top(B, vb, 30 if it else 200)
    if best is None or a * b < best[0]:
        best = (a * b, a, b, x.copy())
    hist.append(a * b)
    g = ub ** 2 - va ** 2
    x = x - eta * g
    if it % 80 == 79:
        eta *= 0.7
prod, a, b, x = best
w = np.exp(x)
A_ = np.abs(U / w[None, :])
B_ = np.abs(L * w[:, None])
schurA = math.sqrt(A_.sum(axis=0).max() * A_.sum(axis=1).max())
schurB = math.sqrt(B_.sum(axis=0).max() * B_.sum(axis=1).max())
k = np.arange(n)
N = n - 1
# Gesamtgewicht relativ zur Balancierung in v38 (U, L schon mit s_k skaliert): log W_k = log s_k + x_k; hier nur x_k (Feinanteil) berichten
sel = lambda lo, hi: (k >= lo * N) & (k <= hi * N)
fein_spanne = x.max() - x.min()
print(f"{f.name}: N = {N}, ‖Ĵ⁻¹‖ = {nX:.5f}; min_W ‖U W⁻¹‖·‖W L‖ = {prod:.4f} ({a:.3f} · {b:.3f}); Schur-Summen {schurA:.3f} / {schurB:.3f}; "
      f"Feinanteil log w: Spanne {fein_spanne:.3f}, Steigung (0.1–0.9) {np.polyfit(k[sel(0.1, 0.9)], x[sel(0.1, 0.9)], 1)[0]:.5f}; "
      f"Verlauf {hist[0]:.3f} → {hist[len(hist) // 2]:.3f} → {hist[-1]:.3f}", flush=True)
json.dump(dict(N=N, prod=prod, a=a, b=b, schur=[schurA, schurB], x=[float(t) for t in x], normX=nX),
          open(f.with_name(f.stem + "_opt.json"), "w"))
