"""D07 — |J_N(p_N)^{-1}| in gewichteten ℓ1- und ℓ2-Normen (gleiche Gewichte w_k = |L|^k (1+k)^s in Urbild und Bild).

ℓ1: max_j Σ_k |G̃_kj| w̃_k / w̃_j ;  ℓ2: grösster Singulärwert von D G̃ D^{-1};  G̃ = |L|-umskalierte Inverse.
Zusätzlich die Diagonal-nahe Summe Σ_{|d|≤j/10} |G̃_{j+d,j}| (untere Schranke für jede ℓ1-Norm mit langsam veränderlichem
Gewicht), die laut D06 wie 0.34·ln j wächst (Sprung des Symbols bei L, L̄).
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [20, 40, 60, 80, 100, 120, 140, 160, 180, 194]
SS = [0, 0.5, 1, 2, 3, 4, 5, 5.4]
p = K.loesungen(max(NS), bits=1200)
L = 1.3745570107437075
aus = {}
t0 = time.time()
for N in NS:
    ctx.prec = 500 + 3 * N
    a = [arb(x.mid()) for x in p[N]]
    _, J = K.FJ(a, N)
    n = N + 1
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]
    z = {}
    for s in SS:
        w = (1.0 + k) ** s
        A = G * w[:, None] / w[None, :]
        col = np.abs(A).sum(axis=0)
        z[str(s)] = {"l1": float(col.max()), "l1_j": int(col.argmax()), "l2": float(np.linalg.norm(A, 2))}
    nah = []
    for j in range(10, n):
        D = max(1, j // 10)
        nah.append(sum(abs(G[j + d, j]) for d in range(-D, D + 1) if 0 <= j + d <= N))
    z["nah_max"] = float(max(nah))
    aus[N] = z
    print(f"N={N:3d} ({time.time()-t0:.0f}s)  nah={z['nah_max']:.3f}  " +
          "  ".join(f"s={s}: l1={z[str(s)]['l1']:.3f}@{z[str(s)]['l1_j']} l2={z[str(s)]['l2']:.3f}" for s in SS), flush=True)
json.dump({str(k): v for k, v in aus.items()}, open(K.HIER + "/d07_normen.json", "w"), indent=1)
