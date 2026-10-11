"""D08 — Die Eckschicht von J_N^{-1}: Spalten j ≥ N − B (B = 30). Norm (ℓ1, ℓ2) in Gewichten |L|^k (1+k)^s gegen N,
Profil der stärksten Eckspalte, und Vergleich mit dem Rest (Spalten j < N − B).
"""
import sys

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [60, 80, 100, 120, 140, 160, 180, 194]
B = 30
p = K.loesungen(max(NS), bits=1200)
L = 1.3745570107437075
for N in NS:
    ctx.prec = 500 + 3 * N
    a = [arb(x.mid()) for x in p[N]]
    _, J = K.FJ(a, N)
    n = N + 1
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]
    zeile = []
    for s in (0, 1, 2, 3):
        w = (1.0 + k) ** s
        A = G * w[:, None] / w[None, :]
        ecke, rest = A[:, N - B:], A[:, :N - B]
        zeile.append(f"s={s}: Ecke l1={np.abs(ecke).sum(axis=0).max():.3f} l2={np.linalg.norm(ecke, 2):.3f} | "
                     f"Rest l1={np.abs(rest).sum(axis=0).max():.3f} l2={np.linalg.norm(rest, 2):.3f}")
    j = N - B + int(np.argmax(np.abs(G[:, N - B:]).sum(axis=0)))
    prof = np.abs(G[:, j])
    teile = [prof[int(t0 * j):int(t1 * j)].sum() for (t0, t1) in ((0, .25), (.25, .5), (.5, .75), (.75, .9), (.9, 1.0))]
    print(f"N={N}: stärkste Eckspalte j={j}; ℓ1-Masse (s=0) in k/j ∈ [0,.25),[.25,.5),[.5,.75),[.75,.9),[.9,1): "
          + ", ".join(f"{x:.3f}" for x in teile) + f"; k ≥ j: {prof[j:].sum():.3f}")
    for z in zeile:
        print("    " + z, flush=True)
