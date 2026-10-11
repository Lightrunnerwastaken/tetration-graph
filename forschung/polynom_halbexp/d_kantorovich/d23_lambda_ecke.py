"""D23 — Hängt die Eckmode am lokalen Exponenten bei L?  Strang Bs λ-Deformation des linearen Operators:
J_λ δ = δ∘p + λ (p'∘p) δ  (p = Π_N f_K).  Lokal bei L verschiebt das den Konkurrenz-Exponenten auf
s_c(λ) = 1 − γ + 2 log λ / L  (γ = 2πi/L; Strang B 20:27).  Gemessen: σ1, σ2 der |L|-umskalierten Inversen (letzte B
Spalten) für λ ∈ {0.5, 0.75, 1, 1.5, 2} und mehrere N.
"""
import json
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
LAMS = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.5, 0.75, 1.0, 1.5, 2.0]
B = 100
PI = json.load(open(K.GEMEINSAM + "/lin_PiN_fK_bis_1000.json"))
Lr = 1.3745570107437075
Lc = complex(0.31813150520476413, 1.3372357014306895)
gam = 2j * np.pi / Lc
t0 = time.time()
for N in NS:
    n = N + 1
    ctx.prec = int(0.85 * N) + 300
    a = [arb(s) for s in PI[str(N)]["b"][:n]]
    P = K.potenzen(a, N)
    fpf = arb_poly([0])
    for mm in range(1, n):
        fpf += mm * a[mm] * P[mm - 1]
    for lam in LAMS:
        lam_a = arb(lam)
        J = arb_mat(n, n)
        for k in range(n):
            for j in range(n):
                J[k, j] = P[j][k] + (lam_a * fpf[k - j] if j <= k else 0)
        E = arb_mat(n, B)
        for c in range(B):
            E[N - B + 1 + c, c] = 1
        X = J.solve(E, algorithm="approx")
        kk = np.arange(n)
        cols = np.arange(N - B + 1, n)
        G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * (Lr ** kk)[:, None] / (Lr ** cols)[None, :]
        S = np.linalg.svd(G, compute_uv=False)
        sc = 1 - gam + 2 * np.log(lam) / Lc
        print(f"N={N} λ={lam}: s_c = {sc.real:+.3f}{sc.imag:+.3f}i | σ1={S[0]:.4f} σ2={S[1]:.4f} σ3={S[2]:.4f}  ({time.time()-t0:.0f}s)",
              flush=True)
