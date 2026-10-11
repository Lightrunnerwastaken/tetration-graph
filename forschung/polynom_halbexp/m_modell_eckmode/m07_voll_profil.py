"""M07 — Volle Inverse (arb) im normierten radialen Quadrat-Modell, SVD, Randprofil der adjungierten Eckmoden V_i nahe L
mit dem VOLLEN Vektor (m05 benutzte nur die letzten B Spalten — nahe L zu ungenau).
Aufruf: python m07_voll_profil.py eps theta_grad "N1,N2,..."
Speichert m07_N{N}.npz (S[:8], U[:, :6], Vt[:6]) und druckt σ_i |V_i(θ+τ)|.
"""
import sys
import time

import numpy as np
from flint import arb_mat, ctx

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
eps = float(sys.argv[1])
thg = float(sys.argv[2])
NS = [int(x) for x in sys.argv[3].split(",")]
th = np.deg2rad(thg)
coef = [eps, 1 - 2 * eps * np.cos(th), eps]
taus = np.array([-4, -2, -1, 0, 0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 30, 45, 60, 90])
for N in NS:
    t0 = time.time()
    ctx.prec = int(0.9 * N) + 300
    J = K.jmat(coef, N, reell=True)
    n = N + 1
    E = arb_mat.eye(n) if hasattr(arb_mat, "eye") else None
    if E is None:
        E = arb_mat(n, n)
        for i in range(n):
            E[i, i] = 1
    X = J.solve(E, algorithm="approx")
    G = K.mat_np(X, True)
    t1 = time.time()
    U, S, Vt = np.linalg.svd(G)
    np.savez(f"m07_N{N}_e{eps}_t{thg}.npz", S=S[:12], U=U[:, :8], Vt=Vt[:8])
    k = np.arange(n)
    print(f"N={N} (Inverse {t1-t0:.0f}s): σ = {', '.join(f'{s:.5g}' for s in S[:8])}", flush=True)
    for i in range(4):
        if S[i] < 1.6:
            break
        vals = np.array([abs(np.sum(Vt[i] * np.exp(1j * k * (th + np.deg2rad(t))))) for t in taus])
        print(f"  Paar-Mode {i+1}: log(σ|V(θ+τ)|):", flush=True)
        print("   τ°  " + " ".join(f"{t:7.1f}" for t in taus))
        print("   log " + " ".join(f"{np.log(S[i]*v):7.3f}" for v in vals), flush=True)
