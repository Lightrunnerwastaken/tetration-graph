"""M05 — Randprofil der adjungierten Eckmode im normierten radialen Quadrat-Modell (Bedingung (E_L)).
λ_N(τ) := (1/N) log( σ₁ |V₁(θ+τ)| ),  V₁(φ) = Σ_k v_k e^{ikφ}  (v₁ = rechter Singulärvektor der Inversen, letzte B Koordinaten).
Aufruf: python m05_adjungiert_profil.py eps theta_grad "N1,N2,..." [B]
"""
import json
import sys
import time

import numpy as np

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
eps = float(sys.argv[1])
thg = float(sys.argv[2])
NS = [int(x) for x in sys.argv[3].split(",")]
B = int(sys.argv[4]) if len(sys.argv) > 4 else 150
th = np.deg2rad(thg)
coef = [eps, 1 - 2 * eps * np.cos(th), eps]
taus = np.array([-8, -4, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 20, 40, 60, 80, 103.38])
out = {"eps": eps, "theta": thg, "taus": taus.tolist(), "N": {}}
t0 = time.time()
for N in NS:
    S, U, Vt, cols = K.ecksingulaer(coef, N, rho=1.0, B=B, k=2)
    v = Vt[0, :]
    # V(φ) = Σ_c v_c e^{i k_c φ}, k_c = cols
    vals = []
    for tg in taus:
        p = th + np.deg2rad(tg)
        vals.append(abs(np.sum(v * np.exp(1j * cols * p))))
    vals = np.array(vals)
    lam = np.log(S[0] * vals) / N
    out["N"][N] = {"sigma1": float(S[0]), "V": vals.tolist(), "lam": lam.tolist()}
    print(f"N={N} ({time.time()-t0:.0f}s): σ₁ = {S[0]:.5g}; σ₁|V(θ)| = {S[0]*vals[taus.tolist().index(0)]:.4g}", flush=True)
    print("   τ[°]:  " + " ".join(f"{t:8.2f}" for t in taus))
    print("   λ_N :  " + " ".join(f"{x:8.5f}" for x in lam))
with open(f"m05_profil_eps{eps}_th{thg}.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh)
