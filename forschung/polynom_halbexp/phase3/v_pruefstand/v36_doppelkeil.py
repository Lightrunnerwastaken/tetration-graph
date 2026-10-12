"""Strang V — Idee VE2 (Doppelkeil): Linse über Γin (Ks Profil, Höhe S) und Einbuchtung über Γout mit f(γ_out) ⊂ D̄ (Profil −A(1 − e^{−cτ/A}),
A = 0.7, c = 1.0: max|F| auf dem eingebuchteten Bogen < 1, geprüft in Gleitkomma). Gemessen ‖Ĵ_N‖ (Fehlskalierung?) und ‖Ĵ_N⁻¹‖.
Aufruf: python v36_doppelkeil.py MODELL N1,N2,... S"""
import os
import sys
import numpy as np
os.environ.setdefault("V12_EIN", "0.7,1.0")
M_, Ns_, S_ = sys.argv[1], sys.argv[2], sys.argv[3]
sys.argv = ["v12", M_, Ns_, S_, "12", "1.3", "doppel"]
import v12_konform_S as V12  # noqa: E402
sys.stdout.reconfigure(encoding="utf-8")
F, it, d = V12.theodorsen()
res, kmin, kmax, mono, distL = V12.pruefe_abbildung(F)
print(f"[{M_}, Doppelkeil, S = {S_}, Einbuchtung {os.environ['V12_EIN']}] Theodorsen {it} It.; Randfehler {res:.1e}, min|k|/ρ = {kmin:.4f}, "
      f"max|k|/ρ = {kmax:.4f}, Δarg_min = {mono:.1e}", flush=True)
for N in [int(x) for x in Ns_.split(",")]:
    bits = int(1.2 * N) + 400
    J = V12.J_hat(N, F, bits)
    Jf = np.array([[float(J[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
    sv, r = V12.norm_inverse(J)
    print(f"N = {N}: ‖Ĵ_N‖ ≈ {np.linalg.norm(Jf, 2):.3e}, ‖Ĵ_N⁻¹‖ = {sv[0]:.6f} (weitere {sv[1]:.4f}, {sv[2]:.4f}), Residuum {r:.1e}", flush=True)
