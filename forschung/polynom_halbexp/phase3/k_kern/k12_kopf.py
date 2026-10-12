"""k12 — Kopf-Zerlegung (LU) der Abschnittsinversen und Grad-Antwortprofil (A2c Schritt 1).

Exakt (bewiesen, BERICHT §11): J_N⁻¹ = Σ_{n≤N} c_n ℓ_nᵀ mit c_n = J_n⁻¹e_n (Spalte n von U⁻¹) und ℓ_nᵀ = e_nᵀJ_n⁻¹/(J_n⁻¹)_nn (Zeile n von L⁻¹);
Ĵ_N⁻¹ = Σ_{n≤N} c_n (ℓ_nᵀ C_k⁻¹). Der Kopfterm bei n = N ist X_N − (X_{N−1} ⊕ 0) = c_N ℓ_Nᵀ, aus EINER Inversen X_N ablesbar:
  c_N = X_N e_N,  ℓ_Nᵀ = e_NᵀX_N/(X_N)_NN.
Gemessen je N: ‖Kopfterm‖ in x (ℓ²_ρ → ℓ²_ρ) und in w (ℓ²(w) → ℓ²_ρ), ‖c_N‖, ‖ℓ_N‖, ‖ℓ_NᵀC⁻¹‖, Pivot 1/(X_N)_NN,
Grad-Antwortprofil R_N(k) = ‖X_N e_k‖ (Rate (1/N)ln R_N(ξN)) und G_N(L) = ‖X_N ẽ_L‖.
Aufruf: python k12_kopf.py kubisch 300,400,500,600
"""
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, Lc, RHO  # noqa: E402
import k04_konform as K4  # noqa: E402
from k05_widom import inv  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    name = sys.argv[1]
    Ns = [int(s) for s in sys.argv[2].split(",")]
    F, d = K4.theodorsen()
    for N in Ns:
        t0 = time.time()
        bits = int(1.2 * N) + 300
        J = J_tilde(name, N, bits=bits)
        ctx.prec = bits
        X = K4.float_mat(inv(J))
        C = K4.C_matrix(F, N)
        Ci = K4.float_mat(inv(C))
        n = N
        cN = X[:, n]
        rowN = X[n, :]
        piv = X[n, n]
        lN = rowN / piv
        lw = lN @ Ci
        kopf_x = np.linalg.norm(cN) * np.linalg.norm(lN)
        kopf_w = np.linalg.norm(cN) * np.linalg.norm(lw)
        R = np.linalg.norm(X, axis=0)
        eL = (RHO / Lc) ** np.arange(N + 1)
        GL = np.linalg.norm(X @ eL)
        sX = np.linalg.svd(X, compute_uv=False)[0]
        print("N = %4d  ‖J̃⁻¹‖ = %.4g | Kopf: ‖c_N‖ = %.3e, ‖ℓ_N‖_x = %.3e, ‖ℓ_NᵀC⁻¹‖_w = %.3e, |1/(X)_NN| = %.3e; ‖Kopfterm‖ x = %.4g, w = %.4g | G_N(L) = %.4g  (%.0f s)"
              % (N, sX, np.linalg.norm(cN), np.linalg.norm(lN), np.linalg.norm(lw), abs(1 / piv), kopf_x, kopf_w, GL, time.time() - t0), flush=True)
        xs = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
        print("          Grad-Antwort (1/N)ln R_N(ξN):  " + "  ".join("ξ=%.2f: %+.4f" % (x, np.log(R[min(int(round(x * N)), N)]) / N) for x in xs), flush=True)


if __name__ == "__main__":
    main()
