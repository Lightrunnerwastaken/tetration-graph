"""Strang V — Idee VE4: faktorisierter Schur-Test nach Ks Satz K5(5). Ĵ_N = L̂_N U_N (LU ohne Pivot, L̂ = C_kL unten mit Einsen... hier L̂ unten, U oben),
Ĵ_N⁻¹ = U_N⁻¹ L̂_N⁻¹ = (U⁻¹W⁻¹)(WL̂⁻¹) für jede Diagonalmatrix W. Gemessen: min über W (iterativ ausgeglichen) von ‖U⁻¹W⁻¹‖·‖WL̂⁻¹‖ gegen
‖Ĵ_N⁻¹‖, und die Schur-Summen (max Zeilen-/Spaltensumme der Beträge) beider Faktoren. Kubisch, Linse wie v12 (Ks Profil oder Keil).
Aufruf: python v37_lu_schur.py N1,N2,... S [tau1] [p]
"""
import math
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

ARGS = sys.argv[1:]
Ns = [int(x) for x in ARGS[0].split(",")]
S = ARGS[1]
TAU1 = ARGS[2] if len(ARGS) > 2 else "12"
PEXP = ARGS[3] if len(ARGS) > 3 else "1.3"
sys.argv = ["v12", "kubisch", "100", S, TAU1, PEXP]
import v12_konform_S as V12  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def lu_ohne_pivot(J, n):
    """Doolittle in arb: J = L U, L unten mit Einsen."""
    A = [[arb(J[i, j].mid()) for j in range(n)] for i in range(n)]      # Mittelpunkt-Arithmetik (Kugelradien würden in der Elimination explodieren)
    for k in range(n):
        piv = A[k][k]
        for i in range(k + 1, n):
            f = arb((A[i][k] / piv).mid())
            A[i][k] = f
            rowk = A[k]
            rowi = A[i]
            for j in range(k + 1, n):
                rowi[j] = arb((rowi[j] - f * rowk[j]).mid())
    L = arb_mat(n, n)
    U = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            if i > j:
                L[i, j] = A[i][j]
            elif i == j:
                L[i, j] = 1
                U[i, j] = A[i][j]
            else:
                U[i, j] = A[i][j]
    return L, U


def inv(M, n):
    Id = arb_mat(n, n)
    for i in range(n):
        Id[i, i] = 1
    return M.solve(Id, algorithm="approx")


def zu_np(M, n):
    return np.array([[float(M[i, j].mid()) for j in range(n)] for i in range(n)])


def main():
    F, it, d = V12.theodorsen()
    print(f"[kubisch, Linse S = {S}, τ₁ = {TAU1}°, p = {PEXP}]", flush=True)
    for N in Ns:
        t0 = time.time()
        n = N + 1
        bits = int(3 * N) + 800              # LU ohne Pivot verliert viel Präzision
        J = V12.J_hat(N, F, bits)
        ctx.prec = bits
        L, U = lu_ohne_pivot(J, n)
        Jia = inv(J, n)
        Jm = arb_mat([[arb(Jia[i, j].mid()) for j in range(n)] for i in range(n)])
        Ui = zu_np(Jm * L, n)           # U⁻¹ = Ĵ⁻¹L
        Li = zu_np(U * Jm, n)           # L⁻¹ = UĴ⁻¹
        Ji = zu_np(Jia, n)
        nJ = np.linalg.norm(Ji, 2)
        print(f"   Diagnose: max|U⁻¹| = {np.nanmax(np.abs(Ui)):.2e}, max|L⁻¹| = {np.nanmax(np.abs(Li)):.2e}, endlich: {np.isfinite(Ui).all()} {np.isfinite(Li).all()}; "
              f"‖Ĵ⁻¹ − U⁻¹L⁻¹‖ = {np.linalg.norm(Ji - Ui @ Li):.1e}", flush=True)
        # Ausgleich: w_k so, dass Spalte k von U⁻¹W⁻¹ und Zeile k von WL⁻¹ gleiche Norm haben (Iteration auf ‖·‖₂ der Faktoren)
        cu = np.linalg.norm(Ui, axis=0)
        rl = np.linalg.norm(Li, axis=1)
        w = np.sqrt(cu / rl)
        best = None
        for _ in range(30):
            A_ = Ui / w[None, :]
            B_ = Li * w[:, None]
            a, b = np.linalg.norm(A_, 2), np.linalg.norm(B_, 2)
            if best is None or a * b < best[0]:
                best = (a * b, a, b, w.copy())
            # Feinausgleich über die führenden Singulärvektoren (gedämpft, begrenzt)
            try:
                ua = np.abs(np.linalg.svd(A_)[2][0])
                vb = np.abs(np.linalg.svd(B_)[0][:, 0])
            except np.linalg.LinAlgError:
                break
            schritt = np.clip(np.log((ua + 1e-30) / (vb + 1e-30)), -1, 1)
            w = w * np.exp(0.15 * schritt)
        prod, a, b, w = best
        A_ = np.abs(Ui / w[None, :])
        B_ = np.abs(Li * w[:, None])
        schurA = math.sqrt(A_.sum(axis=0).max() * A_.sum(axis=1).max())
        schurB = math.sqrt(B_.sum(axis=0).max() * B_.sum(axis=1).max())
        print(f"N = {N}: ‖Ĵ⁻¹‖ = {nJ:.4f}; min_W ‖U⁻¹W⁻¹‖·‖WL̂⁻¹‖ ≈ {prod:.4f} ({a:.3f} · {b:.3f}); Schur-Summen {schurA:.3f} / {schurB:.3f}; "
              f"ohne Gewicht ‖U⁻¹‖·‖L̂⁻¹‖ = {np.linalg.norm(Ui, 2) * np.linalg.norm(Li, 2):.2e}  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
