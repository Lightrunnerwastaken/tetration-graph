"""k10 — Affines Eckmodell: f_aff(x) = L + μ(x − L) (nur Fixpunkt L, Multiplikator μ = √L, m = μ konstant).

In der x-Basis ist J_aff = C_φ + μI OBERE Dreiecksmatrix (Komposition mit einer affinen Abbildung senkt den Grad nicht), C_k ist UNTERE
Dreiecksmatrix. Daher exakt (bewiesen, Zeile): Ĵ_N⁻¹ = (J̃⁻¹)_N (C_k⁻¹)_N und die Widom-Differenz
    D_N = Ĵ_N⁻¹ − P_N 𝒥⁻¹ P_N = − P_N J̃⁻¹ Q_N C_k⁻¹ P_N      (Hankel-artiges Produkt zweier Eckblöcke).
Test: trifft D_N(affin) die Eck-Korrektur des kubischen/Kneser-Modells (k09: 0.5123/0.5057, 0.2345/0.2232, 0.0917/0.0875)?
Gleiche Kurve γ (Kreis + Linse über arg ∈ (arg L, 2π − arg L), S = 0.0166) wie k04/k09.
Aufruf: python k10_affin.py 600 200
"""
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import Lacb, L_RE, L_IM  # noqa: E402
import k04_konform as K4  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def J_affin(N, bits):
    """J̃_aff (ρ-skaliert) als acb_mat: [x^k](μx + b)^j ρ^{k−j} + μ δ_kj, b = L(1 − μ)."""
    ctx.prec = bits
    L = Lacb()
    mu = L.sqrt()
    b = L * (1 - mu)
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    J = acb_mat(N + 1, N + 1)
    # Spalte j: Σ_k C(j,k) μ^k b^{j−k} ρ^{k−j}  =  Σ_k C(j,k) (μ)^k (b/ρ)^{j−k}
    bq = b / rho
    for j in range(N + 1):
        c = acb(1)  # C(j,k) μ^k (b/ρ)^{j−k} für k = 0 : (b/ρ)^j
        terms = [acb(0)] * (j + 1)
        # rekursiv über k: t_k = C(j,k) μ^k bq^{j−k}
        t = bq ** j
        for k in range(j + 1):
            terms[k] = t
            if k < j:
                t = t * (j - k) / (k + 1) * mu / bq
        for k in range(j + 1):
            J[k, j] = terms[k]
        J[j, j] = J[j, j] + mu
    return J


def to_acb(A):
    n, m = A.nrows(), A.ncols()
    B = acb_mat(n, m)
    for i in range(n):
        for j in range(m):
            B[i, j] = acb(A[i, j])
    return B


def inv(A):
    n = A.nrows()
    I = acb_mat(n, n)
    for i in range(n):
        I[i, i] = acb(1)
    return A.solve(I, algorithm="approx")


def block(A, n):
    B = acb_mat(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = A[i, j]
    return B


def cmat(X):
    n, m = X.nrows(), X.ncols()
    out = np.zeros((n, m), dtype=complex)
    for i in range(n):
        for j in range(m):
            z = X[i, j]
            out[i, j] = complex(float(z.real.mid()), float(z.imag.mid()))
    return out


def main():
    Nr = int(sys.argv[1])
    N = int(sys.argv[2])
    bits = int(1.2 * Nr) + 300
    t0 = time.time()
    F, d = K4.theodorsen()
    J = J_affin(Nr, bits)
    ctx.prec = bits
    C = to_acb(K4.C_matrix(F, Nr))
    Jh = C * J
    Xr = cmat(inv(Jh))
    Yr = cmat(inv(J))
    sJ = np.linalg.svd(Yr, compute_uv=False)
    n = N + 1
    XN = cmat(inv(block(Jh, n)))
    D = XN - Xr[:n, :n]
    U, sD, Vh = np.linalg.svd(D)
    # exakte Formel: D_N = − (J̃⁻¹)_{≤N, >N} (C⁻¹)_{>N, ≤N}
    Ci = cmat(inv(C))
    Dform = -Yr[:n, n:] @ Ci[n:, :n]
    print("affin  σ₁(J̃_%d⁻¹) = %.5g | ‖Ĵ_%d⁻¹‖ = %.6f, ‖Ĵ_%d⁻¹‖ = %.6f | σ(D_%d) = %s | ‖D − D_formel‖ = %.2e  (%.0f s)"
          % (Nr, sJ[0], Nr, np.linalg.norm(Xr, 2), N, np.linalg.norm(XN, 2), N,
             " ".join("%.5f" % s for s in sD[:6]), np.linalg.norm(D - Dform, 2), time.time() - t0), flush=True)
    # Lage der führenden Vektoren: Indexschwerpunkt (Zeilen = x-Grad, Spalten = w-Grad)
    for i in range(3):
        u = np.abs(U[:, i]) ** 2
        v = np.abs(Vh[i, :]) ** 2
        ks = np.arange(n) / N
        print("  Mode %d: σ = %.5f, x-Seite Schwerpunkt ξ = %.3f (Max bei %.3f), w-Seite Schwerpunkt ξ = %.3f (Max bei %.3f)"
              % (i + 1, sD[i], (u * ks).sum(), ks[np.argmax(u)], (v * ks).sum(), ks[np.argmax(v)]), flush=True)


if __name__ == "__main__":
    main()
