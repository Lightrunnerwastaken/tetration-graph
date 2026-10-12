"""k05 — Widom-Struktur der konformen Abschnitte: Ĵ_N ist der linke obere Block der unendlichen Matrix Ĵ_∞ = C_k J̃_∞ (geschachtelt).
Vergleich Ĵ_N^{-1} mit dem Block P_N Ĵ_{N'}^{-1} P_N (N' gross ≈ Ersatz für die Inverse 𝒥^{-1} des unendlichen Operators):
  ‖Ĵ_N^{-1} − P_N Ĵ_{N'}^{-1} P_N‖ und wo die Differenz sitzt (Zeilen = Unbekannte x, Spalten = w-Gleichungen, nahe der Kante N?).
Bei Toeplitz-Abschnitten wäre die Differenz eine reine Kantenkorrektur (Widom: T_N(a)^{-1} ≈ P_N T(a^{-1}) P_N + W_N P_N T(ã^{-1}) P_N W_N).
Aufruf: python k05_widom.py kubisch 600 200,300,400
"""
import math
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde  # noqa: E402
import k04_konform as K4  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def inv(A):
    n = A.nrows()
    I = arb_mat(n, n)
    for i in range(n):
        I[i, i] = arb(1)
    return A.solve(I, algorithm="approx")


def block(A, n):
    B = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = A[i, j]
    return B


def main():
    name = sys.argv[1]
    Nr = int(sys.argv[2])
    Ns = [int(x) for x in sys.argv[3].split(",")]
    F, d = K4.theodorsen()
    bits = int(1.2 * Nr) + 300
    t0 = time.time()
    J = J_tilde(name, Nr, bits=bits)
    ctx.prec = bits
    C = K4.C_matrix(F, Nr)
    Jh = C * J
    Xr = K4.float_mat(inv(Jh))
    print("Referenz N' = %d fertig (%.0f s), ‖Ĵ_{N'}⁻¹‖ = %.4f" % (Nr, time.time() - t0, np.linalg.norm(Xr, 2)), flush=True)
    for N in Ns:
        n = N + 1
        XN = K4.float_mat(inv(block(Jh, n)))
        D = XN - Xr[:n, :n]
        nD = np.linalg.norm(D, 2)
        # Lage: Norm der Differenz auf Kanten-Teilblöcken (letzte m Zeilen bzw. Spalten) gegen den Rest
        out = []
        for m in (10, 30, 100):
            kante = np.linalg.norm(D[:, n - m:], 2)           # Spalten nahe der Kante (w-Gleichungen hohen Grades)
            innen = np.linalg.norm(D[:, :n - m], 2)
            out.append("m=%d: Kantenspalten %.2e / Rest %.2e" % (m, kante, innen))
        zeilen = [np.linalg.norm(D[n - m:, :], 2) for m in (10, 30, 100)]
        print("N = %d: ‖Ĵ_N⁻¹‖ = %.4f, ‖Ĵ_N⁻¹ − P_NĴ_{N'}⁻¹P_N‖ = %.3e;  %s;  Kantenzeilen (m = 10/30/100): %.2e / %.2e / %.2e"
              % (N, np.linalg.norm(XN, 2), nD, "; ".join(out), *zeilen), flush=True)


if __name__ == "__main__":
    main()
