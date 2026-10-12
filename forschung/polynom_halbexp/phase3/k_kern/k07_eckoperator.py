"""k07 — Der Eckoperator der Abschnittsmethode: E_N := P_N^w Ĵ Q_N 𝒥^{-1} P_N^w  (Q_N = I − P_N, 𝒥 = Ĵ_∞ = C_k J̃_∞).

Exakt (bewiesen, siehe BERICHT §8): Ĵ_N · (P_N 𝒥^{-1} P_N^w) = I − E_N auf P_N(w); ist I − E_N invertierbar, so
  Ĵ_N^{-1} = P_N 𝒥^{-1} P_N^w (I − E_N)^{-1},  also ‖Ĵ_N^{-1}‖ ≤ ‖𝒥^{-1}‖ · ‖(I − E_N)^{-1}‖.
Numerisch wird 𝒥^{-1} durch Ĵ_{N'}^{-1} (N' ≫ N) ersetzt: E_N ≈ I − Ĵ_N · [Ĵ_{N'}^{-1}]_{N×N}.
Gemessen: Singulärwerte von E_N und von (I − E_N)^{-1}, Lage (Winkel, Index) der führenden Singulärvektoren von E_N.
Aufruf: python k07_eckoperator.py kubisch 600 100,150,200
"""
import math
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, Lc  # noqa: E402
import k04_konform as K4  # noqa: E402
from k05_widom import inv, block  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


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
    Xr = inv(Jh)
    print("Referenz N' = %d (%.0f s)" % (Nr, time.time() - t0), flush=True)
    th = np.linspace(0, math.pi, 721)
    kw = np.exp(1j * th) * np.exp(np.polynomial.polynomial.polyval(np.exp(1j * th), F[:4096]))
    argk = np.degrees(np.angle(kw))
    for N in Ns:
        n = N + 1
        JN = block(Jh, n)
        XB = block(Xr, n)
        E = arb_mat(n, n)
        P = JN * XB
        for i in range(n):
            for j in range(n):
                E[i, j] = (arb(1) if i == j else arb(0)) - P[i, j]
        Ef = K4.float_mat(E)
        U, S, Vt = np.linalg.svd(Ef)
        sIE = np.linalg.svd(np.eye(n) - Ef, compute_uv=False)
        print("\nN = %d: σ(E_N) = %s ;  ‖(I − E_N)⁻¹‖ = %.4f, σ_min(I − E_N) = %.4f, Anzahl σ(E_N) > 0.01: %d"
              % (N, " ".join("%.3g" % s for s in S[:8]), 1 / sIE[-1], sIE[-1], int(np.sum(S > 0.01))))
        for r in range(2):
            v = Vt[r, :]
            u = U[:, r]
            mv = [np.sum(v[int(a * n / 10):int((a + 1) * n / 10)] ** 2) for a in range(10)]
            mu = [np.sum(u[int(a * n / 10):int((a + 1) * n / 10)] ** 2) for a in range(10)]
            wv = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * th), v))
            wu = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * th), u))
            print("  r=%d σ=%.3g  rechts (w-Eingang): Masse je Zehntel %s, Winkel %.1f° | links (w-Ausgang): %s, Winkel %.1f°"
                  % (r, S[r], " ".join("%.2f" % m for m in mv), argk[np.argmax(wv)], " ".join("%.2f" % m for m in mu), argk[np.argmax(wu)]))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
