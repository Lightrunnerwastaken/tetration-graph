"""k08 — Jacobi-Reduktion: Ĵ_N^{-1} = T11 − T12 T22^{-1} T21 mit T = Ĵ_{N'}^{-1} (Blöcke: Grade ≤ N | N < Grad ≤ N').

Im Grenzfall N' → ∞ ist T = 𝒥^{-1} = J_∞^{-1}C_k^{-1} BESCHRÄNKT, und (S_γ) folgt aus sup_N ‖T22^{-1}‖ < ∞ (Schwanzabschnitte von 𝒥^{-1}).
Gemessen: σ_min(T22), ‖T12 T22^{-1} T21‖ (= Widom-Korrektur D_N), Probe der Identität, Lage der kleinsten Singulärvektoren von T22.
Aufruf: python k08_schwanz.py kubisch 600 100,150,200,300
"""
import math
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, ctx  # noqa: E402
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
    T = K4.float_mat(inv(Jh))
    print("Referenz N' = %d (%.0f s), ‖T‖ = %.4f" % (Nr, time.time() - t0, np.linalg.norm(T, 2)), flush=True)
    th = np.linspace(0, math.pi, 721)
    kw = np.exp(1j * th) * np.exp(np.polynomial.polynomial.polyval(np.exp(1j * th), F[:4096]))
    argk = np.degrees(np.angle(kw))
    for N in Ns:
        n = N + 1
        T11, T12, T21, T22 = T[:n, :n], T[:n, n:], T[n:, :n], T[n:, n:]
        U, S, Vt = np.linalg.svd(T22)
        D = T12 @ np.linalg.solve(T22, T21)
        XN = K4.float_mat(inv(block(Jh, n)))
        probe = np.linalg.norm(XN - (T11 - D), 2)
        u = U[:, -1]                     # x-Seite (Grade N+1..N'), kleinster Singulärwert
        v = Vt[-1, :]                    # w-Seite
        m = len(u)
        mu = [np.sum(u[int(a * m / 10):int((a + 1) * m / 10)] ** 2) for a in range(10)]
        phx = np.linspace(0, math.pi, 721)
        wu = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * phx), np.concatenate([np.zeros(n), u])))
        wv = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * th), np.concatenate([np.zeros(n), v])))
        print("N = %d: σ_min(T22) = %.4f (nächste %.4f, %.4f), ‖T22⁻¹‖ = %.3f, ‖T12T22⁻¹T21‖ = %.3f, Probe Jacobi %.1e;"
              " kleinster Vektor: x-Winkel %.1f°, w-Bildwinkel %.1f°, Masse je Zehntel der Schwanzgrade %s"
              % (N, S[-1], S[-2], S[-3], 1 / S[-1], np.linalg.norm(D, 2), probe, math.degrees(phx[np.argmax(wu)]),
                 argk[np.argmax(wv)], " ".join("%.2f" % x for x in mu)), flush=True)


if __name__ == "__main__":
    main()
