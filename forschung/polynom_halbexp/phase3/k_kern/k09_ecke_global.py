"""k09 — Ist das Eck-Grenzobjekt (Widom-Korrektur D_N an L) lokal bestimmt? Test mit einer Störung, die in L, L̄ von 4. Ordnung verschwindet.

Modelle: kubisch, jet4_ε (f + (ε/ρ⁷)(x² − 2Re L x + ρ²)⁴; L, μ, f'', f''' in L unverändert, global anders).
Gemessen: σ₁(J̃_{N'}⁻¹) (globale Eckrate), ‖Ĵ_N⁻¹‖, Singulärwerte von D_N = Ĵ_N⁻¹ − P_NĴ_{N'}⁻¹P_N (N' ≥ 3N), Lage der Vektoren.
Lokal ⇔ D_N ändert sich (bis auf o(ε)) nicht, während σ₁ sich in 1. Ordnung ändert.
Aufruf: python k09_ecke_global.py 600 200 kubisch,jet4_1e-4,jet4_-1e-4
"""
import math
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import ctx  # noqa: E402
from k_kern import J_tilde  # noqa: E402
import k04_konform as K4  # noqa: E402
from k05_widom import inv, block  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    Nr = int(sys.argv[1])
    N = int(sys.argv[2])
    modelle = sys.argv[3].split(",")
    F, d = K4.theodorsen()
    bits = int(1.2 * Nr) + 300
    for name in modelle:
        t0 = time.time()
        J = J_tilde(name, Nr, bits=bits)
        ctx.prec = bits
        C = K4.C_matrix(F, Nr)
        Jh = C * J
        Xr = K4.float_mat(inv(Jh))
        Yr = K4.float_mat(inv(J))
        sJ = np.linalg.svd(Yr, compute_uv=False)
        n = N + 1
        XN = K4.float_mat(inv(block(Jh, n)))
        D = XN - Xr[:n, :n]
        sD = np.linalg.svd(D, compute_uv=False)
        print("%-12s σ₁(J̃_%d⁻¹) = %.5g, σ₂ = %.5g | ‖Ĵ_%d⁻¹‖ = %.6f, ‖Ĵ_%d⁻¹‖ = %.6f | σ(D_%d) = %s  (%.0f s)"
              % (name, Nr, sJ[0], sJ[1], Nr, np.linalg.norm(Xr, 2), N, np.linalg.norm(XN, 2), N,
                 " ".join("%.5f" % s for s in sD[:6]), time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
