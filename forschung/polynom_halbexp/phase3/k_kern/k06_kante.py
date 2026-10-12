"""k06 — Gestalt der Kantenkorrektur D_N = Ĵ_N^{-1} − P_N Ĵ_{N'}^{-1} P_N (Widom-Differenz) in skalierten Koordinaten.

Für die Grenzoperator-Frage (Leitung 22:00): Hat D_N in ξ = j/N (Unbekannte, x-Grad) und η = i/N (Gleichung, w-Grad) eine N-unabhängige Gestalt?
Wo im Winkel sitzen die führenden Singulärvektoren (Ecke L oder ganzer Γin-Bogen)? Wie viele Singulärwerte sind gross (Rang)?
Ausgabe je N: Singulärwerte von D_N (die ersten 8), Indexprofil (Masse je Zehntel von N) und Winkelmaximum der ersten Singulärvektoren:
  Unbekannten-Seite u (Polynom in x auf |x| = ρ), Gleichungs-Seite v (Polynom in w auf |w| = 1, Winkel des Bildpunkts k(w)).
Aufruf: python k06_kante.py kubisch 600 150,200,300
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
argL = math.degrees(math.atan2(Lc.imag, Lc.real))


def winkelprofil_x(vec, n=721):
    phi = np.linspace(0, math.pi, n)
    v = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * phi), vec))
    return math.degrees(phi[np.argmax(v)]), v


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
    print("Referenz N' = %d (%.0f s)" % (Nr, time.time() - t0), flush=True)
    # Winkel auf |w| = 1 → Bildwinkel arg k(w)
    th = np.linspace(0, math.pi, 721)
    kw = np.exp(1j * th) * np.exp(np.polynomial.polynomial.polyval(np.exp(1j * th), F[:4096]))
    argk = np.degrees(np.angle(kw))
    for N in Ns:
        n = N + 1
        XN = K4.float_mat(inv(block(Jh, n)))
        D = XN - Xr[:n, :n]
        U, S, Vt = np.linalg.svd(D)
        print("\nN = %d: σ(D_N) = %s" % (N, " ".join("%.3g" % s for s in S[:8])))
        for r in range(3):
            u = U[:, r]
            v = Vt[r, :]
            mu = np.array([np.sum(u[int(a * n / 10):int((a + 1) * n / 10)] ** 2) for a in range(10)])
            mv = np.array([np.sum(v[int(a * n / 10):int((a + 1) * n / 10)] ** 2) for a in range(10)])
            phx, _ = winkelprofil_x(u)
            vw = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * th), v))
            phw = argk[np.argmax(vw)]
            print("  r=%d σ=%.3g  u (x-Seite): Masse je Zehntel %s  Max-Winkel %.1f°" % (r, S[r], " ".join("%.2f" % m for m in mu), phx))
            print("            v (w-Seite): Masse je Zehntel %s  Max-Winkel (Bild) %.1f°" % (" ".join("%.2f" % m for m in mv), phw))
        sys.stdout.flush()
    print("\n(arg L = %.2f°)" % argL)


if __name__ == "__main__":
    main()
