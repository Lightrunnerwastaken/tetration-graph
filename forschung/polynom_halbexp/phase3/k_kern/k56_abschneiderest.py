"""k56 — Hebel (i): Ist der Abschneide-Rest R_n = (I − T_n) J u_n auf dem Γout-Bogen polynomial beschränkt?
u_n = J_n⁻¹ T_n k_L (Eck-Pol-Antwort, k_L(x) = 1/(L − x)), J u = u∘f + m·u, m = f'∘f. Kubisch.
Da T_n J u_n = T_n k_L (Galerkin), gilt exakt R_n(x) = u_n(f(x)) + m(x)u_n(x) − T_n k_L(x).
Ausgabe je Punkt y = ρe^{iφ} auf dem Γout-Bogen: |f(y)|/ρ, |R_n|, |u_n(f(y))|, |k_L(y) − m u_n(y)| (≈ s(f(y))), |m u_n(y)|, |T_n k_L(y)|.
Aufruf: python k56_abschneiderest.py n   (eine Zeile pro Punkt; Raten bildet k56_raten im Bericht aus mehreren n)
"""
import math
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

# Winkel auf dem Γout-Bogen: Abstand δ (Grad) unterhalb arg L, dann über 0° bis nahe −arg L
DELTAS = [0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 76.62, 96.62, 116.62, 136.62, 146.62, 151.62, 152.74]


def main():
    n = int(sys.argv[1])
    t0 = time.time()
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    J = acb_mat(J_tilde("kubisch", n, bits=bits)); ctx.prec = bits
    bv = acb_mat(n + 1, 1)
    for k in range(n + 1):
        v = rho ** k / L ** (k + 1)
        bv[k, 0] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
    X = J.solve(bv, algorithm="approx")
    u = [mp.mpc(mp.mpf(X[k, 0].real.mid().str(nd, radius=False)), mp.mpf(X[k, 0].imag.mid().str(nd, radius=False))) for k in range(n + 1)]

    def un(x):
        s_ = mp.mpc(0); w = x / rho
        for k in range(n, -1, -1):
            s_ = s_ * w + u[k]
        return s_

    ck = [mp.mpf(c.mid().str(nd, radius=False)) for c in kubik_koeff()]

    def f(x): return ck[0] + x * (ck[1] + x * (ck[2] + x * ck[3]))

    def fp(x): return ck[1] + x * (2 * ck[2] + x * 3 * ck[3])

    argL = float(mp.arg(L))
    print("n = %d (bits %d, Löser %.0f s); arg L = %.3f°; |u_n(L)|·|d_n| zum Vergleich: |d_n| = ρ/n = %.3e" % (n, bits, time.time() - t0, math.degrees(argL), float(rho) / n))
    print("  %8s %8s %8s %12s %12s %12s %12s %12s" % ("φ/Grad", "|y−L|", "|f|/ρ", "|R_n|", "|u_n(f)|", "|s(f)|≈", "|m u_n|", "|T_n k_L|"))
    for dlt in DELTAS:
        ph = argL - math.radians(dlt)
        y = rho * mp.expjpi(ph / math.pi)
        fy = f(y); m = fp(fy)
        uy = un(y); ufy = un(fy)
        kl = 1 / (L - y)
        tk = (1 - (y / L) ** (n + 1)) / (L - y)
        R = ufy + m * uy - tk
        sf = kl - m * uy
        print("  %8.2f %8.4f %8.4f %12s %12s %12s %12s %12s" % (math.degrees(ph), float(abs(y - L)), float(abs(fy) / rho), mp.nstr(abs(R), 4), mp.nstr(abs(ufy), 4), mp.nstr(abs(sf), 4), mp.nstr(abs(m * uy), 4), mp.nstr(abs(tk), 4)))
    print("fertig n = %d nach %.0f s" % (n, time.time() - t0))


if __name__ == "__main__":
    main()
