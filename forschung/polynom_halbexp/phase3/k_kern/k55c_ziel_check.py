"""k55c — Ist die zusammengesetzte Funktion c_n = s_reg + (C/d_n)V((x−L)/d_n) überhaupt das richtige Ziel? Vergleich mit der exakten Eck-Pol-Antwort u_n
(J̃_n⁻¹ b, b = ρ-Koeff. von 1/(L − x)) an Punkten längs Σ₁ (u_A = 80°), Σ₂ (u_A = −110°) und auf f(Γout) nahe L. Kubisch.
Aufruf: python k55c_ziel_check.py 120
"""
import math
import os
import sys

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def main():
    n = int(sys.argv[1])
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    al = L / (4 * mp.pi); a = 1 + 1j * al; mu = mp.sqrt(L)
    bP = 0.5 + 1j * al
    lamP = mp.mpc("0.603424", "1.165020"); lamPt = -lamP * rho ** 2 / L ** 2
    Ck = -lamPt * mp.gamma(0.5 + 1j * al) / ((mu + 1 / mu) * mp.gamma(0.5 + 2j * al))
    cs = -1 / (mu + 1 / mu)
    dn = rho * mp.power(n, -a)
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
    gb = mp.gamma(bP)
    cV = [mp.gamma(bP - 1j * al * k) / (mp.factorial(k) * gb) / (k + 1) for k in range(1500)]
    def V(z):
        w = -lamPt * z; s_ = mp.mpc(0); p = mp.mpc(1)
        for c in cV:
            s_ += c * p; p *= w
        return s_
    ck = [complex(float(x.mid())) for x in kubik_koeff()]
    def f(x): return ck[0] + x * (ck[1] + x * (ck[2] + x * ck[3]))
    print("n = %d: Vergleich u_n (exakt) mit (C/d_n)V + c_s/(x−L)-freiem Teil; Spalten: Punkt, |x−L|, |ζ|, u_n·d_n, (C)V(ζ), Verhältnis" % n)
    pts = []
    for uA, nm in ((80.0, "Σ1"), (-110.0, "Σ2")):
        for r in (0.005, 0.02, 0.05, 0.1):
            th = math.radians(uA + 71.46) + math.radians(1.623) * math.log(r / float(rho))
            pts.append(("%s r=%.3f" % (nm, r), complex(L) + r * complex(math.cos(th), math.sin(th))))
    for r in (0.02, 0.05, 0.1, 0.2):
        # Punkt auf dem Γout-Bogen im Abstand r von L, dann f-Bild
        ph = float(mp.arg(L)) - 2 * math.asin(r / (2 * float(rho)))
        y = float(rho) * complex(math.cos(ph), math.sin(ph))
        pts.append(("f(Γout) r_y=%.3f" % r, f(y)))
    for name, x in pts:
        xm = mp.mpc(x)
        z = (xm - L) / dn
        a1 = un(xm) * dn
        a2 = Ck * V(z)
        print("  %-18s |x−L| = %.3f  |ζ| = %7.2f  d_n·u_n = %-24s  C·V = %-24s  |Verh.| = %s" % (name, abs(x - complex(L)), float(abs(z)), mp.nstr(a1, 6), mp.nstr(a2, 6), mp.nstr(abs(a1 / a2), 4)))


if __name__ == "__main__":
    main()
