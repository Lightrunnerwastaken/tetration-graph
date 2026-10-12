"""k40 — Direkter Test der Null-Spirale am echten adjungierten Kopf: |Λ_n(L + r e^{iθ})|/|Λ_n(L)| längs Strahlen θ = θ_c − φ (θ_c = Γin-Tangente,
φ = Winkel nach aussen) für r = A·ρn^{−Re a}, gegen das Eckprofil F_A(u) (mpmath-α; u = λ_A ζ_G mit exaktem w = ρ/x).
Vorhersage (k39): Null-Spirale bei φ_null(r) ≈ 6.7° + 3.74°·log10(0.01/r); für φ < φ_null wächst der Kopf (relativ zum Eckwert) nach aussen.
Aufruf: python k40_strahlen_echt.py 600,1000
"""
import math
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    mp.mp.dps = 60
    Lm = mp.mpc(L_RE, L_IM)
    rho = abs(Lm)
    al = Lm / (4 * mp.pi)
    a = 1 + 1j * al
    ba = 0.5 - 2j * al
    lam = mp.mpc(-0.203297, -1.406008)
    cY = [mp.gamma(ba - 1j * al * k) / (mp.factorial(k) * mp.gamma(ba)) for k in range(700)]

    def F(u):
        s = mp.mpc(0); p = mp.mpc(1)
        for c in cY:
            s += c * p; p *= -u
        return s
    thc = mp.arg(Lm) + mp.pi / 2
    phis = [-10, -5, 0, 2.5, 5, 7.5, 10, 15, 20, 30]
    As = [1, 2, 4, 8, 12, 16, 24, 32, 48]
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
        l = J.transpose().solve(e, algorithm="approx")
        lf = [acb(l[k, 0]) for k in range(n + 1)]
        rhoa = acb((arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt())

        def Lam(x):
            q = rhoa / acb(x.real, x.imag)
            s = acb(0)
            for k in range(n, -1, -1):
                s = s * q + lf[k]
            return complex(s.mid())
        L0 = complex(Lm)
        lL = abs(Lam(L0))
        scale = float(rho) * n ** (-0.8935861640191741)
        dn = mp.power(n, -a)
        wL = rho / Lm
        print("n = %d (%.0f s), Kegelskala %.3e; Einträge log10(|Λ|/|Λ(L)|) Kopf / Modell" % (n, time.time() - t0, scale), flush=True)
        for phi in phis:
            th = float(thc) - math.radians(phi)
            row = []
            for A in As:
                r = A * scale
                x = L0 + r * complex(math.cos(th), math.sin(th))
                vk = math.log10(abs(Lam(x)) / lL)
                u = lam * (rho / mp.mpc(x) - wL) / dn
                vm = float(mp.log10(abs(F(u))))
                row.append("A%g:%+.2f/%+.2f" % (A, vk, vm))
            print("  φ = %5.1f°: %s" % (phi, " ".join(row)), flush=True)


if __name__ == "__main__":
    main()
