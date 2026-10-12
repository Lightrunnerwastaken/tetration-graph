"""k47 — Innenprofil der Antwort auf den ECK-POL: u_n = J_n⁻¹T_nk_L (k_L = 1/(L − x)), auf der Kegelskala x = L + d_nζ, d_n = ρn^{−a} (komplex).
Vorhersage (BERICHT §25): d_n·u_n(L + d_nζ) → C·V(ζ), V(ζ) = ∫₀¹F_P(λ̃_P sζ)ds, λ̃_P = −λ_Pρ²/L², und
C = −λ̃_P Γ(½ + iα)/((μ + 1/μ)Γ(½ + 2iα)) aus der Anpassung von V(ζ) ~ Γ(½+2iα)/(Γ(½+iα)λ̃_Pζ) an den äusseren Pol −1/((μ+1/μ)(x − L)).
Aufruf: python k47_polantwort_profil.py 200,400,800
"""
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    mp.mp.dps = 40
    Lm = mp.mpc(L_RE, L_IM)
    rho = abs(Lm)
    al = Lm / (4 * mp.pi)
    a = 1 + 1j * al
    mu = mp.sqrt(Lm)
    bp = 0.5 + 1j * al
    lamP = mp.mpc("0.603424", "1.165020")
    lamPt = -lamP * rho ** 2 / Lm ** 2
    C = -lamPt * mp.gamma(0.5 + 1j * al) / ((mu + 1 / mu) * mp.gamma(0.5 + 2j * al))
    cP = [mp.gamma(bp - 1j * al * k) / (mp.factorial(k) * mp.gamma(bp)) for k in range(300)]

    def V(z):  # ∫₀¹ F_P(λ̃ s z) ds = Σ c_k (−λ̃ z)^k/(k+1)
        s = mp.mpc(0); p = mp.mpc(1); w = -lamPt * z
        for k, c in enumerate(cP):
            s += c * p / (k + 1); p *= w
        return s
    print("C (Vorhersage) = %s" % mp.nstr(C, 8), flush=True)
    zs = [0, 0.5, 1, 2, 4, 1j, 2j, -1, -2, -2j, 2 + 2j, -2 + 2j]
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        Jc = acb_mat(J)
        rhoa = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        La = acb(arb(L_RE), arb(L_IM))
        b = acb_mat(n + 1, 1)
        q = acb(rhoa) / La
        p = acb(1) / La
        for k in range(n + 1):     # ρ-skaliert: [x^k]k_L · ρ^k = ρ^k/L^{k+1}
            b[k, 0] = p
            p = p * q
        X = Jc.solve(b, algorithm="approx")
        dn = rho * mp.power(n, -a)
        row = []
        for z in zs:
            x = Lm + dn * z
            xx = acb(arb(str(mp.re(x))), arb(str(mp.im(x)))) / acb(rhoa)
            s = acb(0)
            for k in range(n, -1, -1):
                s = s * xx + X[k, 0]
            u = mp.mpc(float(s.real.mid()), float(s.imag.mid()))
            val = dn * u
            mod = C * V(z)
            row.append("ζ=%s: %s / %s" % (mp.nstr(mp.mpc(z), 3), mp.nstr(val, 5), mp.nstr(mod, 5)))
        print("n = %d (%.0f s):\n   " % (n, time.time() - t0) + "\n   ".join(row), flush=True)


if __name__ == "__main__":
    main()
