"""k35 — Wo sitzen die Normen der Köpfe?  Primal: ‖c_n‖² = (1/2π)∫|δ(ρe^{iφ})|²dφ (Kreis).  Adjungiert (w-Norm mit Linse, K1″(b)):
‖r̂_n‖² = (1/2π)∫|Λ(γ(α(θ)))|²dθ, γ = Linsenkurve (k04: S = 0.0166, τ₁ = 12°, p = 1.3), α(θ) = Theodorsen-Randkorrespondenz.
Ausgabe: kumulierte Masse als Funktion des Abstands |x − L| in Kegeleinheiten A = |x − L|/(ρ n^{−Re a}) (beide Ecken L, L̄ zusammen),
Vergleich der Gesamtnormen mit k12, und die Profile |Λ(γ)|/|Λ(L)| längs γ gegen das Modell |Y(u)| (u = λ_A ζ_G, Gs Koordinate).
Aufruf: python k35_normverteilung.py 300,600
"""
import math
import os
import sys
import time

import mpmath as mp
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, Lc, RHO  # noqa: E402
import k04_konform as K4  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
RE_A = 0.8935861640191741
IM_A = 0.02531610080330162
LAM_A = complex(-0.203297, -1.406008)
ALPHA = None  # mpc, in main gesetzt (Doppelgenauigkeit verfälscht die Reihe bei |u| > 25, Fehler in der ersten Fassung)


def theodorsen_alpha(M):
    th = 2 * math.pi * np.arange(M) / M
    n = np.fft.fftfreq(M, d=1.0 / M)
    alpha = th.copy()
    for _ in range(80):
        g = K4.s_profil(alpha)
        Hg = np.real(np.fft.ifft(-1j * np.sign(n) * np.fft.fft(g)))
        an = th + Hg
        if np.max(np.abs(an - alpha)) < 1e-15:
            alpha = an
            break
        alpha = an
    return th, alpha


def horner(coef, q):
    s = acb(0)
    for k in range(len(coef) - 1, -1, -1):
        s = s * q + coef[k]
    return s


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    M = 2 ** 13
    th, alpha = theodorsen_alpha(M)
    gam = RHO * np.exp(K4.s_profil(alpha)) * np.exp(1j * alpha)
    argL = math.atan2(Lc.imag, Lc.real)
    mp.mp.dps = 60
    global ALPHA
    ALPHA = mp.mpc(L_RE, L_IM) / (4 * mp.pi)
    cY = [mp.gamma(0.5 - 2j * ALPHA - 1j * ALPHA * k) / (mp.factorial(k) * mp.gamma(0.5 - 2j * ALPHA)) for k in range(400)]

    def Ymod(u):
        s = mp.mpc(0); p = mp.mpc(1)
        for c in cY:
            s += c * p; p *= -u
        return complex(s)

    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
        c = J.solve(e, algorithm="approx")
        l = J.transpose().solve(e, algorithm="approx")
        rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        cf = [acb(c[k, 0]) for k in range(n + 1)]
        lf = [acb(l[k, 0]) for k in range(n + 1)]
        scale = RHO * n ** (-RE_A)
        # primal auf dem Kreis
        phi = 2 * math.pi * np.arange(M) / M
        dv = np.array([complex(horner(cf, acb(math.cos(p), math.sin(p))).mid()) for p in phi])  # δ(ρe^{iφ}) = Σ c_k e^{ikφ}
        lv = np.array([complex(horner(lf, acb(rho) / acb(z.real, z.imag)).mid()) for z in gam])
        L = complex(Lc)
        dL = complex(horner(cf, acb(L.real, L.imag) / acb(rho)).mid())
        lL = complex(horner(lf, acb(rho) / acb(L.real, L.imag)).mid())
        nc = math.sqrt(np.mean(np.abs(dv) ** 2))
        nl = math.sqrt(np.mean(np.abs(lv) ** 2))
        print("n = %d (%.0f s): ‖c‖ = %.4e (Koeff.-Norm %.4e), ‖r̂‖ ≈ %.4e;  |δ(L)| = %.4e, |Λ(L)| = %.4e; τ ≈ %.4f; pivot c_nn = %.4e"
              % (n, time.time() - t0, nc, float(sum(abs(x) ** 2 for x in cf).sqrt().mid()), nl, abs(dL), abs(lL), nc * nl / abs(complex(cf[n].mid())),
                 abs(complex(cf[n].mid()))), flush=True)
        # Abstand in Kegeleinheiten (beide Ecken)
        Ac = np.minimum(np.abs(RHO * np.exp(1j * phi) - L), np.abs(RHO * np.exp(1j * phi) - L.conjugate())) / scale
        Ag = np.minimum(np.abs(gam - L), np.abs(gam - L.conjugate())) / scale
        for name, A_, v, tot in [("primal (Kreis)", Ac, dv, nc), ("adjung. (Linse)", Ag, lv, nl)]:
            row = []
            for Am in [1, 2, 4, 8, 16, 32, 64, 128, 1e9]:
                frac = math.sqrt(np.sum(np.abs(v[A_ <= Am]) ** 2) / np.sum(np.abs(v) ** 2))
                row.append("A≤%s: %.3f" % (("%g" % Am) if Am < 1e8 else "∞", frac))
            imax = np.argmax(np.abs(v))
            print("   %s: Normanteil sqrt(Masse) %s | Maximum bei A = %.1f (|·|/Ecke = %.3g)" % (name, " ".join(row), A_[imax], np.abs(v[imax]) / (abs(dL) if v is dv else abs(lL))))
        # Profil längs γ nahe L (Γin-Seite) gegen Modell
        wL = RHO / L
        sel = [j for j in range(M) if abs(gam[j] - L) / scale <= 40 and gam[j].imag > 0]
        sel = sorted(sel, key=lambda j: alpha[j])
        out = []
        for j in sel[:: max(1, len(sel) // 14)]:
            w = RHO / gam[j]
            zG = (w - wL) / (n ** (-(RE_A + 1j * IM_A)))
            u = LAM_A * zG
            out.append("A=%+.1f: %.3g/%.3g" % (np.sign(alpha[j] - argL) * abs(gam[j] - L) / scale, abs(lv[j]) / abs(lL), abs(Ymod(u))))
        print("   adjung. |Λ(γ)|/|Λ(L)| gegen |Y(u)|: " + " ".join(out), flush=True)


if __name__ == "__main__":
    main()
