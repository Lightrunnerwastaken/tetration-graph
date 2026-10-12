"""k43 — Test der linsenunabhängigen primalen Grenze (BERICHT §21.5) an einem Modell mit schnellerer Spirale: kubisch, |μ| = 2, arg μ = arg √L
(Im a = 0.110 statt 0.025). Primaler Kopf (arb) längs des KREISES |x| = ρ auf beiden Seiten von L, log10|δ(x)/δ(L)| für A = |x − L|/(ρ|n^{−a}|),
gegen das Modell F_{β_P}(λ_P ζ_x) (λ_P aus dem Taylor-Koeffizienten c_1 beim grössten n). Vorhersage: die Γout-Hälfte wandert mit 6.3°/e-Faktor in den
Wachstumssektor von F_P — die Kreisnorm des primalen Kopfes wird nicht-lokal.
Aufruf: python k43_mu2_primal_kreis.py 200,400,800
"""
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, Lacb, L_RE, L_IM  # noqa: E402
from k30_primal_profil import taylor_bei  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    mp.mp.dps = 60
    Lm = mp.mpc(L_RE, L_IM)
    rho = abs(Lm)
    mu_m = 2.0 * mp.expjpi(mp.arg(mp.sqrt(Lm)) / mp.pi)
    al = mp.log(mu_m) / (2 * mp.pi)
    a = 1 + 1j * al
    bp = 0.5 + 1j * al
    cP = [mp.gamma(bp - 1j * al * k) / (mp.factorial(k) * mp.gamma(bp)) for k in range(900)]

    def P(u):
        s = mp.mpc(0); p = mp.mpc(1)
        for c in cP:
            s += c * p; p *= -u
        return s
    res = {}
    lamP = None
    for n in sorted(Ns, reverse=True):
        t0 = time.time()
        bits = 3 * n + 500
        mu = acb(arb(str(mp.re(mu_m))), arb(str(mp.im(mu_m))))
        J = J_tilde("kubisch", n, bits=bits, mu=mu)
        ctx.prec = bits
        e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
        c = J.solve(e, algorithm="approx")
        rhoa = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        coef = [acb(c[k, 0]) / acb(rhoa) ** k for k in range(n + 1)]
        L = Lacb()
        dn = mp.power(n, -a)
        dx = -(Lm ** 2 / rho) * dn
        if lamP is None:
            tp = taylor_bei(coef, L, acb(complex(dx).real, complex(dx).imag), 3)
            c1 = complex(float(tp[1].real.mid()), float(tp[1].imag.mid())) / complex(float(tp[0].real.mid()), float(tp[0].imag.mid()))
            lamP = -mp.mpc(c1) * mp.gamma(bp) / mp.gamma(bp - 1j * al)
            print("λ_P (aus c_1 bei n = %d) = %s" % (n, mp.nstr(lamP, 6)), flush=True)

        def delta(x):
            s = acb(0)
            xx = acb(complex(x).real, complex(x).imag)
            for k in range(n, -1, -1):
                s = s * xx + coef[k]
            return complex(float(s.real.mid()), float(s.imag.mid()))
        dL = abs(delta(complex(Lm)))
        scale = float(rho * abs(dn))
        rows = []
        for side, sg in [("Γin ", 1), ("Γout", -1)]:
            vals = []
            for A in [2, 4, 8, 16, 32, 64]:
                tau = sg * A * scale / float(rho)
                x = rho * mp.expjpi((mp.arg(Lm) + tau) / mp.pi)
                vk = mp.log10(abs(delta(x)) / dL)
                vm = mp.log10(abs(P(lamP * (x - Lm) / dx)))
                vals.append("A%d:%+.2f/%+.2f" % (A, float(vk), float(vm)))
            rows.append("%s %s" % (side, " ".join(vals)))
        print("n = %d (%.0f s): log10|δ/δ(L)| längs des Kreises, Kopf/Modell\n   %s\n   %s" % (n, time.time() - t0, rows[0], rows[1]), flush=True)


if __name__ == "__main__":
    main()
