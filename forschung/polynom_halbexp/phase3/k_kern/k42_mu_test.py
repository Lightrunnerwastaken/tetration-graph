"""k42 — Unabhängige Prüfung der Eckprofil-Formel und der Drehung an kubischen Modellen mit anderem Multiplikator μ (L fest, arg μ = arg √L,
|μ| ∈ {1.6, 2.0}): α = log μ/(2π), a = 1 + iα (Im a = ln|μ|/2π: 0.075 bzw. 0.110 statt 0.025 — die Drehung ist 3–4× schneller).
Für jedes n: adjungierter und primaler Kopf (arb), Taylor bei w_L (adj., Skala n^{−a}) bzw. L (primal, Skala −(L²/ρ)n^{−a}); q_k-Test gegen
β_A = ½ − 2iα, β_P = ½ + iα; Stationarität des Profils in komplexer gegen reelle Skala (Wiener-Abstand zwischen n und 2n).
Aufruf: python k42_mu_test.py 1.6 100,200,400   (Präzision 2n + 400 Bit; Kontrolle mit 1.5n + 300 bei kleinstem n)
"""
import os
import sys
import time

import mpmath as mp
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx  # noqa: E402
from k_kern import J_tilde, Lacb, L_RE, L_IM  # noqa: E402
from k30_primal_profil import taylor_bei, q_test  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
Lm = mp.mpc(L_RE, L_IM)
KT = 24


def q_modell(beta, al, K):
    g = lambda j: mp.loggamma(beta - 1j * al * j)
    return [complex(mp.exp(g(k) + g(k - 2) - 2 * g(k - 1))) for k in range(2, K)]


def koepfe(mu_abs, n, bits):
    mu_m = mu_abs * mp.expjpi(mp.arg(mp.sqrt(Lm)) / mp.pi)
    mu = acb(arb(str(mp.re(mu_m))), arb(str(mp.im(mu_m))))
    J = J_tilde("kubisch", n, bits=bits, mu=mu)
    ctx.prec = bits
    e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
    c = J.solve(e, algorithm="approx")
    l = J.transpose().solve(e, algorithm="approx")
    return c, l, mu_m


def main():
    mu_abs = float(sys.argv[1])
    Ns = [int(s) for s in sys.argv[2].split(",")]
    mu_m = mu_abs * mp.expjpi(mp.arg(mp.sqrt(Lm)) / mp.pi)
    al = mp.log(mu_m) / (2 * mp.pi)
    a = 1 + 1j * al
    print("|μ| = %.2f: α = %s, a = %s, Drehung %.2f°/e-Faktor" % (mu_abs, mp.nstr(al, 6), mp.nstr(a, 6), float(mp.degrees(mp.im(a)))), flush=True)
    profs = {}
    for n in Ns:
        for bits in ([int(1.5 * n) + 300, 2 * n + 400] if n == Ns[0] else [2 * n + 400]):
            t0 = time.time()
            c, l, _ = koepfe(mu_abs, n, bits)
            rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
            L = Lacb()
            na = complex(mp.power(n, -a))
            dA = acb(na.real, na.imag)
            dP = -(L * L / acb(rho)) * dA
            ta = taylor_bei([acb(l[k, 0]) for k in range(n + 1)], acb(rho) / L, dA, KT)
            tp = taylor_bei([acb(c[k, 0]) / acb(rho) ** k for k in range(n + 1)], L, dP, KT)
            ta = np.array([complex(x.real.mid(), x.imag.mid()) for x in ta]); ta /= ta[0]
            tp = np.array([complex(x.real.mid(), x.imag.mid()) for x in tp]); tp /= tp[0]
            print("n = %d, %d Bit (%.0f s): adj c1..c3 %s | prim c1..c3 %s" % (n, bits, time.time() - t0,
                  " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in ta[1:4]), " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in tp[1:4])), flush=True)
        profs[n] = (ta, tp)
        qa, qp = q_test(ta), q_test(tp)
        print("   q_k adj (k=2..7):  " + " ".join("%.4f%+.4fi" % (complex(x).real, complex(x).imag) for x in qa[:6]))
        print("   q_k prim (k=2..7): " + " ".join("%.4f%+.4fi" % (complex(x).real, complex(x).imag) for x in qp[:6]), flush=True)
    print("Modell adj  β_A = ½ − 2iα: " + " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in q_modell(0.5 - 2j * al, al, 8)))
    print("Modell prim β_P = ½ + iα:  " + " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in q_modell(0.5 + 1j * al, al, 8)))
    # Stationarität: Profil bei 2n gegen n, komplexe Skala (wie berechnet) vs reelle Skala (umrechnen: Faktor (n2/n1)^{i Im a} pro Potenz)
    ns = sorted(profs)
    for i in range(len(ns) - 1):
        n1, n2 = ns[i], ns[i + 1]
        for idx, name in [(0, "adj"), (1, "prim")]:
            s1, s2 = profs[n1][idx], profs[n2][idx]
            w = 1.05 ** np.arange(KT)
            dk = np.sum(np.abs(s2 - s1) * w)
            rot = complex(mp.power(mp.mpf(n2) / n1, 1j * mp.im(a)))  # Umrechnung auf reelle Skala n^{−Re a}
            s2r = s2 * rot ** np.arange(KT)
            dr = np.sum(np.abs(s2r - s1) * w)
            print("   %s %d→%d: Wiener-Abstand komplexe Skala %.4f, reelle Skala %.4f" % (name, n1, n2, dk, dr), flush=True)


if __name__ == "__main__":
    main()
