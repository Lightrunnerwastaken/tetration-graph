"""k30 — Eckprofil des PRIMALEN Kopfes auf der Kegelskala: Taylor-Koeffizienten bei L, Konvergenz in n, Vergleich mit dem Kandidaten
P(u) = Σ_k Γ(β_P − iαk)(−u)^k/k!  mit β_P = ½ + iα (Pole von Γ(½ + iα(s+1)) = rechte Nullstellen von 1 + μ^{−s−1}, primale Kopfgleichung
P(μz) + μP(z) = Quelle + Überlauf).  Zur Kontrolle der Pipeline auch der adjungierte Kopf in w = ρ/x (muss Gs S_n reproduzieren).

Köpfe in arb: J̃_n c = e_n,  J̃_nᵀ ℓ = e_n.  Primal δ(x) = Σ c_k (x/ρ)^k, Taylor in ζ mit x = L + d·ζ, d = −(L²/ρ)·n^{−a}
(so dass ζ dieselbe Koordinate ist wie Gs ζ_G = (w − w_L)/n^{−a} in linearer Näherung).  Adjungiert 𝒜(w) = Σ ℓ_k w^k, w = w_L + n^{−a}ζ.
Aufruf: python k30_primal_profil.py 200,300,400
"""
import os
import sys
import time

import mpmath as mp
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_poly, ctx  # noqa: E402
from k_kern import J_tilde, Lacb, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
Lm = mp.mpc(L_RE, L_IM)
ALPHA = Lm / (4 * mp.pi)
A_EXP = 1 / (2j * mp.pi / (2j * mp.pi - Lm / 2))
KT = 30  # Zahl der Taylor-Koeffizienten


def taylor_bei(poly_coeffs, x0, d, K):
    """Taylor-Koeffizienten von p(x0 + d·ζ) bis ζ^{K−1} (acb, exakt über Horner-Shift)."""
    n = len(poly_coeffs) - 1
    # Horner-Taylor-Shift: b = Koeffizienten von p(x0 + y) in y
    b = list(poly_coeffs)
    for i in range(n):
        for j in range(n - 1, i - 1, -1):
            b[j] = b[j] + x0 * b[j + 1]
    out = []
    dp = acb(1)
    for k in range(K):
        out.append(b[k] * dp)
        dp = dp * d
    return out


def q_test(c):
    r = [k * c[k] / c[k - 1] for k in range(1, len(c))]
    return [r[i] / r[i - 1] for i in range(1, len(r))]


def q_modell(beta, K):
    g = lambda j: mp.loggamma(beta - 1j * ALPHA * j)
    return [complex(mp.exp(g(k) + g(k - 2) - 2 * g(k - 1))) for k in range(2, K)]


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    res = {}
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        e = arb_mat(n + 1, 1)
        e[n, 0] = arb(1)
        c = J.solve(e, algorithm="approx")
        l = J.transpose().solve(e, algorithm="approx")
        rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        L = Lacb()
        na = complex(mp.power(n, -A_EXP))
        d_x = -(L * L / acb(rho)) * acb(na.real, na.imag)
        # primal: Koeffizienten in x
        pc = [acb(c[k, 0]) / acb(rho) ** k for k in range(n + 1)]
        tp = taylor_bei(pc, L, d_x, KT)
        # adjungiert in w
        wl = acb(rho) / L
        ac = [acb(l[k, 0]) for k in range(n + 1)]
        ta = taylor_bei(ac, wl, acb(na.real, na.imag), KT)
        tp = np.array([complex(x.real.mid(), x.imag.mid()) for x in tp])
        ta = np.array([complex(x.real.mid(), x.imag.mid()) for x in ta])
        res[n] = (tp / tp[0], ta / ta[0])
        print("n = %d (%.0f s): |δ(L)| = %.3e, |𝒜(w_L)| = %.3e" % (n, time.time() - t0, abs(tp[0]), abs(ta[0])), flush=True)
        print("   primal  c_1..c_5: " + " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in res[n][0][1:6]))
        print("   adjung. c_1..c_5: " + " ".join("%.5f%+.5fi" % (x.real, x.imag) for x in res[n][1][1:6]))
        qp = q_test(res[n][0])
        qa = q_test(res[n][1])
        print("   primal  q_2..q_9: " + " ".join("%.4f%+.4fi" % (complex(x).real, complex(x).imag) for x in qp[:8]))
        print("   adjung. q_2..q_9: " + " ".join("%.4f%+.4fi" % (complex(x).real, complex(x).imag) for x in qa[:8]))
        np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "k30_profile_%d.npz" % n), primal=res[n][0], adj=res[n][1])
    for name, beta in [("β_P = ½ + iα", 0.5 + 1j * ALPHA), ("β_A = ½ − 2iα", 0.5 - 2j * ALPHA)]:
        print("Modell %s: q_2..q_9: %s" % (name, " ".join("%.4f%+.4fi" % (x.real, x.imag) for x in q_modell(beta, 10))))


if __name__ == "__main__":
    main()
