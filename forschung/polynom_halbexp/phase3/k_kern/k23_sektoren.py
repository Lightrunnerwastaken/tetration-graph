"""k23 — Schritt a1′ (Vergleichsproblem): in welchen Richtungen sind die Köpfe nahe L homogene Eckmoden, auf welcher Skala?

Vorhersage (BERICHT §17): primaler Kopf homogen in S_in = {θ : x, f(x) innen} = (arg L + 90°, arg L + 90° + πα′) = (166.6°, 308.3°),
adjungierter Kopf in S_out = {θ : ζ, g(ζ) aussen} = (166.6° − πα′, 166.6°) = (24.9°, 166.6°); Öffnung je πα′ = 141.7° (α′ = 1 − arg μ/π);
Lücke (−51.7°, 24.9°) um die Γout-Tangente (−13.4°) der Öffnung arg L. Auf der Kegelskala r ~ ρn^{−Re a} keiner homogen (Quelle/Überlauf).
Gemessen: Defekte |δ∘f + mδ|/|mδ| und |ℒΛ + mΛ|/|mΛ| auf Kreisen |x − L| = r·ρ·n^{−Re a}·A für A ∈ {1, 4, 16} und 24 Richtungen.
Aufruf: python k23_sektoren.py 300,600
"""
import cmath
import math
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_poly, ctx  # noqa: E402
from k_kern import J_tilde, kubik_koeff, Lacb, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
RE_A = 0.8935861640191741


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        e = arb_mat(n + 1, 1)
        e[n, 0] = arb(1)
        c = J.solve(e, algorithm="approx")
        l = J.transpose().solve(e, algorithm="approx")
        cf = [acb(c[k, 0]) for k in range(n + 1)]
        lf = [acb(l[k, 0]) for k in range(n + 1)]
        f = acb_poly([acb(x) for x in kubik_koeff()])
        fp = f.derivative()
        L = Lacb()
        rho = acb((arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt())
        dpoly = acb_poly([cf[k] / rho ** k for k in range(n + 1)])
        mu = fp(L)

        def Lam(t):
            q = rho / t
            s = acb(0)
            for k in range(n, -1, -1):
                s = s * q + lf[k]
            return s

        def m(x):
            return fp(f(x))

        def g(y):
            x = L + (y - L) / mu
            for _ in range(60):
                x = x - (f(x) - y) / fp(x)
            return x

        print("n = %d (%.0f s), Kegelskala ρn^{-Re a} = %.4f" % (n, time.time() - t0, float(abs(complex(rho.mid()))) * n ** (-RE_A)), flush=True)
        for A in [1.0, 4.0, 16.0]:
            r = float(abs(complex(rho.mid()))) * n ** (-RE_A) * A
            rowp, rowa = [], []
            for th in range(-60, 330, 15):
                y = L + acb(r) * acb(cmath.exp(1j * math.radians(th)))
                dy = dpoly(y)
                prim = abs(complex(((dpoly(f(y)) + m(y) * dy) / (m(y) * dy)).mid()))
                Ly = Lam(y)
                gy = g(y)
                LL = y * Lam(gy) * (1 / fp(gy)) / gy
                adj = abs(complex(((LL + m(y) * Ly) / (m(y) * Ly)).mid()))
                rowp.append((th, prim)); rowa.append((th, adj))
            def fmt(row):
                return " ".join("%d:%s" % (th, ("%.0e" % v) if v < 0.1 else ("%.1f" % v if v < 10 else ">10")) for th, v in row)
            print("  A = %4.1f (r = %.4f)  primal: %s" % (A, r, fmt(rowp)), flush=True)
            print("                          adjung.: %s" % fmt(rowa), flush=True)


if __name__ == "__main__":
    main()
