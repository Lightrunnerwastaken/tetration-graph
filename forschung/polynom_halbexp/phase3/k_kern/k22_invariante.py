"""k22 — Erster Test zu H-A (Abel/Koenigs-Torus): ist die Paarungsform ω = m·Λ·δ·dy/y der Köpfe unter der lokalen Dynamik f invariant?

Herleitung (BERICHT §16): Erfüllen δ (primal) und Λ (adjungiert, Cauchy-Paarung ∮Λδ dx/x) nahe L die homogenen Eckgleichungen
δ∘f = −m·δ und (ℒΛ) + mΛ = 0 (ℒ = lokaler Transferoperator, Zweig g = f⁻¹ bei L), so ist ω(y) := m(y)Λ(y)δ(y)dy/y unter y ↦ f(y) invariant,
in der Abelkoordinate s also eine ½-periodische Dichte W(s) (Funktion auf dem Torus). Koordinatenfrei prüfbar ohne Koenigs-Funktion:
  Q(y) := [m(f(y))·Λ(f(y))·δ(f(y))·f'(y)/f(y)] / [m(y)·Λ(y)·δ(y)/y]  ≈ 1 ?
Köpfe: δ_n(x) = Σ_k c_{n,k}(x/ρ)^k mit J̃_n c_n = e_n, Λ_n(t) = Σ_k ℓ_{n,k}(ρ/t)^k mit J̃_nᵀℓ_n = e_n (in arb, volle Präzision).
Zusätzlich die lokalen Homogenitätsdefekte: primal |δ∘f + mδ|/|mδ|, adjungiert |ℒΛ + mΛ|/|mΛ| (ℒΛ(y) = yΛ(g(y))g'(y)/g(y)).
Aufruf: python k22_invariante.py 200,300,400
"""
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ctx  # noqa: E402
from k_kern import J_tilde, kubik_koeff, Lacb, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def solve_vec(A, b):
    return A.solve(b, algorithm="approx")


THETAS = [200, 240, 280]


def main():
    global THETAS
    if len(sys.argv) > 2:
        THETAS = [int(v) for v in sys.argv[2].split(",")]
    Ns = [int(s) for s in sys.argv[1].split(",")]
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        e = arb_mat(n + 1, 1)
        e[n, 0] = arb(1)
        c = solve_vec(J, e)
        l = solve_vec(J.transpose(), e)
        cf = [acb(c[k, 0]) for k in range(n + 1)]
        lf = [acb(l[k, 0]) for k in range(n + 1)]
        co = [acb(x) for x in kubik_koeff()]
        f = acb_poly(co)
        fp = f.derivative()
        L = Lacb()
        rho = acb((arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt())
        dpoly = acb_poly([cf[k] / rho ** k for k in range(n + 1)])          # δ(x)
        def Lam(t):
            q = rho / t
            s = acb(0)
            for k in range(n, -1, -1):
                s = s * q + lf[k]
            return s
        def m(x):
            return fp(f(x))
        def g(y):   # lokaler Umkehrzweig von f bei L (Newton vom Startwert L + (y − L)/μ)
            mu = fp(L)
            x = L + (y - L) / mu
            for _ in range(60):
                x = x - (f(x) - y) / fp(x)
            return x
        print("n = %d (%.0f s): ‖c‖-Grösse %s, ‖ℓ‖-Grösse %s" % (n, time.time() - t0,
              "%.2e" % max(abs(complex(x.mid())) for x in cf), "%.2e" % max(abs(complex(x.mid())) for x in lf)), flush=True)
        import cmath
        for r in [0.03, 0.06, 0.12]:
            for th in THETAS:
                y = L + acb(r) * acb(cmath.exp(1j * th * cmath.pi / 180))
                fy = f(y)
                Ly, Lfy = Lam(y), Lam(fy)
                dy, dfy = dpoly(y), dpoly(fy)
                num = m(fy) * Lfy * dfy * fp(y) / fy
                den = m(y) * Ly * dy / y
                Qv = complex((num / den).mid())
                prim = abs(complex(((dfy + m(y) * dy) / (m(y) * dy)).mid()))
                gy = g(y)
                LL = y * Lam(gy) * (1 / fp(gy)) / gy           # ℒΛ(y) = yΛ(g(y))g'(y)/g(y), g' = 1/f'(g)
                adj = abs(complex(((LL + m(y) * Ly) / (m(y) * Ly)).mid()))
                print("   y = L + %.2f·e^{i%d°}: Q = %+.4f%+.4fi (|Q| = %.4f)   Defekt primal %.3e, adjungiert %.3e   |δ(y)| = %.2e, |Λ(y)| = %.2e"
                      % (r, th, Qv.real, Qv.imag, abs(Qv), prim, adj, abs(complex(dy.mid())), abs(complex(Ly.mid()))), flush=True)


if __name__ == "__main__":
    main()
