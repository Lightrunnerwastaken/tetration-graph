"""k24 — Aussenverhalten der Köpfe in ihren homogenen Keilen: welche Mellin-Exponenten s (Lösungen von μ^{s−1} = −1, s_k = 1 + (2k+1)πi/log μ)?
Lokale Steigung d log δ / d log r längs eines Strahls im homogenen Keil (primal θ = 240°, adjungiert θ = 90°), r = A·ρn^{−Re a}, A = 8 … 256.
Vorhersage: Re s ∈ {…, −3.45, 5.45, …}; Im-Teil zeigt sich als Phasendrehung d arg/d log r.
Aufruf: python k24_potenzen.py 300,600"""
import cmath, math, sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, acb, arb_mat, acb_poly, ctx
from k_kern import J_tilde, kubik_koeff, Lacb, L_RE, L_IM
sys.stdout.reconfigure(encoding="utf-8")
RE_A = 0.8935861640191741
def main():
    for n in [int(s) for s in sys.argv[1].split(",")]:
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits); ctx.prec = bits
        e = arb_mat(n + 1, 1); e[n, 0] = arb(1)
        c = J.solve(e, algorithm="approx"); l = J.transpose().solve(e, algorithm="approx")
        L = Lacb(); rho = acb((arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt())
        dpoly = acb_poly([acb(c[k, 0]) / rho ** k for k in range(n + 1)])
        lf = [acb(l[k, 0]) for k in range(n + 1)]
        def Lam(t):
            q = rho / t; s = acb(0)
            for k in range(n, -1, -1): s = s * q + lf[k]
            return s
        sc = float(abs(complex(rho.mid()))) * n ** (-RE_A)
        print("n = %d, Kegelskala %.4f" % (n, sc), flush=True)
        for name, fun, th in [("primal", dpoly, 240.0), ("adjung.", Lam, 90.0)]:
            As = [8, 11.3, 16, 22.6, 32, 45.3, 64, 90.5, 128, 181, 256]
            vals = []
            for A in As:
                r = sc * A
                if r > 0.5: break
                v = complex(fun(L + acb(r) * acb(cmath.exp(1j * math.radians(th)))).mid())
                vals.append((A, r, v))
            out = []
            for (A1, r1, v1), (A2, r2, v2) in zip(vals, vals[1:]):
                sl = (math.log(abs(v2)) - math.log(abs(v1))) / (math.log(r2) - math.log(r1))
                ph = (cmath.phase(v2 / v1)) / (math.log(r2) - math.log(r1))
                out.append("A≈%.0f: Re s = %+.2f, Im s ≈ %+.2f" % (math.sqrt(A1 * A2), sl, ph))
            print("  %s θ = %.0f°: " % (name, th) + " | ".join(out), flush=True)
main()
