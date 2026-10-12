"""k52 — Adjungierte Rate auf Kreisen |t| = ρ(1+s) (ausserhalb D) und Exponenten-Budget (BERICHT §30):
  h_n(t) := (1/n)·ln|Λ_n(t)/Λ_n(L)|   für s ∈ {0, 0.005, 0.0166, 0.03}, alle Winkel, n aus der Liste.
Ausgabe je s: max_θ h_n (Lage) ausserhalb der Winkel-Umgebung ±w° von L, L̄ (w = 3, 10), und max über alle θ.
Aufruf: python k52_adjungiert_aussen.py 300,600
"""
import math
import os
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, acb_poly, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def main():
    Ns = [int(s) for s in sys.argv[1].split(",")]
    Lf = complex(float(L_RE[:20]), float(L_IM[:20]))
    argL = math.degrees(math.atan2(Lf.imag, Lf.real))
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        Jc = acb_mat(J)
        e = acb_mat(n + 1, 1); e[n, 0] = acb(1)
        l = Jc.transpose().solve(e, algorithm="approx")
        rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        L = acb(arb(L_RE), arb(L_IM))
        Lam = acb_poly([l[k, 0] for k in range(n + 1)])
        LL = Lam(rho / L)
        rhof = float(rho.mid())
        print("n = %d (%.0f s)" % (n, time.time() - t0), flush=True)
        for s in (0.0, 0.005, 0.0166, 0.03):
            vals = []
            for j in range(1440):
                th = 360.0 * j / 1440
                r = rhof * (1 + s)
                t = acb(r * math.cos(math.radians(th)), r * math.sin(math.radians(th)))
                v = abs(Lam(rho / t) / LL)
                h = float((v.log() / n).mid())
                d = min(abs(th - argL), abs(th - (360 - argL)))
                vals.append((h, th, d))
            row = []
            for w in (3, 10, 0):
                cand = [x for x in vals if x[2] >= w] if w else vals
                m = max(cand)
                row.append("ohne ±%d°: max h = %+.5f bei %.2f°" % (w, m[0], m[1]) if w else "alle: max h = %+.5f bei %.2f°" % (m[0], m[1]))
            gout = max(x for x in vals if (x[1] < argL - 10 or x[1] > 360 - argL + 10))
            gin = max(x for x in vals if (argL + 10 < x[1] < 360 - argL - 10))
            print("  s = %.4f: %s | Γout-Bogen (±10° frei): max %+.5f bei %.1f° | Γin-Bogen: max %+.5f bei %.1f°"
                  % (s, " | ".join(row), gout[0], gout[1], gin[0], gin[1]), flush=True)


if __name__ == "__main__":
    main()
