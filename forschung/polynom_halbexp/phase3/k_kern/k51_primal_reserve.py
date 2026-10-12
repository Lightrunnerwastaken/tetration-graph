"""k51 — Exponentielle Reserve der primalen Köpfe auf K_{r0}? (BERICHT §30)
Für feste Punkte z (Abstand r0 von L, verschiedene Richtungen, dazu innere Punkte) wird R(n, z) = |δ_n(z)/δ_n(L)| für mehrere n gemessen und geprüft,
ob R exponentiell (konstante Rate (1/Δn)·ln) oder nur wie eine Potenz fällt (Kegeltheorie im Potenzsektor: R ~ (r0/r_n)^{−3.45} ~ n^{−3.08}).
Aufruf: python k51_primal_reserve.py 200,300,400,600
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
    pts = []
    for r0 in (0.1, 0.05):
        for th in (150, 168, 190, 220, 256, 290, 320, 345):
            z = Lf + r0 * complex(math.cos(math.radians(th)), math.sin(math.radians(th)))
            if abs(z) <= abs(Lf) + 1e-12:
                pts.append(("r0=%.2f θ=%d°" % (r0, th), z))
    for name, s in (("0.9L", 0.9), ("0.7L", 0.7), ("0.5L", 0.5), ("0", 0.0)):
        pts.append((name, Lf * s))
    pts.append(("ρ·e^{i·0}", complex(abs(Lf), 0)))
    pts.append(("−ρ", complex(-abs(Lf), 0)))
    res = {p[0]: [] for p in pts}
    for n in Ns:
        t0 = time.time()
        bits = int(1.3 * n) + 400
        J = J_tilde("kubisch", n, bits=bits)
        ctx.prec = bits
        Jc = acb_mat(J)
        e = acb_mat(n + 1, 1); e[n, 0] = acb(1)
        c = Jc.solve(e, algorithm="approx")
        rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
        L = acb(arb(L_RE), arb(L_IM))
        Del = acb_poly([c[k, 0] for k in range(n + 1)])
        DL = Del(L / rho)
        for name, z in pts:
            v = abs(Del(acb(z.real, z.imag) / rho) / DL)
            res[name].append(float(v.log().mid()))
        print("n = %d fertig (%.0f s)" % (n, time.time() - t0), flush=True)
    print("ln|δ_n(z)/δ_n(L)| je n = %s; Raten je Schritt (Δln/Δn) und lokaler Potenz-Exponent (Δln/Δln n):" % Ns)
    for name, _ in pts:
        v = res[name]
        rates = ["%+.5f" % ((v[i + 1] - v[i]) / (Ns[i + 1] - Ns[i])) for i in range(len(Ns) - 1)]
        pows = ["%+.2f" % ((v[i + 1] - v[i]) / math.log(Ns[i + 1] / Ns[i])) for i in range(len(Ns) - 1)]
        print("  %-16s ln R = %s | Rate %s | Potenz %s" % (name, " ".join("%.2f" % x for x in v), " ".join(rates), " ".join(pows)), flush=True)


if __name__ == "__main__":
    main()
