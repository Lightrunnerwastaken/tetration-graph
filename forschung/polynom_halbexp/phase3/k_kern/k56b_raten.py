"""k56b — Auswertung k56: Raten (1/Δn)·Δln|R_n| je Punkt des Γout-Bogens gegen die Bernstein-Walsh-Schranke ln(|f(y)|/ρ),
dazu das Koeffizientenprofil (1/n)ln|ũ_k| der Eck-Pol-Antwort (ρ-skaliert) bei n = 200, 400.
Aufruf: python k56b_raten.py   (liest k56_100.log … k56_400.log)
"""
import math
import os
import sys

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
HIER = os.path.dirname(os.path.abspath(__file__))


def lies(n):
    zeilen = open(os.path.join(HIER, "k56_%d.log" % n), encoding="utf-8").read().splitlines()[2:]
    out = []
    for z in zeilen:
        t = z.split()
        if len(t) < 8 or t[0] == "fertig":
            continue
        out.append((float(t[0]), float(t[1]), float(t[2]), float(t[3])))
    return out


def koeff(n):
    bits = int(1.3 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); rho = abs(L)
    J = acb_mat(J_tilde("kubisch", n, bits=bits)); ctx.prec = bits
    bv = acb_mat(n + 1, 1)
    for k in range(n + 1):
        v = rho ** k / L ** (k + 1)
        bv[k, 0] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
    X = J.solve(bv, algorithm="approx")
    return [float(mp.log(abs(mp.mpc(mp.mpf(X[k, 0].real.mid().str(nd, radius=False)), mp.mpf(X[k, 0].imag.mid().str(nd, radius=False)))) + mp.mpf(10) ** (-nd + 5))) for k in range(n + 1)]


def main():
    ns = [100, 200, 300, 400]
    d = {n: lies(n) for n in ns}
    print("Raten je Grad auf dem Γout-Bogen: Steigung von ln|R_n| zwischen n = 100→200, 200→300, 300→400; BW = ln(|f|/ρ); Verhältnis = Rate(300→400)/BW")
    for i, (ph, dy, fr, _) in enumerate(d[100]):
        lr = [math.log(d[n][i][3]) for n in ns]
        sl = [(lr[j + 1] - lr[j]) / 100 for j in range(3)]
        bw = math.log(fr)
        print("  φ = %7.2f°  |y−L| = %.4f  Raten %+.4f %+.4f %+.4f   BW %.4f   Verh. %s" % (ph, dy, sl[0], sl[1], sl[2], bw, ("%.2f" % (sl[2] / bw)) if bw > 1e-3 else "—"))
    for n in (200, 400):
        lk = koeff(n)
        prof = ["%.2f:%+.3f" % (t, lk[int(round(t * n))] / n) for t in (0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0)]
        print("n = %d: (1/n)ln|ũ_k| bei k/n = " % n + "  ".join(prof))


if __name__ == "__main__":
    main()
