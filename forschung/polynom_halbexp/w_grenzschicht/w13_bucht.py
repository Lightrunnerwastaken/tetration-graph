"""W13 — Rand der „Bucht“ (Bereich innerhalb |t|<1 nahe Γout, in dem die adjungierte Eckmode auf dem Sockel bleibt), gegen N.
Genaue Auswertung (arb). Vergleich: gespiegelter A-Bogen (W̄): 1/(ã₀+R̃₀) = 0.743 bei 0°; K': 1/F̃(1) = 0.616.
Aufruf: python w13_bucht.py 400,600,800
"""
import sys, time
import numpy as np
from flint import acb, arb, arb_mat, ctx
import w_kern as W
sys.stdout.reconfigure(encoding="utf-8")
for N in [int(x) for x in sys.argv[1].split(",")]:
    n = N + 1; ctx.prec = int(1.5 * N) + 400
    coef, _ = W.kubik_koeff(); rho = W.RHO
    d = np.load(W.HIER + f"/w01_N{N}.npz"); u1 = d["U"][:, 0]
    t0 = time.time()
    J = W.J_mat(W.kubik_a(N, coef), N); rr = arb(rho); Dk = [rr ** k for k in range(n)]
    Z = J.transpose().solve(arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)]), algorithm="approx")
    y = [Z[k, 0] / Dk[k] for k in range(n)]
    def Y(t):
        w = 1 / acb(t.real, t.imag); s = acb(0)
        for v in y[::-1]: s = s * w + v
        return s
    def pot(t):
        m = abs(Y(t)).mid(); return float(m.log()) / N
    floor = pot(3.0 + 0j)
    out = [f"N={N} ({time.time()-t0:.0f}s) Sockel(1/N ln|y(3)|)={floor:.4f}"]
    for phi in (0, 30, 50):
        rs = np.linspace(0.6, 1.0, 81)
        ps = [pot(r * np.exp(1j * np.deg2rad(phi))) for r in rs]
        # Rand: grösstes r mit pot > floor + 0.01
        idx = [i for i, p in enumerate(ps) if p > floor + 0.01]
        rb = rs[max(idx)] if idx else float("nan")
        out.append(f"φ={phi}: Buchtrand r≈{rb:.3f}, pot(0.7)={ps[20]:.4f}")
    print(" | ".join(out), flush=True)
