"""W53 — Nullstellen der exakten Interpolanten Q_N (Q_N∘f + m·Q_N = O(x^N), monisch) und A_N (adjungiert), im Vergleich mit |x| = ρ, L und ∂Ω_y."""
import json
import sys

import numpy as np
from flint import acb, acb_poly, arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho, L = W.RHO, W.L
out = {}
for N in [int(v) for v in sys.argv[1].split(",")]:
    n = N + 1
    ctx.prec = int(2.0 * N) + 400
    a_ = W.kubik_a(N, coef)
    J = W.J_mat(a_, N)
    A1 = arb_mat(N, N, [J[i, j] for i in range(N) for j in range(N)])
    qp = A1.solve(arb_mat(N, 1, [-J[i, N] for i in range(N)]), algorithm="approx")
    ap = A1.transpose().solve(arb_mat(N, 1, [-J[N, j] for j in range(N)]), algorithm="approx")
    Q = acb_poly([acb(qp[i, 0].mid()) for i in range(N)] + [acb(1)])            # in x
    A = acb_poly([acb(ap[i, 0].mid()) for i in range(N)] + [acb(1)])            # in w = 1/x
    zq = np.array([complex(z) for z in Q.roots(maxprec=4 * ctx.prec)]) / rho
    za = np.array([1 / complex(z) for z in A.roots(maxprec=4 * ctx.prec)]) / rho
    for nm, z in (("Q_N (vorwärts)", zq), ("A_N (adjungiert)", za)):
        r = np.abs(z)
        print(f"N={N} {nm}: |t| Quantile 5/25/50/75/95 % = " + " ".join(f"{np.quantile(r, p):.3f}" for p in (0.05, 0.25, 0.5, 0.75, 0.95))
              + f"; Anteil |t|<1: {np.mean(r < 1):.3f}; Anteil |arg t| < arg L: {np.mean(np.abs(np.angle(z)) < np.angle(L)):.3f}", flush=True)
    out[N] = {"q_re": zq.real.tolist(), "q_im": zq.imag.tolist(), "a_re": za.real.tolist(), "a_im": za.imag.tolist()}
    json.dump(out, open(W.HIER + "/w53_nullstellen.json", "w"))
