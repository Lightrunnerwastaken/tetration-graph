"""W52 — Raten der biorthogonalen Paare: ln‖Q_n‖, ln‖A_n‖ (ρ-skaliert), ln|h_n| für n = 100…600. Wer trägt c?"""
import sys

import numpy as np
from flint import arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho = W.RHO
NS = [int(v) for v in sys.argv[1].split(",")]
Nmax = max(NS)
ctx.prec = int(1.5 * Nmax) + 300
a_ = W.kubik_a(Nmax, coef)
J = W.J_mat(a_, Nmax)
print("  n    ln‖Q_n‖   ln‖A_n‖   ln|h_n|   ln(‖Q‖‖A‖/|h|)")
for N in NS:
    n = N + 1
    A1 = arb_mat(N, N, [J[i, j] for i in range(N) for j in range(N)])
    qp = A1.solve(arb_mat(N, 1, [-J[i, N] for i in range(N)]), algorithm="approx")
    ap = A1.transpose().solve(arb_mat(N, 1, [-J[N, j] for j in range(N)]), algorithm="approx")
    q = [qp[i, 0] for i in range(N)] + [arb(1)]
    a = [ap[i, 0] for i in range(N)] + [arb(1)]
    h = sum((a[i] * sum((J[i, j] * q[j] for j in range(n)), arb(0)) for i in range(n)), arb(0))
    lq = 0.5 * np.log(sum(float((q[k] * arb(rho) ** k).mid()) ** 2 for k in range(n)))
    la = 0.5 * np.log(sum(float((a[k] * arb(rho) ** (-k)).mid()) ** 2 for k in range(n)))
    lh = float(abs(h).mid().log())
    # Normierung: Q monisch im unskalierten Sinn; skaliert: Leitkoeffizient ρ^N bzw. ρ^{−N} — beide Normen relativ zum Leitkoeffizienten:
    print(f"{N:4d}  {lq - N*np.log(rho):8.3f}  {la + N*np.log(rho):8.3f}  {lh:8.3f}   {lq + la - lh:8.3f}", flush=True)
