"""X10 — „Überlauf“ der Eckmoden: wohin bildet J_∞ die Fast-Kernvektoren u_i von J_N ab? (kubisches Modell)

Für die führenden Singulärvektoren u_i von J_N⁻¹ (Ausgänge, ℓ²_ρ; hochgenau per Block-Inversiteration in arb) werden berechnet:
  σ_i,  |J_N u_i| = 1/σ_i (quadratischer Rest),  das Überlaufprofil |[J_∞u_i]_{N+j}| für j = 1…jmax (Zeilen jenseits der Abschneidung)
und die kumulierte Überlaufnorm ‖T_{(N,N+d]}J_∞u_i‖ für d = 1…jmax. Erklärt, wie viele Zusatzzeilen d(N) der Rechteck-Abschnitt braucht
(Lemma 3.1): u_i wird von R_{N,d} „gesehen“, sobald ‖T_{(N,N+d]}J_∞u_i‖ ≳ min_{Γin}|m|.
Aufruf: python x10_ueberlauf_eckmoden.py N [k Moden] [jmax]
"""
import os
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
K = int(sys.argv[2]) if len(sys.argv) > 2 else 4
JMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 40
ctx.prec = int(1.6 * N) + 400
HIER = os.path.dirname(os.path.abspath(__file__))
t0 = time.time()
Lr = arb("0.31813150520476413531265425158766451720351761387139986692237")
Li = arb("1.33723570143068940890116214319371061253950213846051241687")
R2 = Lr * Lr + Li * Li
rho = R2.sqrt()
mu_r = ((rho + Lr) / 2).sqrt()
mu_i = ((rho - Lr) / 2).sqrt()
rL_re = mu_i / (2 * Li)
rL_im = -(mu_r - 1) / (2 * Li)
beta = rL_im / Li
alpha = rL_re - beta * Lr
c = [R2 * alpha, 1 + R2 * beta - 2 * Lr * alpha, alpha - 2 * Lr * beta, beta]
f = arb_poly(c)
m = f.derivative()(f)
M = N + JMAX
rk = [arb(1)]
for k in range(1, M + 1):
    rk.append(rk[-1] * rho)
A = arb_mat(M + 1, N + 1)
pot = arb_poly([1])
for j in range(N + 1):
    col = pot + m * arb_poly([0] * j + [1])
    for k in range(min(col.degree(), M) + 1):
        A[k, j] = col[k] * rk[k] / rk[j]
    pot = (pot * f).truncate(M + 1)
JN = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        JN[i, j] = A[i, j]
Jinv = JN.inv()
JinvT = Jinv.transpose()
print(f"N={N}: Matrix ({time.time() - t0:.0f}s), prec {ctx.prec}", flush=True)

# Block-Inversiteration: U ← J⁻¹ J⁻ᵀ U (Eigenvektoren von J⁻¹J⁻ᵀ = Σσ²uuᵀ), Start: Einheitsvektoren an den oberen Graden (u lebt dort)
rng = np.random.default_rng(1)
U = arb_mat(N + 1, K)
for i in range(K):
    for j in range(N + 1):
        U[j, i] = arb(float(rng.standard_normal()) * (1.0 if j > 0.3 * N else 1e-3))


def orthonorm(Mx):
    """Gram–Schmidt in arb (Spalten)."""
    n, k = Mx.nrows(), Mx.ncols()
    Q = arb_mat(n, k)
    for i in range(k):
        v = [Mx[r, i] for r in range(n)]
        for p in range(i):
            s = sum((Q[r, p] * v[r] for r in range(n)), arb(0))
            v = [v[r] - s * Q[r, p] for r in range(n)]
        nv = sum((x * x for x in v), arb(0)).sqrt()
        for r in range(n):
            Q[r, i] = v[r] / nv
    return Q


for it in range(7):
    U = orthonorm(Jinv * (JinvT * U))
print(f"   Iteration fertig ({time.time() - t0:.0f}s)", flush=True)
# Rayleigh–Ritz im Unterraum: σ² aus U^T (J⁻¹J⁻ᵀ) U
W = Jinv * (JinvT * U)
H = U.transpose() * W
Hf = np.array([[float(H[i, j].mid()) for j in range(K)] for i in range(K)])
ew, ev = np.linalg.eigh((Hf + Hf.T) / 2)
order = np.argsort(ew)[::-1]
print("   σ_i (Ritz): " + ", ".join(f"{np.sqrt(ew[o]):.5g}" for o in order), flush=True)
AU = A * U                                            # J_∞ u (Zeilen 0…N+JMAX)
for idx, o in enumerate(order):
    coeff = ev[:, o]
    col = [sum((AU[r, i] * arb(float(coeff[i])) for i in range(K)), arb(0)) for r in range(M + 1)]
    sq = float(sum((x * x for x in col[: N + 1]), arb(0)).sqrt().mid())
    over = [abs(float(col[N + j].mid())) for j in range(1, JMAX + 1)]
    cum = np.sqrt(np.cumsum(np.array(over) ** 2))
    dsee = next((j + 1 for j in range(JMAX) if cum[j] > 0.5), None)
    print(f"   Mode {idx + 1}: σ = {np.sqrt(ew[o]):.5g}, |J_N u| = {sq:.3e}; Überlauf |[J_∞u]_N+j| j=1..10: "
          + " ".join(f"{x:.2e}" for x in over[:10]), flush=True)
    print(f"            kumuliert ‖T_(N,N+d]J_∞u‖ für d = 1,2,3,5,8,12,20,{JMAX}: "
          + " ".join(f"{cum[j - 1]:.3f}" for j in (1, 2, 3, 5, 8, 12, 20, JMAX) if j <= JMAX)
          + f";  erstes d mit > 0.5: {dsee}", flush=True)
print(f"   [{time.time() - t0:.0f}s]", flush=True)
