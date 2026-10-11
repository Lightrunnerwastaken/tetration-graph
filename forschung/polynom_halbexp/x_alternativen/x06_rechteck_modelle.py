"""X06 — Rechteckige Abschnitte T_{N+d}J_∞|_{P_N} für polynomiale Modelle (kubisch wie D/O/W, radiales Quadrat wie M).

σ_min(T_{N+d}J_∞|_{P_N}) in ℓ²_ρ für eine Liste von d (inkrementell: G_d = G_{d−1} + b_dᵀb_d; nur bei den gewünschten d invertieren).
Aufruf: python x06_rechteck_modelle.py modell N d1,d2,…   (modell = kubisch | radial)
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
modell, N = sys.argv[1], int(sys.argv[2])
DL = sorted(int(x) for x in sys.argv[3].split(","))
ctx.prec = int(2.6 * N) + 500
t0 = time.time()

if modell == "kubisch":
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
elif modell == "radial":
    # f(z) = z + ε(z − e^{iθ})(z − e^{−iθ}) = z + ε(z² − 2cosθ z + 1), Exp-Lage (Strang M): θ = 76.618°, ε = 0.31382
    th = arb("76.618") * arb.pi() / 180
    eps = arb("0.31382")
    c = [eps, 1 - 2 * eps * th.cos(), eps]
    rho = arb(1)
else:
    raise SystemExit("modell?")
f = arb_poly(c)
fp = f.derivative()
m = fp(f)
deg = len(c) - 1
M = N + max(DL)
rk = [arb(1)]
for k in range(1, M + 1):
    rk.append(rk[-1] * rho)
A = arb_mat(M + 1, N + 1)
pot = arb_poly([1])
for j in range(N + 1):
    col = pot + m * arb_poly([0] * j + [1])
    inv_rj = 1 / rk[j]
    for k in range(min(col.degree(), M) + 1):
        A[k, j] = col[k] * rk[k] * inv_rj
    pot = (pot * f).truncate(M + 1)
print(f"{modell} N={N}: f = {[float(x.mid()) for x in c]}, Matrix ({time.time() - t0:.0f}s), prec {ctx.prec}", flush=True)
G = arb_mat(N + 1, N + 1)
JN = arb_mat(N + 1, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        JN[i, j] = A[i, j]
G = JN.transpose() * JN
erg = []
for d in range(0, max(DL) + 1):
    if d > 0:
        row = arb_mat(1, N + 1)
        for j in range(N + 1):
            row[0, j] = A[N + d, j]
        G = G + row.transpose() * row
    if d not in DL:
        continue
    Gi = G.inv()
    Gf = np.array([[float(Gi[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
    ew = np.sort(np.linalg.eigvalsh((Gf + Gf.T) / 2))[::-1]
    s = 1 / np.sqrt(ew[:4])
    print(f"   d={d}: σ_min = " + ", ".join(f"{x:.5f}" for x in s) + f"   [{time.time() - t0:.0f}s]", flush=True)
