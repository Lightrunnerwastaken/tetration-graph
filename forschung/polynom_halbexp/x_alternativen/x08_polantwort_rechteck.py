"""X08 — Pol-Antworten quadratisch gegen rechteckig (kubisches Modell): hat die Kleinste-Quadrate-Projektion Π̂ = R⁺T_{N+d}J_∞ die Spitze?

y_η := Π k_η − T_N k_η, k_η = 1/(η − x), R(η) = Σ_j ρ^j |y_j| (Cs Grösse R_N, ℓ¹_ρ).  Π = J_N⁻¹T_N J_∞ (quadratisch) bzw. Π̂ = R⁺T_{N+d}J_∞.
J_∞k_η = 1/(η − f(x)) + m(x)/(η − x), m = f'∘f.  Punkte: −ρ (ferner D-Bogen, Ds/Cs Spitze), D-Bogen t = 0.1 und 0.5, L (Ecke),
Schlitz L(1 + N^{−0.883}), Scheitel a₀+R₀ von D(f(0), |L − f(0)|).
Aufruf: python x08_polantwort_rechteck.py N d
"""
import sys
import time

import mpmath as mp
from flint import acb, acb_mat, acb_poly, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
N, d = int(sys.argv[1]), int(sys.argv[2])
ctx.prec = int(2.6 * N) + 500
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
M = N + d
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
JNi = JN.inv()
Gi = (A.transpose() * A).inv()
Rp = Gi * A.transpose()                     # (N+1)×(M+1)
print(f"N={N}, d={d}: Inverse fertig ({time.time() - t0:.0f}s), prec {ctx.prec}", flush=True)

L = acb(Lr, Li)
a0 = c[0]
R0 = (L - a0).abs_upper()
fa = acb_poly([acb(x) for x in c])
ma = acb_poly([acb(m[k]) for k in range(m.degree() + 1)])
punkte = {
    "−ρ": acb(-rho),
    "D t=0.5": acb(rho) * (acb(0, 1) * (arb.atan2(Li, Lr) + (arb.pi() - arb.atan2(Li, Lr)) * arb("0.5"))).exp(),
    "D t=0.1": acb(rho) * (acb(0, 1) * (arb.atan2(Li, Lr) + (arb.pi() - arb.atan2(Li, Lr)) * arb("0.1"))).exp(),
    "Ecke L": L,
    "Schlitz": L * (1 + arb(N) ** arb("-0.883")),
    "A-Scheitel": acb(a0 + (L - a0).abs_lower()),
}


def taylor_inv(poly, eta, n):
    """Taylor-Koeffizienten von 1/(eta − poly(x)) bis Grad n."""
    q = acb_poly([eta]) - poly
    inv = [acb(0)] * (n + 1)
    inv[0] = 1 / q[0]
    deg = q.degree()
    for k in range(1, n + 1):
        s = acb(0)
        for i in range(1, min(k, deg) + 1):
            s += q[i] * inv[k - i]
        inv[k] = -s / q[0]
    return inv


for name, eta in punkte.items():
    a1 = taylor_inv(fa, eta, M)                              # 1/(η − f(x))
    kz = [1 / eta ** (k + 1) for k in range(M + 1)]         # 1/(η − x)
    mk = (ma * acb_poly(kz)).truncate(M + 1)
    rhs = acb_mat(M + 1, 1)
    for k in range(M + 1):
        rhs[k, 0] = (a1[k] + mk[k]) * rk[k]
    rq = acb_mat(N + 1, 1)
    for k in range(N + 1):
        rq[k, 0] = rhs[k, 0]
    xq = acb_mat(JNi) * rq
    xr = acb_mat(Rp) * rhs
    Rq = sum(abs(complex(((xq[j, 0]) - kz[j] * rk[j]).mid())) for j in range(N + 1))
    Rr = sum(abs(complex(((xr[j, 0]) - kz[j] * rk[j]).mid())) for j in range(N + 1))
    print(f"   {name:10s}: R_quadr = {Rq:.5g}   R_rechteck = {Rr:.5g}   [{time.time() - t0:.0f}s]", flush=True)
