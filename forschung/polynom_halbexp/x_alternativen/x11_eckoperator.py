"""X11 — Der kanonische Eckoperator 𝒞_N (Rang ≤ d) und seine Singulärdaten (kubisches Modell), Daten für W.

Mit Z := B J_N⁻¹ (d×(N+1)) ist für eine Eingabe w = T_{N+d}J_∞g = (w_top, w_over) ∈ ℂ^{N+1+d} (ℓ²_ρ-skaliert)
  c(g) = Π_N g − Π̂ g = 𝒞_N w,   𝒞_N := J_N⁻¹ Zᵀ (I_d + Z Zᵀ)⁻¹ [Z | −I_d]    (Rang ≤ d; BERICHT §3.5/3.6).
Ausgabe: Singulärwerte von 𝒞_N (= Eckverstärkung in der natürlichen Rechteck-Normierung), Vergleich mit σ_i(J_N⁻¹);
Lage der rechten Singulärvektoren (Eingangsseite, Zeilenindex k bzw. Überlaufzeilen) und der linken (Ausgänge, als Polynom auf |x| = ρ).
npz: Uc (Ausgänge), Vc (Eingänge), sc.
Aufruf: python x11_eckoperator.py N d
"""
import os
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
N, d = int(sys.argv[1]), int(sys.argv[2])
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
B = arb_mat(d, N + 1)
for i in range(d):
    for j in range(N + 1):
        B[i, j] = A[N + 1 + i, j]
Zt = JN.transpose().solve(B.transpose())          # (N+1)×d
Z = Zt.transpose()
Y = JN.solve(Zt)                                  # J_N⁻¹Zᵀ
G = Z * Zt
for i in range(d):
    G[i, i] = G[i, i] + 1
Gi = G.inv()
# rechter Faktor F := Gi·[Z | −I]  (d × (N+1+d))
F = arb_mat(d, N + 1 + d)
GZ = Gi * Z
for i in range(d):
    for j in range(N + 1):
        F[i, j] = GZ[i, j]
    for j in range(d):
        F[i, N + 1 + j] = -Gi[i, j]
# robuste Singulärwerte: σ² = Eigenwerte von (YᵀY)(FFᵀ) (d×d), in arb/mpmath
import mpmath as mp
mp.mp.dps = int(ctx.prec / 3.4)
YY = Y.transpose() * Y
FF = F * F.transpose()
P_ = YY * FF
Pm = mp.matrix([[mp.mpf(P_[i, j].mid().str(int(mp.mp.dps) - 5, radius=False)) for j in range(d)] for i in range(d)])
ewp = sorted([abs(mp.sqrt(abs(e))) for e in mp.eig(Pm, left=False, right=False)], reverse=True)
print("   σ(𝒞_N) robust (arb/mpmath): " + ", ".join(mp.nstr(e, 6) for e in ewp), flush=True)
Yf = np.array([[float(Y[i, j].mid()) for j in range(d)] for i in range(N + 1)])
Ff = np.array([[float(F[i, j].mid()) for j in range(N + 1 + d)] for i in range(d)])
# Rang-d-SVD: Y = Q1 R1, Fᵀ = Q2 R2  ⇒ 𝒞 = Q1 (R1 R2ᵀ) Q2ᵀ
Q1, R1 = np.linalg.qr(Yf)
Q2, R2m = np.linalg.qr(Ff.T)
Uc_s, sc, Vc_sT = np.linalg.svd(R1 @ R2m.T)
Uc = Q1 @ Uc_s
Vc = Q2 @ Vc_sT.T
print(f"N={N} d={d}: σ(𝒞_N) = " + ", ".join(f"{s:.5g}" for s in sc) + f"   [{time.time() - t0:.0f}s]", flush=True)
argL = np.angle(complex(float(Lr.mid()), float(Li.mid())))
phis = np.deg2rad(np.array([0, 40, 60, 70, 74, 76.6, 78, 80, 90, 120, 150, 180]))
for i in range(min(d, 4)):
    u, v = Uc[:, i], Vc[:, i]
    k = np.arange(N + 1)
    vin = v[: N + 1]
    over = np.linalg.norm(v[N + 1:])
    schw = np.sum(k * vin ** 2) / max(np.sum(vin ** 2), 1e-300)
    w = np.abs(vin) ** 2
    q = np.cumsum(w) / max(w.sum(), 1e-300)
    k05, k50, k95 = (int(np.searchsorted(q, t)) for t in (0.05, 0.5, 0.95))
    uz = np.array([abs(np.polyval(u[::-1], np.exp(1j * p))) for p in phis])
    print(f"   Mode {i + 1}: σ = {sc[i]:.5g}; Eingang: Anteil Überlaufzeilen {over:.3f}, Zeilen-Quantile 5/50/95 % bei k/N = "
          f"{k05 / N:.3f}/{k50 / N:.3f}/{k95 / N:.3f}; Ausgang |u(ρe^{{iφ}})| bei φ = 0/40/60/70/74/76.6/78/80/90/120/150/180°: "
          + " ".join(f"{x:.2g}" for x in uz), flush=True)
np.savez_compressed(os.path.join(HIER, f"x11_N{N}_d{d}.npz"), Uc=Uc, Vc=Vc, sc=sc)
