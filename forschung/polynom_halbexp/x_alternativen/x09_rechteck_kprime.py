"""X09 — Eckfreie Projektion Π̂ und Eckteil c auf ∂K' (K' = D̄ ∪ f(D̄)), kubisches Modell, grosse N; Daten der d-dim. Eckmatrix für W.

Billige Formel (bewiesen, BERICHT §3.5): mit Z := B J_N⁻¹ (d×(N+1), B = Zeilen N+1…N+d von J_∞|_{P_N}) ist
  Π̂g = Π_N g − c(g),   c(g) = J_N⁻¹ Zᵀ (I_d + Z Zᵀ)⁻¹ o(g),   o(g) := [J_∞(Π_N g − g)]_{N+1…N+d} = BΠ_N g − v(g)
(Kleinste-Quadrate mit d Zusatzzeilen; nur EINE LU-Zerlegung von J_N, alles reell ⇒ Zᵀ statt Z*). Alles ℓ²_ρ-skaliert.
Gemessen für η auf ∂K' (Γin und f(Γout), obere Hälfte, dicht bei L), am Schlitz und an Fernpunkten:
  R(η) = |Π_N k_η − T_N k_η|, R̂(η) = |Π̂k_η − T_N k_η|, Rc(η) = |c(k_η)|   (ℓ¹_ρ, wie C).
Ausserdem: Singulärwerte von Z (Eckwachstum), npz mit Z, J_N⁻¹Zᵀ für W.
Aufruf: python x09_rechteck_kprime.py N d
"""
import json
import os
import sys
import time

import numpy as np
from flint import acb, acb_mat, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
N, d = int(sys.argv[1]), int(sys.argv[2])
ctx.prec = int(float(sys.argv[3]) * N) + 400 if len(sys.argv) > 3 else int(1.6 * N) + 400
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
print(f"N={N} d={d} prec={ctx.prec}: Matrix ({time.time() - t0:.0f}s)", flush=True)

# --- Punkte auf ∂K'
fz = lambda z: ((c[3] * z + c[2]) * z + c[1]) * z + c[0]
argL = float(arb.atan2(Li, Lr).mid())
rho_f = float(rho.mid())
punkte = []      # (Name, η als acb, Lagecode)
taus = [1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 0.8, 1.2, 1.8]
for t in taus:   # Γin-Seite: φ = arg L + t (rad)
    ph = arb(str(argL + t))
    punkte.append((f"Γin τ=+{t}", acb(rho) * (acb(0, 1) * ph).exp(), "in"))
punkte.append(("Γin −ρ", acb(-rho), "in"))
for t in taus[:-3] + [0.4, 0.8, 1.2]:   # f(Γout)-Seite: x = ρe^{i(arg L − t)}, η = f(x)
    if argL - t <= 0:
        continue
    ph = arb(str(argL - t))
    x = acb(rho) * (acb(0, 1) * ph).exp()
    punkte.append((f"f(Γout) τ=−{t}", fz(x), "out"))
punkte.append(("f(Γout) f(ρ)", acb(fz(rho)), "out"))
alpha_s = arb("0.787172")
L = acb(Lr, Li)
for v in (0.3, 1.0, 3.0):
    punkte.append((f"Schlitz v={v}", L * (1 + arb(v) * arb(N) ** (-alpha_s)), "schlitz"))
punkte.append(("Ecke L", L, "ecke"))

# --- rechte Seiten: T_{N+d} J_∞ k_η (skaliert), J_∞k = 1/(η − f) + m/(η − x)
fa = [acb(x) for x in c]
ma = [acb(m[k]) for k in range(m.degree() + 1)]


def taylor_inv_f(eta, n):
    q = [eta - fa[0]] + [-x for x in fa[1:]]
    inv = [acb(0)] * (n + 1)
    inv[0] = 1 / q[0]
    for k in range(1, n + 1):
        s = acb(0)
        for i in range(1, min(k, 3) + 1):
            s += q[i] * inv[k - i]
        inv[k] = -s / q[0]
    return inv


npk = len(punkte)
W = acb_mat(N + 1, npk)
V = acb_mat(d, npk)
KZ = []
for p, (name, eta, _) in enumerate(punkte):
    a1 = taylor_inv_f(eta, M)
    kz = [1 / eta ** (k + 1) for k in range(M + 1)]
    for k in range(M + 1):
        mk = acb(0)
        for i in range(min(k, len(ma) - 1) + 1):
            mk += ma[i] * kz[k - i]
        val = (a1[k] + mk) * rk[k]
        if k <= N:
            W[k, p] = val
        else:
            V[k - N - 1, p] = val
    KZ.append([kz[k] * rk[k] for k in range(N + 1)])
JNc = acb_mat(JN)
X = JNc.solve(W)                                     # Π_N k_η (skaliert)
Zt = acb_mat(JN.transpose()).solve(acb_mat(B.transpose()))    # Zᵀ = J_N^{−T} Bᵀ, (N+1)×d
Z = Zt.transpose()
Y = JNc.solve(Zt)                                    # J_N⁻¹ Zᵀ
G = acb_mat(d, d)
ZZ = Z * Zt
for i in range(d):
    for j in range(d):
        G[i, j] = ZZ[i, j] + (1 if i == j else 0)
O = acb_mat(B) * X - V                               # Überlauf-Residuum der quadratischen Lösung
Cc = Y * G.solve(O)                                  # Eckteil c
print(f"   Lösungen ({time.time() - t0:.0f}s)", flush=True)


def l1(col):
    return sum(abs(complex(v.mid())) for v in col)


zs = np.linalg.svd(np.array([[complex(Z[i, j].mid()) for j in range(N + 1)] for i in range(d)]), compute_uv=False)
print("   σ(Z = B J_N⁻¹): " + ", ".join(f"{s:.4g}" for s in zs), flush=True)
erg = []
for p, (name, eta, lage) in enumerate(punkte):
    y = [X[k, p] - KZ[p][k] for k in range(N + 1)]
    yc = [Cc[k, p] for k in range(N + 1)]
    yh = [y[k] - yc[k] for k in range(N + 1)]
    Rq, Rh, Rc = l1(y), l1(yh), l1(yc)
    erg.append({"name": name, "lage": lage, "eta": [float(eta.real.mid()), float(eta.imag.mid())], "R": Rq, "Rhat": Rh, "Rc": Rc})
    print(f"   {name:18s} η={complex(eta.mid()):.5f}: R = {Rq:.5g}  R̂ = {Rh:.5g}  Rc = {Rc:.5g}", flush=True)
rad = max(float(X[k, p].rad()) for k in range(0, N + 1, 37) for p in range(npk))
print(f"   max Radius (Stichprobe) {rad:.2e}  [{time.time() - t0:.0f}s]", flush=True)
json.dump({"N": N, "d": d, "sigma_Z": [float(s) for s in zs], "punkte": erg},
          open(os.path.join(HIER, f"x09_N{N}_d{d}.json"), "w"), indent=1)
print("fertig", flush=True)
np.savez_compressed(os.path.join(HIER, f"x09_N{N}_d{d}.npz"),
                    Z=np.array([[complex(Z[i, j].mid()) for j in range(N + 1)] for i in range(d)]),
                    Y=np.array([[complex(Y[i, j].mid()) for j in range(d)] for i in range(N + 1)]))
