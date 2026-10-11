"""X20 — Kanten-Aussage für (a″) mit KLEINEM d (Kneser, ℓ²_ρ):
  (1) Eckmoden von J_N: σ_i(J_N⁻¹) > 1.5, linke/rechte Singulärvektoren u_i, v_i (Blockpotenz auf J⁻ᵀJ⁻¹ in arb);
  (2) Überlauf-Profil der Eckmoden: |[J_∞v_i]_{N+j}|ρ^{N+j}, j = 1…D, und ‖T_{>N}J_∞v_i‖ (bis Zeile N+D);
  (3) σ_min(R_d), R_d = T_{N+d}J_∞|_{P_N}, für d = 0, 1, 2, …, D (über (R_dᵀR_d)⁻¹ = J⁻¹(I + Z_dᵀZ_d)⁻¹J⁻ᵀ, Z_d = B_dJ⁻¹, Woodbury);
  (4) σ_min(B_dV_r) (Überlauf der Eckmoden in den ersten d Zusatzzeilen) und der Winkel zwischen B_dV_r und B_d(Bulk).
Aufruf: python x20_kante.py N [D=24]
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
N = int(sys.argv[1])
D = int(sys.argv[2]) if len(sys.argv) > 2 else 24
M = N + D
ctx.prec = int(1.8 * N) + 400
mp.mp.dps = 60
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = arb(mp.nstr(abs(Lm), 55))
rk = [arb(1)]
for k in range(1, M + 1):
    rk.append(rk[-1] * rho)


def mullow(a, b, n):
    return (a * b).truncate(n)


def inv_series(a, n):
    g = arb_poly([1 / a[0]])
    k = 1
    while k < n:
        k = min(2 * k, n)
        g = mullow(g, arb_poly([2]) - mullow(a.truncate(k), g, k), k)
    return g


def mid(x):
    return mp.mpf(x.mid().str(50, radius=False))


t0 = time.time()
F = arb_poly([arb(s) for s in A_str[: M + 1]])
ex = arb_poly([1 / arb.fac_ui(k) for k in range(M + 1)])
m = mullow(ex, inv_series(F.derivative(), M + 1), M + 1)
A = arb_mat(M + 1, N + 1)
pw = arb_poly([1])
for j in range(N + 1):
    for k in range(M + 1):
        A[k, j] = (pw[k] + (m[k - j] if k >= j else 0)) * rk[k] / rk[j]
    pw = mullow(pw, F, M + 1)
J = arb_mat(N + 1, N + 1)
B = arb_mat(D, N + 1)
for i in range(N + 1):
    for j in range(N + 1):
        J[i, j] = A[i, j]
for i in range(D):
    for j in range(N + 1):
        B[i, j] = A[N + 1 + i, j]
Ji = J.inv()
G = Ji.transpose() * Ji
print(f"N = {N}, D = {D}: Matrizen ({time.time() - t0:.0f}s)", flush=True)

# (1) Eckmoden: Blockpotenz auf G = (J Jᵀ)⁻¹ → linke Singulärvektoren u_i von J
b = 6
rng = np.random.default_rng(3)
X = arb_mat([[arb(float(x)) for x in row] for row in rng.standard_normal((N + 1, b))])
for _ in range(50):
    Y = G * X
    S = Y.transpose() * Y
    Sm = mp.matrix([[mid(S[i, j]) for j in range(b)] for i in range(b)])
    Li = mp.cholesky(Sm) ** -1
    X = Y * arb_mat([[arb(mp.nstr(Li[j, i], 50)) for j in range(b)] for i in range(b)])
H = X.transpose() * G * X
Hm = mp.matrix([[mid(H[i, j]) for j in range(b)] for i in range(b)])
ew, ev = mp.eigsy(Hm)
order = sorted(range(b), key=lambda i: -ew[i])
sig = [mp.sqrt(ew[i]) for i in order]
U = X * arb_mat([[arb(mp.nstr(ev[r, i], 50)) for i in order] for r in range(b)])   # Spalten: u_i (orthonormal)
r = sum(1 for s in sig if s > 1.5)
print("σ_i(J_N⁻¹) = " + ", ".join(mp.nstr(s, 8) for s in sig) + f";  Eckmoden r = {r}", flush=True)
V = Ji * U                                                     # Spalten: σ_i v_i
for i in range(b):
    for k in range(N + 1):
        V[k, i] = V[k, i] / arb(mp.nstr(sig[i], 50))
BV = B * V
for i in range(max(r, 1) + 1):
    prof = [float(abs(BV[j, i]).mid()) for j in range(D)]
    tot = float(np.sqrt(sum(p * p for p in prof)))
    kum = np.sqrt(np.cumsum(np.array(prof) ** 2)) / tot
    print(f"  Mode {i + 1} (σ = {mp.nstr(sig[i], 6)}): ‖Überlauf bis N+{D}‖ = {tot:.4f}; |[J_∞v]_{{N+j}}|, j = 1…8: "
          + ", ".join(f"{p:.3f}" for p in prof[:8]) + f"; Anteil der Überlauf-Norm in j ≤ 1/2/3/6/12: "
          + ", ".join(f"{kum[j - 1]:.3f}" for j in (1, 2, 3, 6, 12) if j <= D), flush=True)

# (3) σ_min(R_d) genau: (R_dᵀR_d)⁻¹ = J⁻¹J⁻ᵀ − (J⁻¹Z_dᵀ)(I + Z_dZ_dᵀ)⁻¹(Z_dJ⁻ᵀ), dann float-eigvalsh
# (5) Winkel: Q_d = Orthonormalbasis der Zeilen von B_d (Gram-Schmidt in arb), C = Q_dᵀV_r; σ_min(C) > 0 ⇔ Eckraum wird von den
#     normierten Überlauf-Funktionalen gesehen (Kanten-Aussage in skalenfreier Form)
Z = B * Ji
JJ = Ji * Ji.transpose()
Vr = [[V[k, i] for k in range(N + 1)] for i in range(r)]
Q = []
res = {}
for d in range(0, D + 1):
    if d >= 1:
        q = [B[d - 1, k] for k in range(N + 1)]
        for _ in range(2):
            for qq in Q:
                c = sum((a * b for a, b in zip(q, qq)), arb(0))
                q = [a - c * b for a, b in zip(q, qq)]
        nq = sum((a * a for a in q), arb(0)).sqrt()
        Q.append([a / nq for a in q])
    if d not in (0, 1, 2, 3, 4, 6, 8, 12, 16, 24) or d > D:
        continue
    if d == 0:
        K = JJ
    else:
        Zd = arb_mat(d, N + 1)
        for a in range(d):
            for k in range(N + 1):
                Zd[a, k] = Z[a, k]
        Wd = (arb_mat([[arb(1) if a == c else arb(0) for c in range(d)] for a in range(d)]) + Zd * Zd.transpose()).inv()
        JZ = Ji * Zd.transpose()
        K = JJ - JZ * Wd * JZ.transpose()
    Kf = np.array([[float(K[a, c].mid()) for c in range(N + 1)] for a in range(N + 1)])
    lam = np.linalg.eigvalsh((Kf + Kf.T) / 2)[-1]
    res[d] = float(1 / np.sqrt(lam))
    if d >= 1 and r >= 1:
        C = np.array([[float(sum((a * b for a, b in zip(Q[x], Vr[y])), arb(0)).mid()) for y in range(r)] for x in range(d)])
        sv = np.linalg.svd(C, compute_uv=False)
        ang = sv.min() if d >= r else 0.0
    else:
        ang = 0.0
    print(f"  d = {d:2d}: σ_min(R_d) = {res[d]:.5f};  σ_min(Q_dᵀV_r) = {ang:.4f}  ({time.time() - t0:.0f}s)", flush=True)
(HIER / f"x20_N{N}.json").write_text(json.dumps({"N": N, "D": D, "sigma": [str(s) for s in sig], "r": r, "smin_Rd": res}), encoding="utf-8")
