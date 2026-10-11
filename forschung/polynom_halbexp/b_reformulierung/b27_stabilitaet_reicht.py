"""b27 (Runde 2): Satz B-Stab numerisch nachvollziehen.
Π_N f_K = J_N⁻¹ T_N[eˣ + g·f_K] exakt (g = f_K'∘f_K = eˣ/f_K').  Norm A_s: |c|_s² = Σ |c_k|² |L|^{2k} (1+k)^{2s}.
(i)  S_N(s) := |J_N⁻¹|_{A_s→A_s} (Spektralnorm von D_w J⁻¹ D_w⁻¹, D_w = diag(|L|^k (1+k)^s));
(ii) |T_N J_∞ f_K|_s  (beschränkt, da eˣ + g f_K ∈ A_s für s < 4.95);  (iii) |Π_N f_K|_s ≤ S_N · |T_N J_∞ f_K|_s."""
import json, sys
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
aL = abs(L)


def norm_s(v, s):
    return mp.sqrt(mp.fsum(abs(v[k]) ** 2 * aL ** (2 * k) * (1 + k) ** (2 * s) for k in range(len(v))))


# volle Reihe eˣ + g f_K bis Grad 400 (Kneser-Koeffizienten bis 1000), für die Norm der unabgeschnittenen Funktion
M_ = 400
Kc = [mp.mpf(c) for c in K[:M_ + 2]]
ex_ser = [mp.mpf(1)]
for k in range(1, M_ + 1):
    ex_ser.append(ex_ser[-1] / k)
fp = [Kc[k] * k for k in range(1, M_ + 2)]
gk = []
for n in range(M_ + 1):
    gk.append((ex_ser[n] - mp.fsum(fp[k] * gk[n - k] for k in range(1, n + 1))) / fp[0])
v_full = [ex_ser[n] + mp.fsum(gk[i] * Kc[n - i] for i in range(n + 1)) for n in range(M_ + 1)]
for s in (1, 2, 3):
    print(f"s={s}: |eˣ + g f_K|_s (Grad ≤ 400) = {mp.nstr(norm_s(v_full, s), 6)},  |f_K|_s = {mp.nstr(norm_s(Kc[:M_+1], s), 6)},  |g|_s = {mp.nstr(norm_s(gk, s), 6)}")

for N in (40, 80, 120, 160):
    pw, g, Ji = baue(N, J=N)
    Jinv = [[mp.mpf(Ji[i, j].mid().str(60, radius=False)) for j in range(N + 1)] for i in range(N + 1)]
    rhs = v_full[:N + 1]
    PiF = [mp.fsum(Jinv[i][j] * rhs[j] for j in range(N + 1)) for i in range(N + 1)]
    for s in (1, 2, 3):
        w = [aL ** k * (1 + k) ** s for k in range(N + 1)]
        A = np.array([[float(w[i] * Jinv[i][j] / w[j]) for j in range(N + 1)] for i in range(N + 1)])
        S = np.linalg.norm(A, 2)
        r = norm_s(rhs, s); p = norm_s(PiF, s)
        print(f"N={N} s={s}: S_N = {S:.4f}   |T_N J_∞f_K|_s = {mp.nstr(r, 6)}   |Π_N f_K|_s = {mp.nstr(p, 6)}   Quotient = {mp.nstr(p / r, 4)}", flush=True)
