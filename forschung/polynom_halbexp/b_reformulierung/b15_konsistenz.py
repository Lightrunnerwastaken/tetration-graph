"""b15: Lax-Schema für das lineare Teilproblem (LP): P_N − T_N f_K = J_N⁻¹ r_N mit r_N = T_N[Schwanz_N(f_K)∘f_K].
Norm l¹ mit Gewichten |L|^k (1+k)^s (Strang C: |J_N⁻¹| beschränkt für s = 4). Ist r_N klein (Konsistenz)?"""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
aL = abs(L)
def norm(v, s):
    return mp.fsum(abs(v[k]) * aL ** k * (1 + k) ** s for k in range(len(v)))
for N in (30, 60, 90, 120, 150):
    pw, g, Ji = baue(N, J=700)
    r = [arb(0)] * (N + 1)
    for j in range(N + 1, 701):
        cj = arb(K[j])
        for k in range(N + 1):
            r[k] += cj * pw[j][k]
    rv = [mp.mpf(x.mid().str(40, radius=False)) for x in r]
    x = Ji * arb_mat([[v] for v in r])
    xv = [mp.mpf(x[k, 0].mid().str(40, radius=False)) for k in range(N + 1)]
    for s in (0, 4):
        print(f"N={N} s={s}: |r_N| = {mp.nstr(norm(rv, s), 3)},  |J⁻¹ r_N| = {mp.nstr(norm(xv, s), 3)},  Quotient = {mp.nstr(norm(xv, s)/norm(rv, s), 3)}", flush=True)
