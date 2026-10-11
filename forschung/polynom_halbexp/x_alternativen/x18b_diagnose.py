"""X18b — Diagnose zu x18: Grösse der Störungsteile in ℓ²_ρ (N klein):
‖η‖, ‖T_N[η'∘f]‖, ‖T_N[(m'/f')η]‖ für die Zufallsvarianten; für p_N: ‖p_N − T_N f‖, ‖T_N[p_N'∘p_N] − m‖; für ft: ‖T_N[(T_Nf)'∘(T_Nf)] − m‖.
Dazu max_{|x|=ρ}|f_K(x)|/ρ. Aufruf: python x18b_diagnose.py N
"""
import sys

sys.argv = [sys.argv[0], sys.argv[1], "f"]
import x18_eck_familie as X  # noqa: E402  (führt nur Variante f aus)
from flint import arb_poly  # noqa: E402
import numpy as np  # noqa: E402

N = X.N


def norm(p):
    return float(np.sqrt(sum(float(abs(p[k] * X.rk[k]).mid()) ** 2 for k in range(min(p.degree(), N) + 1))))


def comp_T(g, F):
    out = arb_poly([0])
    for k in range(g.degree(), -1, -1):
        out = X.mullow(out, F, N + 1) + arb_poly([g[k]])
    return out


FN = arb_poly(X.F_K)
m = X.M_EX.truncate(N + 1)
print(f"‖m‖ = {norm(m):.3e}")
for v in ("zt", "zh", "zv"):
    e = arb_poly(X.stoerung(v))
    t1 = comp_T(e.derivative(), FN)
    t2 = X.mullow(X.MP_F, e, N + 1)
    print(f"{v}: ‖η‖ = {norm(e):.3e}, ‖T_N[η'∘f]‖ = {norm(t1):.3e}, ‖T_N[(m'/f')η]‖ = {norm(t2):.3e}")
P = arb_poly(X.F_P)
print(f"pN: ‖p_N − T_N f‖ = {norm(P - FN):.3e}, ‖T_N[p_N'∘p_N] − m‖ = {norm(comp_T(P.derivative(), P) - m):.3e}")
print(f"ft: ‖T_N[(T_Nf)'∘(T_Nf)] − m‖ = {norm(comp_T(FN.derivative(), FN) - m):.3e}")
