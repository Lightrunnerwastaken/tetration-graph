"""b14: Lineares Teilproblem (LP): Π_N f_K = J_N⁻¹ T_N[e^x + (f_K'∘f_K) f_K] (J_N an f_K). Vergleich mit p_N und f_K."""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
Kc = [mp.mpf(c) for c in K]
kreis = [L * mp.expjpi(2 * mp.mpf(j) / 96) for j in range(96)]   # |z| = |L|, j = 0: z = L
for N in (30, 60, 90, 120, 150, 180):

    pw, g, Ji = baue(N, J=N)
    # rechte Seite: T_N[e^x + g·f_K]
    ex_ser = [arb(1)]
    for k in range(1, N + 1):
        ex_ser.append(ex_ser[-1] / k)
    gf = (g * arb_poly([arb(c) for c in K[:N + 1]])).truncate(N + 1)
    b = arb_mat([[ex_ser[k] + gf[k]] for k in range(N + 1)])
    x = Ji * b
    PiF = [mp.mpf(x[k, 0].mid().str(60, radius=False)) for k in range(N + 1)]
    d_lin = max(abs(horner(PiF, z) - horner(Kc, z)) for z in kreis)
    d_nl = max(abs(horner(P[N], z) - horner(Kc, z)) for z in kreis)
    d_lp = max(abs(horner(PiF, z) - horner(P[N], z)) for z in kreis)
    d1 = max(abs(horner(PiF, z) - horner(Kc, z)) for z in [mp.expjpi(2 * mp.mpf(j) / 48) for j in range(48)])
    print(f"N={N}: max_|z|=|L| |Π_N f_K − f_K| = {mp.nstr(d_lin,3)}, |p_N − f_K| = {mp.nstr(d_nl,3)}, |Π_N f_K − p_N| = {mp.nstr(d_lp,3)};"
          f"  max_|z|=1 |Π_N f_K − f_K| = {mp.nstr(d1,3)}", flush=True)
