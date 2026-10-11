"""b12: Spektrum von DΨ_N = −S⁻¹Γ (Linearisierung der Jet-Fixpunktform) an p_N: Eigenwerte nahe −1? Spektralradius?"""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 400
mp.mp.dps = 120
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
RAW = json.load(open(G + "polynome_bis_194.json"))
for N in (10, 20, 30, 40):
    a = [arb(s) for s in RAW[str(N)]]
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(N + 1))
    fpf = arb_poly([0])
    for m in range(1, N + 1):
        fpf += m * a[m] * P[m - 1]
    h = arb_poly([0] + a[1:])
    Hp = [arb_poly([1])]
    for _ in range(N):
        Hp.append((Hp[-1] * h).truncate(N + 1))
    H = mp.matrix(N + 1, N + 1); S = mp.matrix(N + 1, N + 1); Mm = mp.matrix(N + 1, N + 1)
    a0 = mp.mpf(RAW[str(N)][0])
    for k in range(N + 1):
        for i in range(N + 1):
            H[k, i] = mp.mpf(Hp[i][k].mid().str(110, radius=False))
            S[k, i] = mp.binomial(i, k) * a0 ** (i - k) if i >= k else 0
            Mm[k, i] = mp.mpf(fpf[k - i].mid().str(110, radius=False)) if i <= k else 0
    Gam = H ** -1 * Mm
    D = -(S ** -1) * Gam
    ev = mp.eig(D, left=False, right=False)
    ev = sorted(ev, key=lambda x: abs(x + 1))
    rho = max(abs(x) for x in ev)
    print(f"N={N}: Spektralradius von DΨ = {mp.nstr(rho, 5)}; Eigenwerte nächst −1: {[mp.nstr(x, 6) for x in ev[:4]]}")
