"""D01 — Induktion N -> N+1: Rest r_{N+1}, Schur-Komplement s, Newton-Vorhersage, Schritt c_N = p_{N+1} - pad(p_N).

Setzt man pad(p_N) (a_{N+1} = 0) in das System N+1 ein, so ist F_{N+1} = r_{N+1} e_{N+1}.
Newton-Schritt: delta = -r J_{N+1}^{-1} e_{N+1} = -r * w,  w = (-J_N^{-1} b, 1)/s,  s = Schur-Komplement.
"""
import json
import sys
import time

import mpmath as mp
from flint import acb, arb, arb_mat, arb_poly, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 160
BITS = 1200
p = K.loesungen(NMAX + 1, bits=BITS)
ctx.prec = BITS
L = K.L_wert()
absL = abs(L)
kn = [arb(s) for s in json.load(open(K.GEMEINSAM + "/kneser_halbexp_koeff_80.json"))]
RHOS = ["0.5", "1", "1.3", "1.37455701074370750324"]
NPKT = 256
kreis = [(acb(0, 1) * arb.pi() * 2 * j / NPKT).exp() for j in range(NPKT)]


def sup_kreis(v, R):
    pol = arb_poly(v)
    from flint import acb_poly
    q = acb_poly([acb(x) for x in v])
    return max(abs(q(R * z)).mid() for z in kreis)


aus = {}
t0 = time.time()
for N in range(1, NMAX + 1):
    a = p[N]
    a1 = p[N + 1]
    ap = a + [arb(0)]
    F, J = K.FJ(ap, N + 1)
    r = F[N + 1]
    resmax = max(abs(F[k].mid()) for k in range(N + 1))
    e = [arb(0)] * (N + 1) + [arb(1)]
    w = K.loese(J, e)
    s = 1 / w[N + 1]
    delta = [-r * x for x in w]
    c = [a1[k] - ap[k] for k in range(N + 2)]
    fehl = [delta[k] - c[k] for k in range(N + 2)]
    zeile = {"r": r, "s": s, "t_newton": -r / s, "a_top": a1[N + 1],
             "resmax_lo": resmax}
    for rho in RHOS:
        zeile["c_" + rho[:4]] = K.gew_norm(c, rho)
        zeile["w_" + rho[:4]] = K.gew_norm(w, rho)
        zeile["relfehl_" + rho[:4]] = K.gew_norm(fehl, rho) / K.gew_norm(c, rho)
    zeile["c_sup_L"] = sup_kreis(c, absL)
    zeile["w_sup_L"] = sup_kreis(w, absL)
    zeile["w_sup_1"] = sup_kreis(w, arb(1))
    if N + 1 <= 80:
        zeile["kn_top"] = kn[N + 1]
    aus[N] = {k: float(v.mid()) if hasattr(v, "mid") else float(v) for k, v in zeile.items()}
    if N % 10 == 0 or N < 6:
        z = aus[N]
        print(f"N={N:3d} r={z['r']:+.3e} s={z['s']:+.3e} t_N={z['t_newton']:+.3e} a_top={z['a_top']:+.3e} "
              f"c1.3={z['c_1.3']:.2e} c|L|={z['c_1.37']:.2e} csupL={z['c_sup_L']:.2e} w1.3={z['w_1.3']:.2e} "
              f"rel={z['relfehl_1.3']:.1e}  ({time.time()-t0:.0f}s)", flush=True)
json.dump(aus, open(K.HIER + "/d01_induktion.json", "w"), indent=0)
