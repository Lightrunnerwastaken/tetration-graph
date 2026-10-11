"""W42 — Nullstellenverteilung der Eckmode y₁ (= σ₁V₁) gegen das Leitordnungsmodell (Test der globalen Potentialstruktur, Programm §22.6(A)).
Vorhersage: Φ_y = H auf {H > 0} (Ω₁ ∪ Linse), = 0 (Sockel) sonst ⇒ alle N Nullstellen liegen asymptotisch auf {H = 0} = Buchtrand B ∪ Linsenrand,
mit Dichte (N/2π)|∇H| (Fluss von H aus {H > 0}; Gesamtmasse 1 wegen H ~ −ln|t| bei 0).
Aufruf: python w42_nullstellen.py N1,N2
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, arb, arb_mat, ctx

import w_kern as W
from w40_lade import lade

sys.stdout.reconfigure(encoding="utf-8")
rho, L = W.RHO, W.L
H, Lam, rB, pg = lade()
coef, _ = W.kubik_koeff()
res = {}
for N in [int(v) for v in sys.argv[1].split(",")]:
    n = N + 1
    ctx.prec = int(1.5 * N) + 400
    ctx.cap = n + 10
    t0 = time.time()
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    U = d["U"]
    a = W.kubik_a(N, coef)
    P = W.potenzen(a, N)
    J = W.J_mat(a, N, P)
    Dk = [arb(rho) ** k for k in range(n)]
    Z = J.transpose().solve(arb_mat(n, 1, [Dk[k] * arb(float(U[k, 0])) for k in range(n)]), algorithm="approx")
    p = acb_poly([acb(Z[k, 0].mid()) for k in range(n)])          # Y(x) = Σ Z_k x^{−k} = p(1/x)
    ws = p.roots(maxprec=4 * ctx.prec)
    xs = np.array([complex(1 / w) for w in ws]) / rho             # Nullstellen in t = x/ρ
    Hz = H(xs)
    inD = np.abs(xs) < 1
    print(f"N={N}: {len(xs)} Nullstellen ({time.time()-t0:.0f}s); innen {inD.sum()}, aussen {(~inD).sum()}; "
          f"|H| an den Nullstellen: Median {np.median(np.abs(Hz)):.4f}, 90%-Quantil {np.quantile(np.abs(Hz), 0.9):.4f}", flush=True)
    res[N] = {"re": xs.real.tolist(), "im": xs.imag.tolist()}
    json.dump(res, open(W.HIER + "/w42_nullstellen.json", "w"))
