"""W38 — Test der lokalen Ecktheorie: H (= Φ_T aussen, = Φ_y innen; asymptotische Raten) auf Kreisen um L.
Vorhersage (Verklebung W_A ~ h(W_A), Deckverschiebung 2πi − log μ in s = log(t − L)):  H ≈ Re(c₁ ζ^{κ₁}), κ₁ = 2πi/(2πi − log μ) = 1.1182 − 0.0317i,
Knotenlinien im Winkelabstand π/Re κ₁ = 161°, positiv auf der Γin-Seite. Daten: genaue Raten aus N = 400…1000 (kubisch).
Punkte: t = L + r e^{i(arg L + ψ)}, ψ ∈ [−180°, 180°] (ψ = 0 äussere Normale, +90° Γin-Tangente, −90° Γout-Tangente); innen (|t| < ρ) wird Φ_y
benutzt (ausser in der Bucht), aussen Φ_T.
Aufruf: python w38_eckmode_lokal.py N1,N2,...
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho, L = W.RHO, W.L
argL = np.angle(L)
RADII = [0.03, 0.06, 0.12]
PSIS = list(range(-180, 181, 15))
res = {}
for N in [int(x) for x in sys.argv[1].split(",")]:
    n = N + 1
    ctx.prec = int(1.5 * N) + 400
    ctx.cap = n + 10
    t0 = time.time()
    d = np.load(W.HIER + f"/w01_N{N}.npz")
    U, S = d["U"], d["S"]
    r = max(1, sum(1 for s in S[:U.shape[1]] if s > 1.40))
    a = W.kubik_a(N, coef)
    P = W.potenzen(a, N)
    J = W.J_mat(a, N, P)
    rr_ = arb(rho)
    Dk = [rr_ ** k for k in range(n)]
    Zs = J.transpose().solve(arb_mat(n, r, [Dk[k] * arb(float(U[k, i])) for k in range(n) for i in range(r)]), algorithm="approx")
    x0 = acb_poly([acb(c) for c in a])
    mser = acb_poly([0])
    for j in range(1, n):
        mser += j * acb(a[j]) * acb_poly(P[j - 1].coeffs())
    mser = mser.truncate(n)
    uD = [[arb(float(U[k, i])) * Dk[k] for k in range(n)] for i in range(r)]

    def Yabs(x):
        w = 1 / acb(x.real, x.imag)
        tot = 0.0
        for i in range(r):
            s = acb(0)
            for k in range(N, -1, -1):
                s = s * w + Zs[k, i]
            m = abs(s).mid()
            tot += float(m * m) if m > 0 else 0.0
        return 0.5 * np.log(tot)

    def Tabs(eta):
        et = acb(eta.real, eta.imag)
        s1 = list(acb_series((acb_poly([et]) - x0).coeffs(), prec=n).inv().coeffs())
        s1 += [acb(0)] * (n - len(s1))
        kz = [et ** (-(j + 1)) for j in range(n)]
        s2 = list((mser * acb_poly(kz)).truncate(n).coeffs())
        s2 += [acb(0)] * (n - len(s2))
        tot = 0.0
        for i in range(r):
            s = acb(0)
            for j in range(n):
                s += Zs[j, i] * (s1[j] + s2[j]) - uD[i][j] * kz[j]
            m = abs(s).mid()
            tot += float(m * m) if m > 0 else 0.0
        return 0.5 * np.log(tot)
    out = {"L_T": Tabs(L), "L_Y": Yabs(L)}
    for rad in RADII:
        for ps in PSIS:
            t = L + rad * np.exp(1j * (argL + np.deg2rad(ps)))
            out[f"{rad}_{ps}"] = Tabs(t) if abs(t) > rho else Yabs(t)
    res[N] = out
    json.dump(res, open(W.HIER + "/w38_lokal.json", "w"))
    print(f"N={N} ({time.time()-t0:.0f}s) ok", flush=True)
