"""W41 — Satz W1'' mit MODELL-definierter Linse: Λ := Γin-Profil des Leitordnungsmodells (w40), Linse ℬ_x = {|t| ≥ ρ: P_ext[Λ](t) > x·ln(|t|/ρ)}.
Teil 1 (ohne N): Linsenbreite, M_ℬ = max g_W̄ auf ℬ_x, ε_Γ-Schranke 1 − M_ℬ/d, für x = 0.7 und 0.9.
Teil 2 (arb, je N): Hypothese (H*) auf γ* = ∂ℬ_x ∪ h(A-Bogen): ln(Σ_i|Y_i|²)^{1/2} auf ∂ℬ_x (Eckmoden σ_i > 1.40), dazu T auf ∂ℬ_x.
Aufruf: python w41_hstern_modell.py [N1,N2,...]   (ohne Argument nur Teil 1)
"""
import json
import sys
import time

import numpy as np

import w_kern as W
import w24_rand_genau as M
from w40_lade import lade

sys.stdout.reconfigure(encoding="utf-8")
rho, L = W.RHO, W.L
argL = np.angle(L)
H, Lam, _, _ = lade()
PH = np.linspace(-np.pi, np.pi, 8001)[:-1]
tt = np.degrees(np.abs(PH)) - np.degrees(argL)
LV = np.zeros_like(PH)
LV[tt > 0] = np.maximum(Lam(tt[tt > 0]), 0.0)


def P_ext(t):
    r = abs(t) / rho
    ker = (r ** 2 - 1) / (r ** 2 - 2 * r * np.cos(np.angle(t) - PH) + 1)
    return np.mean(LV * ker)


def s_lens(ph, x):
    lo, hi = 0.0, 0.3
    if P_ext(rho * np.exp(1e-7 + 1j * ph)) <= x * 1e-7:
        return 0.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if P_ext(rho * np.exp(mid + 1j * ph)) > x * mid:
            lo = mid
        else:
            hi = mid
    return hi


XS = [0.7, 0.9]
TAUS = [0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 103.38]
pts = {}
for x in XS:
    Mb = 0.0
    for t in np.linspace(0.1, 180 - np.degrees(argL), 240):
        ph = argL + np.deg2rad(t)
        Mb = max(Mb, M.green(rho * np.exp(s_lens(ph, x) + 1j * ph)))
    print(f"x = {x}: Linsenbreite s(τ) = " + " ".join(f"{t}°:{s_lens(argL+np.deg2rad(t), x):.5f}" for t in TAUS)
          + f";  M_ℬ = {Mb:.5f} ⇒ ε_Γ ≥ {1-Mb/0.05:.3f} (d = 0.05), {1-Mb/0.1:.3f} (d = 0.1)")
    for t in TAUS:
        ph = argL + np.deg2rad(t)
        pts[f"x{x}_Z{t}"] = rho * np.exp(s_lens(ph, x) + 1j * ph)
print("Schlitz: P_ext[Λ](L(1+u)) gegen x ln|t|/ρ: " + " ".join(f"u={u}: {P_ext(L*(1+u)):.2e}/{0.9*np.log(1+u):.2e}" for u in (1e-3, 1e-2, 0.05)))
if len(sys.argv) < 2:
    sys.exit()

from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

coef, _ = W.kubik_koeff()
res = {}
for N in [int(v) for v in sys.argv[1].split(",")]:
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
    Dk = [arb(rho) ** k for k in range(n)]
    Zs = J.transpose().solve(arb_mat(n, r, [Dk[k] * arb(float(U[k, i])) for k in range(n) for i in range(r)]), algorithm="approx")
    x0 = acb_poly([acb(c) for c in a])
    mser = acb_poly([0])
    for j in range(1, n):
        mser += j * acb(a[j]) * acb_poly(P[j - 1].coeffs())
    mser = mser.truncate(n)
    uD = [[arb(float(U[k, i])) * Dk[k] for k in range(n)] for i in range(r)]

    def Yabs(xx):
        w = 1 / acb(xx.real, xx.imag)
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
    out = {}
    for nm, z in pts.items():
        out["Y_" + nm] = Yabs(z)
        out["T_" + nm] = Tabs(z)
    res[N] = out
    json.dump(res, open(W.HIER + "/w41_hstern.json", "w"), indent=1)
    print(f"N={N} (r={r}, {time.time()-t0:.0f}s): max ln Y auf ∂ℬ_0.9 = {max(v for k, v in out.items() if k.startswith('Y_x0.9')):.2f}", flush=True)
