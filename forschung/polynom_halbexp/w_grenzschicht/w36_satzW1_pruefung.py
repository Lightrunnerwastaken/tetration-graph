"""W36 — Genaue Prüfung der Hypothesen von Satz W1 (kubisch, arb mit ctx.cap):
(V3) Y_i(h(η)) für η auf dem A-Bogen (Bucht, innerhalb D) — alle Eckmoden;
(Rand) Eckanteil T_i(η) = ⟨u_i, y_η⟩ auf dem Linsenrand Z⁺ = {ρe^{s(φ)+iφ}} (s aus w35, x₀ = 0.7) und auf dem A-Bogen;
(G1) Lage der Urbilder h(η), η ∈ Z⁺, relativ zur Linse (s-Wert gegen Linsenbreite).
Ausgabe je N: ln-Werte; danach (mit w37) asymptotische Raten. Aufruf: python w36_satzW1_pruefung.py N1,N2,...
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, ctx

import w_kern as W
import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
rho, L = W.RHO, W.L
argL = np.angle(L)
a0 = coef[0]
R0 = abs(L - a0)
psiL = np.angle(L - a0)
Fp = np.polynomial.Polynomial(coef)
P35 = np.load(W.HIER + "/w35_profil.npz")
taus35, Lam = P35["taus"], P35["Lam_env"]
X0 = 0.7
PH = np.linspace(-np.pi, np.pi, 8001)[:-1]


def Lam_of_phi(phi):
    t = np.degrees(np.abs(np.angle(np.exp(1j * phi)))) - np.degrees(argL)
    return np.where(t > 0, np.interp(t, taus35, Lam, left=0.0, right=Lam[-1]), 0.0)


LV = Lam_of_phi(PH)


def P_ext(t):
    r = abs(t) / rho
    ker = (r ** 2 - 1) / (r ** 2 - 2 * r * np.cos(np.angle(t) - PH) + 1)
    return np.mean(LV * ker)


def s_lens(ph):
    lo, hi = 0.0, 0.3
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if P_ext(rho * np.exp(mid + 1j * ph)) > X0 * mid:
            lo = mid
        else:
            hi = mid
    return hi


ZT = [0.25, 0.5, 1, 2, 4, 8, 15, 30, 60, 103.38]
ZP = [(f"Z{t}", rho * np.exp(s_lens(argL + np.deg2rad(t)) + 1j * (argL + np.deg2rad(t)))) for t in ZT]
AD = [0.005, 0.02, 0.05, 0.1, 0.2, 0.38, 0.6, 1.0, 1.5, 1.9]
AP = [(f"A{dd}", M.a_bogen(dd)) for dd in AD]
# (G1): Urbilder der Z⁺-Punkte
print("(G1) Urbild h(η) für η ∈ Z⁺: |h|/ρ − 1 gegen Linsenbreite e^{s(arg h)} − 1 dort")
for nm, z in ZP:
    rr = (Fp - z).roots()
    hh = rr[np.argmin(abs(rr - z))]
    sh = s_lens(np.angle(hh)) if abs(np.angle(hh)) > argL else 0.0
    print(f"   {nm}: |h|/ρ − 1 = {abs(hh)/rho - 1:+.5f}, Linse dort {np.exp(sh) - 1:.5f}, weitere Urbilder |x|/ρ = "
          + ", ".join(f"{abs(x)/rho:.2f}" for x in rr if abs(x - hh) > 1e-9))
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
    out = {}

    def Yabs(x):        # (Σ_i |Y_i(x)|²)^{1/2}, Y_i(x) = Σ_k Z_{k,i} x^{−k}
        w = 1 / acb(x.real, x.imag)
        tot = 0.0
        for i in range(r):
            s = acb(0)
            for k in range(N, -1, -1):
                s = s * w + Zs[k, i]
            m = abs(s).mid()
            tot += float(m * m) if m > 0 else 0.0
        return np.sqrt(tot)

    def Tabs(eta):      # (Σ_i |T_i(η)|²)^{1/2}
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
        return np.sqrt(tot)
    for nm, z in AP:
        rr = (Fp - z).roots()
        hh = rr[np.argmin(abs(rr))] if float(nm[1:]) > 0.3 else rr[np.argmin(abs(rr - L))]
        out["V3_" + nm] = float(np.log(Yabs(hh)))
        out["T_" + nm] = float(np.log(Tabs(z)))
    for nm, z in ZP:
        out["T_" + nm] = float(np.log(Tabs(z)))
        out["Y_" + nm] = float(np.log(Yabs(z)))
    out["T_L"] = float(np.log(Tabs(L)))
    res[N] = out
    json.dump(res, open(W.HIER + "/w36_pruefung.json", "w"), indent=1)
    print(f"N={N} (r={r}, {time.time()-t0:.0f}s): ln T(L) = {out['T_L']:.2f} | max ln T auf Z⁺ = {max(out['T_'+nm] for nm, _ in ZP):.2f} | "
          f"max ln T auf A = {max(out['T_'+nm] for nm, _ in AP):.2f} | max ln Y(h(A)) = {max(out['V3_'+nm] for nm, _ in AP):.2f}", flush=True)
