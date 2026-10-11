"""W26 — Wie w01 + w24, aber für KNESER: A_f = P_N J∞(f_K) P_N, f_K-Koeffizienten (50 Stellen) aus gemeinsam/kneser_halbexp_koeff_4000.json,
Multiplikator m = f'∘f = eˣ/f'(x) (Reihe bis Grad N). Volle umskalierte Inverse → SVD (Eckmoden u_i), dann Randdaten R_N, Eckteile K_i an denselben
85 Punkten wie w24 (W̄ von Kneser: a₀ = f_K(0)). Ausgabe w26_N{N}.npz (S, U[:, :6]) und w26_N{N}.json (wie w24).
Aufruf: python w26_kneser_rand.py N1,N2,...
"""
import json
import sys
import time

import numpy as np
from flint import acb, acb_poly, acb_series, arb, arb_mat, arb_poly, arb_series, ctx

import w_kern as W
import w24_rand_genau as M

sys.stdout.reconfigure(encoding="utf-8")
KO = json.load(open(W.GEMEINSAM + "/kneser_halbexp_koeff_4000.json"))["a"]
# Geometrie für Kneser: a₀, R₀ (überschreibt die kubischen Werte in M für Punkte/Green)
M.a0 = float(KO[0])
M.R0 = abs(M.L - M.a0)
M.psiL = np.angle(M.L - M.a0)
PTS = M.punkte()
for N in [int(x) for x in sys.argv[1].split(",")]:
    n = N + 1
    ctx.prec = int(0.9 * N) + 400
    ctx.cap = n + 10
    t0 = time.time()
    a = [arb(KO[k]) for k in range(n)]
    P = W.potenzen(a, N)
    fp = arb_poly(a).derivative()
    ex = arb_series([arb(1) / arb.fac_ui(k) for k in range(n)], prec=n)
    mser_arb = (ex / arb_series(list(fp.coeffs()) + [arb(0)] * (n - len(fp.coeffs())), prec=n)).coeffs()
    mser_arb = list(mser_arb) + [arb(0)] * (n - len(mser_arb))
    J = W.J_mat_m(a, mser_arb, N, P)
    G = W.inverse_umskaliert(J, N)
    U, S, Vt = np.linalg.svd(G)
    np.savez(W.HIER + f"/w26_N{N}.npz", S=S, U=U[:, :6], Vt=Vt[:6])
    print(f"N={N} ({time.time()-t0:.0f}s): σ = {S[0]:.5g} {S[1]:.5g} {S[2]:.5g} {S[3]:.5g}", flush=True)
    x0 = acb_poly([acb(c) for c in a])
    mser = acb_poly([acb(c) for c in mser_arb])
    E = arb_mat(n, 2 * len(PTS))
    KZ = []
    for i, (nm, z) in enumerate(PTS):
        zz = acb(z.real, z.imag)
        s1 = list(acb_series((acb_poly([zz]) - x0).coeffs(), prec=n).inv().coeffs())
        s1 += [acb(0)] * (n - len(s1))
        gl = [zz ** (-(k + 1)) for k in range(n)]
        s2 = list((mser * acb_poly(gl)).truncate(n).coeffs())
        s2 += [acb(0)] * (n - len(s2))
        for k in range(n):
            w = s1[k] + s2[k]
            E[k, 2 * i] = w.real
            E[k, 2 * i + 1] = w.imag
        KZ.append(gl)
    X = J.solve(E, algorithm="approx")
    rr = arb(W.RHO)
    Dk = [rr ** k for k in range(n)]
    uarb = [[arb(float(U[k, c])) for k in range(n)] for c in range(6)]
    out = {"S": [float(s) for s in S[:8]], "pts": {}}
    for i, (nm, z) in enumerate(PTS):
        yD = [(acb(X[k, 2 * i], X[k, 2 * i + 1]) - KZ[i][k]) * Dk[k] for k in range(n)]
        Rn = sum((abs(v) for v in yD), arb(0))
        Ks = []
        for c in range(6):
            s = acb(0)
            for k in range(n):
                s += uarb[c][k] * yD[k]
            Ks.append(float(abs(s).mid().log()) if abs(s).mid() > 0 else -1e9)
        out["pts"][nm] = {"z": [z.real, z.imag], "lnR": float(Rn.mid().log()), "lnK": Ks}
    json.dump(out, open(W.HIER + f"/w26_N{N}.json", "w"))
    print(f"   Randdaten ({time.time()-t0:.0f}s): lnR(L) = {out['pts']['L']['lnR']:.3f}, lnK1(L) = {out['pts']['L']['lnK'][0]:.3f}", flush=True)
