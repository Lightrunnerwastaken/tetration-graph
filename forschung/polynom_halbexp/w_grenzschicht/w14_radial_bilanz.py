"""W14 — Radiales Modell (Strang M): f(z) = z + ε(z − e^{iθ})(z − e^{−iθ}), auf |z| = 1: f(e^{iφ}) = s(φ)e^{iφ}.
Frage: welche Winkelbereiche der adjungierten Eckmode tragen die Bilanz (J^T y)_k = 0 auf dem Steg (k/N ≈ 0.95)?
(J^T y)_k = (1/2π)∫ Y(φ)[s(φ)^k + m(e^{iφ})] e^{ikφ} dφ,  Y(φ) = Σ_j y_j e^{−ijφ},  y = J^{-T} u₁ (hochgenau).
Zerlegung mit glatter Zerlegung der Eins: Steg (Γin ohne Beine), Γout, Beine bei ±θ.
Aufruf: python w14_radial_bilanz.py N
"""
import sys
import time

import mpmath as mp
import numpy as np
from flint import acb, arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
n = N + 1
TH = np.deg2rad(76.618)
EPS = 0.31382
ctx.prec = int(2.5 * N) + 300
mp.mp.prec = ctx.prec
# Koeffizienten: f = z + ε(z² − 2cosθ z + 1) = ε + (1 − 2ε cosθ) z + ε z²
co = [EPS, 1 - 2 * EPS * np.cos(TH), EPS]
a = [arb(c) for c in co] + [arb(0)] * (n - 3)
J = W.J_mat(a, N)
t0 = time.time()
G = W.inverse_umskaliert(J, N, rho=1.0)
U, S, Vt = np.linalg.svd(G)
print(f"N={N}: σ = {S[:4]}  ({time.time()-t0:.0f}s)", flush=True)
u1 = U[:, 0]
Z = J.transpose().solve(arb_mat(n, 1, [arb(float(x)) for x in u1]), algorithm="approx")
y = [mp.mpf(Z[k, 0].mid().str(int(ctx.prec * 0.3), radius=False)) for k in range(n)]
M = 8 * N
phis = [2 * mp.pi * i / M for i in range(M)]
# Y(φ) auf dem Gitter
Yv = []
for ph in phis:
    w = mp.expjpi(-ph / mp.pi)          # e^{−iφ}
    s_ = mp.mpc(0)
    for c in y[::-1]:
        s_ = s_ * w + c
    Yv.append(s_)
print(f"Y ausgewertet ({time.time()-t0:.0f}s)", flush=True)


def sfun(ph):
    return 1 + 2 * EPS * (mp.cos(ph) - mp.cos(TH))


def mfun(ph):
    z = mp.expjpi(ph / mp.pi)
    fz = EPS + (1 - 2 * EPS * mp.cos(TH)) * z + EPS * z * z
    return 1 + EPS * (2 * fz - 2 * mp.cos(TH))


def bump(x, a_, b_, w_):
    """glatte Indikatorfunktion von [a,b] (periodisch), Übergangsbreite w."""
    def st(t):
        return 0.5 * (1 + mp.tanh(t / w_))
    return st(x - a_) * st(b_ - x)


dl = mp.mpf(np.deg2rad(6.0))
wt = mp.mpf(np.deg2rad(1.5))
thm = mp.mpf(TH)


def regionen(ph):
    p = ph if ph <= mp.pi else ph - 2 * mp.pi          # in (−π, π]
    beinL = bump(p, thm - dl, thm + dl, wt)
    beinLq = bump(p, -thm - dl, -thm + dl, wt)
    steg = bump(abs(p), thm + dl, mp.pi + 1, wt)
    rest = 1 - beinL - beinLq - steg
    return {"Steg": steg, "Gout": rest, "BeinL": beinL, "BeinLq": beinLq}


Rv = [regionen(ph) for ph in phis]
sv = [sfun(ph) for ph in phis]
mv = [mfun(ph) for ph in phis]
names = ["Steg", "Gout", "BeinL", "BeinLq"]
print("k/N | |u_k| | m-Term je Region (log10|·|) | C-Term je Region | Summe gesamt")
for xk in (0.3, 0.6, 0.8, 0.88, 0.92, 0.95, 0.97):
    k = int(xk * N)
    ek = [mp.expjpi(k * ph / mp.pi) for ph in phis]
    mt, ct = {}, {}
    tot = mp.mpc(0)
    for nm in names:
        mt[nm] = sum(Yv[i] * Rv[i][nm] * mv[i] * ek[i] for i in range(M)) / M
        ct[nm] = sum(Yv[i] * Rv[i][nm] * sv[i] ** k * ek[i] for i in range(M)) / M
        tot += mt[nm] + ct[nm]
    def lg(z):
        return float(mp.log10(abs(z))) if abs(z) > 0 else -999
    print(f"{xk:.2f} | {np.log10(abs(u1[k])):6.2f} | " + " ".join(f"{nm}:{lg(mt[nm]):6.2f}" for nm in names) + " | "
          + " ".join(f"{nm}:{lg(ct[nm]):6.2f}" for nm in names) + f" | Σ={lg(tot):6.2f}", flush=True)
