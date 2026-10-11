"""W11 — Test: (1/N) log|y(s)| ≈ g_{Ĉ∖W̄}(1/s̄) (Green-Funktion des Äusseren von W̄, Pol ∞, an der Spiegelung am Einheitskreis)
für die adjungierte Eckmode y = J̃^{-T}u₁ (umskaliert, kubisch). Genaue Auswertung in arb.
Aufruf: python w11_green_vergleich.py N
"""
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
n = N + 1
ctx.prec = int(1.5 * N) + 400
coef, _ = W.kubik_koeff()
rho = W.RHO
L = W.L
d = np.load(W.HIER + f"/w01_N{N}.npz")
u1 = d["U"][:, 0]
S1 = d["S"][0]
a = W.kubik_a(N, coef)
J = W.J_mat(a, N)
rr = arb(rho)
Dk = [rr ** k for k in range(n)]
rhs = arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)])
Z = J.transpose().solve(rhs, algorithm="approx")
y = [Z[k, 0] / Dk[k] for k in range(n)]
a0 = coef[0]
R0 = abs(L - a0)


def Y(t):
    w = 1 / acb(t.real, t.imag)
    s = acb(0)
    for v in y[::-1]:
        s = s * w + v
    return s


def lnabs(x):
    m = abs(x).mid()
    return float(m.log()) if m > 0 else -1e9


def green_W(zeta):
    """g von Ĉ∖W̄ (unskaliert), Pol ∞."""
    m = (zeta - L) / (zeta - np.conj(L))
    th1 = np.angle((-rho - L) / (-rho - np.conj(L)))
    th2 = np.angle((a0 + R0 - L) / (a0 + R0 - np.conj(L)))

    def inside(lo, hi, t):
        return ((t - lo) % (2 * np.pi)) < ((hi - lo) % (2 * np.pi))
    if inside(th2, th1, 0.0):
        lo, wid = th2, (th1 - th2) % (2 * np.pi)
    else:
        lo, wid = th1, (th2 - th1) % (2 * np.pi)
    ang = (np.angle(m) - lo) % (2 * np.pi)
    if ang > wid:
        return -1.0        # innerhalb W̄
    tau = (abs(m) ** (np.pi / wid)) * np.exp(1j * ang * np.pi / wid)
    tau0 = np.exp(1j * ((0.0 - lo) % (2 * np.pi)) * np.pi / wid)
    return float(np.log(abs((tau - np.conj(tau0)) / (tau - tau0))))


print(f"N={N}, σ1={S1:.3f}.  Spalten: (1/N)ln|y|  |  g(ρ/s̄)·  (unskaliert: Punkt ρ²/s̄)")
phis = [0, 30, 60, 70, 76.6, 80, 90, 120, 150, 180]
for r in [0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]:
    row1, row2 = [], []
    for p in phis:
        s = r * np.exp(1j * np.deg2rad(p))
        row1.append(lnabs(Y(s)) / N)
        row2.append(green_W(rho / np.conj(s)))      # Spiegelung: rescaled 1/s̄ ↔ unskaliert ρ/s̄
    print(f"r={r:4.2f} y: " + " ".join(f"{x:7.4f}" for x in row1))
    print(f"        g: " + " ".join(f"{x:7.4f}" for x in row2))
