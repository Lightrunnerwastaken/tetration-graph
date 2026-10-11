"""W17 — Urbild-Term der Eckformel D6 auf dem A-Bogen: σ₁|V₁(h(η))| = |y(h(η))|, h(η) ∈ D.
Vergleich: float-Horner (wie aus den gespeicherten Singulärvektoren) gegen genaue arb-Auswertung von y = J̃^{-T}u₁.
Ausgabe (1/N) ln|y(h(η))| für η auf dem A-Bogen in Abständen |η − L| (unskaliert), kubisches Modell.
Aufruf: python w17_abogen_urbild.py N
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
u1, S1, v1 = d["U"][:, 0], d["S"][0], d["Vt"][0]
t0 = time.time()
J = W.J_mat(W.kubik_a(N, coef), N)
rr = arb(rho)
Dk = [rr ** k for k in range(n)]
Z = J.transpose().solve(arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)]), algorithm="approx")
y = [Z[k, 0] / Dk[k] for k in range(n)]
print(f"N={N}: gelöst ({time.time()-t0:.0f}s), σ1 = {S1:.4f}", flush=True)
Fp = np.polynomial.Polynomial(coef)            # unskaliert
a0 = coef[0]
R0 = abs(L - a0)
psiL = np.angle(L - a0)


def y_genau(x):
    t = x / rho
    w = 1 / acb(t.real, t.imag)
    s = acb(0)
    for v in y[::-1]:
        s = s * w + v
    m = abs(s).mid()
    return float(m.log()) / N if m > 0 else -1e9


def y_float(x):
    return np.log(abs(S1 * W.V_auswerten(v1, x))) / N


print(" |η−L|   h(η): |h|/ρ    (1/N)ln|y(h)| float   genau     (1/N)ln|y(η)| genau")
for dist in (0.02, 0.05, 0.1, 0.19, 0.38, 0.6, 1.0):
    # Punkt auf dem A-Bogen im Abstand dist von L (Richtung zum rechten Scheitel)
    dpsi = 2 * np.arcsin(dist / (2 * R0))
    eta = a0 + R0 * np.exp(1j * (psiL - dpsi))
    roots = (Fp - eta).roots()
    h = roots[np.argmin(abs(roots - L))] if dist < 0.3 else roots[np.argmin(abs(roots))]
    print(f" {dist:5.2f}   {abs(h)/rho:.4f}        {y_float(h):+.4f}          {y_genau(h):+.4f}       {y_genau(eta):+.4f}", flush=True)
