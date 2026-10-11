"""W10 — Adjungierte Eckmode y = J̃^{-T}u₁ (arb) GENAU ausgewertet: Gitter innen/aussen, Urbild-Terme am A-Bogen,
𝒴 = 𝓛y direkt (Summe über Urbilder) gegen die c_m-Darstellung. Kubisches Modell.
Aufruf: python w10_innen_genau.py N
"""
import sys
import time

import numpy as np
from flint import acb, arb, arb_mat, arb_poly, ctx

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
n = N + 1
ctx.prec = int(1.5 * N) + 400
coef, _ = W.kubik_koeff()
rho = W.RHO
d = np.load(W.HIER + f"/w01_N{N}.npz")
u1 = d["U"][:, 0]
t0 = time.time()
a = W.kubik_a(N, coef)
J = W.J_mat(a, N)
rr = arb(rho)
Dk = [rr ** k for k in range(n)]
rhs = arb_mat(n, 1, [Dk[k] * arb(float(u1[k])) for k in range(n)])
Z = J.transpose().solve(rhs, algorithm="approx")
y = [Z[k, 0] / Dk[k] for k in range(n)]
print(f"N={N}: gelöst {time.time()-t0:.0f}s", flush=True)
Ft = [arb(float(coef[j])) * rr ** (j - 1) for j in range(4)]
Fp = np.polynomial.Polynomial([float(v.mid()) for v in Ft])
a0 = float(Ft[0].mid())
Lt = W.L / rho


def Y(t):
    """y(t) = Σ y_k t^{−k} in acb."""
    w = 1 / acb(t.real, t.imag)
    s = acb(0)
    for v in y[::-1]:
        s = s * w + v
    return s


def l10(x):
    m = abs(x).mid()
    return float(m.log()) / 2.302585092994046 if m > 0 else -999


phis = [0, 20, 40, 60, 70, 76.6, 80, 90, 120, 180]
print("log10|y(r e^{iφ})| genau; Zeilen r")
print("  r   " + " ".join(f"{p:7.1f}" for p in phis))
for r in [0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1, 1.3]:
    print(f"{r:5.2f} " + " ".join(f"{l10(Y(r*np.exp(1j*np.deg2rad(p)))):7.2f}" for p in phis), flush=True)
R0 = abs(Lt - a0)
print("\nA-Bogen: ψ, Urbilder z_i (|z|, arg), log10|y(z_i)/(z_i F̃'(z_i))|, log10|Σ| = log10|𝒴(s)|")
for psi in [0, 30, 60, 80, 90, 95, 97]:
    s = a0 + R0 * np.exp(1j * np.deg2rad(psi))
    roots = (Fp - s).roots()
    terms = [Y(z) / acb(complex(z * Fp.deriv()(z)).real, complex(z * Fp.deriv()(z)).imag) for z in roots]
    tot = sum(terms, acb(0))
    print(f" ψ={psi:5.1f}: " + "  ".join(f"({abs(z):.3f},{np.degrees(np.angle(z)):6.1f}):{l10(tt):6.2f}" for z, tt in zip(roots, terms))
          + f"  | Σ: {l10(tot):6.2f}")
