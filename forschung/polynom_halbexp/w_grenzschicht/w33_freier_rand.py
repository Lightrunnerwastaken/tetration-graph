"""W33 — Selbstkonsistentes Leitordnungsproblem mit CAUCHY-Bedingung auf Γout (statt Niveaumengen-Iteration, w32 instabil).
Unbekannt: Buchtrand B (gerade Funktion r_B(φ), K Parameter) und v (Ladungen in der Bucht). Bedingungen:
  H := −log|t| + v;  H = 0 auf B;  H(t) = H(h(t)) UND ∂_r H(t) = ∂_r[H∘h](t) für t ∈ Γout (Gleichheit auf offener Menge ⇔ Cauchy-Daten);
  v beschränkt bei ∞ (Σq = 0).
Äussere Minimierung über die B-Parameter (scipy least_squares), innere lineare Kleinste Quadrate für die Ladungen.
Ausgabe: c = H(−1), Profil, Buchtrand(0°), Restnorm. Aufruf: python w33_freier_rand.py [K]
"""
import sys

import numpy as np
from scipy.optimize import least_squares

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
K = int(sys.argv[1]) if len(sys.argv) > 1 else 4
coef, Lf = W.kubik_koeff()
rho = W.RHO
Ft = np.polynomial.Polynomial([coef[j] * rho ** (j - 1) for j in range(4)])
Lt = Lf / rho
phiL = np.angle(Lt)
a0 = Ft(0).real
R0 = abs(Lt - a0)


def h(t):
    r = (Ft - t).roots()
    return r[np.argmin(np.abs(r))]


hv = np.vectorize(h)
psis = np.linspace(np.angle(Lt - a0), -np.angle(Lt - a0), 2001)
B0 = hv(a0 + R0 * np.exp(1j * psis))
o = np.argsort(np.angle(B0))
phi_grid = np.linspace(-phiL, phiL, 301)
r0 = np.interp(phi_grid, np.angle(B0)[o], np.abs(B0)[o])
r0[0] = r0[-1] = 1.0
Gphi = np.linspace(-phiL, phiL, 161)[4:-4]
Gc = np.exp(1j * Gphi)
dlt = 1e-4
hG = hv(Gc)
hGp = hv(Gc * (1 + dlt))
hGm = hv(Gc * (1 - dlt))


def rB_of(a):
    x = phi_grid / phiL
    w = np.exp(sum(a[k] * np.cos(k * np.pi * x / 2) ** 2 * (1 - x ** 2) for k in range(len(a))))   # glatt, =1 an den Ecken
    return 1 - w * (1 - r0)


def baue(a):
    rB = rB_of(a)
    Bc = (rB * np.exp(1j * phi_grid))[3:-3:2]
    srcs = [((rB + f * (1 - rB)) * np.exp(1j * phi_grid))[4:-4:2] for f in (0.2, 0.5, 0.8)]
    src = np.concatenate(srcs)

    def Lm(p):
        return np.log(np.abs(p[:, None] - src[None, :]))
    one = lambda p: np.ones((len(p), 1))
    zero = lambda p: np.zeros((len(p), 1))
    # H = −log|t| + Σq log|t−s| + C
    MB = np.hstack([Lm(Bc), one(Bc)]);  bB = np.log(np.abs(Bc))
    MG = np.hstack([Lm(Gc) - Lm(hG), zero(Gc)]);  bG = np.log(np.abs(Gc)) - np.log(np.abs(hG))
    # radiale Ableitungen per Differenz
    Dp = (Lm(Gc * (1 + dlt)) - Lm(Gc * (1 - dlt))) / (2 * dlt) - (Lm(hGp) - Lm(hGm)) / (2 * dlt)
    bD = (np.log(np.abs(Gc * (1 + dlt))) - np.log(np.abs(Gc * (1 - dlt)))) / (2 * dlt) - (np.log(np.abs(hGp)) - np.log(np.abs(hGm))) / (2 * dlt)
    MD = np.hstack([Dp, zero(Gc)])
    MS = np.hstack([np.ones((1, len(src))), np.zeros((1, 1))])
    Mat = np.vstack([MB, MG, 0.3 * MD, 1e3 * MS])
    rhs = np.concatenate([bB, bG, 0.3 * bD, [0.0]])
    sol, *_ = np.linalg.lstsq(Mat, rhs, rcond=1e-12)
    res = Mat @ sol - rhs
    return sol, src, res


def resid(a):
    return baue(a)[2]


a_init = np.zeros(K)
r_init = resid(a_init)
fit = least_squares(resid, a_init, diff_step=1e-3, max_nfev=200)
sol, src, res = baue(fit.x)
q, C = sol[:-1], sol[-1]


def H(t):
    t = np.atleast_1d(np.asarray(t, complex))
    return -np.log(np.abs(t)) + np.log(np.abs(t[:, None] - src[None, :])) @ q + C


rB = rB_of(fit.x)
print(f"K = {K}: Rest anfangs {np.linalg.norm(r_init):.3e} → {np.linalg.norm(res):.3e}; Parameter {np.round(fit.x, 4)}")
print(f"Buchtrand(0°) = {np.interp(0, phi_grid, rB):.4f} (Start h(A) {np.interp(0, phi_grid, r0):.4f}; gemessen 0.815)")
print(f"c = H(−1) = {H(-1.0)[0]:+.5f}   (gemessen 0.0120)")
print("Profil: " + " ".join(f"{t}:{H(np.exp(1j*(phiL+np.deg2rad(t))))[0]:+.4f}" for t in (1, 2, 4, 8, 15, 30, 60, 103.38)))
hh = 1e-5
print("Stegfrequenz: " + " ".join(f"{p}°:{-(H((1+hh)*np.exp(1j*np.deg2rad(p)))[0]-H((1-hh)*np.exp(1j*np.deg2rad(p)))[0])/(2*hh):.4f}" for p in (80, 90, 120, 180)))
print("Halbmond-Probe H(t) − H(h(t)) bei t = 1.05/1.2/1.3: " + " ".join(f"{H(x)[0]-H(h(x))[0]:+.4f}" for x in (1.05, 1.2, 1.3)))
print("H(1.05), H(1.2), H(1.3) (gemessen 0.30/0.117/0.010 bei N=600): " + " ".join(f"{H(x)[0]:.4f}" for x in (1.05, 1.2, 1.3)))
