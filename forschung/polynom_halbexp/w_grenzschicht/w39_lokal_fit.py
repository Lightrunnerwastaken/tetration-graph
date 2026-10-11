"""W39 — Auswertung w38: asymptotisches H(r, ψ) um L und Test der lokalen Ecktheorie (Deckverschiebung 2πi − log μ).
Modell: auf der verklebten Fläche (von der Nullkurve Z_T gegen den Uhrzeigersinn bis zum Buchtrand B = h(Z_T)) ist
H = Σ_n Re(c_n w^n), w = (ζ e^{−iψ_Z})^κ, κ₁ = 2πi/(2πi − log μ) = 1.1182 − 0.0317i.  Punkte im ausgeschnittenen Keil (zwischen B und Z_T)
werden mit dem exakten Zweig h auf die Fläche abgebildet (H = H∘h im Halbmond) bzw. ausgelassen (Bucht).
"""
import json
import sys

import numpy as np
from scipy.optimize import least_squares

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
coef, _ = W.kubik_koeff()
L, rho = W.L, W.RHO
argL = np.angle(L)
Fp = np.polynomial.Polynomial(coef)
mu = Fp.deriv()(L)
K1 = 2j * np.pi / (2j * np.pi - np.log(mu))
R = json.load(open(W.HIER + "/w38_lokal.json"))
NS = sorted(int(k) for k in R)
Ns = np.array(NS, float)
A = np.vstack([np.ones_like(Ns), np.log(Ns), Ns]).T


def rate(key):
    y = np.array([R[str(N)][key] for N in NS])
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    c1 = np.linalg.lstsq(A[:-1], y[:-1], rcond=None)[0]
    c2 = np.linalg.lstsq(A[1:], y[1:], rcond=None)[0]
    return c[2], max(abs(c1[2] - c[2]), abs(c2[2] - c[2]), 1e-4)


RADII = [0.03, 0.06, 0.12]
PSIS = list(range(-180, 181, 15))
lamL = rate("L_T")[0]
D = {}
for rad in RADII:
    for ps in PSIS[:-1]:
        D[(rad, ps)] = rate(f"{rad}_{ps}")
# Knoten (lineare Interpolation)
print(f"κ₁ = {K1.real:.4f}{K1.imag:+.4f}i, Knotenabstand π/Re κ₁ = {180/K1.real:.1f}°;  arg μ = {np.degrees(np.angle(mu)):.2f}°, |μ| = {abs(mu):.4f}")
print(f"A-Bogen-Tangente ψ = {np.degrees(np.angle(-1j*(L-coef[0])) - argL):.2f}°, Γout-Tangente −90°, Γin +90°")
nodes = {}
for rad in RADII:
    zs = []
    for a, b in ((-90, -75), (75, 90), (90, 105)):
        la, lb = D[(rad, a)][0] - lamL, D[(rad, b)][0] - lamL
        if la * lb < 0:
            zs.append(a + (b - a) * la / (la - lb))
    nodes[rad] = zs
    print(f"r = {rad}: Knoten ψ = " + ", ".join(f"{z:.1f}°" for z in zs) + (f"  Abstand {zs[1]-zs[0]:.1f}°" if len(zs) > 1 else ""))
psiZ = np.mean([nodes[r][0] for r in RADII])
print(f"Nahtrichtung Z_T: ψ_Z = {psiZ:.1f}°  ⇒ Buchtrand B = h(Z_T): ψ_B = ψ_Z − arg μ = {psiZ - np.degrees(np.angle(mu)):.1f}°")
# Punkte auf der Fläche: Winkel θ = ψ − ψ_Z ∈ [0, 360 − arg μ); Halbmondpunkte (ψ ≤ ψ_Z, ausserhalb D) via exaktes h
pts = []
for (rad, ps), (lam, err) in D.items():
    t = L + rad * np.exp(1j * (argL + np.deg2rad(ps)))
    th = (ps - psiZ) % 360
    if th < 360 - np.degrees(np.angle(mu)) - 3 and not (abs(t) < rho and th > 360 - np.degrees(np.angle(mu)) - 20):
        pts.append((t, lam - lamL, err, "direkt"))
    elif abs(t) > rho:      # Halbmond: H(t) = H(h(t))
        rr = (Fp - t).roots()
        hh = rr[np.argmin(abs(rr - L - (t - L) / mu))]
        pts.append((hh, lam - lamL, err, "h"))
print(f"{len(pts)} Punkte, davon {sum(1 for p in pts if p[3]=='h')} über h abgebildet")


def wcoord(t, kap, psz):
    z = (t - L) * np.exp(-1j * (argL + np.deg2rad(psz)))
    ang = np.angle(z) % (2 * np.pi)
    return np.exp(kap * (np.log(abs(z)) + 1j * ang))


def design(kap, psz, nmax):
    cols = []
    for n in range(1, nmax + 1):
        w = np.array([wcoord(p[0], kap, psz) for p in pts]) ** n
        cols += [w.real, -w.imag]
    return np.array(cols).T


y = np.array([p[1] for p in pts])
wt = 1 / np.array([p[2] for p in pts])
for nmax in (1, 2, 3):
    X = design(K1, psiZ, nmax)
    c, *_ = np.linalg.lstsq(X * wt[:, None], y * wt, rcond=None)
    res = X @ c - y
    print(f"κ = κ₁ fest, n ≤ {nmax}: rms {np.sqrt(np.mean(res**2)):.5f}, max {abs(res).max():.5f}; c_n = "
          + ", ".join(f"{c[2*k]:+.4f}{c[2*k+1]:+.4f}i" for k in range(nmax)))


def resid(p, nmax):
    kap = p[0] + 1j * p[1]
    X = design(kap, p[2], nmax)
    c, *_ = np.linalg.lstsq(X * wt[:, None], y * wt, rcond=None)
    return (X @ c - y) * wt


for nmax in (1, 2):
    sol = least_squares(resid, [K1.real, K1.imag, psiZ], args=(nmax,))
    r_ = sol.fun / wt
    print(f"κ frei, n ≤ {nmax}: κ = {sol.x[0]:.4f}{sol.x[1]:+.4f}i, ψ₀ = {sol.x[2]:.1f}°, rms {np.sqrt(np.mean(r_**2)):.5f}, max {abs(r_).max():.5f}")
# reine Radialexponenten auf festen Strahlen
print("Radialexponent je Strahl (ln|H| über ln r, r = 0.03/0.06/0.12):")
for ps in (-60, -30, 0, 30, 60, 120, 150, 165, 180):
    v = [abs(D[(r, ps if ps < 180 else -180)][0] - lamL) for r in RADII]
    print(f"   ψ = {ps}: {np.polyfit(np.log(RADII), np.log(v), 1)[0]:.3f}")
