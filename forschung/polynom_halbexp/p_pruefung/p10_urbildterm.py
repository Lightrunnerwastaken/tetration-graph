"""Strang P, Prüfung 10: Wie gross ist der Urbild-Term V₁(h(η)) der Eckformel D6 wirklich — gegen die Schranke, die (E_L) allein
liefert (äusseres Poisson-Integral des Randprofils)?  Kubisches Modell (Daten von Strang W: w_grenzschicht/w01_N{N}.npz, nur lesen).
V₁(x) = Σ_k v_k (ρ/x)^k, v = erste Zeile von Vt (Eingabe-/adjungierte Eckmode, |L|-umskaliert)."""
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
WD = Path(__file__).resolve().parent.parent / "w_grenzschicht"
L = complex(0.31813150520476413, 1.3372357014306895)
rho = abs(L)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
d = np.load(WD / f"w01_N{N}.npz")
S, Vt, coef = d["S"], d["Vt"], d["coef"]
v = Vt[0]
s1 = S[0]
f = np.polynomial.Polynomial(coef)
fp = f.deriv()


def V(x):
    w = rho / x
    # Horner in w
    s = 0j
    for c in v[::-1]:
        s = s * w + c
    return s


thL = np.angle(L)
# Randprofil
M = 4096
th = thL + np.linspace(-np.pi, np.pi, M, endpoint=False)
prof = np.array([np.log(abs(V(rho * np.exp(1j * t))) + 1e-300) for t in th])
print(f"N = {N}, σ₁ = {s1:.4g}; Profil log(σ₁|V₁|)/N: max {np.max(prof + np.log(s1)) / N:.5f}, "
      f"bei τ = 0: {(np.log(abs(V(L))) + np.log(s1)) / N:.5f}")


def poisson_aussen(x):
    """exp(∫ P_x log|V| dθ) — Poisson–Jensen-Schranke für |V(x)|, |x| > ρ (V holomorph aussen inkl. ∞)."""
    r = abs(x)
    ker = (r ** 2 - rho ** 2) / np.abs(x - rho * np.exp(1j * th)) ** 2 / M
    return np.exp(np.sum(ker * prof))


def urbild(eta, x0):
    x = x0
    for _ in range(80):
        x = x - (f(x) - eta) / fp(x)
    return x


print(" Bogen   τ(°)   |η−L|    |h|/ρ    log(σ₁|V(η)|)/N   log(σ₁|V(h)|)/N   Poisson-Schranke/N")
for tau_deg in [0.5, 1, 2, 4, 8, 16]:
    eta = rho * np.exp(1j * (thL + np.deg2rad(tau_deg)))        # D-Bogen (hinter L)
    h = urbild(eta, L + (eta - L) / np.sqrt(L))
    pb = poisson_aussen(h) if abs(h) > rho else np.nan
    print(f" D    {tau_deg:6.1f}  {abs(eta - L):.4f}  {abs(h) / rho:.5f}   {(np.log(abs(V(eta))) + np.log(s1)) / N:+.5f}"
          f"          {(np.log(abs(V(h))) + np.log(s1)) / N:+.5f}           {(np.log(pb) + np.log(s1)) / N:+.5f}")
a0 = coef[0]
R0 = abs(L - a0)
psiL = np.angle(L - a0)
for dpsi in [0.5, 1, 2, 4, 8, 16]:
    eta = a0 + R0 * np.exp(1j * (psiL - np.deg2rad(dpsi)))      # A-Bogen
    h = urbild(eta, L + (eta - L) / np.sqrt(L))
    print(f" A    {dpsi:6.1f}  {abs(eta - L):.4f}  {abs(h) / rho:.5f}   {(np.log(abs(V(eta))) + np.log(s1)) / N:+.5f}"
          f"          {(np.log(abs(V(h))) + np.log(s1)) / N:+.5f}           (h in D)")
