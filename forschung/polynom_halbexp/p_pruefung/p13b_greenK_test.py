"""Strang P, Prüfung 13b: Tests der K'-Green-Funktion (Modul kprime) — Randwerte, Harmonizität, g' ≤ g_W, Eckenexponent."""
import sys

import mpmath as mp
import numpy as np

import kprime as K
from reihen import kneser_koeff

sys.stdout.reconfigure(encoding="utf-8")
A = kneser_koeff(4000)
mp.mp.dps = 25
L, rho = K.L, K.rho
a0, R0 = float(A[0]), abs(K.L - float(A[0]))


def fK(z):
    s = mp.mpc(0)
    for c in reversed(A):
        s = s * z + c
    return complex(s)


# g_W explizit (B L7.3)
Lb = np.conj(L)
beta = np.angle((a0 + R0 - L) / (a0 + R0 - Lb))
om = np.angle((-rho - L) / (-rho - Lb)) - beta


def gW(z):
    def Wm(z):
        w = (z - L) / (z - Lb) * np.exp(-1j * beta)
        ang = np.angle(w)
        if ang < -1e-12:
            ang += 2 * np.pi
        return abs(w) ** (np.pi / om) * np.exp(1j * ang * np.pi / om)
    t0 = Wm(1e12)
    t = Wm(z)
    return np.log(abs((t - np.conj(t0)) / (t - t0)))


print("1) Randwerte g' (sollen 0 sein):")
for ph in [0.2, 0.6, 1.0, -0.5]:
    z = fK(mp.mpf(rho) * mp.expj(ph))
    print(f"   f_K(ρe^{{i{ph}}}) = {z:.5f}: g' = {K.g(z):+.2e}")
for th in [1.5, 2.2, 3.14159, -2.0]:
    z = rho * np.exp(1j * th)
    print(f"   ρe^{{i{th}}}: g' = {K.g(z):+.2e}")
print("2) Punkte knapp ausserhalb K' (g' > 0, ≤ g_W):")
for z in [2.2326515 * 1.05, 2.2326515 * 1.01, fK(mp.mpf(rho) * mp.expj(0.6)) * 1.02, -rho * 1.05, 2.6 + 0.3j, 0.4 + 2.0j, -1.8 + 0.9j, 3.0, -3.0]:
    gp, gw = K.g(z), gW(z)
    print(f"   z = {z:.4f}: g' = {gp:.6f}, g_W = {gw:.6f}, g' ≤ g_W: {gp <= gw + 1e-12}")
print("3) Harmonizität (5-Punkt-Laplace, h = 1e-3):")
for z0 in [2.6 + 0.3j, -1.8 + 0.9j, 2.3 + 0.05j, L * 1.05]:
    hh = 1e-3
    lap = (K.g(z0 + hh) + K.g(z0 - hh) + K.g(z0 + 1j * hh) + K.g(z0 - 1j * hh) - 4 * K.g(z0)) / hh ** 2
    print(f"   z = {z0:.4f}: Δg' ≈ {lap:+.2e}")
Z = 1e7
print(f"4) cap(K') = {np.exp(-(K.g(Z) - np.log(Z))):.6f}; τ₀' = {K.tau0:.6f}, κ₀' = {K.kappa0:.6f}")
print("5) Eckenexponent: log-Steigung von g'(L(1+u)) zwischen u = 1e-4 und 1e-3:",
      f"{np.log(K.g(L * (1 + 1e-3)) / K.g(L * (1 + 1e-4))) / np.log(10):.5f}  (1/α' = {1 / K.alpha_p:.5f})")
print("6) Randkoordinaten: Γin bei τ hinter L, f_K(Γout) bei t vor L:")
for tau_deg in [0.1, 0.5, 2, 8, 30, 103.38]:
    eta = rho * np.exp(1j * (np.angle(L) + np.radians(tau_deg)))
    print(f"   Γin τ = {tau_deg:7.2f}°: s = {K.s_Gin(eta):+.6e}")
for tdeg in [0.1, 0.5, 2, 8, 30, 76.6]:
    z = fK(mp.mpf(rho) * mp.expj(np.angle(L) - np.radians(tdeg)))
    print(f"   f(Γout) {tdeg:5.1f}° vor L: |η−L| = {abs(z - L):.4e}, s = {K.s_Gout_x(K.Mw(z).real):+.6e}")
