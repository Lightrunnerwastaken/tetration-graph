"""Strang P, Prüfung 13: Green-Funktion und Halbebenen-Koordinate von E' = Ĉ ∖ K', K' = D̄ ∪ f_K(D̄).

Weg: w = log M(ζ), M = (ζ−L)/(ζ−L̄) bildet E' auf den fast geraden Streifen S' = {φ(x) < y < θ₁} ab (Γin ↦ y = θ₁ = arg L,
f_K(Γout) ↦ y = φ(x), L ↦ x = −∞, L̄ ↦ x = +∞). Konforme Abbildung F: Σ = {0 < Im s < π} → S' mit
  F(s) = (H/π)s + iφ∞ + G(s),  G(s) = −(1/π)∫ψ(t)/(e^{s−t} − 1)dt,  ψ(σ) = φ(Re F(σ)) − φ∞   (Fixpunkt-Iteration),
denn G ist beschränkt analytisch mit Im G = 0 oben und Im G = ψ unten (Poisson-Kern des Streifens). τ' = e^{Ψ(w)}, Ψ = F⁻¹:
E' → obere Halbebene, Γin ↦ ℝ⁻, f_K(Γout) ↦ ℝ⁺, L ↦ 0, L̄ ↦ ∞.  g'(ζ) = log|(τ' − τ̄₀)/(τ' − τ₀)|, τ₀ = τ'(∞).
Ausgabe: p13_greenK.npz (Gitter σ, ψ, x(σ)) und Kennzahlen (Kapazität, Eckenexponent, c'_g, κ₀', u₀'(d)).
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import hilbert

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A = [mp.mpf(s) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]]
mp.mp.dps = 30
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
L = complex(Lm)
rho = abs(L)
argL = np.angle(L)


def fK(z):
    s = mp.mpc(0)
    for c in reversed(A):
        s = s * z + c
    return complex(s)


def Mw(z):
    return np.log((z - L) / (z - np.conj(L)))


theta1 = Mw(-rho).imag                      # = arg L
# ---- Kurve f_K(Γout) im w-Bild; dicht an den Enden (|ζ−L| bis ~1e-4)
ts = argL * np.tanh(np.linspace(-5.3, 5.3, 1400))          # t ∈ (−arg L, arg L), Häufung an den Enden
pts = []
for t in ts:
    z = fK(mp.mpf(rho) * mp.expj(t))
    w = Mw(z)
    pts.append((w.real, w.imag, t))
pts = np.array(pts)
order = np.argsort(pts[:, 0])
xs, ys, tt = pts[order, 0], pts[order, 1], pts[order, 2]
print(f"θ₁ = {np.degrees(theta1):.6f}°, x-Bereich der Kurve [{xs[0]:.2f}, {xs[-1]:.2f}], monoton: {np.all(np.diff(xs) > 0)}")
phi_inf_L = ys[0]
phi_inf_R = ys[-1]
phi_inf = np.angle(np.exp(1j * (argL - np.pi / 2 + np.angle(np.sqrt(L))))) - np.pi / 2
print(f"φ an den Enden: {np.degrees(phi_inf_L):.4f}° / {np.degrees(phi_inf_R):.4f}°, theoretisch φ∞ = {np.degrees(phi_inf):.4f}°;"
      f" min φ = {np.degrees(ys.min()):.4f}°; Symmetrie max|φ(x) − φ(−x)| ≈ "
      f"{np.max(np.abs(np.interp(-xs[xs > 0][::-1], xs, ys) - ys[xs > 0][::-1])):.2e}")
H = theta1 - phi_inf
print(f"Streifenhöhe H = {np.degrees(H):.5f}° (= Aussenwinkel), α' = H/π = {H / np.pi:.6f}")

# φ(x) mit exponentieller Fortsetzung an den Enden
spl = CubicSpline(xs, ys - phi_inf)
xl, xr = xs[2], xs[-3]
al = (ys[2] - phi_inf) / np.exp(xl)
ar = (ys[-3] - phi_inf) / np.exp(-xr)


def psi_of_x(x):
    x = np.asarray(x)
    out = np.where(x < xl, al * np.exp(np.minimum(x, 50)), np.where(x > xr, ar * np.exp(-np.maximum(x, -50)), 0.0))
    mid = (x >= xl) & (x <= xr)
    out = out.astype(float)
    out[mid] = spl(x[mid])
    return out


# ---- Fixpunkt auf dem unteren Rand
S, n = 60.0, 4096
sig = np.linspace(-S, S, n, endpoint=False)
h = sig[1] - sig[0]
k = H / np.pi
U = sig[:, None] - sig[None, :]
with np.errstate(divide="ignore", invalid="ignore"):
    Kreg = np.where(np.abs(U) < 1e-8, -0.5, 1.0 / np.expm1(U) - 1.0 / U)
psi = np.zeros(n)
for it in range(60):
    hil = np.imag(hilbert(psi))                       # (1/π) PV ∫ ψ(t)/(σ − t) dt
    ReG = -hil - (1 / np.pi) * (Kreg @ psi) * h
    x_sig = k * sig + ReG
    new = psi_of_x(x_sig)
    d = np.max(np.abs(new - psi))
    psi = new
    if d < 1e-13:
        break
print(f"Fixpunkt: {it + 1} Iterationen, letzte Änderung {d:.1e}; ∫ψ = {psi.sum() * h:.6f}, max|ψ| = {np.max(np.abs(psi)):.4f}")
hil = np.imag(hilbert(psi))
ReG = -hil - (1 / np.pi) * (Kreg @ psi) * h
x_sig = k * sig + ReG
print(f"x(σ) monoton: {np.all(np.diff(x_sig) > 0)}")


def F(s):
    """F(s) für Im s ∈ (0, π) (nicht zu nah an 0)."""
    s = np.atleast_1d(s)
    E = np.exp(s[:, None] - sig[None, :])
    G = -(1 / np.pi) * ((psi[None, :] / (E - 1)).sum(axis=1)) * h
    return k * s + 1j * phi_inf + G


def dF(s):
    s = np.atleast_1d(s)
    E = np.exp(s[:, None] - sig[None, :])
    return k + (1 / np.pi) * ((psi[None, :] * E / (E - 1) ** 2).sum(axis=1)) * h


def Psi(w, s0=None):
    s = (w - 1j * phi_inf) / k if s0 is None else s0
    s = complex(s)
    for _ in range(60):
        ds = (F(s)[0] - w) / dF(s)[0]
        s -= ds
        if abs(ds) < 1e-13:
            break
    return s


def tau(z):
    return np.exp(Psi(Mw(z)))


tau0 = np.exp(Psi(0j))
kappa0 = 2 * np.pi * tau0.imag / abs(tau0) ** 2


def g(z):
    t = tau(z)
    return np.log(abs((t - np.conj(tau0)) / (t - tau0)))


print(f"τ₀' = {tau0:.8f}, κ₀' = {kappa0:.6f}")
# Proben: g auf dem Rand (nahe), Kapazität, Harmonizität, Vergleich mit g_W
Z = 1e6
cap = np.exp(-(g(Z) - np.log(Z)))
print(f"cap(K') = {cap:.6f}   (cap(W̄) = 1.516889)")
for z0 in [2.6 + 0.3j, -1.8 + 0.9j, 0.4 + 2.0j]:
    hh = 1e-3
    lap = (g(z0 + hh) + g(z0 - hh) + g(z0 + 1j * hh) + g(z0 - 1j * hh) - 4 * g(z0)) / hh ** 2
    print(f"  Laplace g' bei {z0}: {lap:.2e}")
zb = fK(mp.mpf(rho) * mp.expj(0.6))
print(f"  g' auf f_K(Γout) (φ = 0.6): {g(zb * (1 + 1e-9)):.2e};  auf Γin (−ρ·1.000001): {g(-rho * 1.000001):.2e}")
print(f"  g'(2.2327·1.05) = {g(2.2326515 * 1.05):.5f}; g'(−2.0) = {g(-2.0):.5f}")
alpha_p = H / np.pi
print("Schlitz ζ = L(1+u): g'/u^{1/α'}")
cg = []
for u in [1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2, 0.3]:
    gv = g(L * (1 + u))
    cg.append(gv / u ** (1 / alpha_p))
    print(f"  u = {u:<6}: g' = {gv:.6e},  g'/u^(1/α') = {gv / u ** (1 / alpha_p):.6f},  τ' = {tau(L * (1 + u)):.5e}, arg τ' = {np.degrees(np.angle(tau(L * (1 + u)))):.3f}°")
for dlev in [0.05, 0.1]:
    lo, hi = 1e-4, 2.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if g(L * (1 + mid)) < dlev:
            lo = mid
        else:
            hi = mid
    print(f"u₀'({dlev}) = {lo:.5f}  (|ζ−L| = {lo * rho:.4f})")
np.savez(HIER / "p13_greenK.npz", sig=sig, psi=psi, x_sig=x_sig, H=H, phi_inf=phi_inf, theta1=theta1, tau0=tau0,
         xs=xs, ys=ys, tt=tt)
print("gespeichert: p13_greenK.npz")
