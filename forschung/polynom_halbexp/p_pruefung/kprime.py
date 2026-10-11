"""Modul (Strang P): Halbebenen-Koordinate τ' und Green-Funktion g' von E' = Ĉ ∖ K', K' = D̄ ∪ f_K(D̄).
Daten aus p13_greenK.npz (Streifen-Fixpunkt). Auswertung von F(s) = kσ… mit Singularitäts-Abzug nahe dem unteren Rand.

Konventionen: w = log M(ζ), M = (ζ−L)/(ζ−L̄); S' = {φ(x) < y < θ₁}; F: Σ = {0<Im s<π} → S', Ψ = F⁻¹; τ' = e^{Ψ};
Γin ↦ τ' < 0 (s = σ + iπ), f_K(Γout) ↦ τ' > 0 (s = σ). g'(ζ) = log|(τ' − τ̄₀)/(τ' − τ₀)|.
"""
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

HIER = Path(__file__).resolve().parent
_d = np.load(HIER / "p13_greenK.npz")
sig, psi, x_sig = _d["sig"], _d["psi"], _d["x_sig"]
H, phi_inf, theta1, tau0 = float(_d["H"]), float(_d["phi_inf"]), float(_d["theta1"]), complex(_d["tau0"])
xs_curve, ys_curve, tt_curve = _d["xs"], _d["ys"], _d["tt"]
h = sig[1] - sig[0]
k = H / np.pi
alpha_p = H / np.pi
L = complex(0.31813150520476413, 1.3372357014306895)
rho = abs(L)
_spl = CubicSpline(sig, psi)
_dspl = _spl.derivative()
_d2spl = _spl.derivative(2)
dpsi = _dspl(sig)


def _kern_int(vals, s):
    """−(1/π) ∫ v(t)/(e^{s−t} − 1) dt mit v auf dem Gitter (v = ψ oder ψ'), Spline-Werte für den Abzug."""
    u = s - sig
    with np.errstate(over="ignore"):
        E = np.exp(u)
    tau = s.imag
    if tau > 0.4:
        return -(1 / np.pi) * np.sum(vals / (E - 1)) * h
    # regulärer Teil: 1/(e^u − 1) − 1/u
    small = np.abs(u) < 1e-6
    reg = np.where(small, -0.5 + 0j, 1.0 / (E - 1) - 1.0 / np.where(small, 1, u))
    R = np.sum(vals * reg) * h
    # Cauchy-Teil ∫ v(t)/(s − t) dt mit Abzug von (v(σ) + v'(σ)(t−σ))·χ, χ = 1/(1+(t−σ)²)
    sg = s.real
    if vals is psi:
        v0, v1 = float(_spl(sg)), float(_dspl(sg))
    else:
        v0, v1 = float(_dspl(sg)), float(_d2spl(sg))
    tt = sig - sg
    chi = 1 / (1 + tt ** 2)
    num = vals - (v0 + v1 * tt) * chi
    C = np.sum(num / (s - sig)) * h + v0 * (-np.pi * 1j / (tau + 1)) + v1 * (-np.pi / (tau + 1))
    return -(1 / np.pi) * (R + C)


def F(s):
    s = complex(s)
    return k * s + 1j * phi_inf + _kern_int(psi, s)


def dF(s):
    s = complex(s)
    return k + _kern_int(dpsi, s)        # G'(s) = −(1/π)∫ψ'(t)/(e^{s−t}−1)dt


def Psi(w, s0=None):
    s = complex((w - 1j * phi_inf) / k) if s0 is None else complex(s0)
    s = complex(s.real, min(max(s.imag, 1e-3), np.pi - 1e-3))
    for _ in range(80):
        ds = (F(s) - w) / dF(s)
        s -= ds
        s = complex(s.real, min(max(s.imag, 1e-9), np.pi - 1e-9))
        if abs(ds) < 1e-13:
            break
    return s


def Mw(z):
    return np.log((z - L) / (z - np.conj(L)))


def tau(z):
    return np.exp(Psi(Mw(z)))


def g(z):
    t = tau(z)
    return float(np.log(abs((t - np.conj(tau0)) / (t - tau0))))


kappa0 = 2 * np.pi * tau0.imag / abs(tau0) ** 2

# Randkoordinaten
_top_spl = None


def s_Gin(eta):
    """Randpunkt η ∈ Γin (|η| = ρ, links): τ'-Koordinate (negativ)."""
    x = Mw(eta).real
    # Re F(σ + iπ) = kσ + (1/π)∫ψ(t)/(e^{σ−t}+1)dt
    sgm = x / k
    for _ in range(60):
        with np.errstate(over="ignore"):
            E = np.exp(sgm - sig)
        val = k * sgm + (1 / np.pi) * np.sum(psi / (E + 1)) * h
        der = k - (1 / np.pi) * np.sum(psi * E / (E + 1) ** 2) * h
        d = (val - x) / der
        sgm -= d
        if abs(d) < 1e-14:
            break
    return -np.exp(sgm)


_inv = CubicSpline(x_sig, sig)


def s_Gout_x(x):
    """Randpunkt auf f_K(Γout) mit w-Realteil x: τ'-Koordinate (positiv)."""
    return np.exp(_inv(x))


def harm_mass(tz, s_lo, s_hi):
    """ω_τ([s_lo, s_hi]) in der oberen Halbebene."""
    return (np.angle(tz - s_hi) - np.angle(tz - s_lo)) / np.pi
