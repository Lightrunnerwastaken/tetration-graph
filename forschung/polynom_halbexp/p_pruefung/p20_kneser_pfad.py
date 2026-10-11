"""Strang P, Prüfung 20 (Anfrage O, Satz O7b): Pfadbedingungen (S1)–(S4) für f_t = (1−t)f_M + t f_K, t ∈ [0,1], f_M = Möbius-exp
(Drehung um L, L̄ mit f_M∘f_M(0) = 1). f_K aus der 4000er-Reihe auf D̄; ausserhalb D̄ über f_K(w) = exp(f_K(Log w)) (|Log w| < ρ).
Gitter, KEIN Einschluss (wie bei O)."""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A = [complex(mp.mpf(s)) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"][:4001]]
A = np.array(A)
dA = np.arange(1, len(A)) * A[1:]
L = complex(mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337)))
rho = abs(L)
Lb = np.conj(L)


def horner(c, z):
    s = np.zeros_like(z, dtype=complex)
    for x in c[::-1]:
        s = s * z + x
    return s


def fK(z):
    z = np.asarray(z, dtype=complex)
    inside = np.abs(z) <= rho * (1 + 1e-12)
    out = np.empty_like(z)
    out[inside] = horner(A, z[inside])
    w = np.log(z[~inside])
    if np.any(np.abs(w) >= rho):
        raise ValueError(f"|Log w| ≥ ρ: {np.max(np.abs(w)):.4f}")
    out[~inside] = np.exp(horner(A, w))
    return out


def fKp(z):
    z = np.asarray(z, dtype=complex)
    inside = np.abs(z) <= rho * (1 + 1e-12)
    out = np.empty_like(z)
    out[inside] = horner(dA, z[inside])
    w = np.log(z[~inside])
    if np.any(np.abs(w) >= rho):
        raise ValueError("|Log w| ≥ ρ")
    out[~inside] = np.exp(horner(A, w)) * horner(dA, w) / z[~inside]
    return out


def M(x):
    return (x - L) / (x - Lb)


def Minv(w):
    return (L - Lb * w) / (1 - w)


lam = np.sqrt(M(1.0) / M(0.0))
if Minv(lam * M(0.0)).real < 0:
    lam = -lam


def fM(x):
    return Minv(lam * M(x))


def fMp(x, h=1e-7):
    return (fM(x + h) - fM(x - h)) / (2 * h)


print(f"f_M(0) = {fM(0.0):.6f} (E: 0.47801), f_M(f_M(0)) = {fM(fM(0.0)):.6f}")
K = 2400
th = 2 * np.pi * (np.arange(K) + 0.5) / K
x = rho * np.exp(1j * th)
fkx, fkpx = fK(x), fKp(x)
fmx, fmpx = fM(x), fMp(x)
argL = np.angle(L)
gin = np.abs(np.angle(x)) > argL
for t in np.linspace(0, 1, 11):
    ft = (1 - t) * fmx + t * fkx
    ftp = (1 - t) * fmpx + t * fkpx
    # (S1) Univalenz: Windungszahl von f_t' (keine Nullstelle von f_t' in D) und einfache Kurve f_t(∂D)
    wn = np.sum(np.angle(np.roll(ftp, -1) / ftp)) / (2 * np.pi)
    sub = ft[::6]
    dmat = np.abs(sub[:, None] - sub[None, :])
    idx = np.arange(len(sub))
    near = np.minimum(np.abs(idx[:, None] - idx[None, :]), len(sub) - np.abs(idx[:, None] - idx[None, :])) < 4
    dmat[near] = np.inf
    selbst = dmat.min()
    # (S3) Γin ↦ D, Γout ↦ aussen
    s3 = np.all(np.abs(ft[gin]) < rho) and np.all(np.abs(ft[~gin]) > rho)
    # m_t = f_t'∘f_t auf ∂D
    try:
        mt = (1 - t) * fMp(ft) + t * fKp(ft)
        mmin = np.min(np.abs(mt))
    except ValueError as e:
        mmin = float("nan")
    mu_t = (1 - t) * fMp(L) + t * np.exp(L / 2)
    sig = np.linspace(-60, 60, 24001)
    s4 = np.min(np.abs(1 + mu_t ** (-1.5 + 1j * sig)))
    print(f"t = {t:.1f}: min|f_t'| (∂D) = {np.min(np.abs(ftp)):.4f}, Windung f_t' = {wn:+.3f}, Selbstabstand f_t(∂D) ≥ {selbst:.4f}, "
          f"(S3) {s3}, min|m_t| = {mmin:.4f}, |μ_t| = {abs(mu_t):.4f}, arg μ_t = {np.degrees(np.angle(mu_t)):.2f}°, min|1+μ_t^(s−1)| = {s4:.6f}")

# ---- Os Bedingungen P1–P5 (o21) für den Kneser-Pfad, Gitter (nicht zertifiziert: f_K nur aus der Koeffizientendatei)
print("\nOs P1–P5 (o21) für f_t = (1−t)f_M + t f_K, Gitter φ ∈ [0, π] (2400 Punkte), t = 0, 0.05, …, 1:")
pol = Minv(np.inf) if False else (L - Lb * 0) # Platzhalter
# Pol von f_M: M(x)·λ = 1 ⇒ M(x) = 1/λ
pol_M = Minv(1 / lam) if abs(1 / lam - 1) > 0 else np.inf
pol_M = (L - Lb / lam) / (1 - 1 / lam)
phis = np.linspace(0, np.pi, 2401)
xx = rho * np.exp(1j * phis)
fk2, fkp2 = fK(xx), fKp(xx)
fm2, fmp2 = fM(xx), fMp(xx)
thL = np.angle(L)
eta = np.radians(1.0)
worst = {"P1": 9, "P2": 9, "P3": 9, "P4": True, "P5": 9}
for t in np.linspace(0, 1, 21):
    ft = (1 - t) * fm2 + t * fk2
    ftp = (1 - t) * fmp2 + t * fkp2
    mt = (1 - t) * fMp(ft) + t * fKp(ft)
    g = np.abs(ft) ** 2 - rho ** 2
    gp = np.gradient(g, phis)
    p4 = (np.all(g[phis <= thL - eta] > 0) and np.all(g[phis >= thL + eta] < 0)
          and np.all(gp[(phis > thL - eta) & (phis < thL + eta)] < 0))
    mu_t = (1 - t) * fMp(L) + t * np.exp(L / 2)
    lg = np.log(mu_t)
    X = np.pi * lg.imag / abs(lg) ** 2
    worst["P1"] = min(worst["P1"], np.min(ftp.real))
    worst["P2"] = min(worst["P2"], np.min(mt.real))
    worst["P3"] = min(worst["P3"], np.min(pol_M.real - ft.real))
    worst["P4"] = worst["P4"] and p4
    worst["P5"] = min(worst["P5"], X)
print(f"   Pol von f_M bei {pol_M:.5f}")
print(f"   P1 min Re f_t' = {worst['P1']:.4f} (> 0?)   P2 min Re m_t = {worst['P2']:.4f}   P3 min (Re p − Re f_t) = {worst['P3']:.4f}   "
      f"P4 Vorzeichen/Monotonie (η = 1°): {worst['P4']}   P5 min X_t = {worst['P5']:.4f} (> 1.5?)")
