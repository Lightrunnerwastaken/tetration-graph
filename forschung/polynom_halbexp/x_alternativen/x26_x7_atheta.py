"""X26 — Konstanten für Lemma X7 (explizit), die beweisbare Familien-Konstante A_θ und die Schlitzlänge u₀(d) für P-I.
 (1) c_in  = min_{η∈Γin, x∈∂D} |η − f(x)|/|η − x|      (Γin = Teil von ∂D ausserhalb f_K(D̄): arg x − arg L ∈ (0°, 103.38°] und gespiegelt)
     c_out = min_{ξ∈Γout, x∈∂D} |f(ξ) − f(x)|/|ξ − x|   (Γout = Rest von ∂D; ξ ↦ η = f(ξ) ∈ f_K(Γout))
     (Minimumprinzip: die Quotienten sind in x holomorph und nullstellenfrei auf D, stetig auf D̄ ⇒ Minimum auf ∂D.)
     X7-Konstante: |w_η|_{ℓ²_ρ} ≤ C_X7·√M, C_X7 = e^{1/2+1/(4M)}/√(1−1/(4M)) · (1/c_K + max_{∂D}|m|)/ρ, c_K = min(c_in, c_out).
 (2) max_{∂D}|m'| (m = eˣ/f', m' = m(1 − f''/f')) und A_θ(beweisbar) = (1/min|f'|)·[K/ρ + (max|m|/ρ + max|m'|/N)/c₀] + ½M2/(c₀N), K = 5 (Lemma X5).
 (3) g_W̄(L(1+u)) für u bis 0.6 (explizite Green-Funktion wie x24) ⇒ u₀(d) mit g = d auf dem Schlitz.
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
A_str = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]
mp.mp.dps = 40
Lm = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rs = abs(Lm)
L = complex(Lm)
rho = float(rs)
argL = float(mp.arg(Lm))
K = len(A_str)
b0 = np.array([float(mp.mpf(a) * rs ** k) for k, a in enumerate(A_str)])            # f(ρw) = Σ b0[k] w^k
b1 = np.array([float(k * mp.mpf(A_str[k]) * rs ** k) for k in range(1, K)]) / rho     # f'(ρw) = Σ b1[k-1] w^{k-1}
b2 = np.array([float(k * (k - 1) * mp.mpf(A_str[k]) * rs ** k) for k in range(2, K)]) / rho ** 2


def F(w, b):
    return np.polynomial.polynomial.polyval(w, b)


nphi = 6000
phi = np.linspace(-np.pi, np.pi, nphi, endpoint=False) + argL
w = np.exp(1j * phi)
x = rho * w
fx, f1x, f2x = F(w, b0), F(w, b1), F(w, b2)
m = np.exp(x) / f1x
m1 = m * (1 - f2x / f1x)
dl = np.degrees((phi - argL + np.pi) % (2 * np.pi) - np.pi)   # Winkel relativ zu arg L in (−180, 180]
# Γin: Bogen von L gegen den Uhrzeigersinn bis zum Bild von −ρ? Ps Gitter: τ ∈ (0, 103.38°] ab arg L (und Spiegelbild bei L̄)
in_mask = (dl > 0) & (dl <= 103.38)
out_mask = (dl < 0) & (dl >= -76.618)
# (1) c_in: η ∈ Γin (oberer Teil), x ∈ ∂D
eta = x[in_mask][::3]
q_in = np.abs(eta[:, None] - fx[None, :]) / np.maximum(np.abs(eta[:, None] - x[None, :]), 1e-300)
c_in = q_in.min()
# c_out: ξ ∈ Γout (oberer Teil), x ∈ ∂D (Diagonale x = ξ ⇒ |f'(ξ)|)
xi = x[out_mask][::3]
fxi = fx[out_mask][::3]
dd = np.abs(xi[:, None] - x[None, :])
q_out = np.where(dd > 1e-12, np.abs(fxi[:, None] - fx[None, :]) / np.maximum(dd, 1e-300), np.inf)
c_out = min(q_out.min(), np.abs(f1x).min())
print(f"(1) c_in = {c_in:.4f}, c_out = {c_out:.4f} (min|f'| = {np.abs(f1x).min():.4f}), max|m| = {np.abs(m).max():.4f}")
for M in (100, 200, 400, 1000):
    cK = min(c_in, c_out)
    C = np.exp(0.5 + 1 / (4 * M)) / np.sqrt(1 - 1 / (4 * M)) * (1 / cK + np.abs(m).max()) / rho
    print(f"    M = {M}: C_X7 = {C:.3f}  (gemessen max|w_η|/√(M+1) ≈ 1.17)")
# (2) A_θ beweisbar
print(f"(2) max|m'| = {np.abs(m1).max():.4f}, min|f'| = {np.abs(f1x).min():.4f}, M2 aus x21 = 0.8183")
for c0 in (0.7319, 0.54, 0.4249):
    for Kr in (5, 3):
        for N in (150, 193):
            At = (1 / np.abs(f1x).min()) * (Kr / rho + (np.abs(m).max() / rho + np.abs(m1).max() / N) / c0) + 0.5 * 0.8183 / (c0 * N)
            print(f"    c₀ = {c0}, Zeilen bis {Kr}N, N = {N}: A_θ = {At:.2f}")
# (3) Green-Funktion W̄ auf dem Schlitz
a0 = float(mp.mpf(A_str[0]))
R0 = abs(L - a0)
T = lambda z: (z - L) / (z - np.conj(L))
ph = np.linspace(0, 2 * np.pi, 200001)
c1 = rho * np.exp(1j * ph)
c1 = c1[np.abs(c1 - a0) > R0]
c2 = a0 + R0 * np.exp(1j * ph)
c2 = c2[np.abs(c2) > rho]
a1, a2, tinf = np.angle(T(c1[len(c1) // 2])), np.angle(T(c2[len(c2) // 2])), np.angle(T(1e9))
lo, hi = sorted([a1, a2])
if not (lo < tinf < hi):
    lo, hi = hi, lo + 2 * np.pi
alpha = (hi - lo) / np.pi


def g(z):
    t = T(z)
    ang = np.angle(t)
    if ang < lo:
        ang += 2 * np.pi
    s = np.abs(t) ** (1 / alpha) * np.exp(1j * (ang - lo) / alpha)
    p = np.exp(1j * ((tinf - lo) % (2 * np.pi)) / alpha)
    return np.log(np.abs((s - np.conj(p)) / (s - p)))


us = [0.05, 0.1, 0.155, 0.2, 0.3, 0.4, 0.5, 0.6]
print("(3) g_W̄(L(1+u)): " + ", ".join(f"u={u}: {g(L * (1 + u)):.4f}" for u in us))
for d in (0.1, 0.2, 0.3, 0.5):
    lo_u, hi_u = 0.0, 3.0
    for _ in range(60):
        mid = 0.5 * (lo_u + hi_u)
        if g(L * (1 + mid)) < d:
            lo_u = mid
        else:
            hi_u = mid
    print(f"    u₀(d = {d}) = {lo_u:.4f}  (|L(1+u₀)| = {rho * (1 + lo_u):.3f})")
