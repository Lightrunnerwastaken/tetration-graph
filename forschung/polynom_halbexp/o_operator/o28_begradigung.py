"""O28 — Konforme Begradigung von Ω̃ = log M(Ω) − iθ₀ = {0 < Im ζ < h(Re ζ)} (Boden = Γin gerade, Deckel = f(Γout) gewellt)
auf den Streifen S_H = {0 < Im w < H}, H = π + β (Höhe an beiden Enden). Ansatz ζ(w) = w + φ(w), φ reell auf Im w = 0:
  φ(w) = c·tanh(πw/(4H)) + ∫ a(σ) e^{iσw} dσ,  a(σ) = −i ê_r(σ)/sinh(σH)  (ê_r: Fourier von e_r = e − c·sech(πt/(2H)), ∫e_r = 0),
  e(t) := Im φ(t + iH) = h(Re ζ(t + iH)) − H   (Fixpunkt-Iteration über die Randzuordnung).
Ausgabe: Konvergenz, Restfehler am Deckel, Endkonstanten φ(±∞), Abklingen von φ − φ(±∞) an den Enden; speichert o28_<modell>.npz
(σ-Gitter, a(σ), c) zur Auswertung von φ an beliebigen Punkten.
"""
import sys

import numpy as np

import o26_streifen as G

name = sys.argv[1] if len(sys.argv) > 1 else "kneser"
sys.stdout.reconfigure(encoding="utf-8")
H = np.pi + G.ell.imag
# Deckelkurve: ζ_top(s) = s + iπ + ℓ_t(s), s ∈ ℝ (Γout-Parameter); als Graph η = h(ξ)
S = np.linspace(-20, 20, 16001)
top = S + 1j * np.pi + G.daten(S, True)[1]
xi_top = top.real
h_top = top.imag
assert np.all(np.diff(xi_top) > 0), "Deckel kein Graph"


def h(xi):
    return np.interp(xi, xi_top, h_top, left=H, right=H)


T = 120.0
n = 2 ** 15
t = np.linspace(-T, T, n, endpoint=False)
dt = t[1] - t[0]
sig = 2 * np.pi * np.fft.fftfreq(n, d=dt)
KAPPE = np.abs(sig) <= 40.0
c = 0.0
a = np.zeros(n, complex)


def phi_line(y):
    """φ(t + iy) auf dem Gitter t (über FFT) plus c·tanh."""
    A_ = np.zeros(n, complex)
    ok = KAPPE
    A_[ok] = a[ok] * np.exp(-sig[ok] * y) * np.exp(1j * sig[ok] * t[0])
    F = np.fft.ifft(A_) * n * (sig[1] - sig[0])
    # ∫ a(σ) e^{iσ(t+iy)} dσ  ≈ Σ a_k e^{iσ_k t} e^{−σ_k y} Δσ ; ifft liefert Σ_k A_k e^{2πi k j/n}/n mit t_j = t0 + j dt
    return F + c * np.tanh(np.pi * (t + 1j * y) / (4 * H))


for it in range(40):
    zt = t + 1j * H + phi_line(H)
    e = h(zt.real) - H
    E = np.sum(e) * dt
    c_new = E / (2 * H)
    er = e - c_new * (1 / np.cosh(np.pi * t / (2 * H)))
    eh = np.fft.fft(er) * dt * np.exp(-1j * sig * t[0]) / (2 * np.pi)      # ê(σ) = (1/2π)∫ e^{−iσt} e(t) dt
    a_new = np.zeros(n, complex)
    ok = KAPPE & (sig != 0)
    a_new[ok] = -1j * eh[ok] / np.sinh(sig[ok] * H)
    diff = np.max(np.abs(a_new - a)) + abs(c_new - c)
    a, c = a_new, c_new
    if diff < 1e-15:
        break
zt = t + 1j * H + phi_line(H)
res = np.abs(zt.imag - h(zt.real))
zb = t + phi_line(0.0)
print(f"{name}: H = {H:.6f}, Iterationen {it + 1}, c = {c:.8f}, max|Im ζ_Boden| = {np.abs(zb.imag).max():.2e}, "
      f"Deckel-Restfehler max = {res[np.abs(t) < T * 0.8].max():.2e}")
mid = np.abs(t) < 50
print(f"φ(±∞) ≈ {phi_line(0)[np.argmin(np.abs(t + 50))].real:+.6f} / {phi_line(0)[np.argmin(np.abs(t - 50))].real:+.6f};  "
      f"max|φ'| am Boden ≈ {np.max(np.abs(np.diff(zb.real) / dt - 1)[mid[:-1]]):.4f}, am Deckel ≈ {np.max(np.abs(np.diff(zt) / dt - 1)[mid[:-1]]):.4f}")
for tt in (-30, -20, -15, -10, -6, -3, 0):
    j = np.argmin(np.abs(t - tt))
    p0 = phi_line(0)[j]
    print(f"  t = {tt:+4d}: φ(t) = {p0.real:+.3e}, φ(t+iH) = {phi_line(H)[j]:+.3e}, φ(t)−φ(−∞) = {(p0 - phi_line(0)[np.argmin(np.abs(t + 50))]).real:+.3e}")
np.savez(f"o28_{name}.npz", sig=sig, a=a, c=c, H=H, t0=t[0], dt=dt, n=n)
