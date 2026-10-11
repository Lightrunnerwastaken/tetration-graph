"""W09 — Lokale Frequenz x(φ) der adjungierten Eckmode auf dem Bogen Γin gegen die Normalableitung der Green-Funktion
von Ĉ∖W̄ (W̄ = D̄ ∪ D̄(a₀,R₀)) und von Ĉ∖K' (K' = D̄ ∪ f(D̄)) — Vermutung: Q ≈ Bernstein–Walsh-extremal ⇒ x = ∂_n g.
Green-Funktion von Ĉ∖W̄ explizit: M(ζ) = (ζ−L)/(ζ−L̄) bildet das Äussere auf einen Sektor der Öffnung απ ab.
Aufruf: python w09_frequenz_green.py N
"""
import sys

import numpy as np

import w_kern as W

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
d = np.load(W.HIER + f"/w01_N{N}.npz")
S, Vt = d["S"], d["Vt"]
coef = d["coef"]
rho = W.RHO
L = W.L
a0 = coef[0]
R0 = abs(L - a0)
n = N + 1
k = np.arange(n)
# komplexe L-Kombination: V1 + i·a V2 mit a so, dass bei L̄ ausgelöscht? -> einfacher: Fenster in φ auf Ṽ1
M = 1 << 16
phig = 2 * np.pi * np.arange(M) / M
Vt1 = np.fft.ifft(np.concatenate([Vt[0], np.zeros(M - n)])) * M      # Σ_k v_k e^{+ikφ}  (V(ρ e^{-iφ})...)
# Achtung: V(ρe^{iφ}) = Σ v_k e^{-ikφ} = conj(Σ v_k e^{ikφ}) (v reell) -> Betrag gleich, Frequenz gleich.


def lokale_frequenz(phi0_deg, breite_deg=3.0):
    w = np.exp(-0.5 * ((np.angle(np.exp(1j * (phig - np.deg2rad(phi0_deg))))) / np.deg2rad(breite_deg)) ** 2)
    sig = Vt1 * w
    spec = np.abs(np.fft.fft(sig)) / M
    kk = np.fft.fftfreq(M, d=1.0 / M)
    sel = (kk >= 0) & (kk <= 1.2 * N)
    j = np.argmax(spec[sel])
    return kk[sel][j] / N


# Green-Funktion von Ĉ∖W̄: Sektor-Abbildung. M(ζ) = (ζ−L)/(ζ−L̄); Rand: zwei Strahlen durch 0.
def green_W(zeta):
    m = (zeta - L) / (zeta - np.conj(L))
    # Strahl-Richtungen: Bild des D-Bogens (über −ρ) und des A-Bogens (über a0+R0)
    th1 = np.angle((-rho - L) / (-rho - np.conj(L)))
    th2 = np.angle((a0 + R0 - L) / (a0 + R0 - np.conj(L)))
    # Sektor des Äusseren: enthält Bild von ∞ = 1
    # Öffnung: von th2 nach th1 (gegen Uhrzeigersinn) oder umgekehrt — wähle die, die arg(1)=0 enthält
    def inside(lo, hi, t):
        return ((t - lo) % (2 * np.pi)) < ((hi - lo) % (2 * np.pi))
    if inside(th2, th1, 0.0):
        lo, wid = th2, (th1 - th2) % (2 * np.pi)
    else:
        lo, wid = th1, (th2 - th1) % (2 * np.pi)
    ang = (np.angle(m) - lo) % (2 * np.pi)
    tau = (np.abs(m) ** (np.pi / wid)) * np.exp(1j * ang * np.pi / wid)       # obere Halbebene
    tau0 = 1.0 ** (np.pi / wid) * np.exp(1j * ((0.0 - lo) % (2 * np.pi)) * np.pi / wid)
    return np.log(np.abs((tau - np.conj(tau0)) / (tau - tau0))), wid / np.pi


argL = np.degrees(np.angle(L))
g0, alpha = green_W(np.array([3.0 + 0j]))
print(f"N={N}: α(W̄) = {alpha:.4f};  a0 = {a0:.4f}, R0 = {R0:.4f}")
print(" φ(°)   x(φ) gemessen   ρ∂_r g_W (r=ρ)   Verhältnis")
for phi in [78, 80, 82, 85, 90, 95, 100, 110, 120, 135, 150, 165, 180]:
    x = lokale_frequenz(phi, breite_deg=2.0 if phi < 90 else 4.0)
    z1 = rho * np.exp(1j * np.deg2rad(phi))
    h = 1e-6
    g1, _ = green_W(np.array([z1 * (1 + h)]))
    gn = g1[0] / h          # ρ ∂_r g ≈ (g(ρ(1+h)) − 0)/h
    print(f" {phi:5.1f}   {x:8.4f}        {gn:8.4f}        {x/gn:6.3f}")
