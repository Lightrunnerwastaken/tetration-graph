"""M15 — 𝒯 im radialen Modell in Streifenkoordinaten (Os ζ = log M, M = (x−L)/(x−L̄), D = Streifen 0 < Im ζ < π).
Unterkante Im ζ = 0 ↔ Γin, Oberkante Im ζ = π ↔ Γout; L bei Re ζ → −∞, L̄ bei +∞.
L²/Poisson-Fassung: Daten x auf beiden Kanten (Mass dφ/2π = w(t)dt), Innenauswertung über den Streifen-Poissonkern.
  (K x)(t) = W_in(t) · P_S[x](G(t))  auf der Unterkante,  G = ξ⁻¹∘f∘ξ,  W_in = 1/m
  (K x)(t) = W_out(t)· P_S[x](H(t))  auf der Oberkante,   H = ξ⁻¹∘h∘ξ,  W_out = f'
Unitär: Y = x·√w. Nyström auf [−T, T] (Gauss-Legendre-Paneele). Ausgaben: Tiefen, λ_min(I + Re K̃) gesamt und für IMS-Stücke
(Ende L: t ≤ −T1, Mitte: |t| ≤ T2, Ende L̄: t ≥ T1), Kommutatorterm ‖Re(K̃ − Σ χ_j K̃ χ_j)‖.
Aufruf: python m15_streifen.py T paneellaenge knoten T1 T2
"""
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
T = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
hp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
nk = int(sys.argv[3]) if len(sys.argv) > 3 else 10
T1 = float(sys.argv[4]) if len(sys.argv) > 4 else 4.0
T2 = float(sys.argv[5]) if len(sys.argv) > 5 else 12.0

eps, thg = 0.31382, 76.618
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)
L = np.exp(1j * th)
Lb = np.conj(L)
f = lambda z: eps + c0 * z + eps * z * z
fp = lambda z: c0 + 2 * eps * z
h = lambda w: (-c0 + np.sqrt(c0 * c0 - 4 * eps * (eps - w) + 0j)) / (2 * eps)
mu = 1 + 2j * eps * np.sin(th)


def xi(zeta):                      # Streifen → Scheibe
    M = np.exp(1j * th) * np.exp(zeta)
    return (Lb * M - L) / (M - 1)


def xi_inv(x):                     # Scheibe → Streifen, Im in (0, π)
    M = (x - L) / (x - Lb)
    z = np.log(M * np.exp(-1j * th))
    # Zweig: Im z in (0, π) für x in D
    z = np.where(z.imag < -1e-12, z + 2j * np.pi, z)
    return z


def xi_abl(zeta):
    M = np.exp(1j * th) * np.exp(zeta)
    return (L - Lb) * M / (M - 1) ** 2


# Knoten
npan = int(round(2 * T / hp))
edges = np.linspace(-T, T, npan + 1)
gx, gw = np.polynomial.legendre.leggauss(nk)
t = np.concatenate([(a + b) / 2 + (b - a) / 2 * gx for a, b in zip(edges[:-1], edges[1:])])
om = np.concatenate([(b - a) / 2 * gw for a, b in zip(edges[:-1], edges[1:])])
n = len(t)
# Unterkante
x_lo = xi(t + 0j)
G = xi_inv(f(x_lo))
W_in = 1 / fp(f(x_lo))
w_lo = np.abs(xi_abl(t + 0j)) / (2 * np.pi)
# Oberkante
x_up = xi(t + 1j * np.pi)
H = xi_inv(h(x_up))
W_out = fp(x_up)
w_up = np.abs(xi_abl(t + 1j * np.pi)) / (2 * np.pi)
print(f"Kontrolle |x|=1: {np.abs(np.abs(x_lo)-1).max():.1e}, {np.abs(np.abs(x_up)-1).max():.1e}; "
      f"Γin-Seite cos φ < cos θ: {np.all(x_lo.real < np.cos(th))}; Γout: {np.all(x_up.real > np.cos(th))}")
print(f"Tiefen: Im G ∈ [{G.imag.min():.4f}, {G.imag.max():.4f}], Im H ∈ [{H.imag.min():.4f}, {H.imag.max():.4f}]; "
      f"arg μ = {np.angle(mu):.4f}, π − arg μ = {np.pi-np.angle(mu):.4f}")
ell = np.log(mu)
print(f"G − t: links {G[0]-t[0]:.5f}, rechts {G[-1]-t[-1]:.5f} (ℓ = {ell:.5f}, −ℓ̄ = {-np.conj(ell):.5f}); "
      f"H − (t+iπ): links {H[0]-t[0]-1j*np.pi:.5f}, rechts {H[-1]-t[-1]-1j*np.pi:.5f}")


def P0(p, tp):    # Poissonkern Unterkante, Dichte bzgl. dt'
    u, y = p.real[:, None], p.imag[:, None]
    return np.sin(y) / (np.cosh(tp[None, :] - u) - np.cos(y)) / (2 * np.pi)


def P1(p, tp):    # Oberkante
    u, y = p.real[:, None], p.imag[:, None]
    return np.sin(y) / (np.cosh(tp[None, :] - u) + np.cos(y)) / (2 * np.pi)


t0 = time.time()
sq_lo, sq_up = np.sqrt(w_lo), np.sqrt(w_up)
so = np.sqrt(om)
# Blöcke (Zeile: Kante der Auswertung, Spalte: Kante der Daten); Kern bzgl. Daten x (Mass dt'), unitär gemacht
B = np.zeros((2 * n, 2 * n), complex)
B[:n, :n] = (sq_lo * W_in)[:, None] * P0(G, t) / sq_lo[None, :]
B[:n, n:] = (sq_lo * W_in)[:, None] * P1(G, t) / sq_up[None, :]
B[n:, :n] = (sq_up * W_out)[:, None] * P0(H, t) / sq_lo[None, :]
B[n:, n:] = (sq_up * W_out)[:, None] * P1(H, t) / sq_up[None, :]
S = np.concatenate([so, so])
A = S[:, None] * B * S[None, :]
I = np.eye(2 * n)
Ah = 0.5 * (A + A.conj().T)
lam = np.linalg.eigvalsh(I + Ah)
print(f"n = {n} je Kante: λ_min(I + Re K̃) gesamt = {lam[0]:.6f}, nächste {lam[1]:.6f}; ‖K̃‖ = {np.linalg.norm(A, 2):.4f} "
      f"({time.time()-t0:.0f}s)", flush=True)


# IMS-Stücke
def prof(s):      # 1 für s ≤ 0, 0 für s ≥ 1, cos-Profil
    s = np.clip(s, 0, 1)
    return np.cos(np.pi * s / 2)


W_ = T2 - T1
chiL = prof((t + T2) / W_)          # 1 für t ≤ −T2, 0 für t ≥ −T1
chiR = prof((T2 - t) / W_)          # 1 für t ≥ T2, 0 für t ≤ T1
chim = np.sqrt(np.clip(1 - chiL ** 2 - chiR ** 2, 0, 1))
CH = {name: np.concatenate([c, c]) for name, c in (("L", chiL), ("m", chim), ("R", chiR))}
Ssum = sum(CH[k][:, None] * A * CH[k][None, :] for k in CH)
Rh = 0.5 * ((A - Ssum) + (A - Ssum).conj().T)
ev = np.linalg.eigvalsh(Rh)
print(f"IMS T1 = {T1}, T2 = {T2} (W = {W_}): ‖Re R‖ = {max(abs(ev[0]), abs(ev[-1])):.4f} (λ_min {ev[0]:.4f}, λ_max {ev[-1]:.4f})")
for name, c in CH.items():
    idx = np.where(c > 1e-12)[0]
    sub = Ah[np.ix_(idx, idx)]
    lm = np.linalg.eigvalsh(np.eye(len(idx)) + sub)[0]
    print(f"  Stück {name}: Träger {len(idx)} Knoten, λ_min(I + Re K̃ | Träger) = {lm:.6f}")
print(f"({time.time()-t0:.0f}s)")
