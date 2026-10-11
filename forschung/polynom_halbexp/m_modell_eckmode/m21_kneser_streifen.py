"""M21 — 𝒯 für Kneser in Streifenkoordinaten (Gleitkomma, Daten Strang A): Nyström-Budget der IMS-Zerlegung und C_h.
Normierte Variable z = x/ρ, Ecke e^{iθ}. Unterkante (Γin): G = ξ⁻¹∘F∘ξ, W_in = 1/F'(F); Oberkante (Γout): H = ξ⁻¹∘h∘ξ, W_out = F'.
Aufruf: python m21_kneser_streifen.py T hp nk T1 T2
"""
import sys
import time

import numpy as np

from mk_kneser import EITH, MU, TH, F_innen, Fp_innen, h_newton

sys.stdout.reconfigure(encoding="utf-8")
T = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
hp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
nk = int(sys.argv[3]) if len(sys.argv) > 3 else 10
T1 = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
T2 = float(sys.argv[5]) if len(sys.argv) > 5 else 13.0
th = TH
L, Lb = EITH, np.conj(EITH)


def xi(zeta):
    M = np.exp(1j * th) * np.exp(zeta)
    return Lb + (Lb - L) / (M - 1)


def xi_inv(x):
    z = np.log((x - L) / (x - Lb) * np.exp(-1j * th))
    return np.where(z.imag < -1e-12, z + 2j * np.pi, z)


npan = int(round(2 * T / hp))
edges = np.linspace(-T, T, npan + 1)
gx, gw = np.polynomial.legendre.leggauss(nk)
t = np.concatenate([(a + b) / 2 + (b - a) / 2 * gx for a, b in zip(edges[:-1], edges[1:])])
om = np.concatenate([(b - a) / 2 * gw for a, b in zip(edges[:-1], edges[1:])])
n = len(t)
x_lo = xi(t + 0j)
Fx = F_innen(x_lo)
G = xi_inv(Fx)
W_in = 1 / Fp_innen(Fx)
x_up = xi(t + 1j * np.pi)
v_up = h_newton(x_up)
H = xi_inv(v_up)
W_out = Fp_innen(x_up)
w_lo = np.sin(th) / (2 * np.pi * (np.cosh(t) - np.cos(th)))
w_up = np.sin(th) / (2 * np.pi * (np.cosh(t) + np.cos(th)))
ell = np.log(MU)
print(f"Kneser: |F| auf Γin ≤ {np.abs(Fx).max():.6f}, |h| auf Γout ≤ {np.abs(v_up).max():.6f}, Rest h {np.abs(F_innen(v_up)-x_up).max():.1e}")
print(f"Tiefen: Im G ∈ [{G.imag.min():.4f}, {G.imag.max():.4f}], Im H ∈ [{H.imag.min():.4f}, {H.imag.max():.4f}]; arg μ = {np.angle(MU):.4f}")
print(f"G − t: links {G[0]-t[0]:.5f}, rechts {G[-1]-t[-1]:.5f} (ℓ = {ell:.5f}); H − (t+iπ): links {H[0]-t[0]-1j*np.pi:.5f}")
print(f"max|Re G − t| = {np.abs(G.real - t).max():.4f}, max|Re H − t| = {np.abs(H.real - t).max():.4f}; "
      f"|W_in| ≤ {np.abs(W_in).max():.4f}, |W_out| ≤ {np.abs(W_out).max():.4f}")


def P0(p, tp):
    u, y = p.real[:, None], p.imag[:, None]
    return np.sin(y) / (np.cosh(tp[None, :] - u) - np.cos(y)) / (2 * np.pi)


def P1(p, tp):
    u, y = p.real[:, None], p.imag[:, None]
    return np.sin(y) / (np.cosh(tp[None, :] - u) + np.cos(y)) / (2 * np.pi)


t0 = time.time()
sq_lo, sq_up = np.sqrt(w_lo), np.sqrt(w_up)
so = np.sqrt(om)
B = np.zeros((2 * n, 2 * n), complex)
B[:n, :n] = (sq_lo * W_in)[:, None] * P0(G, t) / sq_lo[None, :]
B[:n, n:] = (sq_lo * W_in)[:, None] * P1(G, t) / sq_up[None, :]
B[n:, :n] = (sq_up * W_out)[:, None] * P0(H, t) / sq_lo[None, :]
B[n:, n:] = (sq_up * W_out)[:, None] * P1(H, t) / sq_up[None, :]
S = np.concatenate([so, so])
A = S[:, None] * B * S[None, :]
Ah = 0.5 * (A + A.conj().T)
lam = np.linalg.eigvalsh(np.eye(2 * n) + Ah)
print(f"n = {n} je Kante: λ_min(I + Re K̃) = {lam[0]:.6f}; ‖K̃‖ = {np.linalg.norm(A, 2):.4f} ({time.time()-t0:.0f}s)", flush=True)


def prof(s):
    s = np.clip(s, 0, 1)
    return np.cos(np.pi * s / 2)


W_ = T2 - T1
chiL = prof((t + T2) / W_)
chiR = prof((T2 - t) / W_)
chim = np.sqrt(np.clip(1 - chiL ** 2 - chiR ** 2, 0, 1))
CH = {k: np.concatenate([c, c]) for k, c in (("L", chiL), ("m", chim), ("R", chiR))}
Ssum = sum(CH[k][:, None] * A * CH[k][None, :] for k in CH)
Rh = 0.5 * ((A - Ssum) + (A - Ssum).conj().T)
ev = np.linalg.eigvalsh(Rh)
print(f"IMS T1 = {T1}, T2 = {T2}: ‖Re R‖ = {max(abs(ev[0]), abs(ev[-1])):.4f} (λ_min {ev[0]:.4f})")
for name, c in CH.items():
    idx = np.where(c > 1e-12)[0]
    lm = np.linalg.eigvalsh(np.eye(len(idx)) + Ah[np.ix_(idx, idx)])[0]
    print(f"  Stück {name}: λ_min = {lm:.6f}")
# Carleson-Grössen (Oberkante ohne W_out)
Bh = np.zeros((n, 2 * n), complex)
Bh[:, :n] = sq_up[:, None] * P0(H, t) / sq_lo[None, :]
Bh[:, n:] = sq_up[:, None] * P1(H, t) / sq_up[None, :]
Ch = np.linalg.norm(so[:, None] * Bh * S[None, :], 2)
phs = np.linspace(th, 2 * np.pi - th, 20001)
m_in = np.abs(Fp_innen(F_innen(np.exp(1j * phs)))).min()
CE = np.sqrt(1 / m_in ** 2 + Ch ** 2)
print(f"C_h ≈ {Ch:.4f}, m_in = {m_in:.6f}, C_E ≈ {CE:.4f}; Wahrer Wert λ ≈ {lam[0]:.4f} ⇒ σ_min ≥ λ/C_E ≈ {lam[0]/CE:.4f}")
# Verzerrung r und Re H' (für m20-Analogon)
dt = 1e-5
Hp = (xi_inv(h_newton(xi(t + dt + 1j * np.pi))) - xi_inv(h_newton(xi(t - dt + 1j * np.pi)))) / (2 * dt)


def absxi_abl(zeta):
    M = np.exp(1j * th) * np.exp(zeta)
    return 2 * np.sin(th) * np.exp(zeta.real) / np.abs(M - 1) ** 2


r = absxi_abl(t + 1j * np.pi) / absxi_abl(H)
print(f"Re H' ∈ [{Hp.real.min():.4f}, {Hp.real.max():.4f}], r ≤ {r.max():.4f}")
