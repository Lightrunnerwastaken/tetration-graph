"""O04 — Phasenraum-Struktur von Y_N = (I + K_N)^{-1} (kubisches Modell).
(1) Wellenpakete g_{θ0,k0} (Gauss in θ um θ0, Frequenz k0) : ‖Y g‖/‖g‖ und ⟨Yg, g⟩/‖g‖² — „Symbol“ von Y im Phasenraum (θ, k/N).
(2) Singulärvektoren von Y − T_N(χ) − (Eckmoden) mit σ ≈ 1: Lage in k und θ.
Aufruf: python o04_phasenraum.py N
"""
import sys
import time

import numpy as np

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
name = sys.argv[2] if len(sys.argv) > 2 else "kubisch"
n = N + 1
t0 = time.time()
Y = K.Y_matrix(name, N, bits=int(1.0 * N) + 250)
print(f"Y berechnet ({time.time()-t0:.0f}s)", flush=True)
f, fp = K.f_numpy(name)
R = K.RHO
Mq = 1 << 14
phi = 2 * np.pi * np.arange(Mq) / Mq
xx = R * np.exp(1j * phi)
inn = (np.abs(f(xx)) < R).astype(float)
chi = np.fft.fft(inn) / Mq
T = np.array([[chi[(i - j) % Mq] for j in range(n)] for i in range(n)])
T = T.real if np.allclose(T.imag, 0) else T
thL = np.angle(K.Lc)


def paket(th0, k0, breite):
    """Koeffizienten (Orthonormalbasis) eines Wellenpakets: c_k ∝ exp(−(k−k0)²/(2 b²)) e^{−ik θ0} (lokalisiert bei θ0, Frequenz k0)."""
    k = np.arange(n)
    c = np.exp(-(k - k0) ** 2 / (2 * breite ** 2)) * np.exp(-1j * k * th0)
    return c / np.linalg.norm(c)


print("(1) Symbol von Y im Phasenraum: Einträge ⟨Yg,g⟩ (Betrag ‖Yg‖) für Pakete, Breite in k = √N")
ths = [thL + d for d in (-0.6, -0.2, -0.05, 0.05, 0.2, 0.6, 1.5, np.pi - thL)]
ths_lab = ["L−0.6", "L−0.2", "L−.05", "L+.05", "L+0.2", "L+0.6", "L+1.5", "π"]
b = np.sqrt(N)
for kfrac in (0.1, 0.3, 0.5, 0.7, 0.85, 0.95, 1.0):
    k0 = kfrac * N
    zeile = []
    for th in ths:
        g = paket(th, k0, b)
        Yg = Y @ g
        zeile.append(f"{(np.vdot(g, Yg)).real:+.2f}|{np.linalg.norm(Yg):.2f}")
    print(f"  k/N={kfrac:4.2f}: " + "  ".join(f"{lab}:{z}" for lab, z in zip(ths_lab, zeile)), flush=True)
print("  (θ in Γin für θ > arg L bis 2π − arg L; Γout für |θ| < arg L)")

print("(2) Y − T_N(χ): Singulärwerte und Lage der führenden Singulärvektoren")
U, S, Vh = np.linalg.svd(Y)
r = int((S > 1.001).sum())
Yr = Y - (U[:, :r] * S[:r]) @ Vh[:r]
Dm = Yr - T
u, s, vh = np.linalg.svd(Dm)
print(f"  r = {r} Eckmoden abgezogen; σ(Y − Ecke − T_N(χ)) = {np.array2string(s[:12], precision=4)}")
k = np.arange(n)
for i in range(6):
    for nm, vec in (("Ausgang", u[:, i]), ("Eingang", vh[i].conj())):
        w = np.abs(vec) ** 2
        kmean = (w * k).sum()
        top = w[k > N - 30].sum()
        # Lage auf dem Kreis: Polynom auswerten
        pol = np.polynomial.polynomial.polyval(np.exp(1j * phi[::16]), vec)
        th_max = np.degrees(phi[::16][np.argmax(np.abs(pol))])
        print(f"   i={i} σ={s[i]:.4f} {nm}: mittleres k/N = {kmean/N:.3f}, Masse k>N−30: {top:.2f}, Max auf Kreis bei {th_max:.1f}°")
