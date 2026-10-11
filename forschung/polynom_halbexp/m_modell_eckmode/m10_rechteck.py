"""M10 — Rechteck-Abschnitte im normierten radialen Quadrat-Modell: σ_min(R_{N,d}), R_{N,d} = T_{N+d} J_∞|_{P_N} (ℓ², Gewicht 1).
σ_min² = λ_min(RᵀR); berechnet als 1/λ_max((RᵀR)⁻¹) (arb, dann float), da R riesige Einträge hat (~1.48^N).
Aufruf: python m10_rechteck.py eps theta_grad "N1,N2" "d1,d2,..."  (d = -1 bedeutet volle Ausgabe d = N+2)
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
eps = float(sys.argv[1])
thg = float(sys.argv[2])
NS = [int(x) for x in sys.argv[3].split(",")]
DS = [int(x) for x in sys.argv[4].split(",")]
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)
phi = np.linspace(0, 2 * np.pi, 200001)
s = 1 + 2 * eps * (np.cos(phi) - np.cos(th))
m = c0 + 2 * eps * s * np.exp(1j * phi)
gin = np.cos(phi) < np.cos(th)
print(f"ε={eps}, θ={thg}°: min_Γin|m| = {np.abs(m[gin]).min():.6f} bei {np.rad2deg(phi[gin][np.argmin(np.abs(m[gin]))]):.2f}°, "
      f"min_∂D|m| = {np.abs(m).min():.6f}, max_Γin|1/m| = {1/np.abs(m[gin]).min():.6f}", flush=True)


def rmat(N, d):
    """R = T_{N+d} J_∞|_{P_N} als arb_mat ((N+d+1) x (N+1))."""
    M = N + d
    f = arb_poly([arb(eps), arb(c0), arb(eps)])
    mp_ = arb_poly([arb(c0)]) + arb(2 * eps) * f
    P = [arb_poly([1])]
    for _ in range(N):
        P.append((P[-1] * f).truncate(M + 1))
    R = arb_mat(M + 1, N + 1)
    for k in range(N + 1):
        col = P[k]
        mk = (mp_ * arb_poly([0] * k + [1])).truncate(M + 1)
        for j in range(M + 1):
            R[j, k] = col[j] + mk[j]
    return R


for N in NS:
    ctx.prec = int(1.3 * N) + 300
    for d in DS:
        dd = N + 2 if d < 0 else d
        t0 = time.time()
        R = rmat(N, dd)
        G = R.transpose() * R
        n = N + 1
        I = arb_mat(n, n)
        for i in range(n):
            I[i, i] = 1
        Gi = G.solve(I, algorithm="approx")
        Gf = np.array([[float(Gi[i, j].mid()) for j in range(n)] for i in range(n)])
        Gf = 0.5 * (Gf + Gf.T)
        w, V = np.linalg.eigh(Gf)
        smin = 1 / np.sqrt(w[-1])
        v = V[:, -1]
        F = np.abs(np.fft.ifft(np.concatenate([v, np.zeros(8192 - n)])) * 8192)
        kk = np.arange(n)
        kschwer = np.sum(kk * v ** 2)
        print(f"N={N} d={dd}: σ_min(R) = {smin:.6f}, nächste {1/np.sqrt(w[-2]):.6f}, {1/np.sqrt(w[-3]):.6f}; "
              f"Vektor: Schwerpunkt k = {kschwer:.1f}, Kreis-Max bei {360*np.argmax(F)/8192:.2f}° ({time.time()-t0:.0f}s)", flush=True)
