"""M02 — Gestalt der Eckmode im radialen Quadrat-Modell: Koeffizientenprofile und Kreiswerte von u1 (rechts) und v1 (links).
Aufruf: python m02_eckvektoren.py N [alpha]
"""
import sys
import time

import numpy as np

import m_kern as K

sys.stdout.reconfigure(encoding="utf-8")
N = int(sys.argv[1])
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.2283
coef = K.quadratisch(alpha)
rho = K.RHO
t0 = time.time()
S, U, Vt, cols = K.ecksingulaer(coef, N, B=None, k=4)
print(f"N={N}, α={alpha}: σ = {S} ({time.time()-t0:.0f}s)")
u1, v1 = U[:, 0], Vt[0, :]
k = np.arange(N + 1)
for name, vec in (("u1", u1), ("v1", v1)):
    a = np.abs(vec)
    print(f"\n{name}: Koeffizientenprofil |c_k| (Maximum {a.max():.3e} bei k = {int(a.argmax())})")
    for kk in list(range(0, N + 1, max(1, N // 20))) + [N - 5, N - 2, N]:
        print(f"  k={kk:5d} ({kk/N:.3f}): {a[kk]:.3e}")
    w = np.cumsum(a[::-1] ** 2)[::-1]
    for q in (0.5, 0.9, 0.99):
        kk = int(np.argmax(w < (1 - q) * w[0]))
        print(f"  {int(q*100)} % der ℓ²-Masse in k ≥ {kk} ({kk/N:.3f})")
    F = K.auf_kreis(vec)
    M = len(F)
    phi = 360 * np.arange(M) / M
    aF = np.abs(F)
    print(f"  Kreiswerte: max {aF.max():.3e} bei {phi[aF.argmax()]:.2f}°")
    for p in (0, 30, 60, 70, 74, 75, 76, 76.6, 77, 78, 80, 85, 90, 120, 150, 180, 210, 240, 270, 283.4, 290, 300, 330):
        j = int(round(p / 360 * M)) % M
        print(f"   φ={p:6.1f}°: {aF[j]:.3e}")
np.savez(f"m02_vek_N{N}_a{alpha}.npz", S=S, U=U, Vt=Vt)
