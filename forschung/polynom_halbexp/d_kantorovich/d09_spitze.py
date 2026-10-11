"""D09 — Spitze der Eckspalten von J_N^{-1} am Sprungpunkt L: H_N = max_{j ≥ N−40} max_φ |(J^{-1} x^j)(|L|e^{iφ})| / |L|^j,
Lage der Spitze, Breite, und |(J^{-1}x^j)(L)|. Dazu dieselbe Grösse für die Mitte (j = N/2) zum Vergleich (≈ ‖σ‖_∞ = 1.366).
"""
import json
import sys

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else list(range(40, 195, 10)) + [194]
p = K.loesungen(max(NS), bits=1200)
L = 1.3745570107437075
phiL = np.degrees(np.angle(complex(0.31813150520476413, 1.3372357014306895)))
aus = {}
for N in NS:
    ctx.prec = 500 + 3 * N
    a = [arb(x.mid()) for x in p[N]]
    _, J = K.FJ(a, N)
    n = N + 1
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)]) * (L ** k)[:, None] / (L ** k)[None, :]
    M = 8192
    phi = np.degrees(2 * np.pi * np.arange(M) / M)
    best = (0, None, None)
    for j in range(max(0, N - 40), n):
        F = np.abs(np.fft.ifft(np.concatenate([G[:, j], np.zeros(M - n)])) * M)
        i = int(np.argmax(F))
        if F[i] > best[0]:
            best = (float(F[i]), j, float(phi[i]))
    Fm = np.abs(np.fft.ifft(np.concatenate([G[:, N // 2], np.zeros(M - n)])) * M)
    # Breite der Spitze: Winkelbereich mit |δ| > H/2 um die Spitze
    j = best[1]
    F = np.abs(np.fft.ifft(np.concatenate([G[:, j], np.zeros(M - n)])) * M)
    i0 = int(np.argmax(F))
    lo = i0
    while F[lo % M] > best[0] / 2:
        lo -= 1
    hi = i0
    while F[hi % M] > best[0] / 2:
        hi += 1
    breite = (hi - lo) * 360 / M
    aus[N] = {"H": best[0], "j": j, "phi": best[2], "breite_grad": breite, "mitte": float(Fm.max())}
    print(f"N={N:3d}: H={best[0]:.4f} (Spalte {j}, φ={best[2]:.2f}°, Halbwertsbreite {breite:.2f}°)   Mitte j=N/2: {Fm.max():.4f}",
          flush=True)
json.dump({str(k): v for k, v in aus.items()}, open(K.HIER + "/d09_spitze.json", "w"), indent=1)
print("arg L =", round(float(phiL), 3))
