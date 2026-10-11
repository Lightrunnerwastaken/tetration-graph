"""O01 — Singulärwerte von Y_N = (I + K_N)^{-1} in ℓ²_ρ für vier Modelle (kubisch, Möbius elliptisch, Möbius loxodromisch, affin).
Aufruf: python o01_spektrum.py modell N1,N2,... [bits_faktor]
Ausgabe: obere Singulärwerte, Anzahl σ > 1 + 1e-3, > 1 + 1e-5, Anzahl σ > 0.5 (Γin-Anteil), Vergleich mit |Γin|/2π.
"""
import sys
import time

import numpy as np

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
NS = [int(s) for s in sys.argv[2].split(",")]
fak = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0

f, fp = K.f_numpy(name)
Mq = 1 << 14
phi = 2 * np.pi * np.arange(Mq) / Mq
x = K.RHO * np.exp(1j * phi)
inn = np.abs(f(x)) < K.RHO
print(f"Modell {name}: |Γin|/2π = {inn.mean():.4f}", flush=True)
for N in NS:
    t0 = time.time()
    Y = K.Y_matrix(name, N, bits=int(fak * N) + 250)
    s = np.linalg.svd(Y, compute_uv=False)
    n = N + 1
    print(f"N={N} ({time.time()-t0:.0f}s): σ_1..8 = {np.array2string(s[:8], precision=6)}; "
          f"#σ>1.001: {(s > 1.001).sum()}, #σ>1+1e-5: {(s > 1 + 1e-5).sum()}, #σ>0.5: {(s > 0.5).sum()} "
          f"(= {(s > 0.5).sum()/n:.3f}·(N+1)), #σ∈(0.05,0.95): {((s > 0.05) & (s < 0.95)).sum()}, σ_min = {s[-1]:.2e}",
          flush=True)
