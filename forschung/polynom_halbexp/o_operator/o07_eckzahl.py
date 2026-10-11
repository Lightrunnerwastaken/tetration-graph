"""O07 — Anzahl der Eckmoden r(N) und Rest σ_{r+1}(Y_N) − 1 im kubischen Modell bis N groß (volle Matrix, arb).
Aufruf: python o07_eckzahl.py modell N1,N2,...   Ausgabe: obere 10 Singulärwerte, #σ>1.01, #σ>1.001, σ_{r+1}−1 für r = #σ>1.01."""
import sys, time, json
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
out = {}
for N in [int(s) for s in sys.argv[2].split(",")]:
    t0 = time.time()
    Y = K.Y_matrix(name, N, bits=int(0.85 * N) + 200)
    s = np.linalg.svd(Y, compute_uv=False)
    r = int((s > 1.01).sum())
    print(f"{name} N={N} ({time.time()-t0:.0f}s): σ_1..10 = {np.array2string(s[:10], precision=5)}; #σ>1.01 = {r}, #σ>1.001 = {(s>1.001).sum()}, "
          f"σ_(r+1) − 1 = {s[r]-1:.3e}, σ_(r+2) − 1 = {s[r+1]-1:.3e}", flush=True)
    out[N] = [float(x) for x in s[:20]]
    json.dump(out, open(f"o07_{name}.json", "w"))
