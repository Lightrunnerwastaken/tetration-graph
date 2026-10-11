"""O07b — wie O07, zusätzlich die oberen 8 Singulärvektoren speichern (o07b_N.npz) für die Strukturanalyse der Eckmoden-Paare."""
import sys, time, json
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
for N in [int(s) for s in sys.argv[2].split(",")]:
    t0 = time.time()
    Y = K.Y_matrix(name, N, bits=int(0.85 * N) + 200)
    U, s, Vh = np.linalg.svd(Y)
    r = int((s > 1.01).sum())
    print(f"{name} N={N} ({time.time()-t0:.0f}s): σ_1..12 = {np.array2string(s[:12], precision=5)}; #σ>1.01 = {r}, #σ>1.001 = {(s>1.001).sum()}, "
          f"σ_(r+1) − 1 = {s[r]-1:.3e}", flush=True)
    np.savez_compressed(f"o07b_{name}_{N}.npz", s=s, U=U[:, :8], V=Vh[:8].conj().T)
