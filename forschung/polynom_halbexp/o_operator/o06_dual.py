"""O06 — Dualbild bei a₀: Singulärwerte von Y'_N = S(I − Y_N)S^{-1} in H²(D(a₀,R₀)) (Strang-O-Satz: J_N ≅ T^{a₀}_N(I + K^{-1}))."""
import sys, time
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
for N in [int(s) for s in sys.argv[2].split(",")]:
    t0 = time.time()
    Yp = K.Ystrich_matrix(name, N)
    s = np.linalg.svd(Yp, compute_uv=False)
    print(f"{name} N={N} ({time.time()-t0:.0f}s): σ(Y') − 1 = {np.array2string(s[:6]-1, precision=3)}, #σ>0.5: {(s>0.5).sum()} ({(s>0.5).sum()/(N+1):.3f}(N+1))", flush=True)
