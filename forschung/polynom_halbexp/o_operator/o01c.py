import sys, time
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
mu = float(sys.argv[3]) if len(sys.argv) > 3 else None
for N in [int(s) for s in sys.argv[2].split(",")]:
    t0 = time.time()
    Y = K.Y_matrix(name, N, bits=int(1.2*N) + 250, mu=mu)
    s = np.linalg.svd(Y, compute_uv=False)
    print(f"{name} {mu} N={N} ({time.time()-t0:.0f}s): σ−1 = {np.array2string(s[:8]-1, precision=3)} #σ>1+1e-9: {(s>1+1e-9).sum()}", flush=True)
