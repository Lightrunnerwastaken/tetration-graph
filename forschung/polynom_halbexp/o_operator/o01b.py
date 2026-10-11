import sys, time
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
for N in [int(s) for s in sys.argv[2].split(",")]:
    t0 = time.time()
    Y = K.Y_matrix(name, N, bits=int(N) + 250)
    s = np.linalg.svd(Y, compute_uv=False)
    print(f"{name} N={N} ({time.time()-t0:.0f}s): σ−1 = {np.array2string(s[:10]-1, precision=3)}", flush=True)
