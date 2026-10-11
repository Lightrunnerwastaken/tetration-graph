"""O03 — symmetrische elliptische Möbius-Drehung: max(σ(Y_N)) − 1 als Funktion von (c, N). Wo endet die exakte Kontraktion?"""
import sys, time
import numpy as np
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
cs = [float(s) for s in sys.argv[1].split(",")]
Ns = [int(s) for s in sys.argv[2].split(",")]
for c in cs:
    row = []
    for N in Ns:
        Y = K.Y_matrix("moeb_sym", N, bits=int(1.3 * N) + 200, mu=c)
        s = np.linalg.svd(Y, compute_uv=False)
        row.append(f"{s[0]-1:+.2e}/{s[1]-1:+.1e}")
    print(f"c={c:.3f} (1.5c={1.5*c:.3f}, π/2={np.pi/2:.3f}): " + "  ".join(f"N={N}: {r}" for N, r in zip(Ns, row)), flush=True)
