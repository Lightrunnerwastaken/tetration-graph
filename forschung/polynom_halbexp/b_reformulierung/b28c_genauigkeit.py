"""b28c: Ist der wachsende Eck-Singulärwert echt? Vergleich prec 700 vs 1100 und Radien der arb-Inversen (gewichtet)."""
import sys
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
src = open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0]
exec(src)
aL = float(abs(L))
for prec in (700, 1100):
    ctx.prec = prec
    for N in (240, 300):
        pw, g, Ji = baue(N, J=N)
        k = np.arange(N + 1); w = aL ** k
        Jf = np.array([[float(Ji[i, j].mid().str(30, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
        Rf = np.array([[float(Ji[i, j].rad().str(5, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
        A = (w[:, None] * Jf) / w[None, :]; R = (w[:, None] * Rf) / w[None, :]
        sv = np.linalg.svd(A, compute_uv=False)
        print(f"prec={prec} N={N}: σ1 = {sv[0]:.6f}, σ2 = {sv[1]:.6f}, max gewichteter Radius = {R.max():.2e}", flush=True)
