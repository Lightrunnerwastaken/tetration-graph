"""b28: Ist sup_N |J_N⁻¹|_{A_s} wirklich endlich für s > ½?  Theorie-Bedenken: Toeplitz mit Sprungsymbol ist auf ℓ²((1+k)^{2s})
nur für |s| < ½ beschränkt (Hunt–Muckenhoupt–Wheeden).  Wir verfolgen die obersten Singulärwerte von D_w J_N⁻¹ D_w⁻¹ in N
und wo der zugehörige Vektor sitzt (Schwerpunkt des Index), für s = 0, 0.25, 0.5, 0.75, 1, 2."""
import sys
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
aL = float(abs(L))
for N in (60, 120, 180, 240):
    pw, g, Ji = baue(N, J=N)
    Jf = np.array([[float(Ji[i, j].mid().str(30, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    zeile = []
    for s in (0.0, 0.25, 0.5, 0.75, 1.0, 2.0):
        k = np.arange(N + 1)
        w = aL ** k * (1.0 + k) ** s
        A = (w[:, None] * Jf) / w[None, :]
        U, sv, Vt = np.linalg.svd(A)
        v = np.abs(Vt[0]) ** 2; u = np.abs(U[:, 0]) ** 2
        zeile.append(f"s={s}: σ1={sv[0]:.4f} σ2={sv[1]:.4f} (Schwerpunkt rechts {np.dot(k, v):.1f}, links {np.dot(k, u):.1f})")
    print(f"N={N}: " + " | ".join(zeile), flush=True)
