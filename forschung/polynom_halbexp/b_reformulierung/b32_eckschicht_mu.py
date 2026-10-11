"""b32: Liegt die Eckschicht (Spitzen-Vektor v1 und die adjungierte Moden-Funktion y = J⁻ᵀℓ) im Indexband k/N ∈ [1/|μ|, 1]?
Vermutung: Komposition mit f wirkt bei L im Koeffizientenprofil als Streckung x ↦ x/μ (x = k/N); die duale Gleichung koppelt x
an μx und stösst für x > 1/|μ| = 0.853 an die Abschneidung x = 1. Schwerpunkt dort (1 + 1/|μ|)/2 = 0.927 (gemessen 0.93 für v1)."""
import sys
import numpy as np
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
aL = float(abs(L)); mu = complex(mp.exp(L / 2)); print("1/|μ| =", 1 / abs(mu))
for N in (120, 200, 300):
    pw, g, Ji = baue(N, J=N)
    Jf = np.array([[float(Ji[i, j].mid().str(30, radius=False)) for j in range(N + 1)] for i in range(N + 1)])
    k = np.arange(N + 1); w = aL ** k
    A = (w[:, None] * Jf) / w[None, :]
    U, sv, Vt = np.linalg.svd(A)
    def quant(p):
        p = p / p.sum(); c = np.cumsum(p)
        return [k[np.searchsorted(c, q)] / N for q in (0.05, 0.25, 0.5, 0.75, 0.95)]
    print(f"N={N}: σ1={sv[0]:.3f}; v1² Quantile (5/25/50/75/95 %) in k/N: {np.round(quant(Vt[0]**2), 3)}; u1²: {np.round(quant(U[:,0]**2), 3)}")
    # adjungierte Moden-Funktion: ℓ_j = Modus-(+1)-Koeffizient von θ für P = z^j;  y = J⁻ᵀ ℓ  (gewichtet: ỹ = y / w)
    ell = []
    for j in range(N + 1):
        th = [z ** j * gw for (wq, z), gw in zip(pts, gewicht)]
        ell.append(complex(mp.fsum(t * mp.exp(-2j * mp.pi * wq) for t, (wq, z) in zip(th, pts)) / M))
    ell = np.array(ell)
    y = Jf.T @ ell            # Gewichtung: Paarung y·r mit r_k ~ |L|^{-k}  ⇒ relevante Grösse |y_k| |L|^{-k}
    yt = np.abs(y) / w
    print(f"        |y_k|·|L|^-k Quantile in k/N: {np.round(quant(yt**2), 3)};  max bei k/N = {np.argmax(yt)/N:.3f}")
