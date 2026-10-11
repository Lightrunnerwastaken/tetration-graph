"""O18 — Basispunkt für die Index-Homotopie: ist das kubische Modell mit |μ| = 1 (elliptische Fixpunkte L, L̄) akkretiv (Re K_N ≥ 0)?
und ‖Y_N‖? Hochgenau: Hermitescher Teil in arb (Präzision 1.2N + 300), Eigenwerte in mpmath nach Rundung (N ≤ 40) bzw. float (Kontrolle)."""
import sys
import numpy as np
import mpmath as mp
from flint import acb, arb, ctx
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
for b in [float(x) for x in sys.argv[1].split(",")]:
    ctx.prec = 300
    mu = acb(arb(repr(b)).cos(), arb(repr(b)).sin())
    res = []
    for N in (10, 20, 30, 40):
        J, M, kp = K.matrizen("kubisch", N, bits=int(1.3 * N) + 300, mu=mu)
        IK = M.solve(J, algorithm="approx")             # I + K_N (umskaliert)
        n = N + 1
        mp.mp.dps = 60
        H = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                a = (IK[i, j] - (1 if i == j else 0))
                bb = (IK[j, i] - (1 if i == j else 0))
                H[i, j] = (mp.mpf(a.mid().str(50, radius=False)) + mp.mpf(bb.mid().str(50, radius=False))) / 2
        ev = mp.eigsy(H, eigvals_only=True)
        res.append(mp.nstr(min(ev), 4))
    Y = K.Y_matrix("kubisch", 200, bits=500, mu=mu)
    s = np.linalg.svd(Y, compute_uv=False)
    print(f"kubisch μ = e^(i·{b}): min eig Re K_N (N = 10/20/30/40): {res};  ‖Y_200‖ − 1 = {s[0]-1:.2e}, σ₂ − 1 = {s[1]-1:.2e}", flush=True)
