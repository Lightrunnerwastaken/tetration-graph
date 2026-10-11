"""O14 — Aufgabe 3: Lässt sich O5 (|μ| = 1) auf |μ| > 1 schieben? Möbius-Familie f = M⁻¹(μM), μ = r·e^{i·arg√L}, r = 1 … 1.1717.
Für jedes r: σ₁..σ₄(Y_N) und Rest σ_{r(N)+1} − 1 bei N = 100, 200, 300, 400. Erwartung (O3): für r > 1 Eckmoden mit Rate c(r) > 0."""
import sys, time
import numpy as np
from flint import acb, arb
import o_kern as K
sys.stdout.reconfigure(encoding="utf-8")
argmu = 0.6686178507153447
rs = [float(x) for x in sys.argv[1].split(",")]
Ns = [int(x) for x in sys.argv[2].split(",")]
for r in rs:
    zeile = []
    for N in Ns:
        mu = acb(arb(repr(r)) * arb(repr(argmu)).cos(), arb(repr(r)) * arb(repr(argmu)).sin())
        Y = K.Y_matrix("moeb_lox", N, bits=int(1.2 * N) + 250, mu=mu)
        s = np.linalg.svd(Y, compute_uv=False)
        rr = int((s > 1.01).sum())
        zeile.append(f"N={N}: σ₁={s[0]:.4g}, σ₂={s[1]:.4g}, r={rr}, Rest−1={s[rr]-1:.2e}")
    print(f"|μ| = {r}: " + " | ".join(zeile), flush=True)
