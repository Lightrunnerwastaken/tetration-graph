"""D25 — Reelles Modell mit ZWEI abstossenden Randfixpunkten L, L̄ und Multiplikatoren √L, √L̄ (wie Kneser):
f(x) = x + (x² − 2 Re(L) x + |L|²)(α + β x),  α + β L = (μ − 1)/(L − L̄),  μ = f'(L).
Kubisches Polynom ⇒ J_N(f) exakt. Frage: Eckmode wie bei exp (Rate ≈ 0.012) oder wie im komplexen Möbius-Modell (0.077)?
Aufruf: python d25_kubik_modell.py N1,N2,... [mu als Python-Ausdruck, Standard sqrt(L)]
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
MUEXPR = sys.argv[2] if len(sys.argv) > 2 else "np.sqrt(L)"
L = complex(0.31813150520476413, 1.3372357014306895)
PHI = float(sys.argv[3]) if len(sys.argv) > 3 else None   # Fixpunktwinkel in Grad (Standard arg L)
mu = eval(MUEXPR, {"np": np, "L": L})
if PHI is not None:
    L = abs(L) * np.exp(1j * np.deg2rad(PHI))
rL = (mu - 1) / (L - L.conjugate())
beta = rL.imag / L.imag
alpha = rL.real - beta * L.real
Lr = abs(L)
coef = [Lr ** 2 * alpha, 1 + Lr ** 2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]   # x + q(x)(α+βx)
pf = np.polynomial.Polynomial(coef)
print(f"μ = {mu:.4f} (|μ| = {abs(mu):.4f}), α = {alpha:.5f}, β = {beta:.5f}, f(0) = {coef[0]:.4f}; Kontrolle f(L) − L = {abs(pf(L) - L):.1e}, "
      f"f'(L) = {pf.deriv()(L):.4f}", flush=True)
phi = np.linspace(0, 2 * np.pi, 4001)
x = Lr * np.exp(1j * phi)
inn = np.abs(pf(x)) < Lr
print(f"Anteil Γin auf |x| = |L|: {inn.mean():.3f}; max |f| auf dem Kreis {np.abs(pf(x)).max():.3f}", flush=True)
B = 100
t0 = time.time()
for N in NS:
    n = N + 1
    ctx.prec = int(0.9 * N) + 300
    a = [arb(float(c)) for c in coef] + [arb(0)] * (n - 4)
    _, J = K.FJ(a, N)
    E = arb_mat(n, B)
    for c in range(B):
        E[N - B + 1 + c, c] = 1
    X = J.solve(E, algorithm="approx")
    k = np.arange(n)
    cols = np.arange(N - B + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * (Lr ** k)[:, None] / (Lr ** cols)[None, :]
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    M = 8192
    Fu = np.abs(np.fft.ifft(np.concatenate([U[:, 0], np.zeros(M - n)])) * M)
    print(f"N={N} ({time.time()-t0:.0f}s): σ1={S[0]:.4f} σ2={S[1]:.4f} σ3={S[2]:.4f} | u1 max bei {360*int(np.argmax(Fu))/M:.2f}°",
          flush=True)
