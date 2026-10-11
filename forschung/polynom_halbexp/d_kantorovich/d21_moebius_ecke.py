"""D21 — Hat Strang Es Möbius-exp-Modell (Fixpunkte L, L̄, g(0) = 1, elliptisch) dieselbe Eckmode?
f_M(x) = A + B (x − A + B t)/(B + tA − t x), A + iB = L, t = tan(c/2), c = 0.70510 (Strang E).
J_N an T_N f_M (≈ p_N bis auf 0.388^N), letzte BL Spalten von J_N^{-1}, |L|-umskaliert, σ1, σ2, σ3, Lage der Spitze.
Zum Vergleich: Symbol-Sprung bei L wie bei exp (f_M(L) = L), aber |f_M'(L)| = 1 (elliptisch) statt |√L| = 1.172.
"""
import sys
import time

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx

import d_kern as K

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
BL0 = int(sys.argv[2]) if len(sys.argv) > 2 else 100
CC = mp.mpf(sys.argv[3]) if len(sys.argv) > 3 else mp.mpf("0.70510")
mp.mp.dps = 400
A = mp.mpf("0.31813150520476413531265425158766451720351761387139986692237")
Bb = mp.mpf("1.33723570143068940890116214319371061253950213846051241687")
t = mp.tan(CC / 2)
beta0 = Bb + t * A
kappa = Bb * t - A
r = t / beta0
Lr = float(abs(mp.mpc(A, Bb)))


def koeff(N):
    c = [A + Bb / beta0 * kappa]
    for k in range(1, N + 1):
        c.append(Bb / beta0 * (r ** (k - 1) + kappa * r ** k))
    return c


# Kontrolle: f_M(f_M(0)) gegen g(0) = A + B tan(arctan(−A/B) + c)
cf = koeff(60)
f0 = cf[0]
ff0 = A + Bb * (f0 - A + Bb * t) / (Bb + t * A - t * f0)
print("Kontrolle f(f(0)) =", mp.nstr(ff0, 12), " g(0) =", mp.nstr(A + Bb * mp.tan(mp.atan(-A / Bb) + CC), 12),
      " f'(L) =", mp.nstr(abs(Bb ** 2 * (1 + t ** 2) / (Bb + t * A - t * mp.mpc(A, Bb)) ** 2), 8), flush=True)
t0 = time.time()
for N in NS:
    n = N + 1
    BL = min(BL0, n)
    ctx.prec = int(0.9 * N) + 300
    a = [arb(mp.nstr(x, int(ctx.prec * 0.31) + 10)) for x in koeff(N)]
    F, J = K.FJ(a, N)
    res = max(abs(float(F[k].mid())) * Lr ** k for k in range(n))
    E = arb_mat(n, BL)
    for c in range(BL):
        E[N - BL + 1 + c, c] = 1
    X = J.solve(E, algorithm="approx")
    k = np.arange(n)
    cols = np.arange(N - BL + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(BL)] for i in range(n)]) * (Lr ** k)[:, None] / (Lr ** cols)[None, :]
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    M = 8192
    Fu = np.abs(np.fft.ifft(np.concatenate([U[:, 0], np.zeros(M - n)])) * M)
    print(f"N={N} ({time.time()-t0:.0f}s) Rest |F(T_N f_M)| (|L|-gew.) {res:.1e} | σ1={S[0]:.4f} σ2={S[1]:.4f} σ3={S[2]:.4f} | "
          f"u1 max {Fu.max():.3f} bei {360*int(np.argmax(Fu))/M:.2f}°", flush=True)
