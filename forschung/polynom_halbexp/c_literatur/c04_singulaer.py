"""Strang C, c04: kleinste Singulärwerte von D_Y J_N D_X^{-1} (l^2-gewichtet, D = diag(r^k)) für mehrere
Gewichtspaare (rX, rY), um die README-Angabe (N = 30: 1.48e-12, 5.4e-12, 1.81e-11, ...) einzuordnen. J an p_N.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, matrix, svd_r

mp.dps = 90
HIER = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HIER, "..", "gemeinsam", "polynome_bis_194.json"), encoding="utf-8"))


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def jacobi(a):
    n = len(a) - 1
    pot = [[mpf(1)] + [mpf(0)] * n]
    for _ in range(1, n + 1):
        pot.append(mul(pot[-1], a, n))
    dp = [k * a[k] for k in range(1, n + 1)]
    dpp = [mpf(0)] * (n + 1)
    for j, c in enumerate(dp):
        for k in range(n + 1):
            dpp[k] += c * pot[j][k]
    return [[pot[j][k] + (dpp[k - j] if j <= k else mpf(0)) for j in range(n + 1)] for k in range(n + 1)]


L_abs = mpf("1.3745570107437075")
N = 30
a = [mpf(x) for x in P[str(N)]]
J = jacobi(a)
for rX, rY in [(1, 1), (mpf("1.3"), mpf("0.5")), (mpf("1.3"), mpf("1.3")), (L_abs, L_abs), (mpf("0.5"), mpf("0.5"))]:
    M = matrix(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            M[k, j] = J[k][j] * mpf(rY) ** k / mpf(rX) ** j
    s = svd_r(M, compute_uv=False)
    s = sorted([float(x) for x in s])
    print(f"N={N} rX={float(rX):.3f} rY={float(rY):.3f}: kleinste SW {', '.join(f'{x:.2e}' for x in s[:5])}; groesster {s[-1]:.2e}")
