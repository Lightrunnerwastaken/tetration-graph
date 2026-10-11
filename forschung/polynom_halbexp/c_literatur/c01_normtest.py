"""Strang C, Prüfrechnung zur Literaturwahl: Wächst ||J_N^{-1}|| auch bei GLEICHEN Gewichten in X und Y exponentiell?

J_N = Jacobi-Matrix des Systems [x^k](p∘p) - 1/k!, k ≤ N, an der Stelle p = p_N (aus gemeinsam/polynome_bis_194.json).
J[k][j] = [x^k] p^j + ([x^{k-j}] p'∘p falls j ≤ k).
Norm: l^1 mit Gewichten w_k = r^k (1+k)^s in X und Y (gleiches r, gleiches s).
Operatornorm von M: max_j (sum_i |M_ij| w_i) / w_j.
Ausgabe: ||J^{-1}||_{r->r} und ||J||_{r->r} für mehrere r, N.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, matrix, lu_solve, fabs, inverse

mp.dps = 120
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
    # Potenzen p^j bis Grad n
    pot = [[mpf(1)] + [mpf(0)] * n]
    for j in range(1, n + 1):
        pot.append(mul(pot[-1], a, n))
    # p' ∘ p bis Grad n
    dp = [k * a[k] for k in range(1, n + 1)] + [mpf(0)]
    dpp = [mpf(0)] * (n + 1)
    for j, c in enumerate(dp):
        if c == 0:
            continue
        for k in range(n + 1):
            dpp[k] += c * pot[j][k]
    J = matrix(n + 1, n + 1)
    for k in range(n + 1):
        for j in range(n + 1):
            v = pot[j][k]
            if j <= k:
                v += dpp[k - j]
            J[k, j] = v
    return J


def opnorm(M, w):
    n = M.rows
    best = mpf(0)
    for j in range(n):
        s = sum(fabs(M[i, j]) * w[i] for i in range(n))
        best = max(best, s / w[j])
    return best


L_abs = mpf("1.3745570107437075")
for N in (10, 20, 30, 40, 50, 60):
    a = [mpf(x) for x in P[str(N)]]
    J = jacobi(a)
    Ji = inverse(J)
    zeile = []
    for r in (mpf("0.5"), mpf("0.8"), mpf("1.0"), mpf("1.3"), L_abs):
        for s in (0, 4):
            w = [r ** k * (1 + k) ** s for k in range(N + 1)]
            zeile.append(f"r={float(r):.2f},s={s}: |J^-1|={float(opnorm(Ji, w)):.2e} |J|={float(opnorm(J, w)):.2e}")
    print(f"N={N}")
    for z in zeile:
        print("   ", z)
    sys.stdout.flush()
