"""Strang C, c08: (1) |J_N(p_N)^{-1}| in l^1 mit w_k = r^k (1+k)^s (Urbild = Bild) für grosse N — welche Gewichte machen
die Inverse gleichmässig beschränkt?  (2) Paulsens Julia-System korrekt (Zeilen i = 0..N-1 wie bei Paulsen 2016, S. 87).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, factorial
from flint import arb_mat, arb, ctx

mp.dps = 160
ctx.prec = 700
HIER = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HIER, "..", "gemeinsam", "polynome_bis_194.json"), encoding="utf-8"))
L_abs = mpf("1.3745570107437075")


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


def inv_mid(M):
    A = arb_mat([[arb(mp.nstr(x, 150)) for x in z] for z in M])
    Ai = A.inv()
    n = len(M)
    return [[abs(float(Ai[i, j].mid())) for j in range(n)] for i in range(n)]


def norm(Ai, w_in, w_out):
    n = len(Ai)
    return max(sum(Ai[i][j] * w_in[i] for i in range(n)) / w_out[j] for j in range(n))


print("(1) |J_N^-1| in l1 mit r^k (1+k)^s:")
for N in (60, 100, 150, 194):
    Ai = inv_mid(jacobi([mpf(x) for x in P[str(N)]]))
    teile = []
    for r in (float(L_abs), 1.40, 1.45, 1.50):
        for s in (0, 2, 4):
            w = [r ** k * (1 + k) ** s for k in range(N + 1)]
            teile.append(f"r={r:.3f},s={s}:{norm(Ai, w, w):.2f}")
    print(f"   N={N:3d}: " + "  ".join(teile))
    sys.stdout.flush()

print("(2) Julia-System (Paulsen): Zeilen i=0..N-1, Unbekannte c_1..c_N; Abel-System (Walker) zum Vergleich")
for N in (20, 40, 60, 80, 120):
    Mj = [[mpf(j) ** i / factorial(i) - (mpf(1) / factorial(i - j) if j <= i else 0) for j in range(1, N + 1)]
          for i in range(0, N)]
    Ma = [[mpf(j) ** i / factorial(i) - (1 if i == j else 0) for j in range(1, N + 1)] for i in range(0, N)]
    Aj, Aa = inv_mid(Mj), inv_mid(Ma)
    out = []
    for r in (float(L_abs), 1.45):
        w_in = [r ** j for j in range(1, N + 1)]
        w_out = [r ** i for i in range(0, N)]
        out.append(f"r={r:.3f}: Julia {norm(Aj, w_in, w_out):.3e} Abel {norm(Aa, w_in, w_out):.3e}")
    print(f"   N={N:3d}: " + " | ".join(out))
    sys.stdout.flush()
