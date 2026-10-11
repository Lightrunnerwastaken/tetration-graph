"""Strang C, Runde 2, c18: Pol-Antwort für grosse N (250, 300) — sieht R_N die exponentiell wachsende Eckspitze von Strang D?

Druckt R_N(L(1+v/N^0.883))/N^0.883 für v = 0.01, 0.3, 1, 4, 8, ferner R_N an zwei Konturpunkten (Abstand 0.1 von W) und
α_N = |J_N(q_N)^{-1}F_N(q_N)|_0. Kneser-Koeffizienten aus gemeinsam/kneser_halbexp_koeff_1000.json.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc
from flint import arb_mat, arb, ctx

mp.dps = 140
ctx.prec = 1000
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
L = mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)
ALPHA = mpf("0.883")


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def cmul(a, b, n):
    c = [mpc(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def inv_reihe(c, n):
    r = [mpc(0)] * (n + 1)
    r[0] = 1 / c[0]
    for k in range(1, n + 1):
        r[k] = -sum(c[j] * r[k - j] for j in range(1, k + 1)) / c[0]
    return r


for N in [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "250,300")]:
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 135)) for j in range(N + 1)] for k in range(N + 1)])
    Ai = A.inv()
    Aim = [[mpf(Ai[i, j].mid().str(135, radius=False)) for j in range(N + 1)] for i in range(N + 1)]
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(N + 1)]
    d1 = [-sum(Aim[i][j] * F[j] for j in range(N + 1)) for i in range(N + 1)]
    alphaN = sum(abs(d1[i]) * RHO ** i for i in range(N + 1))

    def R(z):
        nenner = [z - q[0]] + [-q[k] for k in range(1, N + 1)]
        v = cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N)
        w = [sum(Aim[i][j] * v[j] for j in range(N + 1)) for i in range(N + 1)]
        return sum(abs(w[i]) * RHO ** i for i in range(N + 1)) / abs(z) ** (N + 1)

    NA = mpf(N) ** ALPHA
    zeile = " ".join(f"v={v}:{mp.nstr(R(L * (1 + v / NA)) / NA, 3)}" for v in (0.01, 0.3, 1, 4, 8))
    kontur = " ".join(f"ζ={mp.nstr(z, 3)}:{mp.nstr(R(z), 3)}" for z in (mpc(-0.01045, 1.473), mpc(1.948, 0), mpc(-1.4746, 0)))
    print(f"N={N}: α_N={mp.nstr(alphaN, 4)} | R_N(L(1+v/N^α))/N^α: {zeile} | Kontur: {kontur}")
    sys.stdout.flush()
