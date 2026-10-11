"""Strang C, Runde 3, c20: die Randgrösse S_N = max_{∂W} R_N für grosse N.

Numerisch liegt das Maximum in der Ecke ζ = L (c19, N ≤ 200). Hier: R_N(L) und R_N an weiteren Randpunkten von W für
N = 250 … 400 (eine lineare Lösung je ζ, kein volles Inverses). Frage: wächst S_N polynomiell (≈ N^{1.1}) oder sieht man
Strang Ds Eckspitze (≈ 1.012^N)?
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, pi, exp, arg
from flint import arb_mat, arb, ctx

mp.dps = 160
ctx.prec = 1300
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
L = mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)
a0 = KN[0]
R0 = abs(L - a0)


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


aL, bL = arg(L), arg(L - a0)
for N in [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "250,300,350,400")]:
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 155)) for j in range(N + 1)] for k in range(N + 1)])

    punkte = {"L": L, "rho*e^(i(argL+0.01))": RHO * exp(1j * (aL + 0.01)), "a0+R0*e^(i(arg(L-a0)-0.01))": a0 + R0 * exp(1j * (bL - 0.01)),
              "-rho": -RHO, "a0+R0": a0 + R0, "rho*e^(2.4i)": RHO * exp(2.4j)}
    vs = []
    for z in punkte.values():
        nenner = [z - q[0]] + [-q[k] for k in range(1, N + 1)]
        vs.append(cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N))
    rhs = arb_mat([[arb(mp.nstr(x, 155)) for v in vs for x in (v[j].real, v[j].imag)] for j in range(N + 1)])
    sol = A.solve(rhs)
    werte = []
    for c, z in enumerate(punkte.values()):
        s = sum(abs(mpc(mpf(sol[j, 2 * c].mid().str(150, radius=False)), mpf(sol[j, 2 * c + 1].mid().str(150, radius=False))))
                * RHO ** j for j in range(N + 1))
        werte.append(s / abs(z) ** (N + 1))
    zeile = " | ".join(f"{k}: {mp.nstr(w, 5)}" for k, w in zip(punkte.keys(), werte))
    print(f"N={N}: {zeile}")
    sys.stdout.flush()
