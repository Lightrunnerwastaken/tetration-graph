"""Strang C, Runde 3, c25: Konsistenz δ1 = −J_N(q_N)^{-1}F_N(q_N) für N = 400, 500 — sieht δ1 die Eckspitze?

Druckt α_N = |δ1|_ρ, max_{|x|=ρ}|δ1(x)| und |δ1| bei arg x = 76.6°…77° (Lage der Spitze nach Strang D/c23), sowie bei 0° und 180°.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, pi, exp, polyval
from flint import arb_mat, arb, ctx

N_LISTE = [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "400,500")]
mp.dps = 190
ctx.prec = 1600
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
RHO = abs(mpc("0.31813150520476413", "1.3372357014306895"))


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


for N in N_LISTE:
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 185)) for j in range(N + 1)] for k in range(N + 1)])
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(N + 1)]
    s = A.solve(arb_mat([[arb(mp.nstr(-x, 185))] for x in F]))
    d = [mpf(s[j, 0].mid().str(180, radius=False)) for j in range(N + 1)]
    alpha = sum(abs(d[j]) * RHO ** j for j in range(N + 1))
    werte = {g: abs(polyval(d[::-1], RHO * exp(1j * pi * g / 180))) for g in (0, 60, 70, 75, 76, 76.6, 77, 78, 80, 90, 120, 180)}
    gmax = max(range(0, 181), key=lambda g: abs(polyval(d[::-1], RHO * exp(1j * pi * g / 180))))
    print(f"N={N}: α_N={mp.nstr(alpha, 5)}; max|δ1| auf |x|=ρ bei {gmax}°: {mp.nstr(abs(polyval(d[::-1], RHO * exp(1j * pi * gmax / 180))), 4)}; "
          + " ".join(f"{g}°:{mp.nstr(w, 3)}" for g, w in werte.items()))
    sys.stdout.flush()
