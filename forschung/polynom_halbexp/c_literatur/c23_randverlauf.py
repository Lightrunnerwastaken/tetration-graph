"""Strang C, Runde 3, c23: Verlauf der Pol-Antwort längs ∂W für N = 200, 300, 400 — wo wächst sie exponentiell?

Punkte: Bogen von D(0,ρ) (Winkel arg L … π, symmetrisch) und Bogen von D(a0,R0) (Winkel um a0 von arg(L−a0) bis 0).
Ausserdem: wo auf |x| = ρ ist die Antwort y_{−ρ}(x) gross (Lage der Spitze)?
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, pi, exp, arg, polyval
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
aL, bL = arg(L), arg(L - a0)


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


punkte = []
for t in [0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.0]:      # Anteil des Wegs von arg L bis π auf dem D-Bogen
    punkte.append(("D", t, RHO * exp(1j * (aL + t * (pi - aL)))))
for t in [0.02, 0.1, 0.3, 0.6, 1.0]:                                # D(a0,R0)-Bogen von arg(L−a0) bis 0
    punkte.append(("A", t, a0 + R0 * exp(1j * (bL * (1 - t)))))

for N in [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "200,300,400")]:
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 155)) for j in range(N + 1)] for k in range(N + 1)])
    vs = []
    for _, _, z in punkte:
        nenner = [z - q[0]] + [-q[k] for k in range(1, N + 1)]
        vs.append(cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N))
    rhs = arb_mat([[arb(mp.nstr(x, 155)) for v in vs for x in (v[j].real, v[j].imag)] for j in range(N + 1)])
    sol = A.solve(rhs)
    zeile = []
    ys = []
    for c, (bog, t, z) in enumerate(punkte):
        y = [mpc(mpf(sol[j, 2 * c].mid().str(150, radius=False)), mpf(sol[j, 2 * c + 1].mid().str(150, radius=False))) / z ** (N + 1)
             for j in range(N + 1)]
        ys.append(y)
        R = sum(abs(y[j]) * RHO ** j for j in range(N + 1))
        zeile.append(f"{bog}{t}:{mp.nstr(R, 4)}")
    print(f"N={N}: " + "  ".join(zeile))
    # Lage des Maximums von |y_η(x)| auf |x| = ρ für η = −ρ (Index des Punktes D1.0)
    idx = [i for i, p in enumerate(punkte) if p[0] == "D" and p[1] == 1.0][0]
    y = ys[idx]
    best = max(((abs(polyval(y[::-1], RHO * exp(1j * pi * g / 180))), g) for g in range(0, 181, 2)))
    werte = " ".join(f"{g}°:{mp.nstr(abs(polyval(y[::-1], RHO * exp(1j * pi * g / 180))), 3)}" for g in (0, 40, 70, 76, 77, 80, 90, 120, 180))
    print(f"   |y_(−ρ)(ρe^(iφ))|: max {mp.nstr(best[0], 4)} bei φ = {best[1]}°; {werte}")
    sys.stdout.flush()
