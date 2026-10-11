"""Strang C, Runde 2, c17: Test des Wachstumsgesetzes (Vermutung R) der Pol-Antwort auf einem Gitter.

Vermutung R: R_N(ζ)^{1/N} → 1/r_W(ζ),  r_W(ζ) = min(|ζ|/ρ, |ζ−a0|/R0),  ρ = |L|, R0 = |L−a0|, a0 = f_K(0).
Beobachtet: rate_obs = (R_160/R_80)^{1/80}. Ausgabe je Gitterpunkt und Statistik von log(rate_obs) − log(1/r_W).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, log
from flint import arb_mat, arb, ctx

mp.dps = 100
ctx.prec = 760
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


def baue(N):
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 95)) for j in range(N + 1)] for k in range(N + 1)])
    Ai = A.inv()
    Aim = [[mpf(Ai[i, j].mid().str(95, radius=False)) for j in range(N + 1)] for i in range(N + 1)]
    return q, pot, Aim


def R(daten, N, z):
    q, pot, Aim = daten
    nenner = [z - q[0]] + [-q[k] for k in range(1, N + 1)]
    v = cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N)
    w = [sum(Aim[i][j] * v[j] for j in range(N + 1)) for i in range(N + 1)]
    return sum(abs(w[i]) * RHO ** i for i in range(N + 1)) / abs(z) ** (N + 1)


N1, N2 = 80, 160
d1, d2 = baue(N1), baue(N2)
abw = []
print(f"a0 = {mp.nstr(a0, 8)}, R0 = {mp.nstr(R0, 8)}, rho = {mp.nstr(RHO, 8)}")
for xr in [-2.2, -1.6, -1.0, -0.4, 0.2, 0.8, 1.4, 1.9, 2.4]:
    for yi in [0.0, 0.6, 1.2, 1.8, 2.4]:
        z = mpc(xr, yi)
        if abs(z - a0) < 0.3 or abs(z - L) < 0.15:
            continue
        rW = min(abs(z) / RHO, abs(z - a0) / R0)
        r1, r2 = R(d1, N1, z), R(d2, N2, z)
        obs = (r2 / r1) ** (mpf(1) / (N2 - N1))
        pred = 1 / rW
        abw.append(float(log(obs) - log(pred)))
        print(f"ζ = {xr:+.1f}{yi:+.1f}i: in W: {'ja ' if rW < 1 else 'nein'} 1/r_W = {mp.nstr(pred, 4):>7}  "
              f"beob. (R_{N2}/R_{N1})^(1/{N2-N1}) = {mp.nstr(obs, 4):>7}  R_{N2} = {mp.nstr(r2, 3)}")
        sys.stdout.flush()
import statistics
print(f"log(beob) − log(vorh): Mittel {statistics.mean(abw):+.4f}, Streuung {statistics.pstdev(abw):.4f}, "
      f"min {min(abw):+.4f}, max {max(abw):+.4f}, n = {len(abw)}")
