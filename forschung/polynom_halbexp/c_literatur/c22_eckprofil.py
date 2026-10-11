"""Strang C, Runde 3, c22: Profil der Pol-Antwort in der Ecke, y_L = A^{-1}B_N[(L−x)^{-1}]_{>N}, gegen den Kopf
v_lo = [(L−x)^{-1}]_{≤N} (Koeffizienten L^{-j-1}, Betrag ρ^{-j-1}).

Druckt für N = 100, 200, 300 das Verhältnis (y_L)_j / (v_lo)_j = (y_L)_j L^{j+1} in Bändern von j/N, und die Anteile der
Normen |y_L|_ρ aus den Bändern. Frage (Strang B): sitzt die Antwort im Band j/N ∈ [1/|√L|, 1] = [0.853, 1], und ist sie dort
≈ β·v_lo mit festem β?
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, sqrt
from flint import arb_mat, arb, ctx

mp.dps = 140
ctx.prec = 1100
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
L = mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)


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


print("1/|sqrt L| =", mp.nstr(1 / abs(sqrt(L)), 6))
for N in [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "100,200,300")]:
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 135)) for j in range(N + 1)] for k in range(N + 1)])
    nenner = [L - q[0]] + [-q[k] for k in range(1, N + 1)]
    v = cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N)
    rhs = arb_mat([[arb(mp.nstr(v[j].real, 135)), arb(mp.nstr(v[j].imag, 135))] for j in range(N + 1)])
    s = A.solve(rhs)
    y = [mpc(mpf(s[j, 0].mid().str(130, radius=False)), mpf(s[j, 1].mid().str(130, radius=False))) / L ** (N + 1) for j in range(N + 1)]
    ratio = [y[j] * L ** (j + 1) for j in range(N + 1)]   # (y_L)_j / (v_lo)_j
    gesamt = sum(abs(y[j]) * RHO ** j for j in range(N + 1))
    print(f"\nN={N}: |y_L|_ρ = {mp.nstr(gesamt, 5)}, (N+1)/ρ = {mp.nstr((N + 1) / RHO, 5)}")
    baender = [(0, 0.5), (0.5, 0.7), (0.7, 0.8), (0.8, 0.853), (0.853, 0.9), (0.9, 0.95), (0.95, 1.0001)]
    for lo, hi in baender:
        js = [j for j in range(N + 1) if lo <= j / N < hi]
        if not js:
            continue
        anteil = sum(abs(y[j]) * RHO ** j for j in js) / gesamt
        mittel = sum(ratio[j] for j in js) / len(js)
        betr = [abs(ratio[j]) for j in js]
        print(f"   j/N∈[{lo},{hi}): Normanteil {mp.nstr(anteil, 3)}; y_j/v_j Mittel {mp.nstr(mittel, 4)}, |·| min {mp.nstr(min(betr), 3)} max {mp.nstr(max(betr), 3)}")
    sys.stdout.flush()
