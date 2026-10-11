"""Strang C, Runde 2, c15: Skalengesetz der Pol-Antwort und Kontur-Schranke für die Konsistenz.

R_N(ζ) := |A^{-1} B_N [(ζ−x)^{-1}]_{>N}|_0 = |A^{-1} v_ζ|_0 / |ζ|^{N+1},  v_ζ = [q^{N+1}/(ζ−q)]_{≤N},  A = J_N(q_N), |·|_0 = l¹_{|L|}.
(1) Radialer Schnitt ζ = L(1+u), u = v/N^α (α = 0.883, Aussenwinkel nach Strang B): Tabelle R_N/N^α gegen v.
(2) Kontur Γ'_d: Rand von W = D(0,|L|) ∪ D(a0,|L−a0|), um d nach aussen geschoben, ohne |ζ−L|,|ζ−L̄| < 0.25: max R_N.
(3) Schranke I_N = 2·(1/2π)∫_0^{u0} |Sprung(u)| R_N(L(1+u)) |L| du  (Leitterm D1, Sprung = |D1|·2|sin πβ|·u^{Re β}),
    verglichen mit α_N = |δ1|_0 und mit |Antwort auf D1·u_β + konj|_0 (aus c14).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, matrix, pi, exp, sin, cos, fabs
from mpmath.calculus.quadrature import GaussLegendre
from flint import arb_mat, arb, ctx

mp.dps = 90
ctx.prec = 720
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
SING = json.load(open(os.path.join(G, "kneser_singularitaet_L.json"), encoding="utf-8"))
D1 = mpc(*SING["D1_form_1_minus_z_durch_L"])
L = mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)
BETA = 1 + 2j * pi / L
ALPHA = mpf("0.883")
SPRUNGFAKTOR = abs(D1) * 2 * abs(sin(pi * BETA))


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
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 85)) for j in range(N + 1)] for k in range(N + 1)])
    Ai = A.inv()
    Aim = [[mpf(Ai[i, j].mid().str(85, radius=False)) for j in range(N + 1)] for i in range(N + 1)]
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(N + 1)]
    d1 = [-sum(Aim[i][j] * F[j] for j in range(N + 1)) for i in range(N + 1)]
    return q, pot, Aim, d1


def R(q, pot, Aim, N, z):
    nenner = [z - q[0]] + [-q[k] for k in range(1, N + 1)]
    v = cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N)
    w = [sum(Aim[i][j] * v[j] for j in range(N + 1)) for i in range(N + 1)]
    return sum(abs(w[i]) * RHO ** i for i in range(N + 1)) / abs(z) ** (N + 1)


a0 = KN[0]
R0 = abs(L - a0)
gq = GaussLegendre(mp)
for N in [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "40,80,120,160")]:
    q, pot, Aim, d1 = baue(N)
    alphaN = sum(abs(d1[i]) * RHO ** i for i in range(N + 1))
    NA = mpf(N) ** ALPHA
    # (1) Skalengesetz
    zeile = []
    for v in (0.01, 0.1, 0.3, 1, 2, 4, 8, 16):
        zeile.append(f"v={v}:{mp.nstr(R(q, pot, Aim, N, L * (1 + v / NA)) / NA, 3)}")
    print(f"N={N}: α_N={mp.nstr(alphaN, 4)} | R_N(L(1+v/N^α))/N^α: " + "  ".join(zeile))
    # (2) Kontur Γ'_d
    for d in (0.05, 0.1, 0.2):
        best = mpf(0)
        bestz = None
        for k in range(72):
            phi = 2 * pi * k / 72
            # Punkt auf ∂W: Strahl von a0/2 aus; nimm den äussersten Schnitt mit beiden Kreisen, dann +d radial
            e = exp(1j * phi)
            # Schnitt mit Kreis um 0 (Radius RHO) und um a0 (Radius R0) entlang Strahl a0/2 + t e
            c0 = a0 / 2
            def t_kreis(m, r):
                # |c0 + t e - m| = r  ->  t^2 + 2t Re((c0-m) conj e) + |c0-m|^2 - r^2 = 0
                b = (c0 - m) * e.conjugate()
                bb = b.real
                disc = bb ** 2 - (abs(c0 - m) ** 2 - r ** 2)
                return -bb + mp.sqrt(disc)
            t = max(t_kreis(0, RHO), t_kreis(a0, R0)) + d
            z = c0 + t * e
            if abs(z - L) < 0.25 or abs(z - L.conjugate()) < 0.25:
                continue
            r = R(q, pot, Aim, N, z)
            if r > best:
                best, bestz = r, z
        print(f"   Γ'_d, d={d}: max R_N = {mp.nstr(best, 3)} bei ζ = {mp.nstr(bestz, 4)}")
    # (3) Schranke über den Schnitt (Leitterm), u in [0, 0.25/|L|]
    u0 = mpf("0.25") / RHO
    knoten = [mpf(0)] + [mpf(x) / NA for x in (0.3, 1, 3, 8, 20)] + [u0]
    knoten = sorted(set([k for k in knoten if k <= u0]))
    I = mpf(0)
    for a, b in zip(knoten[:-1], knoten[1:]):
        for x, wgt in gq.get_nodes(a, b, 3, mp.prec):
            I += wgt * SPRUNGFAKTOR * x ** BETA.real * R(q, pot, Aim, N, L * (1 + x)) * RHO
    I = 2 * I / (2 * pi)
    print(f"   Schranke I_N (Leitterm, Schnitt) = {mp.nstr(I, 4)};  α_N = {mp.nstr(alphaN, 4)};  Verhältnis I/α = {mp.nstr(I / alphaN, 3)}")
    sys.stdout.flush()
