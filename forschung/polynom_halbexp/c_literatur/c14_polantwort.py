"""Strang C, Runde 2, c14: Pol-Antwortfunktion R_N(ζ) und Darstellung der Konsistenz über Cauchy-Integrale.

Exakt: B_N [(ζ−x)^{-1}]_{>N} = ζ^{-(N+1)} v_ζ,  v_ζ := [q^{N+1}/(ζ − q)]_{≤N}   (q = q_N, |ζ| > |q(0)|).
Pol-Antwort:  R_N(ζ) := |A^{-1} B_N [(ζ−x)^{-1}]_{>N}|_s = |A^{-1} v_ζ|_s / |ζ|^{N+1},   A = J_N(q_N).
Modell-Singularität u_β(x) = (1 − x/L)^β (Hauptzweig, Schnitt x = L(1+u), u ≥ 0). Für N+1 > Re β gilt
   [u_β]_{>N}(x) = (1/2πi) ∫_cut Sprung(ζ) · x^{N+1}/(ζ^{N+1}(ζ − x)) dζ,
und damit A^{-1}B_N[u_β]_{>N} = (1/2πi) ∫ Sprung(ζ) ζ^{-(N+1)} A^{-1}v_ζ dζ.  (Vorzeichen/Orientierung wird unten numerisch
an den Koeffizienten k = N+1…N+3 geprüft.)
Vergleich: δ1 = A^{-1}B_N t_N (Knesers Schwanz) mit der Antwort auf D1·u_{1+γ} + konj. (D1 aus gemeinsam/kneser_singularitaet_L.json).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc, matrix, quad, binomial, pi, exp, sin, fabs, log
from flint import arb_mat, arb, ctx

mp.dps = 70
ctx.prec = 600
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
L = mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)
GAM = 2j * pi / L
SING = json.load(open(os.path.join(G, "kneser_singularitaet_L.json"), encoding="utf-8"))


def finde_D1(d):
    for key in ("D1_form_1_minus_z_durch_L", "D1", "D_1"):
        if key in d:
            v = d[key]
            return mpc(v[0], v[1]) if isinstance(v, list) else mpc(v)
    return None


D1 = finde_D1(SING)
print("Schlüssel in kneser_singularitaet_L.json:", list(SING.keys())[:20], "D1 =", D1)


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
    """1/c als Potenzreihe bis Grad n (c[0] ≠ 0)."""
    r = [mpc(0)] * (n + 1)
    r[0] = 1 / c[0]
    for k in range(1, n + 1):
        s = sum(c[j] * r[k - j] for j in range(1, min(k, len(c) - 1) + 1))
        r[k] = -s / c[0]
    return r


def baue(N):
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 65)) for j in range(N + 1)] for k in range(N + 1)])
    Ai = A.inv()
    Aim = matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            Aim[i, j] = mpf(Ai[i, j].mid().str(65, radius=False))
    # δ1 = A^{-1} B_N t  über  B_N t = [t∘q]_{≤N}:  F = qq - exp-Jets,  δ1 = -A^{-1}F
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(N + 1)]
    d1 = [-sum(Aim[i, j] * F[j] for j in range(N + 1)) for i in range(N + 1)]
    return q, pot, Aim, d1


def v_zeta(q, pot, N, z):
    """[q^{N+1}/(ζ − q)]_{≤N}"""
    nenner = [z - q[0]] + [-q[k] for k in range(1, N + 1)]
    return cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N)


def antwort(Aim, v, N):
    return [sum(Aim[i, j] * v[j] for j in range(N + 1)) for i in range(N + 1)]


def norm_s(v, s):
    return sum(abs(v[j]) * RHO ** j * (1 + j) ** s for j in range(len(v)))


for N in (40, 80, 120, 160):
    q, pot, Aim, d1 = baue(N)
    print(f"\n=== N = {N}: |δ1|_0 = {mp.nstr(norm_s(d1, 0), 4)}, |δ1|_4 = {mp.nstr(norm_s(d1, 4), 4)}")
    # (1) R_N auf dem radialen Schnitt und auf Kreisen
    zeile = []
    for u in (1e-3, 1e-2, 3e-2, 0.1, 0.3, 1.0):
        z = L * (1 + u)
        R = norm_s(antwort(Aim, v_zeta(q, pot, N, z), N), 0) / abs(z) ** (N + 1)
        zeile.append(f"u={u:g}:{mp.nstr(R, 3)}")
    print("  R_N(L(1+u)), s=0: " + "  ".join(zeile))
    for Rr in (1.5, 1.8, 2.2, 3.0):
        zeile = []
        for grad in (0, 45, 90, 135, 180):
            z = Rr * exp(1j * pi * grad / 180)
            R = norm_s(antwort(Aim, v_zeta(q, pot, N, z), N), 0) / abs(z) ** (N + 1)
            zeile.append(f"{grad}°:{mp.nstr(R, 3)}")
        print(f"  R_N(|ζ|={Rr}), s=0: " + "  ".join(zeile))
    sys.stdout.flush()
    # (2) Antwort auf D1·u_{1+γ} + konj über das Schnittintegral; Orientierung an Koeffizienten geprüft
    beta = 1 + GAM

    def sprung(u):  # Wert oben minus unten von (1 − ζ/L)^β auf ζ = L(1+u): (−u)^β mit arg ±π
        return (u ** beta) * (exp(1j * pi * beta) - exp(-1j * pi * beta))

    # Prüfung der Darstellung am Koeffizienten k = N+1: [x^{N+1}] u_β = binomial(β, N+1)(−1/L)^{N+1}
    exakt = binomial(beta, N + 1) * (-1 / L) ** (N + 1)
    integ = quad(lambda u: sprung(u) * (L * (1 + u)) ** (-(N + 2)) * L, [0, 1.0 / N, 10.0 / N, 1, mp.inf]) / (2j * pi)
    sig = 1 if abs(integ - exakt) < abs(integ + exakt) else -1
    print(f"  Prüfung [x^(N+1)]: exakt {mp.nstr(exakt, 6)}  Integral {mp.nstr(integ, 6)}  -> Vorzeichen {sig}")

    def integrand_vec(u):
        z = L * (1 + u)
        w = antwort(Aim, v_zeta(q, pot, N, z), N)
        fak = sig * sprung(u) * z ** (-(N + 1)) * L / (2j * pi)
        return [fak * x for x in w]

    # Quadratur in u = v/N (Gauss-Legendre auf Teilintervallen, tanh-sinh wäre teurer)
    knoten = [0, 0.2, 1, 3, 8, 20, 50]
    gl = mp.gauss_legendre if hasattr(mp, "gauss_legendre") else None
    from mpmath.calculus.quadrature import GaussLegendre
    gq = GaussLegendre(mp)
    summe = [mpc(0)] * (N + 1)
    for a, b in zip(knoten[:-1], knoten[1:]):
        nodes = gq.get_nodes(mpf(a) / N, mpf(b) / N, 4, mp.prec)
        for x, wgt in nodes:
            iv = integrand_vec(x)
            summe = [s0 + wgt * y for s0, y in zip(summe, iv)]
    if D1 is not None:
        modell = [2 * (D1 * y).real for y in summe]   # D1·u_β + konj
        diff = [d - m for d, m in zip(d1, modell)]
        print(f"  Antwort auf D1·u_(1+γ)+konj: |·|_0 = {mp.nstr(norm_s(modell, 0), 4)}, |δ1 − Modell|_0 = "
              f"{mp.nstr(norm_s(diff, 0), 4)} (rel {mp.nstr(norm_s(diff, 0) / norm_s(d1, 0), 3)}); "
              f"s=4: {mp.nstr(norm_s(modell, 4), 4)} / {mp.nstr(norm_s(diff, 4), 4)}")
    sys.stdout.flush()
