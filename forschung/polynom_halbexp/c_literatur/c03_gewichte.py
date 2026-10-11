"""Strang C, c03 (Variante von c02 mit Gewichten r^k (1+k)^s und Argmax): Prüfrechnung: affin-kovariante Kantorovich-Grössen (Deuflhard–Heindl 1979 / Kantorovich-Originalform)
an der Stelle x0 = q_N (Knesers Taylor-Polynom vom Grad N), in l^1-Normen mit GLEICHEN Gewichten w_k = r^k in X und Y.

  alpha_N = |J(q_N)^{-1} F_N(q_N)|_w                      (erster Newton-Schritt)
  omega_N = max_{i,j} |J(q_N)^{-1} F''(q_N)[e_i, e_j]|_w / (w_i w_j)   (affin-kovariante Lipschitz-Konstante, lokal)
  beta_N  = |J(q_N)^{-1}|_{w->w},   gamma_N = max_{i,j} |F''(q_N)[e_i,e_j]|_w/(w_i w_j)   (Produktform zum Vergleich)
  h_aff = alpha*omega,  h_prod = alpha*beta*gamma,  d_N = |p_N - q_N|_w (wahrer Abstand).

Hinweis: omega ist hier nur die Ableitung an x0 (nicht das Supremum über eine Kugel) — Diagnose, kein Beweis.
F''(p)[d,e] = d'∘p·e + e'∘p·d + p''∘p·d·e (Jets bis N).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf
from flint import arb_mat, arb, ctx

mp.dps = 110
ctx.prec = 420
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
P = json.load(open(os.path.join(G, "polynome_bis_194.json"), encoding="utf-8"))
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_300.json"), encoding="utf-8"))["a"]]


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def potenzen(a, n):
    pot = [[mpf(1)] + [mpf(0)] * n]
    for _ in range(1, n + 1):
        pot.append(mul(pot[-1], a, n))
    return pot


def kompon(c, pot, n):
    """jet_n( c ∘ p ) für Koeffizientenliste c (Länge ≤ n+1)."""
    r = [mpf(0)] * (n + 1)
    for j, cj in enumerate(c[: n + 1]):
        if cj == 0:
            continue
        for k in range(n + 1):
            r[k] += cj * pot[j][k]
    return r


def shift(v, s, n):
    return ([mpf(0)] * s + v)[: n + 1] if s <= n else [mpf(0)] * (n + 1)


def daten(a):
    n = len(a) - 1
    pot = potenzen(a, n)
    dp = [k * a[k] for k in range(1, n + 1)]
    ddp = [k * (k - 1) * a[k] for k in range(2, n + 1)]
    dpp = kompon(dp, pot, n)
    ddpp = kompon(ddp, pot, n)
    J = [[pot[j][k] + (dpp[k - j] if j <= k else mpf(0)) for j in range(n + 1)] for k in range(n + 1)]
    pp = kompon(a, pot, n)
    F = [pp[k] - mpf(1) / mp.factorial(k) for k in range(n + 1)]
    # F''[e_i, e_j]
    paare, spalten = [], []
    for i in range(n + 1):
        for j in range(i, n + 1):
            v = [mpf(0)] * (n + 1)
            if i >= 1:
                t = shift(pot[i - 1], j, n)
                v = [x + i * y for x, y in zip(v, t)]
            if j >= 1:
                t = shift(pot[j - 1], i, n)
                v = [x + j * y for x, y in zip(v, t)]
            t = shift(ddpp, i + j, n)
            v = [x + y for x, y in zip(v, t)]
            paare.append((i, j))
            spalten.append(v)
    return J, F, paare, spalten


def zu_arb(M):
    return arb_mat([[arb(mp.nstr(x, 105)) for x in zeile] for zeile in M])


def norm_w(vec, w):
    return sum(abs(float(vec[k].mid())) * w[k] for k in range(len(w)))


L_abs = 1.3745570107437075
for N in (30, 50, 70):
    q = KN[: N + 1]
    p = [mpf(x) for x in P[str(N)]]
    J, F, paare, spalten = daten(q)
    Ja = zu_arb(J)
    Ji = Ja.inv()
    dvec = Ji * arb_mat([[arb(mp.nstr(x, 105))] for x in F])
    Bm = arb_mat([[arb(mp.nstr(spalten[c][k], 105)) for c in range(len(spalten))] for k in range(N + 1)])
    JiB = Ji * Bm
    for (r, s) in [(L_abs, 0), (L_abs, 3), (L_abs, 6), (1.45, 0), (1.55, 0), (1.7, 0), (1.55, 3)]:
        w = [r ** k * (1 + k) ** s for k in range(N + 1)]
        alpha = sum(abs(float(dvec[k, 0].mid())) * w[k] for k in range(N + 1))
        dist = sum(abs(float(p[k] - q[k])) * w[k] for k in range(N + 1))
        omega = gamma = 0.0
        for c, (i, j) in enumerate(paare):
            om = sum(abs(float(JiB[k, c].mid())) * w[k] for k in range(N + 1)) / (w[i] * w[j])
            ga = sum(abs(float(spalten[c][k])) * w[k] for k in range(N + 1)) / (w[i] * w[j])
            if om > omega:
                arg = (i, j)
            omega, gamma = max(omega, om), max(gamma, ga)
        beta = max(sum(abs(float(Ji[k, j].mid())) * w[k] for k in range(N + 1)) / w[j] for j in range(N + 1))
        print(f"N={N:3d} r={r:.3f} s={s} argmax={arg}: d={dist:.2e} alpha={alpha:.2e} omega={omega:.2e} h_aff={alpha*omega:.2e} | "
              f"beta={beta:.2e} gamma={gamma:.2e} h_prod={alpha*beta*gamma:.2e}")
    sys.stdout.flush()
