"""Strang C, c09 (Variante von c05 mit Gewichten |L|^k (1+k)^s, s = 0, 3, 4; N bis 194): Wie gross ist die nichtlineare Konstante ENTLANG des tatsächlichen Newton-Zuwachses?

delta = p_N - q_N (wahrer Abstand), Norm l^1 mit w_k = |L|^k.
  omega_eff = |J(q_N)^{-1} F''(q_N)[delta, delta]|_w / |delta|_w^2
  kappa     = |q_N - J(q_N)^{-1}F(q_N) - p_N|_w / |delta|_w^2      (beobachtete quadratische Konstante des 1. Newton-Schritts)
  h_eff     = alpha * omega_eff  (alpha = |J(q_N)^{-1}F(q_N)|_w)
Vergleich mit dem Worst-Case-omega aus c02/c03. Diagnose, kein Beweis.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf
from flint import arb_mat, arb, ctx

mp.dps = 140
ctx.prec = 520
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
    r = [mpf(0)] * (n + 1)
    for j, cj in enumerate(c[: n + 1]):
        if cj == 0:
            continue
        for k in range(n + 1):
            r[k] += cj * pot[j][k]
    return r


def abl(c):
    return [k * c[k] for k in range(1, len(c))] or [mpf(0)]


def arbv(v):
    return arb_mat([[arb(mp.nstr(x, 135))] for x in v])


L_abs = mpf("1.3745570107437075")
for N in (40, 80, 120, 160, 194):
    q = KN[: N + 1]
    p = [mpf(x) for x in P[str(N)]]
    n = N
    pot = potenzen(q, n)
    dq = abl(q)
    dpp = kompon(dq, pot, n)
    ddpp = kompon(abl(dq), pot, n)
    J = [[pot[j][k] + (dpp[k - j] if j <= k else mpf(0)) for j in range(n + 1)] for k in range(n + 1)]
    F = [kompon(q, pot, n)[k] - mpf(1) / mp.factorial(k) for k in range(n + 1)]
    d = [p[k] - q[k] for k in range(n + 1)]
    # F''[d,d] = 2 (d'∘q)·d + (q''∘q)·d^2
    dd = kompon(abl(d), pot, n)
    t1 = mul(dd, d, n)
    t2 = mul(ddpp, mul(d, d, n), n)
    B = [2 * t1[k] + t2[k] for k in range(n + 1)]
    Ja = arb_mat([[arb(mp.nstr(x, 135)) for x in z] for z in J])
    Ji = Ja.inv()
    sF = Ji * arbv(F)
    sB = Ji * arbv(B)
    for s in (0, 3, 4):
      w = [L_abs ** k * (1 + k) ** s for k in range(n + 1)]

      def nrm(v):
          return sum(abs(mpf(v[k])) * w[k] for k in range(n + 1))

      sFl = [mpf(sF[k, 0].mid().str(130, radius=False)) for k in range(n + 1)]
      sBl = [mpf(sB[k, 0].mid().str(130, radius=False)) for k in range(n + 1)]
      nd = nrm(d)
      alpha = nrm(sFl)
      om_eff = nrm(sBl) / nd ** 2
      newton = [q[k] - sFl[k] - p[k] for k in range(n + 1)]
      kappa = nrm(newton) / nd ** 2
      print(f"N={N:3d} s={s}: |delta|={float(nd):.3e} alpha={float(alpha):.3e} omega_eff={float(om_eff):.3e} "
            f"h_eff={float(alpha*om_eff):.3e} kappa={float(kappa):.3e} |q-J^-1F-p|={float(nrm(newton)):.3e}")
      sys.stdout.flush()
