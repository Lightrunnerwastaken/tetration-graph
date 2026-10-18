"""Schichten der Störungsreihe G(y) = e^{G(y-1)} - 1 + u  (Etappe-1-Normierung: G(0)=0, G(-1)=log(1-u)).

H(y) = G(y-1),  u = 2 v^2,  theta = v y,  H = v * h(theta; v),  h = sum_k v^k h_k(theta).
Schichtgleichung (exakt aus H(y+1) = e^{H(y)} - 1 + u):
   sum_{m=1}^{k+1} [h^{(m)}]_{k+1-m}/m!  =  sum_{j=2}^{k+2} [h^j]_{k+2-j}/j!  + 2 delta_{k0}
   h_k(0) = -2^n/n falls k = 2n-1, sonst 0.
Ansatz: h_k ist Polynom in t = tan(theta), L = ln cos(theta), th = theta (algebraisch unabhängig).
Ableitung D = (1+t^2) d/dt - t d/dL + d/dth.  Prüfung gegen die exakten g_n (direkt iteriert).
Aufruf: python s1_schichten.py [KMAX]   Ausgabe: daten/s1_schichten.json
"""
import json, sys
from fractions import Fraction as Fr
from math import factorial
from pathlib import Path
import sympy as sp

HIER = Path(__file__).resolve().parent
(HIER / "daten").mkdir(exist_ok=True)
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
t, L, th = sp.symbols("t L th")


def D(p):
    return sp.expand((1 + t**2) * sp.diff(p, t) - t * sp.diff(p, L) + sp.diff(p, th))


def Dm(p, m):
    for _ in range(m):
        p = D(p)
    return p


h = []
for k in range(KMAX + 1):
    if k == 0:
        h.append(2 * t)
        continue
    hk = sp.Symbol("HK")
    hs = h + [hk]
    lhs_rest = sum(Dm(h[k + 1 - m], m) / sp.factorial(m) for m in range(2, k + 2))
    rhs = 0
    for j in range(2, k + 3):
        w = k + 2 - j
        # [h^j]_w : Summe über Zerlegungen von w in j Teile
        poly_v = sum(h[i] * sp.Symbol("v")**i for i in range(len(h))) + hk * sp.Symbol("v")**k
        rhs += sp.expand(poly_v**j).coeff(sp.Symbol("v"), w) / sp.factorial(j)
    rhs = sp.expand(rhs)
    R = sp.expand(rhs - lhs_rest)          # D(HK) = R  (R enthält HK linear: 2t*HK)
    R0 = sp.expand(R.subs(hk, 0))
    assert sp.expand(R.coeff(hk) - 2 * t) == 0
    # Ansatz
    deg_t, deg_L, deg_th = k + 3, k + 1, k + 1
    cs, ans = [], 0
    for a in range(deg_t + 1):
        for b in range(deg_L + 1):
            for c in range(deg_th + 1):
                s = sp.Symbol(f"c_{a}_{b}_{c}")
                cs.append(s)
                ans += s * t**a * L**b * th**c
    eq = sp.expand(D(ans) - 2 * t * ans - R0)
    eqs = sp.Poly(eq, t, L, th).coeffs()
    init = (-sp.Rational(2**((k + 1) // 2), (k + 1) // 2)) if k % 2 == 1 else 0
    eqs.append(ans.subs({t: 0, L: 0, th: 0}) - init)
    sol = sp.solve(eqs, cs, dict=True)
    assert sol, f"Ansatz für k={k} unlösbar"
    hk_sol = sp.expand(ans.subs(sol[0]).subs({s: 0 for s in cs}))
    h.append(hk_sol)
    print(f"h_{k} =", sp.factor_terms(hk_sol), flush=True)

# ---------------- Gegenprobe gegen exakte g_n(y)
N = (KMAX + 1) // 2 + 1
u = sp.Symbol("u")
# G(k) = f^k(0), f = u + e^x - 1, Reihe in u bis u^N, exakt; g_n Polynom Grad 2n-1 in y aus Stützstellen k = 0..2N
def ser_exp_minus1(a):  # a: Liste Fr, a[0]=0
    res = [Fr(0)] * (N + 1)
    p = [Fr(1)] + [Fr(0)] * N
    for m in range(1, N + 1):
        p = [sum(p[i] * a[j - i] for i in range(j + 1)) for j in range(N + 1)]
        res = [r + q / factorial(m) for r, q in zip(res, p)]
    return res
vals = []
G = [Fr(0)] * (N + 1)
for k in range(2 * N + 1):
    vals.append(G)
    G = ser_exp_minus1(G)
    G[1] += 1
y = sp.Symbol("y")
g = [None]
for n in range(1, N + 1):
    pts = [(k, sp.Rational(vals[k][n].numerator, vals[k][n].denominator)) for k in range(2 * n)]
    g.append(sp.expand(sp.interpolate(pts, y)))
ok = True
for n in range(1, N + 1):
    P = sp.Poly(sp.expand(g[n].subs(y, y - 1)), y)
    for m in range(0, 2 * n):
        k = 2 * n - 1 - m
        if k > KMAX:
            continue
        # [theta^m] h_k
        ser = sp.series(h[k].subs({t: sp.tan(th), L: sp.log(sp.cos(th))}), th, 0, m + 1).removeO()
        pred = sp.Rational(1, 2**n) * sp.expand(ser).coeff(th, m)
        if sp.simplify(pred - P.coeff_monomial(y**m)) != 0:
            ok = False
            print("FEHLER", n, m, pred, P.coeff_monomial(y**m))
print("Gegenprobe gegen exakte g_n (alle Koeffizienten mit Schicht <= KMAX):", ok)
json.dump({"KMAX": KMAX, "h": [str(x) for x in h], "gegenprobe": ok,
           "g_n": [str(x) for x in g[1:]]}, open(HIER / "daten" / "s1_schichten.json", "w"), indent=1)
