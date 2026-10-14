"""Etappe 1, Teil B: exakte Koeffizienten (Brüche) von Iterationen parabolischer Abbildungen und der Störungsreihe.

Ausgabe: daten/b_kombinatorik.json mit allen ganzzahligen Folgen (normiert) für die OEIS-Suche,
plus Selbstprüfungen der bewiesenen Identitäten (Abschnitte mit 'BEWEIS' im Bericht).
"""
import json
from fractions import Fraction as Fr
from math import comb, factorial
from pathlib import Path

import sympy as sp

HIER = Path(__file__).resolve().parent
N = 16                       # Reihen bis x^N
out = {}
pruef = []


def check(name, ok, wert=None):
    pruef.append({"name": name, "ok": bool(ok), "wert": wert})
    print(("OK  " if ok else "FAIL") + " " + name + ("" if wert is None else f"  [{wert}]"), flush=True)


# ---------------------------------------------------------------- Potenzreihen (Liste a[0..N], exakt)
def mul(a, b):
    c = [Fr(0)] * (N + 1)
    for i, ai in enumerate(a):
        if ai:
            for j in range(N + 1 - i):
                c[i + j] += ai * b[j]
    return c


def comp(f, g):
    """f(g(x)), g(0) = 0."""
    assert g[0] == 0
    r = [Fr(0)] * (N + 1)
    p = [Fr(1)] + [Fr(0)] * N
    for k in range(N + 1):
        if f[k]:
            r = [ri + f[k] * pi for ri, pi in zip(r, p)]
        p = mul(p, g)
    return r


def expm1_ser():
    return [Fr(0)] + [Fr(1, factorial(n)) for n in range(1, N + 1)]


def xexp_ser():
    return [Fr(0)] + [Fr(1, factorial(n - 1)) for n in range(1, N + 1)]


X = [Fr(0), Fr(1)] + [Fr(0)] * (N - 1)


def iterate(f, k):
    r = X
    for _ in range(k):
        r = comp(f, r)
    return r


t = sp.Symbol("t")


def frac_iter_polys(f):
    """c_n(t) = [x^n] f^{∘t}: Polynom in t vom Grad ≤ n−1 (Interpolation durch t = 0..N)."""
    its = [iterate(f, k) for k in range(N + 1)]
    polys = []
    for n in range(N + 1):
        pts = [(k, sp.Rational(its[k][n].numerator, its[k][n].denominator)) for k in range(n + 1)]
        polys.append(sp.expand(sp.interpolate(pts, t)) if n >= 1 else sp.Integer(0))
    return its, polys


def julia(polys):
    """Iterativer Logarithmus λ(x) = d/dt f^t |_{t=0}."""
    return [sp.diff(p, t).subs(t, 0) for p in polys]


def ints(xs):
    return [int(v) for v in xs]


for name, f in (("expm1", expm1_ser()), ("xexp", xexp_ser())):
    its, polys = frac_iter_polys(f)
    d = {}
    # ganze Iterationen: n! [x^n] f^k
    for k in (1, 2, 3, 4):
        d[f"iter{k}_nfak"] = [int(its[k][n] * factorial(n)) for n in range(1, N + 1)]
    # Halbiteration und −1-Iteration (Umkehrfunktion)
    half = [p.subs(t, sp.Rational(1, 2)) for p in polys]
    d["halb_koef"] = [str(v) for v in half[1:]]
    d["halb_nenner"] = [int(sp.fraction(v)[1]) for v in half[1:]]
    d["halb_zaehler"] = [int(sp.fraction(v)[0]) for v in half[1:]]
    d["halb_mal_nfak_2hoch"] = [str(half[n] * factorial(n) * 2 ** (n - 1)) for n in range(1, N + 1)]
    # Julia-Funktion
    lam = julia(polys)
    d["julia_koef"] = [str(v) for v in lam[1:]]
    d["julia_zaehler"] = [int(sp.fraction(v)[0]) for v in lam[1:]]
    d["julia_nenner"] = [int(sp.fraction(v)[1]) for v in lam[1:]]
    d["julia_nfak"] = [str(v * factorial(n)) for n, v in enumerate(lam) if n >= 1]
    # Dreieck: n! c_n(t) = Σ_k T(n,k) t^k ; ganzzahlig?
    tri = []
    for n in range(1, N + 1):
        P = sp.Poly(polys[n] * factorial(n), t)
        tri.append([str(P.coeff_monomial(t ** k)) for k in range(n)])
    d["dreieck_nfak_t"] = tri
    # Dreieck in Binomialbasis: n! c_n(t) = Σ_k B(n,k) C(t,k)  (ganzzahlig, weil c_n(k) für ganze k ganzzahlig·/n!)
    bino = []
    for n in range(1, N + 1):
        vals = [its[k][n] * factorial(n) for k in range(n)]
        # Newton-Vorwärtsdifferenzen
        row = []
        cur = list(vals)
        for k in range(n):
            row.append(int(cur[0]))
            cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
        bino.append(row)
    d["dreieck_binomialbasis"] = bino
    d["leit_koef_t^(n-1)_mal_nfak"] = [tri[n - 1][n - 1] for n in range(1, N + 1)]
    out[name] = d
    print(name, "Julia:", d["julia_koef"][:8])
    print(name, "binomial rows:", bino[:6])

# Prüfung: Julia-Gleichung λ(f(x)) = f'(x) λ(x)
f = expm1_ser()
lam = [Fr(int(sp.fraction(v)[0]), int(sp.fraction(v)[1])) for v in [0] + [sp.Rational(s) for s in out["expm1"]["julia_koef"]]]
lam = [Fr(0)] + lam[1:]
lhs = comp(lam, f)
fp = [Fr(n + 1) * f[n + 1] if n + 1 <= N else Fr(0) for n in range(N + 1)]
rhs = mul(fp, lam)
check("Julia-Gleichung λ(eˣ−1) = eˣ λ(x) exakt bis x^%d" % (N - 1), lhs[:N] == rhs[:N])

# Bekanntes: e^{e^x−1}−1 → Bell-Zahlen
bell = [int(sp.bell(n)) for n in range(1, N + 1)]
check("n![xⁿ](eˣ−1)∘(eˣ−1) = Bell(n)", out["expm1"]["iter2_nfak"] == bell)


# ---------------------------------------------------------------- Störungsreihe g_n(y)  (aus e81, unabhängig neu)
y = sp.Symbol("y")
NMAX = 12


QQ = sp.QQ
BERN = {}


def faulhaber(q):
    """S mit S(y) − S(y−1) = q(y−1), S(0) = 0 (wie e81)."""
    s = sp.Poly(0, y, domain=QQ)
    for (m,), c in q.terms():
        if m not in BERN:
            BERN[m] = sp.Poly((sp.bernoulli(m + 1, y) - sp.bernoulli(m + 1, 0)) / (m + 1), y, domain=QQ)
        s += BERN[m] * c
    return s


def verschoben(p, a):
    return sp.Poly(p.as_expr().subs(y, y + a), y, domain=QQ)


def exp_koef(gs, n):
    E = [sp.Poly(1, y, domain=QQ)]
    for k in range(1, n + 1):
        E.append(sum((gs[j] * E[k - j] * j for j in range(1, k + 1) if j < len(gs)), sp.Poly(0, y, domain=QQ)) * QQ(1, k))
    return E


gP = [None, sp.Poly(y, y, domain=QQ)]
for n in range(2, NMAX + 1):
    gm = [None] + [verschoben(gk, -1) for gk in gP[1:]]
    q = verschoben(exp_koef(gm, n)[n], 1)
    gP.append(faulhaber(q))
    print("g_%d fertig" % n, flush=True)
g = [None] + [p.as_expr() for p in gP[1:]]

st = {}
lead = [sp.Poly(g[n], y).LC() for n in range(1, NMAX + 1)]
st["leitkoef"] = [str(v) for v in lead]
# Vermutung/Beweis: Leitkoeffizient = 2^{1−n} T_n/(2n−1)!, T_n Tangentenzahlen
tang = [sp.Integer(abs(sp.euler(2 * n - 1, 0) if False else 0)) for n in range(1, NMAX + 1)]  # Platzhalter


def tangentenzahl(n):   # T_n = (2n−1)! [z^{2n−1}] tan z
    return sp.series(sp.tan(sp.Symbol("z")), sp.Symbol("z"), 0, 2 * n + 1).removeO().coeff(sp.Symbol("z"), 2 * n - 1) * sp.factorial(2 * n - 1)


T = [tangentenzahl(n) for n in range(1, NMAX + 1)]
pred = [sp.Integer(2) ** (1 - n) * T[n - 1] / sp.factorial(2 * n - 1) for n in range(1, NMAX + 1)]
check("Leitkoeffizient von g_n = 2^{1−n}·T_n/(2n−1)! (Tangentenzahlen), n ≤ %d" % NMAX, lead == pred, [str(v) for v in lead[:5]])
# zweiter Koeffizient (y^{2n−2})
sec = [sp.Poly(g[n], y).coeff_monomial(y ** (2 * n - 2)) for n in range(1, NMAX + 1)]
st["zweitkoef"] = [str(v) for v in sec]
st["zweit_durch_leit"] = [str(sp.nsimplify(sec[i] / lead[i])) for i in range(1, NMAX)]
# Werte an ganzen y: n! g_n(y)
for yy in range(1, 6):
    st[f"nfak_g_bei_y={yy}"] = [str(g[n].subs(y, yy) * sp.factorial(n)) for n in range(1, NMAX + 1)]
for yy in (sp.Rational(1, 2), sp.Rational(-1, 2), -1, sp.Rational(3, 2)):
    st[f"g_bei_y={yy}"] = [str(g[n].subs(y, yy)) for n in range(1, NMAX + 1)]
    st[f"g_nenner_bei_y={yy}"] = [int(sp.fraction(sp.nsimplify(g[n].subs(y, yy)))[1]) for n in range(1, NMAX + 1)]
# ganze Werte: G(3) = e^{e^u−1+u}−1+u → n! g_n(3) = B_{n+1} − B_n ... prüfen
B = [sp.bell(n) for n in range(0, NMAX + 3)]
check("n! g_n(2) = 1 für n ≥ 2 (G(2) = e^u − 1 + u)", all(g[n].subs(y, 2) * sp.factorial(n) == 1 for n in range(2, NMAX + 1)))
check("n! g_n(3) = Bell(n+1) für n ≥ 2 (G(3) = e^{e^u−1+u} − 1 + u)", all(g[n].subs(y, 3) * sp.factorial(n) == B[n + 1] for n in range(2, NMAX + 1)),
      st["nfak_g_bei_y=3"][:8])
# Koeffizienten-Dreieck von (2n−1)!·2^{n−1}·g_n(y)? Nenner untersuchen
dreieck = []
for n in range(1, NMAX + 1):
    P = sp.Poly(g[n], y)
    dreieck.append([str(P.coeff_monomial(y ** k)) for k in range(2 * n)])
st["dreieck"] = dreieck
st["nenner_lcm"] = [int(sp.ilcm(*[sp.fraction(sp.Rational(c))[1] for c in row])) for row in dreieck]
# g_n(−1): G(−1) aus G(0)=0 = e^{G(−1)} − 1 + u ⇒ G(−1) = log(1−u)  ⇒ g_n(−1) = −1/n
check("g_n(−1) = −1/n (G(−1) = log(1−u))", all(g[n].subs(y, -1) == sp.Rational(-1, n) for n in range(1, NMAX + 1)))
# Nullstellen: g_n(0) = 0, g_n(1) = 0 (n≥2)
check("g_n(0) = 0 und g_n(1) = 0 für n ≥ 2", all(g[n].subs(y, 0) == 0 and g[n].subs(y, 1) == 0 for n in range(2, NMAX + 1)))
out["stoerungsreihe"] = st
print("Störungsreihe fertig", flush=True)

# ---------------------------------------------------------------- Türme in x um x = 1: ^k(1+h), n![h^n]
def ser_exp(a):          # exp einer Reihe mit a[0] = 0
    r = [Fr(0)] * (N + 1)
    r[0] = Fr(1)
    ap = [Fr(k) * a[k] for k in range(N + 1)]
    for m in range(1, N + 1):
        r[m] = sum(ap[k] * r[m - k] for k in range(1, m + 1)) / m
    return r


ln1h = [Fr(0)] + [Fr((-1) ** (k + 1), k) for k in range(1, N + 1)]
tower = {}
cur = [Fr(1), Fr(1)] + [Fr(0)] * (N - 1)        # ^1(1+h) = 1 + h
for k in range(1, 7):
    tower[f"turm{k}_nfak"] = [int(cur[n] * factorial(n)) for n in range(N + 1)]
    cur = ser_exp(mul(ln1h, cur))             # (1+h)^{cur}
out["tuerme"] = tower
print("tuerme:", {k: v[:9] for k, v in tower.items()}, flush=True)

(HIER / "daten").mkdir(exist_ok=True)
(HIER / "daten" / "b_kombinatorik.json").write_text(json.dumps({"daten": out, "pruefungen": pruef}, indent=1, ensure_ascii=False), encoding="utf-8")
print("geschrieben")
