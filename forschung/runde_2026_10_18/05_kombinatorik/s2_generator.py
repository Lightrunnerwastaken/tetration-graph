"""Erzeuger lambda_u der Abbildung f_u(x) = e^x - 1 + u (formal, Grad: x -> 1, u -> 2).

(a) lambda_u aus exp(lambda d/dx) x = f_u, Grad für Grad (exakte Brüche, dict-Polynome).
(b) Fixpunkte x_pm (e^x - 1 - x = -u) als Reihe in v (u = 2v^2); Probe lambda_u'(x_pm) = x_pm.
(d) Frequenz omega = |x_+|^2/(2 Im x_+), gamma/beta  (Vergleich mit Säkularterm in h_4, s1_schichten).
(e) Bernoulli-Formel für den u-linearen Teil: mu = lambda0 * sum_k B_k^+ (lambda0 d/dx)^k (1/lambda0)/k!.
Aufruf: python s2_generator.py [DMAX]   Ausgabe daten/s2_generator.json
"""
import json, sys
from fractions import Fraction as Fr
from math import factorial
from pathlib import Path
import sympy as sp

HIER = Path(__file__).resolve().parent
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14


def deg(m):
    return m[0] + 2 * m[1]


def mul(p, q):
    r = {}
    for (a, b), c in p.items():
        for (a2, b2), c2 in q.items():
            if a + a2 + 2 * (b + b2) <= DMAX:
                k = (a + a2, b + b2)
                r[k] = r.get(k, 0) + c * c2
    return {k: c for k, c in r.items() if c}


def dx(p):
    return {(a - 1, b): c * a for (a, b), c in p.items() if a > 0}


def add(p, q, s=1):
    r = dict(p)
    for k, c in q.items():
        r[k] = r.get(k, 0) + s * c
    return {k: c for k, c in r.items() if c}


def part(p, d):
    return {k: c for k, c in p.items() if deg(k) == d}


f = {(0, 1): Fr(1)}
for k in range(2, DMAX + 1):
    f[(k, 0)] = Fr(1, factorial(k))
lam = {}
for d in range(2, DMAX + 1):
    acc, term = {}, None
    for m in range(1, d):
        term = dict(lam) if m == 1 else mul(lam, dx(term))
        if m >= 2:
            acc = add(acc, {k: c / factorial(m) for k, c in part(term, d).items()})
    lam = add(lam, add(part(f, d), acc, -1))
# Kontrolle
acc, term = {}, None
for m in range(1, DMAX):
    term = dict(lam) if m == 1 else mul(lam, dx(term))
    acc = add(acc, {k: c / factorial(m) for k, c in term.items()})
print("Lie-Kontrolle exp(lambda d)x - x = f_u - x bis Grad", DMAX, ":", add(acc, f, -1) == {})
x, u, v = sp.symbols("x u v")
for d in range(2, 8):
    print(f"lambda_{d} =", sum(sp.Rational(c.numerator, c.denominator) * x**a * u**b for (a, b), c in part(lam, d).items()))
lam0 = [lam.get((k, 0), Fr(0)) for k in range(DMAX + 1)]
lamu0 = [lam.get((0, b), Fr(0)) for b in range(1, DMAX // 2 + 1)]
lam1 = [lam.get((k, 1), Fr(0)) for k in range(DMAX - 1)]
print("lambda_0 (Julia von e^x-1):", [str(c) for c in lam0[2:]])
print("lambda_u(0):", [str(c) for c in lamu0])
print("mu = d lambda/du |_{u=0}:", [str(c) for c in lam1])

# ---- (b),(d) Fixpunkte
# Gauss-rationale Reihen in v: Liste von (re, im)
K = 12
def gm(a, b):
    c = [(Fr(0), Fr(0))] * (K + 1)
    for i, (ar, ai) in enumerate(a):
        for j in range(K + 1 - i):
            br, bi = b[j]
            cr, ci = c[i + j]
            c[i + j] = (cr + ar * br - ai * bi, ci + ar * bi + ai * br)
    return c
Z = [(Fr(0), Fr(0))] * (K + 1)
xp = list(Z); xp[0] = (Fr(0), Fr(2))
for k in range(1, K + 1):
    # F(xi, v) = sum_{j>=2} v^{j-2} xi^j / j! + 2
    F = list(Z); F[0] = (Fr(2), Fr(0))
    pw = gm(xp, xp)
    for j in range(2, K + 3):
        for i in range(K + 1 - (j - 2)):
            r, im = pw[i]
            fr, fi = F[i + j - 2]
            F[i + j - 2] = (fr + r / factorial(j), fi + im / factorial(j))
        pw = gm(pw, xp)
    fr, fi = F[k]          # c_k = -F_k/(2i) = (i/2) F_k = (-fi/2, fr/2)
    xp[k] = (-fi / 2, fr / 2)
print("xi_+ (re, im) je v^k:", [(str(r), str(i)) for r, i in xp[:8]])
# Probe lambda'(x_+) = x_+ : lambda'(x) = sum a c x^{a-1} u^b, x = v xi, u = 2 v^2
lp = list(Z)
pw = [list(Z) for _ in range(DMAX + 1)]
pw[0] = list(Z); pw[0][0] = (Fr(1), Fr(0))
for a_ in range(1, DMAX + 1):
    pw[a_] = gm(pw[a_ - 1], xp)
for (a_, b_), c in lam.items():
    if a_ == 0:
        continue
    sh = a_ - 1 + 2 * b_
    for i in range(K + 1 - sh):
        r, im = pw[a_ - 1][i]
        lr, li = lp[i + sh]
        lp[i + sh] = (lr + c * a_ * 2**b_ * r, li + c * a_ * 2**b_ * im)
KK = min(K, DMAX - 2)
pr = all(lp[i] == (xp[i - 1] if i >= 1 else (Fr(0), Fr(0))) for i in range(KK))
print("lambda'(x_+) = x_+ bis v^%d:" % (KK - 1), pr)
gam = [r for r, i in xp]; bet = [i for r, i in xp]
vv = sp.Symbol("v")
B = sum(sp.Rational(c.numerator, c.denominator) * vv**k for k, c in enumerate(bet))
G_ = sum(sp.Rational(c.numerator, c.denominator) * vv**k for k, c in enumerate(gam))
print("beta/v =", B); print("gamma/v =", G_)
om = sp.expand(sp.series((B**2 + G_**2) / (2 * B), vv, 0, 11).removeO())
gb = sp.expand(sp.series(G_ / B, vv, 0, 11).removeO())
print("omega/v =", om)
print("gamma/beta =", gb)
# ---- (e) Bernoulli-Formel (u-linear), Laurentreihen in x: Liste mit Exponent-Versatz
N = DMAX - 2
def lmul(a, ea, b, eb, n=N):  # a = sum a_i x^{ea+i}
    c = [Fr(0)] * n
    for i, ai in enumerate(a):
        if ai:
            for j in range(n - i):
                if j < len(b):
                    c[i + j] += ai * b[j]
    return c, ea + eb
l0 = lam0[2:2 + N]                                   # lambda0 = x^2 * l0
inv = [Fr(0)] * N; inv[0] = 1 / l0[0]
for n in range(1, N):
    inv[n] = -sum(l0[k] * inv[n - k] for k in range(1, n + 1)) / l0[0]
term, e = inv, -2                                    # 1/lambda0 = x^{-2} * inv
acc = [Fr(0)] * N                                    # Exponent-Versatz -2
for k in range(0, N + 2):
    Bk = Fr(1, 2) if k == 1 else Fr(sp.bernoulli(k).p, sp.bernoulli(k).q)
    # term hat Versatz e >= -2 ; in acc einsortieren
    for i, c in enumerate(term):
        if e + i + 2 < N:
            acc[e + i + 2] += Bk * c / factorial(k)
    # D term = lambda0 * d/dx term
    dterm = [c * (e + i) for i, c in enumerate(term)]   # Exponent e-1
    term, e = lmul(l0, 2, dterm, e - 1)
    term = term[:N]
mu_f, em = lmul(l0, 2, acc, -2)                       # mu = lambda0 * acc, Versatz 0
mu_f = mu_f[:N - 2]
print("mu (Bernoulli-Formel):", [str(c) for c in mu_f])
okB = all(mu_f[i] == (lam1[i] if i < len(lam1) else 0) for i in range(len(mu_f)))
print("Bernoulli-Formel = mu bis x^%d:" % (len(mu_f) - 1), okB)
json.dump({"DMAX": DMAX, "lambda": {f"{a},{b}": str(c) for (a, b), c in sorted(lam.items())},
           "lambda0": [str(c) for c in lam0], "lambda_u_null": [str(c) for c in lamu0], "mu": [str(c) for c in lam1],
           "xi_plus": str(xp), "omega_durch_v": str(om), "gamma_durch_beta": str(gb), "bernoulli_ok": okB,
           "residuum_ok": pr}, open(HIER / "daten" / "s2_generator.json", "w"), indent=1)
