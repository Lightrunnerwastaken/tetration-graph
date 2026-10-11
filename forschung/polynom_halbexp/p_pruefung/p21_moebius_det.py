"""Strang P, Prüfung 21: E Runde 3 (R3.1): det J_N = (1+τ²)^N·Q_N(τ²)/(1−τ²)^{2N+2}, Q_N irreduzibel vom Grad (N+1)(N+2)/2 + 1 (in τ²),
Q_N(0) = 2^{N+1}. Symmetrische Möbius-Drehung f = (x+τ)/(1−τx), J_N an f (sympy, N = 1…4)."""
import sys

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
x, t = sp.symbols("x tau")
f = (x + t) / (1 - t * x)
fp = sp.diff(f, x)
m = sp.simplify(fp.subs(x, f))
for N in range(1, 5):
    ser = lambda e: sp.series(e, x, 0, N + 1).removeO()
    cols = []
    mser = sp.expand(ser(m))
    fpow = [sp.Integer(1)]
    fser = sp.expand(ser(f))
    for j in range(1, N + 1):
        fpow.append(sp.expand(ser(sp.expand(fpow[-1] * fser))))
    J = sp.zeros(N + 1, N + 1)
    for j in range(N + 1):
        col = sp.expand(ser(fpow[j] + x ** j * mser))
        for k in range(N + 1):
            J[k, j] = sp.simplify(col.coeff(x, k))
    d = sp.factor(sp.simplify(J.det()))
    num, den = sp.fraction(sp.together(d))
    rest = sp.factor(num / (1 + t ** 2) ** N)
    Q = sp.factor(sp.simplify(rest))
    u = sp.symbols("u")
    Qu = sp.Poly(sp.expand(Q).subs(t ** 2, u), u) if sp.expand(Q).has(t) else None
    try:
        Qu = sp.Poly(sp.expand(Q).as_expr().subs(t, sp.sqrt(u)), u)
        deg = Qu.degree()
        irr = Qu.is_irreducible
        q0 = Qu.eval(0)
    except Exception as e:
        deg, irr, q0 = None, None, None
    print(f"N = {N}: Nenner = {sp.factor(den)}; Grad Q_N(u) = {deg} (E: {(N + 1) * (N + 2) // 2 + 1}), irreduzibel: {irr}, Q_N(0) = {q0} (E: {2 ** (N + 1)})")
