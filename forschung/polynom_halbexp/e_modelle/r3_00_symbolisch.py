"""Runde 3, Schritt 0 (Reproduktion): (a) det J_N der Möbius-Drehung symbolisch in t = τ = tan(c/2) faktorisieren;
(b) Vorzeichen von det J_N(c) auf einem Gitter c ∈ (0, π/3) (Resonanzen?).  Aufruf: python r3_00_symbolisch.py NMAX_symb [NMAX_gitter]"""
import sys
import sympy as sp
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 5
t, x = sp.symbols("t x")
f = (x + t) / (1 - t * x)
m = (1 + t ** 2) * (1 - t * x) ** 2 / (1 - t ** 2 - 2 * t * x) ** 2
for N in range(1, NS + 1):
    J = sp.zeros(N + 1, N + 1)
    for j in range(N + 1):
        col = sp.series(f ** j + m * x ** j, x, 0, N + 1).removeO()
        for k in range(N + 1):
            J[k, j] = sp.simplify(col.coeff(x, k))
    d = sp.factor(sp.simplify(J.det()))
    print(f"N = {N}: det J_N = {d}", flush=True)
if len(sys.argv) > 2:
    from moebius_lin import setup
    cs = [mp.mpf(k) / 200 * mp.pi / 3 for k in range(1, 200)]
    for N in range(10, int(sys.argv[2]) + 1, 10):
        s = []
        for c in cs:
            d = setup(c, N, 500)["J"].det()
            s.append(1 if mp.mpf(d.mid().str(30, radius=False)) > 0 else -1)
        print(f"N = {N}: Vorzeichenwechsel bei c = {[mp.nstr(cs[i], 4) for i in range(len(s) - 1) if s[i] != s[i + 1]]}", flush=True)
