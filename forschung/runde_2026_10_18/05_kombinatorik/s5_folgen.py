"""Ganzzahlige Folgen aus den Schichten (für OEIS-Suche) und Kontrolle gegen s1_schichten.json."""
import json
from pathlib import Path
import sympy as sp
HIER = Path(__file__).resolve().parent
d = json.load(open(HIER / "daten" / "s1_schichten.json"))
t, L, th, y = sp.symbols("t L th y")
h = [sp.sympify(s) for s in d["h"]]
g = [sp.sympify(s) for s in d["g_n"]]
NN = 12
def koef(k, m):
    ser = sp.series(h[k].subs({t: sp.tan(th), L: sp.log(sp.cos(th))}), th, 0, m + 1).removeO()
    return sp.expand(ser).coeff(th, m)
out = {}
# dritte Schicht: [y^{2n-3}] g_n(y-1) = 2^{-n} [th^{2n-3}] h_2
s3 = [sp.Rational(1, 2**n) * koef(2, 2 * n - 3) for n in range(2, NN + 1)]
print("[y^(2n-3)] g_n(y-1), n=2..:", s3)
# Kontrolle mit exakten g_n soweit vorhanden
for n in range(2, len(g) + 1):
    P = sp.Poly(sp.expand(g[n - 1].subs(y, y - 1)), y)
    assert P.coeff_monomial(y**(2 * n - 3)) == s3[n - 2], n
norm3 = [s3[n - 2] * sp.factorial(2 * n - 3) * 2**(n - 1) for n in range(2, NN + 1)]
print("  * (2n-3)! 2^(n-1):", norm3)
s2 = [sp.Rational(1, 2**n) * koef(1, 2 * n - 2) for n in range(1, NN + 1)]
norm2 = [s2[n - 1] * sp.factorial(2 * n - 2) * 2**(n - 1) for n in range(1, NN + 1)]
print("zweite Schicht * (2n-2)! 2^(n-1):", norm2)
out["dritte"] = [str(x) for x in s3]; out["dritte_norm"] = [str(x) for x in norm3]; out["zweite_norm"] = [str(x) for x in norm2]
json.dump(out, open(HIER / "daten" / "s5_folgen.json", "w"), indent=1)
