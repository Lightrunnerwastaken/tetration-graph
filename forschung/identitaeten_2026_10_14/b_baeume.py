"""n!·g_n(y) in der Binomialbasis: n! g_n(y) = Σ_j N(n,j)·C(y,j).
Deutung: G(k) = f^{∘k}(0) mit f(x) = u + eˣ − 1 zählt Bäume mit n markierten Blättern (u) und unmarkierten inneren Knoten
(jeder innere Knoten: nichtleere Menge von Kindern, auch genau eins), Höhe ≤ k. N(n,j) = Zahl der Bäume, bei denen genau j
Stufen "echt" sind (Newton-Differenzen). Ausgabe: Dreieck und Zeilensummen für die OEIS-Suche."""
import json
from math import comb, factorial
from pathlib import Path
import sympy as sp

HIER = Path(__file__).resolve().parent
d = json.loads((HIER / "daten" / "b_kombinatorik.json").read_text(encoding="utf-8"))["daten"]["stoerungsreihe"]["dreieck"]
y = sp.Symbol("y")
tri = []
for n, row in enumerate(d, start=1):
    p = sum(sp.Rational(c) * y ** k for k, c in enumerate(row)) * factorial(n)
    vals = [p.subs(y, k) for k in range(2 * n)]
    assert all(v.is_integer for v in vals)
    r, cur = [], list(vals)
    for j in range(2 * n):
        r.append(int(cur[0]))
        cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    tri.append(r)
for r in tri[:6]:
    print(r)
print("Zeilensummen:", [sum(r) for r in tri])
print("letzte Spalte N(n,2n-1):", [r[-1] for r in tri])
(HIER / "daten" / "b_baeume.json").write_text(json.dumps({"newton_dreieck": tri}, indent=1), encoding="utf-8")
