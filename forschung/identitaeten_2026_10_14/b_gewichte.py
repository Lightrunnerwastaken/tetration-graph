"""Beweis-Kontrolle: die zwei obersten Gewichtsschichten der Störungsreihe G(y) = e^{G(y−1)} − 1 + u.

H(y) = G(y−1) = Σ g_n(y−1) uⁿ. Gewicht: y → 1, u → −2. Oberste Schicht (u^n y^{2n−1}) und zweite (u^n y^{2n−2}) erfüllen
  H0' = H0²/2 + u,          H0(0) = 0          ⇒ H0 = √(2u)·tan θ,  θ = y√(u/2)
  H1' = H0·H1 − H0³/12 − u·H0/2,  H1(0) = −u   ⇒ H1 = (u/3)(sec²θ·ln cos θ + 1 − 4 sec²θ)
Geprüft gegen die exakten Polynome g_n (b_kombinatorik.py) für n ≤ 12.
"""
import json
from pathlib import Path
import sympy as sp

HIER = Path(__file__).resolve().parent
d = json.loads((HIER / "daten" / "b_kombinatorik.json").read_text(encoding="utf-8"))["daten"]["stoerungsreihe"]["dreieck"]
y, u, s = sp.symbols("y u s", positive=True)
NM = len(d)
g = [None] + [sum(sp.Rational(c) * y ** k for k, c in enumerate(row)) for row in d]
# H(y) = G(y−1): Koeffizienten von y^{2n−1}, y^{2n−2} in g_n(y−1)
top, sec = [], []
for n in range(1, NM + 1):
    P = sp.Poly(sp.expand(g[n].subs(y, y - 1)), y)
    top.append(P.coeff_monomial(y ** (2 * n - 1)))
    sec.append(P.coeff_monomial(y ** (2 * n - 2)))
# u = s², θ = y s/√2 ; H0 = √2 s tan θ ; H1 = (s²/3)(sec² ln cos + 1 − 4 sec²)
th = y * s / sp.sqrt(2)
H0 = sp.sqrt(2) * s * sp.tan(th)
H1 = s ** 2 / 3 * (sp.log(sp.cos(th)) / sp.cos(th) ** 2 + 1 - 4 / sp.cos(th) ** 2)
ord_ = 2 * NM + 2
S0 = sp.expand(sp.series(H0, s, 0, ord_).removeO())
S1 = sp.expand(sp.series(H1, s, 0, ord_ + 1).removeO())
ok0 = all(sp.simplify(S0.coeff(s, 2 * n).coeff(y, 2 * n - 1) - top[n - 1]) == 0 for n in range(1, NM + 1))
ok1 = all(sp.simplify(S1.coeff(s, 2 * n).coeff(y, 2 * n - 2) - sec[n - 1]) == 0 for n in range(1, NM + 1))
print("oberste Schicht = sqrt(2u) tan(y sqrt(u/2)):", ok0)
print("zweite Schicht  = (u/3)(sec^2 ln cos + 1 - 4 sec^2):", ok1, [str(v) for v in sec[:6]])
(HIER / "daten" / "b_gewichte.json").write_text(json.dumps({"oberste": ok0, "zweite": ok1, "n_max": NM,
    "zweitkoef_g_n(y-1)": [str(v) for v in sec]}, indent=1), encoding="utf-8")
