"""S1: Rigoroser Nachweis (Intervallrechnung, python-flint arb) der Zutaten fuer Satz K2:
  (a) f_K'(x_lo) < 1 < f_K'(x_hi)  => es gibt x* in (x_lo, x_hi) mit f_K'(x*) = 1
  (b) f_K'' > 0 auf [x_lo, e^{x_hi}]
Grundlage: zertifizierte Taylor-Koeffizienten a_0..a_640 (Mittelpunkt +- Fehler) aus
forschung/polynom_halbexp/gemeinsam/kneser_halbexp_koeff_zert_640.json und die bewiesene Schranke
|a_k| rho^k <= 461.36 k^-6 fuer alle k (kn.satz.kneser_e_zertifiziert).
Ausgabe: Mindestwert von f'' auf dem Intervall (untere Schranke) und Lage von x*.
"""
import json
import sys
from pathlib import Path
from flint import arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 200
QUELLE = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam\kneser_halbexp_koeff_zert_640.json")
d = json.loads(QUELLE.read_text(encoding="utf-8"))
K = len(d["a"]) - 1  # 640
A = [arb(m) + arb(0, float(e) * 1.0000001) for m, e, _ in d["a"]]
RHO = arb("1.3745570107437")  # |L| (abgerundet -> Schranke rho^-k etwas groesser, sicher fuer Abschaetzung unten mit r/rho)
RHO_LO = arb("1.37455701074")


def schwanz(j, r):
    """Schranke fuer sum_{k>K} k(k-1)...(k-j+1) |a_k| |x|^{k-j} fuer |x| <= r <= rho, mit |a_k| <= 461.36 k^-6 rho^-k:
    Glied <= 461.36 k^(j-6) rho^-k |x|^(k-j) <= 461.36 k^(j-6) rho^-j; Summe k>K von k^(j-6) <= K^(j-5)/(5-j)."""
    assert j <= 4 and float(arb(r).upper()) <= float(RHO_LO.lower())
    return arb(461.36) / RHO_LO ** j * arb(K) ** (j - 5) / (5 - j)


def ableitung(j, x):
    s = arb(0)
    for k in range(K, j - 1, -1):
        fak = 1
        for i in range(j):
            fak *= (k - i)
        s = s * x + A[k] * fak
    xr = max(float(abs(x).upper()), 0.1)
    return s + arb(0, float(schwanz(j, arb(xr)).upper()))


# (a) x* einschliessen
xlo, xhi = arb("0.2410"), arb("0.2415")
f1lo, f1hi = ableitung(1, xlo), ableitung(1, xhi)
print("f'(x_lo) =", f1lo.str(12), " f'(x_hi) =", f1hi.str(12))
assert f1lo < 1 and f1hi > 1, "x* nicht eingeschlossen"

# (b) f'' > 0 auf [a, b] = [x_lo, e^{x_hi}], Unterteilung: f''(x) >= f''(c) - h * sup|f'''|
a = xlo
b = xhi.exp()
bup = arb(b.upper())
print("Intervall [", a.str(8), ",", bup.str(8), "],  rho =", RHO.str(8))
# globale Schranke fuer |f'''| auf |x| <= bup (alle Koeffizienten-Betraege)
B3 = arb(0)
for k in range(K, 2, -1):
    B3 = B3 * bup + abs(A[k]) * k * (k - 1) * (k - 2)
B3 = B3 + schwanz(3, bup)
print("sup|f'''| auf dem Intervall <=", B3.str(8))
N = 400
minwert = None
lage = None
for i in range(N):
    c = a + (bup - a) * (2 * i + 1) / (2 * N)
    h = (bup - a) / (2 * N)
    unten = ableitung(2, c) - h * B3
    if minwert is None or unten < minwert:
        minwert, lage = unten, c
    assert unten > 0, f"f'' nicht positiv nachweisbar bei {c}"
print("min untere Schranke von f'' =", minwert.str(10), "bei x ~", lage.str(6))
print("f''(0) =", ableitung(2, arb(0)).str(15), " f''(", bup.str(5), ") =", ableitung(2, bup).str(10))
print("ERGEBNIS: f_K'' > 0 auf [x_lo, e^{x_hi}] rigoros (arb, prec 200).")
