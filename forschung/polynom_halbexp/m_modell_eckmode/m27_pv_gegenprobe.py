"""M27 — unabhängige Gegenprobe zu Satz P-V (Strang P): zertifizierte f_K-Koeffizienten a_0…a_640 (gemeinsam/kneser_halbexp_koeff_zert_640.json).
1. Vergleich mit Strang As unabhängiger 50-Stellen-Datei (Motor/DFT + Darboux): max_n |a_n^P − a_n^A|·ρ^n.
2. Funktionalgleichung f_K(f_K(x)) = e^x mit Ps Kugeln (arb), Schwanz n > 640 über Ps Cauchy-Schranke |a_n| ≤ 2.128·1.30^{−n}:
   Residuum auf Kreisen |x| = r — muss die 0 enthalten (Einschluss) und klein sein.
3. Reellwertigkeit/Symmetrie (a_n reell), f_K(0) und f_K'(0) gegen As Werte; f_K(L) ≈ L und f_K'(L) ≈ √L über die Reihe auf |x| = ρ
   (nur mit Schwanzannahme (K2) — angezeigt, nicht als Beweis).
Aufruf: python m27_pv_gegenprobe.py
"""
import json
import os
import sys

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 500
HIER = os.path.dirname(os.path.abspath(__file__))
GEM = os.path.join(HIER, "..", "gemeinsam")
dP = json.load(open(os.path.join(GEM, "kneser_halbexp_koeff_zert_640.json"), encoding="utf-8"))
dA = json.load(open(os.path.join(GEM, "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
L = (-acb(-1).lambertw()).conjugate()
rho = abs(L)
aP = [arb(arb(m).mid(), arb(e).upper()) for m, e, _ in dP["a"]]
# Befund (11.10.): die Dezimal-Mitten haben 140 Stellen; ihre Rundung (bis 5e-141 bei a_0) steckt NICHT in Ps Radien e_n
# (e_0 = 8.8e-143). Variante "korr": Radien um die Dezimalrundung 10^(−Stellen+1)·|a_n| erweitert.
KORR = len(sys.argv) > 1 and sys.argv[1] == "korr"
if KORR:
    neu = []
    for (m, e, _), a in zip(dP["a"], aP):
        ziff = len(m.replace("-", "").replace(".", "").lstrip("0").split("e")[0].split("E")[0])
        exp10 = arb(10) ** (-(ziff - 1))
        neu.append(arb(a.mid(), (a.rad() + abs(a) * exp10 * 10).upper()))
    aP = neu
aA = [arb(x) for x in dA["a"][:641]]
K = len(aP) - 1
# 1. Vergleich
mx, nmx = arb(0), 0
for n in range(K + 1):
    v = abs(aP[n] - aA[n]) * rho ** n
    if v.upper() > mx.upper():
        mx, nmx = arb(v.upper()), n
print(f"1. max_n |a_n^P − a_n^A|·ρ^n = {float(mx):.2e} (bei n = {nmx}); Ps eigene Fehlerschranke max 2.9e-127")
# 2. Funktionalgleichung
M130 = arb("2.1284")          # ≥ 2.1283 ± 4.39e-5 (Ps Datei)
R130 = arb("1.30")


def fK(x):
    """Einschluss von f_K(x) für |x| < 1.30: Σ_{n≤640} a_n x^n (exakter Punkt) + Schwanz M(|x|/1.3)^{641}/(1 − |x|/1.3)."""
    ax = abs(x)
    assert bool(ax < R130)
    s = acb(0)
    for c in reversed(aP):
        s = s * x + c
    q = ax / R130
    sch = M130 * q ** (K + 1) / (1 - q)
    return s + acb(arb(0, sch.upper()), arb(0, sch.upper()))


for r in ("0.3", "0.5", "0.6"):
    worst = arb(0)
    maxy = arb(0)
    enth = True
    for k in range(24):
        ph = 2 * arb.pi() * k / 24
        x = acb(arb(r) * ph.cos(), arb(r) * ph.sin())
        x = acb(arb(x.real.mid()), arb(x.imag.mid()))      # exakter Punkt
        y = fK(x)
        maxy = arb(max(maxy.upper(), abs(y).upper()))
        res = fK(y) - x.exp()
        worst = arb(max(worst.upper(), abs(res).upper()))
        enth = enth and res.contains(0)
    print(f"2. |x| = {r}: max|f_K(x)| ≤ {float(maxy):.4f}, |f_K(f_K(x)) − e^x| ≤ {float(worst):.2e}, Einschluss enthält 0 überall: {enth}")
# 3. Werte
print(f"3. f_K(0) = {aP[0]}, f_K'(0) = {aP[1]}; As Werte a_0 − a_0^P = {float(abs(aA[0]-aP[0]).upper()):.1e}")
# f_K(L), f_K'(L) über die Reihe auf |x| = ρ (mit Schwanzannahme |a_n|ρ^n ≤ 100 n^{−6.4}, n > 640 — NUR Anzeige)
s = acb(0)
for c in reversed(aP):
    s = s * L + c
sp = acb(0)
for n in range(K, 0, -1):
    sp = sp * L + n * aP[n]
sch0 = 100 * arb(640) ** arb("-5.4") / arb("5.4")
sch1 = 100 * arb(640) ** arb("-4.4") / arb("4.4") / rho
print(f"   Reihe bis 640 bei L: |f_K(L) − L| ≈ {float(abs(s - L).mid()):.2e} (Schwanz angenommen ≤ {float(sch0):.1e}), "
      f"|f_K'(L) − √L| ≈ {float(abs(sp - L.sqrt()).mid()):.2e} (≤ {float(sch1):.1e})")
