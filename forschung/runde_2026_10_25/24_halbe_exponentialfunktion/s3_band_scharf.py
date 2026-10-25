"""S3: Scharfe Bandbreite (Satz K4).  D(x) = slog_lin(f(x)) - slog_lin(x) - 1/2 erfuellt D(e^x) = D(x), und
D(R) = { +-(f(u) - u - 1/2) : 0 <= u <= f(0) }.  Also sup|D| = max_{[0,f(0)]} |f(u) - u - 1/2|.
Rigoros fuer f_K mit arb: Gitter + Schranke fuer |f' - 1|.  Gegenprobe: direkte Auswertung von D an Stichproben
mit mpmath (x in (0, 1e6], ueber exp/log-Ruecklauf auf das Grundintervall).
"""
import json
import sys
from pathlib import Path
from flint import arb, ctx
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 200
QUELLE = Path(r"D:\Users\janis\Documents\Matura\tetration-graph\forschung\polynom_halbexp\gemeinsam\kneser_halbexp_koeff_zert_640.json")
d = json.loads(QUELLE.read_text(encoding="utf-8"))
K = len(d["a"]) - 1
A = [arb(m) + arb(0, float(e) * 1.0000001) for m, e, _ in d["a"]]
RHO = arb("1.37455701074")


def schwanz(j):
    return arb(461.36) / RHO ** j * arb(K) ** (j - 5) / (5 - j)


def abl(j, x):
    s = arb(0)
    for k in range(K, j - 1, -1):
        fak = 1
        for i in range(j):
            fak *= (k - i)
        s = s * x + A[k] * fak
    return s + arb(0, float(schwanz(j).upper()))


c = abl(0, arb(0))  # f(0)
# Ableitung von g(u) = u + 1/2 - f(u) ist 1 - f'(u); |1 - f'| <= 1 auf [0, 0.5] (f' in [0.87, 1.0]) -> Lipschitz 0.13
N = 2000
best, lage = None, None
for i in range(N + 1):
    u = c * i / N
    g = u + arb(1) / 2 - abl(0, u)
    if best is None or g.mid() > best.mid():
        best, lage = g, u
L1 = arb("0.13")
h = c / N
print("f(0) =", c.str(20))
print("max_u (u + 1/2 - f(u)) liegt in [", (best - L1 * h).str(8), ",", best.str(12), "] bei u ~", lage.str(8))
# Verfeinern: Nullstelle von f'(u) = 1 im Inneren mit mpmath, dann Wert rigoros
mp.mp.dps = 40
Am = [mp.mpf(m) for m, e, _ in d["a"]]
fm = lambda u: sum(Am[k] * u ** k for k in range(K + 1))
f1m = lambda u: sum(k * Am[k] * u ** (k - 1) for k in range(1, K + 1))
u0 = mp.findroot(lambda u: f1m(u) - 1, mp.mpf(float(lage.mid())))
Mw = u0 + mp.mpf(1) / 2 - fm(u0)
Mr = arb(mp.nstr(u0, 35)) + arb(1) / 2 - abl(0, arb(mp.nstr(u0, 35)))
print("Maximalstelle u0 = x* =", mp.nstr(u0, 25), "(dort f'(u0) = 1)")
print("M = u0 + 1/2 - f(u0) =", mp.nstr(Mw, 25), " rigoros:", Mr.str(15))
print("Rand: f(0) - 1/2 =", (c - arb(1) / 2).str(12))
