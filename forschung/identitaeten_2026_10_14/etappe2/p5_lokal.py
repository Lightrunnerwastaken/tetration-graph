"""(1) Lokale Form direkt prüfen: θ(iy) gegen 2 + h(log(iAy)) und 2 + h(log(iAy) + κ·iy) − iy (Basis e).
(2) Führender Term T1(n) = −h'(c − ln n)/n für riesige n: Betrag·n·ln n und Argument (gibt es einen Grenzwert?).
Aufruf: python p5_lokal.py
"""
import json
from pathlib import Path

import mpmath as mp
from kern import Fix, sexp

mp.mp.dps = 70
fx = Fix(mp.e, K=40)
per = 2j * mp.pi / mp.log(fx.lam)
kor = json.loads(Path("daten/korrektur_e.json").read_text())
A, Bc = mp.mpf(kor["A"]), mp.mpf(kor["B"])
kappa = Bc / A - A / 2
out = {"lokal": [], "gross_n": []}


def red(x, ref):
    m = min(range(-6, 7), key=lambda m: abs(x + m * per - ref))
    return x + m * per


print("(1)  y, |θ − M1|, |θ − M2|   (M1 = 2+h(log(iAy)), M2 = 2 − iy + h(log(iAy) + κ iy))")
ys = [mp.mpf(10) ** (-k) for k in (1, 2, 3, 5, 8, 12, 20, 30, 40)]
ws = sexp("exp(1)", [1j * y for y in ys], 50)
for y, w in zip(ys, ws):
    th = fx.h(w) - 1j * y
    m1 = 2 + fx.h(mp.log(1j * A * y))
    m2 = 2 - 1j * y + fx.h(mp.log(1j * A * y) + kappa * 1j * y)
    d1, d2 = abs(red(th, m1) - m1), abs(red(th, m2) - m2)
    out["lokal"].append([mp.nstr(y, 3), mp.nstr(th, 15), mp.nstr(d1, 4), mp.nstr(d2, 4)])
    print(f"  y=1e{int(mp.log10(y))}: θ={mp.nstr(th,14)}  |θ−M1|={mp.nstr(d1,4)}  |θ−M2|={mp.nstr(d2,4)}")

c = mp.log(1j * A / (2 * mp.pi))
print("(2)  ln n, n·ln n·|T1|, arg T1")
for t in ["8.3", "20", "100", "1e3", "1e5", "1e10", "1e30", "1e100", "1e1000"]:
    t = mp.mpf(t)
    v = c - t
    hp = mp.diff(lambda s: fx.h(s), v)
    T = -hp            # n·T1
    out["gross_n"].append([mp.nstr(t, 5), mp.nstr(abs(T) * t, 10), mp.nstr(mp.arg(T), 10)])
    print(f"  ln n={mp.nstr(t,5):>8}  n ln n |T1| = {mp.nstr(abs(T)*t,10):>14}   arg = {mp.nstr(mp.arg(T),10)}")
Path("daten/lokal_e.json").write_text(json.dumps(out, indent=1))
