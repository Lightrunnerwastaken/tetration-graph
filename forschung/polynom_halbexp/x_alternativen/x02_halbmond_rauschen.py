"""X02 — Ist die „Halbmond-Divergenz“ von Π_N f_K (C c26b, B b42, ARCHITEKTUR 6.4) echt oder ein Rundungsartefakt der 60-stelligen Datei?

Für jeden Punkt x und jedes N: Rauschpegel  ν_N(x) := Σ_k |b_k| |x|^k · 10^{−(s_k)}, s_k = Zahl der gespeicherten signifikanten Stellen
von b_k (aus dem String abgelesen). Ein ehrlich gemessener Wert ist nur dann aussagekräftig, wenn |Π_N f_K(x) − Referenz| ≫ ν_N(x).
Zusätzlich: Differenz zu einer Version, in der jedes b_k zufällig in der letzten gespeicherten Stelle gestört wird (empirisches Rauschen).
"""
import json
import os
import random
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 90
HIER = os.path.dirname(os.path.abspath(__file__))
PI = json.load(open(os.path.join(HIER, "..", "gemeinsam", "lin_PiN_fK_bis_1000.json"), encoding="utf-8"))
pkt = [mp.mpf("1.6"), mp.mpf("1.8"), mp.mpc("1.4", "0.9"), mp.mpf("1.84"), mp.mpf("1.848")]


def stellen(s):
    t = s.strip().lstrip("-").lower()
    mant = t.split("e")[0].replace(".", "").lstrip("0")
    return len(mant)


random.seed(1)
print("N | x | Π_N f_K(x) | Rauschpegel ν | empirisch |Π − Π_gestört|")
for N in (300, 500, 550, 600, 650, 700, 800, 1000):
    S = PI[str(N)]["b"]
    b = [mp.mpf(s) for s in S]
    st = [stellen(s) for s in S]
    bg = []
    for bk, sk, s in zip(b, st, S):
        e = mp.floor(mp.log10(abs(bk))) if bk != 0 else 0
        bg.append(bk + random.choice((-1, 1)) * mp.mpf(10) ** (e - sk + 1) * random.random())
    for x in pkt:
        v = mp.polyval(b[::-1], x)
        nu = sum(abs(bk) * abs(x) ** k * mp.mpf(10) ** (-sk + 1) for k, (bk, sk) in enumerate(zip(b, st)))
        vg = mp.polyval(bg[::-1], x)
        print(f"{N} | {mp.nstr(x, 4)} | {mp.nstr(v, 10)} | {mp.nstr(nu, 3)} | {mp.nstr(abs(v - vg), 3)}   (Stellen min {min(st)})")
    print()
