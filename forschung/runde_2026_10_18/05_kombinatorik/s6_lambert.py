"""Kontrolle: Fixpunkte von e^x-1+u sind x = u-1-W_k(-e^{u-1}) (k=0,-1), gegen die Reihe xi_+(v) aus s2_generator."""
import json
from pathlib import Path
import mpmath as mp
import sympy as sp
mp.mp.dps = 40
HIER = Path(__file__).resolve().parent
d = json.load(open(HIER / "daten" / "s2_generator.json"))
for u in (mp.mpf("0.001"), mp.mpf("0.01")):
    v = mp.sqrt(u / 2)
    for k in (0, -1):
        x = u - 1 - mp.lambertw(-mp.exp(u - 1), k)
        print(f"u={u}, k={k}: x = {mp.nstr(x, 25)},  e^x-1+u-x = {mp.nstr(abs(mp.exp(x) - 1 + u - x), 3)}")
    om = 2 * mp.pi / (2 * mp.pi * 1j / x).real   # nur Kontrolle der Formel omega = |x|^2/(2|Im x|)
    print("   omega = |x|^2/(2|Im x|) =", mp.nstr(abs(x) ** 2 / (2 * abs(x.imag)), 20),
          "  Reihe v(1 - v^4/270 + v^6/42525):", mp.nstr(v * (1 - v**4 / 270 + v**6 / 42525), 20))
