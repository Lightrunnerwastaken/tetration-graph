"""Strang P, Prüfung 8: Strang X (BRETT 23:53) — ist Cs „Halbmond-Divergenz ab N ≈ 600“ ein Rundungsartefakt der 60-Stellen-Datei?
Unabhängig: Π_N f_K(x) aus gemeinsam/lin_PiN_fK_bis_1000.json, Rauschpegel ν = Σ|b_k||x|^k·10^{−59}, und Auswertung mit
zufällig gestörter letzter Stelle. Referenz f_K(x) = exp(f_K(Log x)) (Kneser-Reihe, 4000 Koeff.)."""
import json
import random
import sys
from pathlib import Path

import mpmath as mp

from reihen import horner, kneser_koeff

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 80
random.seed(3)
GEM = Path(__file__).resolve().parent.parent / "gemeinsam"
d = json.loads((GEM / "lin_PiN_fK_bis_1000.json").read_text())
A = kneser_koeff(4000)
for x in [mp.mpf("1.8"), mp.mpc("1.4", "0.9")]:
    ref = mp.exp(horner(A, mp.log(x)))
    print(f"x = {mp.nstr(x, 4)}: f_K(x) = {mp.nstr(ref, 12)}")
    for N in ["400", "500", "600", "650", "700", "800", "1000"]:
        if N not in d:
            continue
        b = [mp.mpf(s) for s in d[N]["b"]]
        val = horner(b, x)
        nu = sum(abs(c) * abs(x) ** k for k, c in enumerate(b)) * mp.mpf(10) ** -59
        bg = [c * (1 + mp.mpf(random.uniform(-1, 1)) * mp.mpf(10) ** -59) for c in b]
        valg = horner(bg, x)
        print(f"   N = {N:>4}: |Π_N f − f_K| = {mp.nstr(abs(val - ref), 3):>9}   Rauschpegel ν = {mp.nstr(nu, 3):>9}   "
              f"|gestört − ungestört| = {mp.nstr(abs(valg - val), 3)}")
