"""X17 — sup_{D̄}|p_N − f_K| und sup_{D̄}|p_N' − f_K'| (auf |x| = ρ, Maximumprinzip) für die echten Lösungen p_N
(gemeinsam/pN_N100_bis_1000.json) — Masstab für die Familien-Störung (R_N bzw. N·R_N) in §3.7.
Aufruf: python x17_abstand_pN.py N1 N2 ...
"""
import json
import sys
from pathlib import Path

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
mp.mp.dps = 40
A = [mp.mpf(s) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]]
PD = json.loads((HIER.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text())["p"]
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(L)


def horner(c, z):
    s = mp.mpc(0)
    ds = mp.mpc(0)
    for a in reversed(c):
        ds = ds * z + s
        s = s * z + a
    return s, ds


phis = [mp.arg(L) + mp.pi * k / 360 for k in range(-720, 720)]
fk = [horner(A, rho * mp.expj(ph)) for ph in phis]
for N in (int(a) for a in sys.argv[1:]):
    p = [mp.mpf(s) for s in PD[str(N)]]
    d0 = d1 = mp.mpf(0)
    arg0 = 0
    for ph, (f, df) in zip(phis, fk):
        v, dv = horner(p, rho * mp.expj(ph))
        if abs(v - f) > d0:
            d0, arg0 = abs(v - f), ph - mp.arg(L)
        d1 = max(d1, abs(dv - df))
    print(f"N = {N}: sup|p_N − f_K| = {mp.nstr(d0, 4)} (bei φ − arg L = {float(mp.degrees(arg0)):+.2f}°), "
          f"sup|p_N' − f_K'| = {mp.nstr(d1, 4)}, Verhältnis/N = {mp.nstr(d1 / d0 / N, 3)}", flush=True)
