"""X19 — Lokale Eckdaten der Störungen aus x18: Fixpunkt L_P von P (P(L_P) = L_P), radiale Verschiebung δr = |L_P| − ρ,
tangentiale δt = ρ·(arg L_P − arg L), Multiplikator μ_P = P'(L_P) − μ, für P = p_N (Polynom) und P = T_N f + η mit den x18-Formen (sup R).
Frage: erklärt die radiale Verschiebung des Fixpunkts (|L_P| ≠ ρ, die Ecke rückt vom Abschneidekreis weg) die Eckempfindlichkeit?
Aufruf: python x19_fixpunkt.py N1 N2 ...
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
mp.mp.dps = 50
A_all = [mp.mpf(s) for s in json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_4000.json").read_text())["a"]]
PD = json.loads((HIER.parent / "gemeinsam" / "pN_N100_bis_1000.json").read_text())["p"]
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho = abs(L)
mu = mp.sqrt(L)


def horner(c, z):
    s = mp.mpc(0)
    ds = mp.mpc(0)
    for a in reversed(c):
        ds = ds * z + s
        s = s * z + a
    return s, ds


def daten(c):
    z = L
    for _ in range(60):
        v, dv = horner(c, z)
        z = z - (v - z) / (dv - 1)
    v, dv = horner(c, z)
    dL = z - L
    return z, abs(z) - rho, rho * (mp.arg(z) - mp.arg(L)), dv - mu


def sup_rand(c):
    zz = [rho * mp.expj(2 * mp.pi * k / 2000) for k in range(2000)]
    return max(abs(horner(c, z)[0]) for z in zz)


for N in (int(a) for a in sys.argv[1:]):
    TNf = A_all[: N + 1]
    pN = [mp.mpf(s) for s in PD[str(N)]]
    eta_p = [a - b for a, b in zip(pN, TNf)]
    R = sup_rand(eta_p)
    zf, drf, dtf, dmf = daten(TNf)
    print(f"N = {N}: R = {mp.nstr(R, 4)};  T_N f selbst: δr = {mp.nstr(drf, 4)}, δt = {mp.nstr(dtf, 4)}, δμ = {mp.nstr(dmf, 4)}  (Abschneide-Referenz)")
    q = [rho ** 2, -2 * mp.re(L), 1]
    formen = {"pN": eta_p, "e0": [mp.mpf(1)], "e1": [0, mp.mpf(1)], "eL": q, "eL2": [sum(q[i] * q[k - i] for i in range(3) if 0 <= k - i < 3) for k in range(5)],
              "e2": [0, 0, mp.mpf(1)], "e3": [0, 0, 0, mp.mpf(1)], "e4": [0, 0, 0, 0, mp.mpf(1)]}
    for name, e in formen.items():
        e = list(e) + [0] * (N + 1 - len(e))
        if name != "pN":
            s = R / sup_rand(e)
            e = [x * s for x in e]
        c = [a + b for a, b in zip(TNf, e)]
        z, dr, dt, dm = daten(c)
        print(f"   {name:4s}: η(L)/R = {mp.nstr(horner(e, L)[0] / R, 4)}, η'(L)/R = {mp.nstr(horner(e, L)[1] / R, 4)};  "
              f"(δr − δr_f)/R = {mp.nstr((dr - drf) / R, 4)}, (δt − δt_f)/R = {mp.nstr((dt - dtf) / R, 4)}, (δμ − δμ_f)/R = {mp.nstr((dm - dmf) / R, 4)}", flush=True)
