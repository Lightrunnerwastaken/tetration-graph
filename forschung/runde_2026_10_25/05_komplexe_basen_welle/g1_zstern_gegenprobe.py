"""Gegenprobe beim Eintragen (25.10.2026) zu K7, K8, K10: Wo liegt die Zusatzsingularitaet z* von theta_U?

(1) Bahn 0 -> 1 -> b -> b^b -> ... : anziehender Zyklus? (2) Loesungen von sexp(z) = sexp(n) nahe 0 fuer die Bahnwerte sexp(n) nahe 1
(Newton mit dem Motor, eigener Zwischenspeicher cache/ ueber kern.py). Vergleich: Prony-Ort z* aus e3_prony.py.
Aufruf: python g1_zstern_gegenprobe.py
"""
import json
from pathlib import Path

import mpmath as mp
from kern import sexp, schliessen

HIER = Path(__file__).resolve().parent
DPS = 30
mp.mp.dps = DPS + 20
FAELLE = [("exp(1.1*(sin(1)+I*cos(1)))", "welle_exp_1.1__sin_1__I_cos_1____30_0.05_512_m0_U.json"),
          ("exp(1.2*(sin(1)+I*cos(1)))", "welle_exp_1.2__sin_1__I_cos_1____30_0.05_512_m0_U.json"),
          ("exp(1.1*exp(I)*exp(-exp(I)))", "welle_exp_1.1_exp_I__exp_-exp_I____30_0.05_512_m0_U.json"),
          ("exp(1.2*exp(I)*exp(-exp(I)))", "welle_exp_1.2_exp_I__exp_-exp_I____30_0.05_512_m0_U.json")]

for basis, fn in FAELLE:
    d = json.loads((HIER / "daten" / fn).read_text(encoding="utf-8"))
    A = mp.mpmathify(d["A"].replace(" ", ""))
    b = sexp(basis, [mp.mpc(1)], DPS)[0]
    lb = mp.log(b)
    # (1) Bahn von 0 und Zyklus
    z = mp.mpc(0); bahn = []
    for _ in range(3000):
        z = mp.exp(lb * z); bahn.append(z)
        if abs(z) > 1e6:
            break
    zyk = None
    if abs(bahn[-1]) < 1e6:
        for p in range(1, 13):
            if abs(bahn[-1 - p] - bahn[-1]) < 1e-15:
                zyk = bahn[-p:]; break
    if zyk:
        mult = mp.fprod([lb * mp.exp(lb * c) for c in zyk])
        print(f"{basis}: A = {mp.nstr(A, 8)}; anziehender {len(zyk)}-Zyklus, |Multiplikator| = {mp.nstr(abs(mult), 5)}, "
              f"Punkte {[mp.nstr(c, 6) for c in zyk]}")
    else:
        print(f"{basis}: A = {mp.nstr(A, 8)}; kein anziehender Zyklus (Periode <= 12) gefunden")
    # (2) Urbilder der Bahnwerte nahe 1 in der Naehe von 0
    E = [mp.mpc(0)]
    for _ in range(14):
        E.append(mp.exp(lb * E[-1]))          # E[k] = sexp(k-1)
    for k in range(2, 15):
        if abs(E[k] - 1) < 0.12:
            z = (E[k] - 1) / A
            h = mp.mpf(10) ** -10
            for _ in range(15):
                a, c, f = sexp(basis, [z + h, z - h, z], DPS)
                st = (f - E[k]) / ((a - c) / (2 * h))
                z -= st
                if abs(st) < mp.mpf(10) ** -18:
                    break
            print(f"   sexp({k - 1}) = {mp.nstr(E[k], 8)}:  sexp(z) = sexp({k - 1}) bei z = {mp.nstr(z, 10)}"
                  f"  ({'oben' if mp.im(z) > 0 else 'unten'}; rho = e^(2 pi Im z) = {mp.nstr(mp.exp(2 * mp.pi * mp.im(z)), 6)})")
schliessen()
