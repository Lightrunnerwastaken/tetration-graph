"""O40 — Tabellen für (KB) aus o38_<modell>.json: alle Grössen minus ½·ln N (Ziel κ_E = 1/2).
Spalten: Ecke (Y_1, T_1, G, G^T), Kreise A = 0.5/1/2/4 in Ws Konvention (innen Y, aussen T): Einzelmode KB_1, Summe über r Moden KB_r,
SVD-freie Majorante max(G innen, G^T aussen) =: KB_G (gilt für JEDEN Einheitsvektor u, also alle Eckmoden zugleich).
Aufruf: python o40_kb_tabelle.py [modell ...]
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
for name in (sys.argv[1:] or ["kubisch"]):
    d = json.load(open(f"o38_{name}.json", encoding="utf-8"))
    print(f"== {name}: Werte − ½ln N (κ_E = 1/2 ⇔ Spalte beschränkt)")
    print("   N  | lnσ₁ | Ecke Y₁   T₁    G    G^T | " + " | ".join(f"A={a}: KB₁  KB_r  KB_G" for a in ["0.5", "1", "2", "4"]))
    for N in sorted(d, key=int):
        o = d[N]
        h = 0.5 * np.log(int(N))
        e = o["ecke"]
        row = f"{int(N):5d} | {np.log(o['s'][0]):4.2f} | {e['Y'][0] - h:5.2f} {e['T'][0] - h:5.2f} {e['G'] - h:5.2f} {e['GT'] - h:5.2f} |"
        for a in ["0.5", "1", "2", "4"]:
            P = o["kreis"][a]
            k1 = max((p["T"][0] if p["aussen"] else p["Y"][0]) for p in P) - h
            kr = max((p["Tr"] if p["aussen"] else p["Yr"]) for p in P) - h
            kg = max((p["GT"] if p["aussen"] else p["G"]) for p in P) - h
            row += f" {k1:5.2f} {kr:5.2f} {kg:5.2f} |"
        print(row)
