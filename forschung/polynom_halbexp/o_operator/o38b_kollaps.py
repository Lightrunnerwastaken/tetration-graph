"""O38b — Kollaps der adjungierten Eckprofile auf der Skala N^{−1/Re κ₁}: ln|Y_1(A, ψ)| − ln|Y_1(L)| (innen) bzw. ln|T_1| − ln|T_1(L)| (aussen)
je ψ (Grad, von der äusseren Normalen) für A = 1, 2, 4 und alle N in o38_<modell>.json.  Aufruf: python o38b_kollaps.py [modell] [feld]
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1] if len(sys.argv) > 1 else "kubisch"
feld = sys.argv[2] if len(sys.argv) > 2 else "Y"
d = json.load(open(f"o38_{name}.json", encoding="utf-8"))
Ns = sorted(d, key=int)
PSIS = list(range(-180, 180, 10))
for A_ in ["1", "2", "4"]:
    print(f"A = {A_}: {feld}-Profil relativ zum Eckwert (Zeilen ψ, Spalten N = {', '.join(Ns)}); a = aussen")
    for j, ps in enumerate(PSIS):
        if ps % 20:
            continue
        row = []
        for N in Ns:
            p = d[N]["kreis"][A_][j]
            e = d[N]["ecke"]
            if feld == "Y":
                row.append(p["Y"][0] - e["Y"][0])
            elif feld == "T":
                row.append(p["T"][0] - e["T"][0])
            else:
                row.append(p[feld] - e[feld])
        print(f"  ψ={ps:+4d}{'a' if d[Ns[-1]]['kreis'][A_][j]['aussen'] else ' '} " + " ".join(f"{x:7.2f}" for x in row))
