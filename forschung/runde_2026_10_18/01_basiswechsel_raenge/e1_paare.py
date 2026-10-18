"""E1: Rang-4- und Rang-5-Phase für Basispaare: Mittel mu, Welle (Spannweite, Fourier-Moden).

Aufruf: python e1_paare.py M b1 d1 b2 d2 ...   (Basen als Motor-Ausdrücke, z. B. exp(1) 2)
Ergebnis: daten/e1_<b>_<d>_M<M>.json mit Phi4, Phi5 auf dem Gitter th_j = j/M und Fourier-Koeffizienten.
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import lib

mp.mp.dps = 120
D = Path(__file__).resolve().parent / "daten"
D.mkdir(exist_ok=True)
M = int(sys.argv[1])
args = sys.argv[2:]
tets = {}


def tet(n):
    if n not in tets:
        tets[n] = lib.Tet(n)
    return tets[n]


def fourier(v):
    M = len(v)
    out = []
    for k in range(M // 2 + 1):
        s = mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(k * j) / M) for j in range(M)) / M
        out.append(s)
    return out


for i in range(0, len(args), 2):
    bb, dd = args[i], args[i + 1]
    t = time.time()
    Tb, Td = tet(bb), tet(dd)
    ths = [mp.mpf(j) / M for j in range(M)]
    p4 = [lib.phi4(Tb, Td, th) for th in ths]
    p5 = [lib.phi5(Tb, Td, th) for th in ths]
    f4, f5 = fourier(p4), fourier(p5)
    res = {"b": bb, "d": dd, "M": M, "L_b": mp.nstr(Tb.L, 40), "L_d": mp.nstr(Td.L, 40),
           "lam_b": mp.nstr(Tb.lam, 40), "lam_d": mp.nstr(Td.lam, 40),
           "mu4": mp.nstr(mp.re(f4[0]), 50), "mu5": mp.nstr(mp.re(f5[0]), 50),
           "pp4": mp.nstr(max(p4) - min(p4), 20), "pp5": mp.nstr(max(p5) - min(p5), 20),
           "phi4": [mp.nstr(v, 50) for v in p4], "phi5": [mp.nstr(v, 50) for v in p5],
           "f4": [[mp.nstr(mp.re(c), 30), mp.nstr(mp.im(c), 30)] for c in f4],
           "f5": [[mp.nstr(mp.re(c), 30), mp.nstr(mp.im(c), 30)] for c in f5]}
    (D / f"e1_{lib._tag(bb)}_{lib._tag(dd)}_M{M}.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"({bb},{dd}) mu4={mp.nstr(mp.re(f4[0]), 20)} mu5={mp.nstr(mp.re(f5[0]), 20)} pp4={mp.nstr(max(p4)-min(p4), 6)} "
          f"pp5={mp.nstr(max(p5)-min(p5), 6)} |a1|: {mp.nstr(abs(f4[1]), 6)} / {mp.nstr(abs(f5[1]), 6)}  [{time.time()-t:.0f}s]", flush=True)
