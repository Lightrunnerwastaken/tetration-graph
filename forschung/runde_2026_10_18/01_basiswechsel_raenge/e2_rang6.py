"""E2: Rang-6-Phase Phi6 für Basispaare (Kneser-Hexation aus dem Motor), Vergleich mit Rang 4 und 5 (aus E1).

Aufruf: python e2_rang6.py M b1 d1 b2 d2 ...    Ergebnis: daten/e2_<b>_<d>_M<M>.json
"""
import json, sys, time
from pathlib import Path
import mpmath as mp
import lib, lib6
mp.mp.dps = 120
D = Path(__file__).resolve().parent / "daten"
M = int(sys.argv[1]); args = sys.argv[2:]
cache = {}
def get(n):
    if n not in cache:
        T = lib.Tet(n); cache[n] = (T, lib6.Hex(T))
    return cache[n]
def fourier(v):
    M = len(v)
    return [mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(k * j) / M) for j in range(M)) / M for k in range(M // 2 + 1)]
for i in range(0, len(args), 2):
    bb, dd = args[i], args[i + 1]
    t = time.time()
    (Tb, Hb), (Td, Hd) = get(bb), get(dd)
    p6 = [lib6.phi6(Tb, Td, Hb, Hd, mp.mpf(j) / M) for j in range(M)]
    f6 = fourier(p6)
    res = {"b": bb, "d": dd, "M": M, "mu6": mp.nstr(mp.re(f6[0]), 40), "pp6": mp.nstr(max(p6) - min(p6), 20),
           "phi6": [mp.nstr(v, 40) for v in p6], "f6": [[mp.nstr(mp.re(c), 25), mp.nstr(mp.im(c), 25)] for c in f6]}
    (D / f"e2_{lib._tag(bb)}_{lib._tag(dd)}_M{M}.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"({bb},{dd}) mu6={mp.nstr(mp.re(f6[0]), 20)} pp6={mp.nstr(max(p6)-min(p6), 6)} |a1|={mp.nstr(abs(f6[1]), 6)} "
          f"|a2|={mp.nstr(abs(f6[2]), 6)} [{time.time()-t:.0f}s]", flush=True)
