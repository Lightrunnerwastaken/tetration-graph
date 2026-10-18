"""E3: Rang-7-Phase Phi7 (hyper(7) aus dem Motor, regulär am linken Fixpunkt der Hexation).

Aufruf: python e3_rang7.py M b1 d1 ...    Ergebnis: daten/e3_<b>_<d>_M<M>.json
"""
import json, sys, time
from pathlib import Path
import mpmath as mp
import lib, lib6, lib7
mp.mp.dps = 120
D = Path(__file__).resolve().parent / "daten"
M = int(sys.argv[1]); args = sys.argv[2:]
cache = {}
def get(n):
    if n not in cache:
        T = lib.Tet(n); H = lib6.Hex(T); cache[n] = (T, H, lib7.Hyp7(H))
    return cache[n]
def fourier(v):
    M = len(v)
    return [mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(k * j) / M) for j in range(M)) / M for k in range(M // 2 + 1)]
for i in range(0, len(args), 2):
    bb, dd = args[i], args[i + 1]
    t = time.time()
    (Tb, Hb, Yb), (Td, Hd, Yd) = get(bb), get(dd)
    p = [lib7.phi7(Tb, Td, Hb, Hd, Yb, Yd, mp.mpf(j) / M) for j in range(M)]
    f = fourier(p)
    res = {"b": bb, "d": dd, "M": M, "mu7": mp.nstr(mp.re(f[0]), 40), "pp7": mp.nstr(max(p) - min(p), 20),
           "phi7": [mp.nstr(v, 40) for v in p], "f7": [[mp.nstr(mp.re(c), 25), mp.nstr(mp.im(c), 25)] for c in f]}
    (D / f"e3_{lib._tag(bb)}_{lib._tag(dd)}_M{M}.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"({bb},{dd}) mu7={mp.nstr(mp.re(f[0]), 20)} pp7={mp.nstr(max(p)-min(p), 6)} |a1|={mp.nstr(abs(f[1]), 6)} "
          f"|a2|={mp.nstr(abs(f[2]), 6)} |a3|={mp.nstr(abs(f[3]), 6)} [{time.time()-t:.0f}s]", flush=True)
