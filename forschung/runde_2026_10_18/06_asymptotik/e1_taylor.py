"""E1: Taylor-Koeffizienten von Kneser-sexp um 0 und um -1 per Cauchy-DFT aus Motorwerten.
Aufruf: python e1_taylor.py <zentrum> <r> <N> <dps>   -> schreibt daten/taylor_<zentrum>_<r>_<N>_<dps>.json
"""
import sys, json, time
from pathlib import Path
import mpmath as mp
from motor import sexp
c0 = mp.mpf(sys.argv[1]); r = mp.mpf(sys.argv[2]); N = int(sys.argv[3]); dps = int(sys.argv[4])
mp.mp.dps = dps + 10
pts = [c0 + r * mp.expjpi(mp.mpf(2 * j) / N) for j in range(N // 2 + 1)]
t = time.time()
v = sexp(pts, dps)
print("Motor", time.time() - t, flush=True)
vals = v + [mp.conj(v[N - j]) for j in range(N // 2 + 1, N)]
co = []
for n in range(N // 2):
    s = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(2 * j * n) / N) for j in range(N)) / N / r ** n
    co.append(mp.re(s))
Path("daten").mkdir(exist_ok=True)
out = Path("daten") / f"taylor_{sys.argv[1]}_{sys.argv[2]}_{N}_{dps}.json"
out.write_text(json.dumps([mp.nstr(c, dps + 5) for c in co]), encoding="utf-8")
for n in range(0, 12):
    print(n, mp.nstr(co[n], 30))
