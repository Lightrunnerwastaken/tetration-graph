"""E7: Koeffizienten d_n von g(z) = sexp(z) - log(z+2) um 0 (Radius 3) per Cauchy-DFT, für den Regimewechsel bei n ~ 500.
Aufruf: python e7_gkoeff.py <r> <N> <dps>  -> daten/gkoeff_<r>_<N>_<dps>.json
"""
import sys, json, time
from pathlib import Path
import mpmath as mp
from motor import sexp
r = mp.mpf(sys.argv[1]); N = int(sys.argv[2]); dps = int(sys.argv[3])
mp.mp.dps = dps + 10
pts = [r * mp.expjpi(mp.mpf(2 * j) / N) for j in range(N // 2 + 1)]
t = time.time(); v = sexp(pts, dps); print("Motor", time.time() - t, flush=True)
v = [a - mp.log(z + 2) for a, z in zip(v, pts)]
vals = v + [mp.conj(v[N - j]) for j in range(N // 2 + 1, N)]
# DFT per Rekursion der Einheitswurzeln
w = [mp.expjpi(-mp.mpf(2 * j) / N) for j in range(N)]
co = []
for n in range(N // 2):
    s = mp.fsum(vals[j] * w[(j * n) % N] for j in range(N)) / N / r ** n
    co.append(mp.re(s))
out = Path("daten") / f"gkoeff_{sys.argv[1]}_{N}_{dps}.json"
out.write_text(json.dumps([mp.nstr(c, 30) for c in co]), encoding="utf-8")
print("fertig", time.time() - t)
