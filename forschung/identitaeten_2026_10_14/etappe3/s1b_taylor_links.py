"""Taylor-Koeffizienten links: Zentren −0.6 … −1.8 (Radius 0.2, bei −1.8 Radius 0.15). Aufruf: python s1b_taylor_links.py DPS M"""
import json, sys, time
import mpmath as mp
from k3 import sexp, HIER
dps = int(sys.argv[1]); M = int(sys.argv[2])
mp.mp.dps = dps + 20
out = {"dps": dps, "M": M, "zentren": {}, "radien": {}}
for k in range(3, 10):
    c = -mp.mpf(k) / 5
    r = mp.mpf("0.15") if k == 9 else mp.mpf("0.2")
    t = time.time()
    zs = [c + r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
    v = sexp(zs, dps)
    co = [mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(j * n) / M) for j in range(M)) / M / r ** n for n in range(M // 2)]
    out["zentren"][mp.nstr(c, 5)] = [mp.nstr(mp.re(x), dps + 5) for x in co]
    out["radien"][mp.nstr(c, 5)] = mp.nstr(r, 5)
    print(f"c={mp.nstr(c,3)} r={mp.nstr(r,3)} {time.time()-t:.1f}s c1={mp.nstr(mp.re(co[1]),20)}", flush=True)
(HIER / "daten" / f"taylor_links_{dps}_{M}.json").write_text(json.dumps(out), encoding="utf-8")
