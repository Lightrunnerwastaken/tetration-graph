"""Taylor-Koeffizienten von Kneser-sexp (Basis e) an Zentren c in [-0.4,0.4] aus Cauchy-Kreisen r=0.2 (|Im| <= 0.2, sicherer Motorbereich).
Aufruf: python s1_taylor.py DPS M   -> daten/taylor_<DPS>_<M>.json
"""
import json, sys, time
import mpmath as mp
from k3 import sexp, HIER

dps = int(sys.argv[1]); M = int(sys.argv[2])
mp.mp.dps = dps + 20
r = mp.mpf("0.2")
zentren = [mp.mpf(k) / 5 for k in (-2, -1, 0, 1, 2)]
out = {"dps": dps, "M": M, "r": "0.2", "zentren": {}}
for c in zentren:
    t = time.time()
    zs = [c + r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
    v = sexp(zs, dps)
    if c == 0 or True:
        pass
    co = []
    for n in range(M // 2):
        s = mp.fsum(v[j] * mp.expjpi(-2 * mp.mpf(j * n) / M) for j in range(M)) / M / r ** n
        co.append(s)
    # reell erzwingen (Symmetrie), Imaginärteil als Fehlermass
    imfehler = max(abs(mp.im(x) * mp.mpf(0.2) ** n) for n, x in enumerate(co))
    out["zentren"][mp.nstr(c, 5)] = [mp.nstr(mp.re(x), dps + 5) for x in co]
    print(f"c={mp.nstr(c,3)}  {time.time()-t:.1f}s  max|Im c_n| r^n = {mp.nstr(imfehler,3)}  c1={mp.nstr(mp.re(co[1]),25)}", flush=True)
(HIER / "daten" / f"taylor_{dps}_{M}.json").write_text(json.dumps(out), encoding="utf-8")
