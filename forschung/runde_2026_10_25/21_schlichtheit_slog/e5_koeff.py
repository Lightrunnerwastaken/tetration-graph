"""E5: Koeffizienten c_n von sexp um 0 gegen log(z+2): d_n = c_n - (-1)^(n+1)/(n 2^n); Wurzeltest fuer den Restradius."""
import json, mpmath as mp
mp.mp.dps = 130
from pathlib import Path
d = json.loads((Path(__file__).resolve().parent.parent.parent / "identitaeten_2026_10_14/etappe3/daten/taylor_130_256.json").read_text())
c = [mp.mpf(x) for x in d["zentren"]["0.0"]]
for n in [10, 20, 40, 60, 80, 100, 120, 127]:
    dn = c[n] - (-1) ** (n + 1) / (n * mp.mpf(2) ** n)
    print(n, mp.nstr(dn, 8), " |d_n|^(-1/n) =", mp.nstr(abs(dn) ** (-mp.mpf(1) / n), 8), " d_n/d_(n-1) =", mp.nstr(dn / (c[n-1] - (-1) ** n / ((n-1) * mp.mpf(2) ** (n-1))), 8))
