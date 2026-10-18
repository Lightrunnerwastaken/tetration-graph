"""Rest d_n = c_n - [(z-c)^n] log(z+2) an allen Zentren c (etappe3-Daten, nur lesend): Abfallrate und Quotienten."""
import json
from pathlib import Path
import mpmath as mp
mp.mp.dps = 140
HIER = Path(__file__).resolve().parent
E3 = HIER.parents[1] / "identitaeten_2026_10_14" / "etappe3" / "daten"
a = json.load(open(E3 / "taylor_130_256.json"))["zentren"]
for key in ["-0.4", "-0.2", "0.0", "0.2", "0.4"]:
    c0 = mp.mpf(key); cc = [mp.mpf(s) for s in a[key]]
    R = c0 + 2
    d = [cc[n] - (-1) ** (n + 1) / (n * R ** n) if n else 0 for n in range(len(cc))]
    rate = [abs(d[n]) ** (mp.mpf(1) / n) for n in (40, 60, 80)]
    q = [d[n + 1] / d[n] for n in (60, 70, 80)]
    print(key, "Abstand zu -3:", mp.nstr(c0 + 3, 4), " |d_n|^(1/n) n=40,60,80:", [mp.nstr(r, 6) for r in rate],
          " d_{n+1}/d_n:", [mp.nstr(x, 6) for x in q])
