"""Etappe 4: Kneser-Tetration zur Basis 2 auf dem Raster x=k/16 (k=-8..32) und an komplexen Hoehen.
Ein einzelner Aufruf von gemeinsam.kneser (40 Stellen)."""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[2] / 'experimente'))
import mpmath as mp
import gemeinsam
mp.mp.dps = 40
pts = [mp.mpf(k) / 16 for k in range(-8, 33)] + [mp.mpc(0.5, y) for y in ('0.25', '0.5', '0.75', '1')] + [mp.mpc('0.25', '0.25'), mp.mpc('0.75', '0.5')]
t = time.time()
w = gemeinsam.kneser('sexp', '2', pts, dps=40)
out = {'dps': 40, 'sekunden': time.time() - t,
       'werte': [{'x': [mp.nstr(mp.re(p), 20), mp.nstr(mp.im(mp.mpc(p)), 20)], 'sexp': [mp.nstr(mp.re(v), 38), mp.nstr(mp.im(v), 38)]} for p, v in zip(pts, w)]}
(HIER / 'data' / 'e4_kneser_raster.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
print('fertig', round(time.time() - t, 1), 's', len(pts), 'Punkte')
