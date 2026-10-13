"""Kneser-Tetration zur Basis 2 an den Testhoehen (ueber experimente/gemeinsam.kneser, unveraendert)."""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[2] / 'experimente'))
import mpmath as mp
import gemeinsam
mp.mp.dps = 40
pts = [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), mp.mpf('1.5'), mp.mpf('1.25'), mp.mpf('1.75'),
       mp.mpc('0.5', '0.5'), mp.mpc('0.5', '1'), mp.mpc('0.25', '0.25'), mp.mpf('-0.5'), mp.mpf(2), mp.mpf(1)]
t = time.time()
w = gemeinsam.kneser('sexp', '2', pts, dps=40)
out = {'basis': '2', 'dps': 40, 'sekunden': time.time() - t,
       'werte': [{'x': [mp.nstr(mp.re(p), 20), mp.nstr(mp.im(mp.mpc(p)), 20)], 'sexp': [mp.nstr(mp.re(v), 38), mp.nstr(mp.im(v), 38)]} for p, v in zip(pts, w)]}
(HIER / 'data' / 'kneser_basis2.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
for p, v in zip(pts, w):
    print(mp.nstr(p, 6), mp.nstr(v, 30))
