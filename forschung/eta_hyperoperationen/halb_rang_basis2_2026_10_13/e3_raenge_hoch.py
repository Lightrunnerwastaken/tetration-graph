"""Etappe 3: viele eps-Ordnungen fuer die ganzzahligen Raenge n=1..Nmax (Borel-Analyse je Rang).
Speichert nach jedem Rang (Zeitbudget)."""
import sys, json, time, argparse
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import exppoly as E
from flint import acb
ap = argparse.ArgumentParser(); ap.add_argument('--J', type=int, default=36); ap.add_argument('--N', type=int, default=3)
ap.add_argument('--prec', type=int, default=2000); ap.add_argument('--budget', type=float, default=1080)
a = ap.parse_args(); E.setprec(a.prec)
P = ['0.25', '0.5', '0.75', '-0.5', '1.5', '0.5+0.5j', '0.5+1j', '3', '0.25+0.25j']
t0 = time.time(); g = [E.d_base()] + [{} for _ in range(a.J - 1)]
out = {'J': a.J, 'prec': a.prec, 'punkte': P, 'raenge': {}}
for n in range(1, a.N + 1):
    g = E.rank_step(g, a.J)
    vals = {}
    for p in P:
        z = complex(p); row = []
        for x in g:
            v = E.evaluate(x, acb(z.real, z.imag)) if x else acb(0)
            row.append([v.real.mid().str(50, radius=False), v.imag.mid().str(50, radius=False)])
        vals[p] = row
    out['raenge'][str(n)] = vals
    (HIER / 'data' / f'e3_raenge_J{a.J}.json').write_text(json.dumps(out), encoding='utf-8')
    print('Rang', n, round(time.time() - t0, 1), 's', flush=True)
    if time.time() - t0 > a.budget / 3 * 1.0 and n < a.N:
        print('Budget knapp, Abbruch'); break
