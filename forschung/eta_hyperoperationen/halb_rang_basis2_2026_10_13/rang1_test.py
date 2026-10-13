"""Rang 1: Jet von T(f_eps), f_eps=x+1+eps(2^x-x-1), Ordnungen j=1..J, Summe bei eps=1 gegen Kneser."""
import sys, json, time, argparse
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import exppoly as E
from flint import acb
ap = argparse.ArgumentParser(); ap.add_argument('--J', type=int, default=20); ap.add_argument('--prec', type=int, default=500)
a = ap.parse_args(); E.setprec(a.prec)
kn = json.loads((HIER / 'data' / 'kneser_basis2.json').read_text())
pts = [(complex(float(w['x'][0]), float(w['x'][1])), complex(float(w['sexp'][0]), float(w['sexp'][1]))) for w in kn['werte']]
t = time.time()
g = [E.d_base()] + [{} for _ in range(a.J - 1)]
e = E.rank_step(g, a.J)
print('Zeit', round(time.time() - t, 1), 's; Grade', [E.degree(x) for x in e])
rows = []
for x, ref in pts:
    s = x + 1; terms = []; hp = []
    for j in range(a.J):
        v = E.evaluate(e[j], acb(x.real, x.imag)) if e[j] else acb(0)
        vc = complex(float(v.real.mid()), float(v.imag.mid()))
        terms.append(vc); s += vc
        hp.append([v.real.mid().str(50, radius=False), v.imag.mid().str(50, radius=False), float(v.rad())])
    rows.append({'x': [x.real, x.imag], 'kneser': [ref.real, ref.imag], 'terme': [[t_.real, t_.imag] for t_ in terms], 'summe': [s.real, s.imag], 'terme_hp': hp})
    print(x, "Summe-Kneser", abs(s - ref))
(HIER / 'data' / f'rang1_J{a.J}.json').write_text(json.dumps({'J': a.J, 'prec': a.prec, 'zeilen': rows}, indent=1), encoding='utf-8')
