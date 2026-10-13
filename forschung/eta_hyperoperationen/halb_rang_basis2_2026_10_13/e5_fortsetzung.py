"""Etappe 5: Rangfortsetzung der vollen Werte an komplexen Hoehen: Newton (Lagrange) in N durch
N=3..Nmax, ausgewertet bei N=3.5 und 4.5; Vergleich mit dem perturbativen Halbrang (Etappe 1)."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
import mpmath as mp
mp.mp.dps = 40
d = json.loads((HIER / 'data' / (sys.argv[1] if len(sys.argv) > 1 else 'e5_werte_N7.json')).read_text())
a = json.loads((HIER / 'data' / 'auswertung_N40J8R7.json').read_text())
def newton(seq, t):
    co = []; s = list(seq)
    while s:
        co.append(s[0]); s = [b - c for c, b in zip(s, s[1:])]
    tot = 0; b = mp.mpf(1); ps = []
    for p, c in enumerate(co):
        if p: b *= (t - p + 1) / p
        tot += b * c; ps.append(tot)
    return ps
pert = {}
for p, v in a['halbrang_jet'].items():
    z = complex(p)
    pert[(z.real, z.imag)] = (z + 1 + sum(complex(float(r), float(i)) for r, i in v['u'][:7]),
                             z + 1 + sum(complex(float(r), float(i)) for r, i in v['u32'][:7]))
out = {}
for key, row in d['werte'].items():
    x, y = key.split('|')
    vals = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in row]
    ok = [abs(v) < 1e90 for v in vals]
    n = 0
    while n < len(vals) and ok[n]: n += 1
    if n < 3: continue
    v35 = newton(vals[:n], mp.mpf('0.5'))
    v45 = newton(vals[1:n], mp.mpf('0.5')) if n >= 4 else None
    z = (float(x), float(y))
    s = 'z=%s+%si  Raenge 3..%d  H3.5: %s (letzte Korrektur %.1e)' % (x, y, 2 + n, mp.nstr(v35[-1], 9), float(abs(v35[-1] - v35[-2])))
    if z in pert:
        s += '  | perturbativ %s, Abstand %.1e' % (mp.nstr(mp.mpc(pert[z][0]), 9), float(abs(v35[-1] - pert[z][0])))
    print(s)
    out[key] = {'n_raenge': n, 'H35': [str(mp.re(v35[-1])), str(mp.im(v35[-1]))], 'korr': float(abs(v35[-1] - v35[-2])),
                'H45': [str(mp.re(v45[-1])), str(mp.im(v45[-1]))] if v45 else None}
(HIER / 'data' / 'e5_fortsetzung.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
