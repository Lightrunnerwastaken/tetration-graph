"""Rangleiter bei eps=1 (x fest): H_3, H_3.5, H_4, H_4.5, H_5 aus den Jets (Ordnungen 1..7, optimaler Abbruch),
und Residuum der Rangrekursion H_4.5(x+1) - H_3.5(H_4.5(x)) bei eps=1 mit den abgeschnittenen Jets."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from halbschritt import lade
from summen import opt_abbruch
mp.mp.dps = 60
d = lade('N40J8R7'); a = json.loads((HIER / 'data' / 'auswertung_N40J8R7.json').read_text())
J0 = 7
def ganzrang(p, n):
    return [mp.mpc(mp.mpf(x), mp.mpf(y)) for x, y in d['werte'][p][n][0]][:J0]
out = {}
print('x      H_3=2^x     H_3.5        H_4          H_4.5        H_5   (je opt. Abbruch, Fehlerskala kleinstes Glied)')
for p in ['0.25', '0.5', '0.75', '1.25', '1.5', '1.75', '-0.5']:
    x = mp.mpc(complex(p)); row = {}
    for name, terme in [('H3', ganzrang(p, 0)), ('H3.5', [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in a['halbrang_jet'][p]['u']][:J0]),
                        ('H4', ganzrang(p, 1)), ('H4.5', [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in a['halbrang_jet'][p]['u32']][:J0]),
                        ('H5', ganzrang(p, 2))]:
        if name == 'H3':
            s, kl = 2 ** x, 0
        else:
            s, kl, _ = opt_abbruch(x, terme)
        row[name] = [float(mp.re(s)), float(kl)]
    out[p] = row
    print(p.ljust(6), '  '.join('%.6f(%.0e)' % tuple(v) for v in row.values()))
# Rangrekursion bei eps=1
from auswertung import wert
print('Residuum H_4.5(x+1) - H_3.5(H_4.5(x)) bei eps=1, Jets bis Ordnung 7:')
res = {}
for x, x1 in [('0.25', '1.25'), ('0.5', '1.5'), ('0.75', '1.75'), ('-0.5', '0.5'), ('0.5+0.5j', '1.5+0.5j')]:
    w = sum(mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in a['halbrang_jet'][x]['u32'][:J0])
    w1 = sum(mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in a['halbrang_jet'][x1]['u32'][:J0])
    # Ableitungen von u (Summe der Ordnungen 1..7) bei x+1
    ud = [sum(wert(d, x1, j, r, mp.mpf('0.5')) for j in range(J0)) for r in range(J0)]
    U = w + sum(ud[r] * w ** r / mp.factorial(r) for r in range(J0))  # U(x+1+w) - (x+2)
    r_ = w1 - U
    res[x] = float(abs(r_))
    print('  x=%s: %.2e' % (x, abs(r_)))
out['rekursion_eps1'] = res
(HIER / 'data' / 'leiter.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
