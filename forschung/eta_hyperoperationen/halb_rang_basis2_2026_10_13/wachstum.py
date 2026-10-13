"""Wachstum der ganzzahligen Rangjets c_n^{(j)}(x) in n: log10|c| und Quotienten."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
from halbschritt import lade, folge
import mpmath as mp
mp.mp.dps = 60
d = lade(sys.argv[1] if len(sys.argv) > 1 else 'N40J12')
out = {}
for p in ['0.5', '1.5', '0.25']:
    print('x =', p)
    out[p] = {}
    for j in range(d['J']):
        f = folge(d, p, j)
        N = len(f) - 1
        lg = [float(mp.log10(abs(v))) if v != 0 else -99 for v in f]
        q = [f[n] / f[n - 1] for n in range(N - 3, N + 1)]
        q2 = mp.sqrt(f[N] / f[N - 2]) if f[N - 2] != 0 else 0
        out[p][j + 1] = {'log10': lg, 'quot_letzte': [mp.nstr(mp.re(z), 6) for z in q]}
        print('  j=%2d' % (j + 1), ' '.join('%5.1f' % lg[n] for n in range(0, N + 1, 5)),
              '  q:', ' '.join(mp.nstr(mp.re(z), 4) for z in q), ' sqrt(c_N/c_N-2):', mp.nstr(mp.re(q2), 4))
(HIER / 'data' / 'wachstum.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
