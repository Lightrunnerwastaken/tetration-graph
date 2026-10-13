"""Etappe 3, Alternative: Rangfortsetzung der vollen (nicht entwickelten) Werte H_N(2,x), N=3..9,
aus der Vorrunde (../half_height_values.json: Kneser fuer N=4, Fortsetzungskette fuer N>=5).
Newtonreihe im Rang ab Basis N0, Levin-u als Beschleunigung, Vergleich mit den perturbativen Werten."""
import json
from pathlib import Path
import mpmath as mp
HIER = Path(__file__).resolve().parent
mp.mp.dps = 40
d = json.loads((HIER.parent / 'half_height_values.json').read_text())
pert = {'half': {'3.5': '1.44302', '4.5': '1.46827'}, 'one_and_half': {'3.5': '2.78231', '4.5': '2.72266'}}
out = {}
def newton(seq, t):
    co = []; s = list(seq)
    while s:
        co.append(s[0]); s = [b - a for a, b in zip(s, s[1:])]
    ps = []; tot = 0; b = mp.mpf(1)
    for p, c in enumerate(co):
        if p: b *= (t - p + 1) / p
        tot += b * c; ps.append(tot)
    return ps
for key in ['half', 'one_and_half']:
    v = [mp.mpf(r[key]) for r in d]  # N=3..9
    out[key] = {}
    for ziel, N0 in [('3.5', 3), ('4.5', 4)]:
        seq = v[N0 - 3:]
        rows = []
        for M in range(4, len(seq) + 1):
            ps = newton(seq[:M], mp.mpf('0.5'))
            try:
                L = mp.levin(method='levin', variant='u'); lv, _ = L.update_psum(ps)
            except Exception:
                lv = mp.nan
            rows.append((M, ps[-1], lv))
        out[key][ziel] = [[M, str(p), str(l)] for M, p, l in rows]
        print('x=%s H_%s:' % ('0.5' if key == 'half' else '1.5', ziel),
              '; '.join('M=%d Newton %s Levin %s' % (M, mp.nstr(p, 8), mp.nstr(l, 8)) for M, p, l in rows),
              ' | perturbativ', pert[key][ziel])
(HIER / 'data' / 'e3_vollwerte.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
