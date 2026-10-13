"""Etappe 6: Rangfolgen H_3..H_Rmax (einheitlich) gegen die Hauptrunde: Differenzen in N, Newton bei 3.5."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import numpy as np, warnings
warnings.filterwarnings('ignore')
from e6_lade import lade_raenge, H
Rmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
R = lade_raenge(Rmax)
Rmax = max(R)
alt = json.loads((HIER / 'data' / 'e5_werte_N9.json').read_text())['werte']
def diffs(v):
    out = [np.array(v)]
    while len(out[-1]) > 1:
        out.append(np.diff(out[-1]))
    return out
def newton(v, t=0.5):
    co = [d[0] for d in diffs(v)]
    tot = 0; b = 1.0; ps = []
    for p, c in enumerate(co):
        if p: b *= (t - p + 1) / p
        tot += b * c; ps.append(tot)
    return ps
out = {}
for key in ['0.25|0', '0.5|0', '0.75|0', '1.5|0', '-0.5|0', '0.5|0.25', '0.5|0.5', '0.25|0.5', '1.5|0.5']:
    x, y = key.split('|'); z = complex(float(x), float(y))
    ein = [H(R, r, z)[0] for r in range(3, Rmax + 1)]
    al = [complex(float(a), float(b)) for a, b in alt[key][:Rmax - 2]]
    alle_alt = [complex(float(a), float(b)) for a, b in alt[key]]
    de = diffs(ein); da = diffs(al)
    ne = newton(ein); na = newton(al); na9 = newton(alle_alt)
    print('z=%-8s hoechste Differenz: einheitlich %.2e, Hauptrunde %.2e | Newton 3.5: einheitlich %s (letzte Korr %.1e), Hauptrunde gleiche Raenge %s, Hauptrunde N<=9 %s' % (
        key, abs(de[-1][0]), abs(da[-1][0]), np.round(ne[-1], 6), abs(ne[-1] - ne[-2]), np.round(na[-1], 6), np.round(na9[-1], 6)))
    out[key] = {'einheitlich': [[v.real, v.imag] for v in ein], 'newton35_einheitlich': [ne[-1].real, ne[-1].imag],
                'newton35_haupt_gleich': [na[-1].real, na[-1].imag], 'newton35_haupt_N9': [na9[-1].real, na9[-1].imag],
                'diffs_einheitlich': [abs(d[0]) for d in de], 'diffs_haupt': [abs(d[0]) for d in da]}
(HIER / 'data' / ('e6_rangfolge_R%d.json' % Rmax)).write_text(json.dumps(out, indent=1), encoding='utf-8')
