"""Etappe 6: einheitlich Kneser-artige Raenge gegen die Raenge der Hauptrunde (e5_werte_N9.json):
D = H_K - H_alt; ist D eine theta-Welle?  Modell D(x) ~ H'(x) * (a sin 2pi x + b (1 - cos 2pi x))."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import numpy as np
from e6_lade import lade_raenge, H
R = lade_raenge(9)
alt = json.loads((HIER / 'data' / 'e5_werte_N9.json').read_text())['werte']
out = {}
for r in sorted(R):
    if r < 4: continue
    print('Stufe H%d' % r)
    D = {}
    for key, row in alt.items():
        x, y = key.split('|'); z = complex(float(x), float(y))
        hk = H(R, r, z)[0]; ha = complex(float(row[r - 3][0]), float(row[r - 3][1]))
        dk = (H(R, r, z + 1e-5)[0] - H(R, r, z - 1e-5)[0]) / 2e-5
        D[key] = ((hk - ha), (hk - ha) / dk)
    q = lambda k: D[k][1]
    a = (q('0.25|0') - q('0.75|0')).real / 2; b = (q('0.25|0') + q('0.75|0')).real / 2
    model = lambda z: a * np.sin(2 * np.pi * z) + b * (1 - np.cos(2 * np.pi * z))
    out[r] = {}
    for key in ['0.25|0', '0.5|0', '0.75|0', '-0.5|0', '1.25|0', '1.5|0', '1.75|0', '0.5|0.25', '0.5|0.5', '0.25|0.5', '1.5|0.5']:
        x, y = key.split('|'); z = complex(float(x), float(y))
        print('   z=%-9s D=%-26s D/H\'=%-26s Modell %s' % (key, '%.3e%+.3ej' % (D[key][0].real, D[key][0].imag),
              '%.3e%+.3ej' % (D[key][1].real, D[key][1].imag), '%.3e%+.3ej' % (model(z).real, model(z).imag)))
        out[r][key] = {'D': [D[key][0].real, D[key][0].imag], 'D_durch_Hstrich': [D[key][1].real, D[key][1].imag], 'modell': [model(z).real, model(z).imag]}
    out[r]['a_b'] = [a, b]
(HIER / 'data' / 'e6_vergleich.json').write_text(json.dumps({str(k): v for k, v in out.items()}, indent=1), encoding='utf-8')
