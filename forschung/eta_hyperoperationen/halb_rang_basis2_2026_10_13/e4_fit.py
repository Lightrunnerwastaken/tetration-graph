"""Etappe 4: Fit D0(eps,x) = C(x) * exp(-B/eps) an die stabilen Hoehen (Pade [24/24] und [27/27] einig),
Form C(x) gegen a*sin(2 pi x)+b*(1-cos 2 pi x), Hochrechnung auf eps=1."""
import json, math
from pathlib import Path
import numpy as np
HIER = Path(__file__).resolve().parent
d = json.loads((HIER / 'data' / 'e4_scan_0.6_0.5_0.45_0.4_0.36_A40_N700.json').read_text())
eps = sorted(float(e) for e in d)
out = {}
for x in ['0.125', '0.25', '0.375', '0.5', '0.625', '0.75', '0.875', '-0.25', '1.25']:
    D = []; ok = True
    for e in eps:
        r = d[str(e)][x]['D0']; a, b = r['24'], r['27']
        if a is None or b is None or abs(a - b) > 0.1 * abs(a):
            ok = False
        D.append((a + b) / 2)
    if not ok or any(v * D[0] <= 0 for v in D):
        print('x=%6s instabil oder Vorzeichenwechsel:' % x, ['%.2e' % v for v in D]); continue
    X = np.array([1 / e for e in eps]); Y = np.log(np.abs(D))
    B, lnC = np.polyfit(X, Y, 1); B = -B
    res = Y - (lnC - B * X)
    C = math.copysign(math.exp(lnC), D[0])
    out[x] = {'B': B, 'C': C, 'D0_eps1_hochgerechnet': C * math.exp(-B), 'max_fitrest_ln': float(np.abs(res).max())}
    print('x=%6s  B=%.3f  C=%+.3f  max Fitrest(ln)=%.2f  -> D0(eps=1) ~ %+.1e' % (x, B, C, np.abs(res).max(), C * math.exp(-B)))
Bm = np.mean([v['B'] for v in out.values()])
print('mittleres B = %.3f   (2*pi*ln2... = %.3f, 2pi/ln2/... )' % (Bm, 2 * math.pi * math.log(2) * 2))
# Form bei festem B
for e in eps:
    D = {x: (d[str(e)][x]['D0']['24'] + d[str(e)][x]['D0']['27']) / 2 for x in ['0.25', '0.75']}
    a = (D['0.25'] - D['0.75']) / 2; b = (D['0.25'] + D['0.75']) / 2
    print('eps=%.2f  a*e^{B/eps}=%+.3f  b*e^{B/eps}=%+.3f' % (e, a * math.exp(Bm / e), b * math.exp(Bm / e)))
out['B_mittel'] = Bm
(HIER / 'data' / 'e4_fit.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
