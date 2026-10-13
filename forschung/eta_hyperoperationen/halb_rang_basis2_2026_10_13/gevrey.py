"""Divergenzanalyse der eps-Reihe von T(f_eps): Gevrey-Fit und Borel-Pade (auch seitlich)."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp, numpy as np
from borel import borel_pade
mp.mp.dps = 80
d = json.loads((HIER / 'data' / (sys.argv[1] if len(sys.argv) > 1 else 'rang1_J52.json')).read_text())
out = {}
for r in d['zeilen'][:10]:
    x = complex(*r['x']); ref = mp.mpc(*r['kneser'])
    c = [mp.mpc(x + 1)] + [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b, _ in r['terme_hp']]
    J = len(c) - 1
    lg = np.array([float(mp.log(abs(c[j]))) for j in range(1, J + 1)])
    j = np.arange(1, J + 1)
    # Einhuellende: lokale Maxima ab j=15
    env = [i for i in range(15, J - 1) if lg[i] >= lg[i - 1] and lg[i] >= lg[i + 1]] + [J - 1]
    # Fit log|c_j| - lgamma(j) = a - j log A  (Gevrey 1, beta=0)
    from math import lgamma
    y = np.array([lg[i] - lgamma(j[i]) for i in env])
    p = np.polyfit(j[env], y, 1)
    A = float(np.exp(-p[0]))
    res = {'x': [x.real, x.imag], 'gevrey1_A': A, 'env_j': [int(j[i]) for i in env]}
    sums = {}
    for (L, M) in [(24, 24), (25, 26), (26, 26)]:
        for ang in [0.0, 0.15, -0.15]:
            try:
                v = borel_pade(c, L, M, richtung=ang)
                sums['%d/%d@%+.2f' % (L, M, ang)] = [float(mp.re(v - ref)), float(mp.im(v - ref))]
            except Exception as ex:
                sums['%d/%d@%+.2f' % (L, M, ang)] = str(ex)
    res['borel_minus_kneser'] = sums
    out[str(x)] = res
    print(x, 'A(Gevrey-1-Fit) = %.3f' % A, ' 2pi/ln2 = %.3f' % (2 * np.pi / np.log(2)))
    for k, v in sums.items():
        print('    ', k, v if isinstance(v, str) else '%+.2e %+.2ei' % tuple(v))
(HIER / 'data' / 'gevrey_rang1.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
