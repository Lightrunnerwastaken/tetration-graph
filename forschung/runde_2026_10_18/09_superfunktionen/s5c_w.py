"""Nullstellen von w = 1 + theta' in {Im z > y0}: Windungszahl der geschlossenen Kurve x -> w(x + i y0)."""
import numpy as np, json, os
from zoo import DATEN
erg = {}
for nm in ['exp', 'q1', 'cosh', 'q4']:
    an = np.array([complex(*v) for v in json.load(open(os.path.join(DATEN, f'theta_{nm}.json')))['a']])
    n = np.arange(len(an)); r = {}
    for y0 in [0.05, 0.02, 0.01]:
        x = np.linspace(0, 1, 20001); q = np.exp(2j*np.pi*(x[:, None] + 1j*y0)*n[None, :])
        w = 1 + (q * (2j*np.pi*n*an)[None, :]).sum(1)
        ph = np.unwrap(np.angle(w)); r[y0] = [(ph[-1]-ph[0])/(2*np.pi), float(abs(w).min())]
    print(nm, r); erg[nm] = r
json.dump({k: {str(a): b for a, b in v.items()} for k, v in erg.items()}, open(os.path.join(DATEN, 'w_windung.json'), 'w'))
