"""Wendepunkte der reellen Kneser-Superfunktion: Vorzeichen von S'' auf (x_s, 1], x_s = Verzweigungspunkt (-2).
Fortpflanzung: S''(x+1) = f''(S)S'^2 + f'(S)S''; f', f'' > 0 fuer S > 0 (alle Zoo-Abb.)."""
import numpy as np, json, os
from zoo import *
from s1_loesen import ABB
def abl(s, x, r, K=32):
    r = np.broadcast_to(np.asarray(r, float), x.shape)[:, None]
    w = np.exp(2j*np.pi*np.arange(K)/K); z = x[:, None] + r*w[None, :]
    v = s.Sweit(z.ravel()).reshape(z.shape); c = np.fft.fft(v, axis=1)/K
    return c[:, 1].real/r[:, 0], 2*c[:, 2].real/r[:, 0]**2
erg = {}
for nm in ['exp', 'q1', 'cosh', 'q4']:
    a = ABB[nm][0](); r = {}
    for e in 'AB':
        s = Loesung.laden(a, f'{nm}_{e}')
        x = np.linspace(-1.99, 1.0, 2991)
        rad = np.minimum(0.2, 0.5*(x + 2))
        d1, d2 = abl(s, x, rad)
        wech = np.where(np.sign(d2[1:]) != np.sign(d2[:-1]))[0]
        xw = []
        for i in wech:   # Bisektion
            lo, hi = x[i], x[i+1]
            for _ in range(50):
                m = 0.5*(lo+hi); dm = abl(s, np.array([m]), min(0.2, 0.5*(m+2)))[1][0]
                if np.sign(dm) == np.sign(d2[i]): lo = m
                else: hi = m
            xw.append(0.5*(lo+hi))
        r[e] = dict(wendepunkte=xw, minS1=float(d1.min()), d2_rechts_min=float(d2[x > (xw[-1] if xw else 0) + 0.05].min()))
    print(nm, r['A'], ' B:', r['B']['wendepunkte'], ' |A-B| =', [abs(p-q) for p, q in zip(r['A']['wendepunkte'], r['B']['wendepunkte'])])
    erg[nm] = r
json.dump(erg, open(os.path.join(DATEN, 'wendepunkt.json'), 'w'), indent=1)
