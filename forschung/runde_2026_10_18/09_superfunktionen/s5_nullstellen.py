"""Argumentprinzip: Nullstellen von S und S' in Rechtecken [x0,x0+1] x [d, Y] der oberen Halbebene."""
import numpy as np, json, os, sys
from zoo import *
from s1_loesen import ABB


def Sp(s, z, r=0.005, K=8):
    w = np.exp(2j*np.pi*np.arange(K)/K)
    zz = z[:, None] + r*w[None, :]
    v = s.Sweit(zz.ravel()).reshape(zz.shape)
    return (v*np.conj(w)[None, :]).mean(1)/r

def rand(x0, d, Y, m):
    t = np.linspace(0, 1, m, endpoint=False)
    unten = x0 + t + 1j*d; rechts = x0+1 + 1j*(d + (Y-d)*t)
    oben = x0+1 - t + 1j*Y; links = x0 + 1j*(Y - (Y-d)*t)
    return np.concatenate([unten, rechts, oben, links, [x0+1j*d]])

def wind(v):
    ph = np.unwrap(np.angle(v)); return (ph[-1]-ph[0])/(2*np.pi)

erg = {}
for nm in (sys.argv[1:] or ['q1', 'cosh', 'q4', 'exp']):
    a = ABB[nm][0](); s = Loesung.laden(a, f'{nm}_A')
    zeilen = []
    for x0 in np.arange(-5.9, 1.2, 0.5):
        print(nm, x0, flush=True) if x0 < -5.8 else None
        z = rand(x0, 0.02, 8.0, 800)
        S = s.Sweit(z); D = Sp(s, z)
        maxsprung = np.max(abs(np.diff(np.unwrap(np.angle(D)))))
        zeilen.append([float(x0), wind(S), wind(D), float(np.min(abs(D))), float(np.min(abs(S))), float(maxsprung)])
    print(nm, ' x0, Windung S, Windung S\', min|S\'|, min|S|, max Phasenschritt S\'')
    for r in zeilen: print('   %+.1f  %+.6f  %+.6f  %.2e  %.2e  %.3f' % tuple(r))
    erg[nm] = zeilen
json.dump(erg, open(os.path.join(DATEN, 'nullstellen.json'), 'w'), indent=1)
