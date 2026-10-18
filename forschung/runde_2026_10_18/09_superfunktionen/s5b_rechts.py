"""Nullstellen von S in H rechts (Re z in [0,3]): Windung mit feiner Abtastung, dann Newton."""
import numpy as np, json, os, sys
from zoo import *
from s1_loesen import ABB
exec(open('s5_nullstellen.py', encoding='utf-8').read().split('erg = {}')[0].split('from s1_loesen import ABB')[1])
erg = {}
for nm in ['q4', 'q1', 'cosh', 'exp']:
    a = ABB[nm][0](); s = Loesung.laden(a, f'{nm}_A'); zeilen = []
    for x0 in [0.0, 0.5, 1.0, 1.5, 2.0]:
        z = rand(x0, 0.02, 6.0, 4000); S = s.Sweit(z)
        st = np.max(abs(np.diff(np.unwrap(np.angle(S)))))
        zeilen.append([x0, wind(S), st])
    print(nm, [('%.1f: %+.3f (Schritt %.2f)' % tuple(r)) for r in zeilen], flush=True)
    # Nullstellen suchen: Gitter
    X, Y = np.meshgrid(np.linspace(0, 2.5, 251), np.linspace(0.02, 6, 300))
    V = s.Sweit((X + 1j*Y).ravel()).reshape(X.shape)
    kand = []
    idx = np.argsort(abs(V).ravel())[:40]
    for i in idx:
        z = complex(X.ravel()[i] + 1j*Y.ravel()[i])
        for _ in range(30):
            v = s.Sweit(np.array([z]))[0]; d = (s.Sweit(np.array([z+1e-6]))[0]-s.Sweit(np.array([z-1e-6]))[0])/2e-6
            z = z - v/d
        v = abs(s.Sweit(np.array([z]))[0])
        if v < 1e-10 and z.imag > 0 and all(abs(z-k) > 1e-6 for k in kand): kand.append(z)
    print('  Nullstellen von S in H:', [f'{k:.12f}' for k in kand])
    for k in kand:
        print('   S(z+1) =', s.Sweit(np.array([k+1]))[0], " S'(z+1) ~", abs(s.Sweit(np.array([k+1+1e-5]))[0]-s.Sweit(np.array([k+1-1e-5]))[0])/2e-5)
    erg[nm] = dict(wind=zeilen, nullst=[[k.real, k.imag] for k in kand])
json.dump(erg, open(os.path.join(DATEN, 'nullstellen_rechts.json'), 'w'), indent=1)
