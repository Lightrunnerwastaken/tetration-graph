"""Loest die Kneser-Superfunktionen des Zoos in zwei Einstellungen und speichert sie."""
import numpy as np, time, sys, json
from zoo import *
ABB = {'exp': (lambda: abb_exp(1.0), 1.0),
       'q1': (lambda: abb_quadrat(1.0, 2), 2.0),       # x^2+1, S(0)=f(c)=2
       'q4': (lambda: abb_quadrat(1.0, 4), 2.0),       # x^4+1
       'cosh': (lambda: abb_cosh(0.0), float(np.cosh(1.0)))}
EINST = {'A': (40.0, 1.0, 16), 'B': (50.0, 0.75, 16)}
if __name__ == '__main__':
    namen = sys.argv[1:] or list(ABB)
    for nm in namen:
        mk, s0 = ABB[nm]; a = mk()
        print(nm, 'L =', a.L, ' lam =', a.lam, ' |lam| =', abs(a.lam), ' arg lam =', np.angle(a.lam), flush=True)
        vals = {}
        for e, (T, P, n) in EINST.items():
            t0 = time.time()
            s = Loesung(a, 0.0, s0, T=T, P=P, n=n); r = s.loese(maxit=60, info=False)
            s.speichern(f'{nm}_{e}')
            z = np.array([0.5, -0.5, 0.3 + 0.2j, 0.1j])
            vals[e] = s.S(z)
            print(f'  {e}: res {r:.1e}  S(.5)={vals[e][0].real:.16f}  S(-.5)={vals[e][1].real:.16f}  ({time.time()-t0:.0f}s)', flush=True)
        print('  max|A-B| =', np.max(abs(vals['A'] - vals['B'])), flush=True)
