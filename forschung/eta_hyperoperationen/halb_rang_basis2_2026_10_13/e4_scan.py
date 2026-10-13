"""Etappe 4: eps-Abtastung von D0(eps,x) = K_eps(x) - S_0(eps,x), S_0 = Borel-Pade-Laplace laengs der
reellen Achse; dazu kleinstes Glied tau.  Aufruf: python e4_scan.py 0.6,0.5 [A N]"""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from e4_kneser_eps import Loeser, fixpunkt
from e4_theta import lade_reihe
mp.mp.dps = 50
J, W = lade_reihe('e4_rang1_J84.json')
EPS = [float(e) for e in sys.argv[1].split(',')]
A = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
N = int(sys.argv[3]) if len(sys.argv) > 3 else 700
XS = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, -0.25, 1.25]
PADE = {}
def S0(x, eps, n):
    key = (x, n)
    if key not in PADE:
        bh = [W[x][j] / mp.factorial(j) for j in range(J)]
        PADE[key] = mp.pade(bh[:2 * n + 1], n, n)
    p, q = PADE[key]
    fn = lambda t: mp.exp(-t) * mp.polyval(p[::-1], eps * t) / mp.polyval(q[::-1], eps * t)
    return x + 1 + eps * mp.quad(fn, [0, 5, 10, 20, 40, 80, mp.inf])
out = {}
L = 0.8247 + 1.5674j
for eps in sorted(EPS, reverse=True):
    t0 = time.time()
    L = fixpunkt(eps, L)
    S = Loeser(eps, A, N, L); res = S.loese(iters=6000, lam=0.3, tol=1e-15); S.normieren()
    out[str(eps)] = {}
    for x in XS:
        K = S.sexp(x)[0]
        terms = [abs(complex(W[x][j])) * eps ** (j + 1) for j in range(J)]
        tau = min(terms[2:])
        d = {}
        for n in (24, 27):
            try:
                d[n] = float(mp.re(K - S0(x, eps, n)))
            except ZeroDivisionError:
                d[n] = None
        out[str(eps)][str(x)] = {'K': [K.real, K.imag], 'D0': d, 'tau': tau}
        print('eps=%.3f x=%6.3f  D0[24]=%+.4e  D0[27]=%+.4e  tau=%.2e' % (eps, x, d[24] or 0, d[27] or 0, tau), flush=True)
    print('   Loeser-Residuum %.1e, %.0f s' % (res, time.time() - t0), flush=True)
(HIER / 'data' / ('e4_scan_%s_A%g_N%d.json' % ('_'.join(sys.argv[1].split(',')), A, N))).write_text(json.dumps(out, indent=1), encoding='utf-8')
