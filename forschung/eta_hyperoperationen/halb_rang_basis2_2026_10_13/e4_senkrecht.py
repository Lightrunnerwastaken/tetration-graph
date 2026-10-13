"""Etappe 4: senkrechtes Wachstum von D0 = K_eps - S_0 (Borel-Pade reell) bei eps=0.4 und 0.5.
Erwartung fuer eine theta-Mode: |D0(x+iy)| ~ e^{2 pi |y|}."""
import sys, json, math
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from e4_kneser_eps import Loeser, fixpunkt
from e4_theta import lade_reihe
mp.mp.dps = 50
J, W = lade_reihe('e4_rang1_J84.json')
PTS = [0.5, 0.5 + 0.25j, 0.5 + 0.5j, 0.5 + 0.75j, 0.5 + 1j, 0.75, 0.75 + 0.5j, 0.25, 0.25 + 0.25j]
def S0(x, eps, n):
    bh = [W[x][j] / mp.factorial(j) for j in range(J)]
    p, q = mp.pade(bh[:2 * n + 1], n, n)
    fn = lambda t: mp.exp(-t) * mp.polyval(p[::-1], eps * t) / mp.polyval(q[::-1], eps * t)
    return x + 1 + eps * mp.quad(fn, [0, 5, 10, 20, 40, 80, mp.inf])
out = {}
L = 0.8247 + 1.5674j
for eps in (0.5, 0.4):
    L = fixpunkt(eps, L)
    S = Loeser(eps, 40.0, 700, L); S.loese(iters=6000, lam=0.3, tol=1e-15); S.normieren()
    out[str(eps)] = {}
    for x in PTS:
        K = S.sexp(x)[0]
        v = [complex(K - S0(x, eps, n)) for n in (24, 27)]
        out[str(eps)][repr(x)] = [[z.real, z.imag] for z in v]
        print('eps=%.1f x=%-10s |D0| [24]=%.2e [27]=%.2e   e^{2pi y}=%.0f' % (eps, x, abs(v[0]), abs(v[1]), math.exp(2 * math.pi * complex(x).imag)), flush=True)
(HIER / 'data' / 'e4_senkrecht.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
