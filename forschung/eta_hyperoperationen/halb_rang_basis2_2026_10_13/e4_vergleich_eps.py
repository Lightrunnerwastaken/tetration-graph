"""Etappe 4: Kneser-artige Superfunktion von f_eps gegen die perturbative eps-Reihe (84 Ordnungen).

D(x) = F_eps(x) - S_opt(x;eps), S_opt = Partialsumme bis vor das kleinste Glied (ab j=3).
tau = kleinstes Glied.  Ist |D| << tau ueber alle eps, ist Kneser asymptotisch gleich der Reihe
(theta jenseits aller Ordnungen); ist |D|/tau gross, traegt Kneser ein nichtperturbatives theta."""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import numpy as np
import mpmath as mp
from e4_kneser_eps import Loeser, fixpunkt
from e4_theta import lade_reihe
J, W = lade_reihe('e4_rang1_J84.json')
xs = [k / 16 for k in range(-8, 25)]
EPS = [float(e) for e in (sys.argv[1].split(',') if len(sys.argv) > 1 else '1,0.8,0.6,0.5,0.4,0.35,0.3'.split(','))]
A = float(sys.argv[2]) if len(sys.argv) > 2 else 40.0
N = int(sys.argv[3]) if len(sys.argv) > 3 else 700
out = {}
L = 0.8247 + 1.5674j
for eps in sorted(EPS, reverse=True):
    t0 = time.time()
    L = fixpunkt(eps, L)
    S = Loeser(eps, A, N, L)
    res = S.loese(iters=6000, lam=0.3, tol=1e-15)
    c = S.normieren()
    rows = []
    for x in xs:
        terms = [complex(W[x][j]) * eps ** (j + 1) for j in range(J)]
        ab = [abs(t) for t in terms]
        js = min(range(2, J), key=lambda j: ab[j])
        Sopt = x + 1 + sum(terms[:js])
        F = S.sexp(x)[0]
        rows.append([x, (F - Sopt).real, ab[js], js + 1, F.real])
    Dmax = max(abs(r[1]) for r in rows); tmax = max(r[2] for r in rows)
    ratio = max(abs(r[1]) / r[2] for r in rows)
    out[str(eps)] = {'L': [L.real, L.imag], 'residuum': float(res), 'c': c, 'zeilen': rows}
    print('eps=%.3g  |D|max=%.2e  tau_max=%.2e  max|D|/tau=%.2f  j*~%d  Loeser-Res %.1e  %.0fs' % (
        eps, Dmax, tmax, ratio, rows[len(rows) // 2][3], res, time.time() - t0), flush=True)
    for r in rows[::4]:
        print('     x=%7.4f  D=%+.3e  tau=%.2e  j*=%d' % (r[0], r[1], r[2], r[3]))
(HIER / 'data' / ('e4_vergleich_eps_%s.json' % '_'.join(str(e) for e in EPS))).write_text(json.dumps(out, indent=1), encoding='utf-8')
