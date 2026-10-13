"""Etappe 4: Laplace-Summen der Rang-1-Reihe in verschiedenen Richtungen phi bei kleinem eps,
gegen die Kneser-artige Loesung von f_eps (e4_kneser_eps).  Borel-Pade [n/n].
B(s) = sum c_j s^j/j!,  S_phi(eps) = int_0^{inf e^{i phi}} e^{-s/eps} B(s) ds/eps (+ x+1).
Ausgabe: Pole der Pade-Approximanten und K - S_phi."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import numpy as np
import mpmath as mp
from e4_kneser_eps import Loeser, fixpunkt
from e4_theta import lade_reihe
mp.mp.dps = 50
J, W = lade_reihe('e4_rang1_J84.json')
EPS = [float(e) for e in sys.argv[1].split(',')]
XS = [0.25, 0.5, 0.75]
PHIS = [0.0, 0.1, 0.2, 0.3, 0.45, 0.6, 0.9, 1.2]
out = {}
L = 0.8247 + 1.5674j
for eps in EPS:
    L = fixpunkt(eps, L)
    S = Loeser(eps, 40.0, 700, L); S.loese(iters=6000, lam=0.3, tol=1e-15); S.normieren()
    out[str(eps)] = {}
    for x in XS:
        K = S.sexp(x)[0]
        b = [W[x][j] / mp.factorial(j + 1) for j in range(J)]   # Borel von sum c_j eps^j (j>=1): c_j -> s^{j-1}/(j-1)!? siehe unten
        # Reihe g(eps)=sum_{j>=1} c_j eps^j = eps * sum_{j>=0} c_{j+1} eps^j ; Borel von h(eps)=sum c_{j+1} eps^j:
        bh = [W[x][j] / mp.factorial(j) for j in range(J)]
        row = {}
        for n in (24, 27, 30, 32):
            try:
                p, q = mp.pade(bh[:2 * n + 1], n, n)
            except ZeroDivisionError:
                continue
            pole = sorted(mp.polyroots(q[::-1], maxsteps=300, extraprec=300), key=abs)[:4]
            vals = {}
            for phi in PHIS:
                w = mp.expj(phi)
                fn = lambda t: mp.exp(-t * w) * mp.polyval(p[::-1], eps * t * w) / mp.polyval(q[::-1], eps * t * w) * w
                h = mp.quad(fn, [0, 5, 10, 20, 40, 80, mp.inf])
                val = x + 1 + eps * h
                vals[phi] = val
            row[n] = vals
            print('eps=%.2f x=%.2f [%d/%d] Pole %s' % (eps, x, n, n, ' '.join(mp.nstr(z, 4) for z in pole)))
            print('     K - S_phi: ' + '  '.join('phi=%.1f:%s' % (phi, mp.nstr(K - v, 4)) for phi, v in vals.items()))
            print('     K - (S_phi+S_-phi)/2 (reell): ' + '  '.join('phi=%.1f:%s' % (phi, mp.nstr(mp.re(K - v), 4)) for phi, v in vals.items()), flush=True)
        out[str(eps)][str(x)] = {str(n): {str(phi): [float(mp.re(K - v)), float(mp.im(K - v))] for phi, v in vals.items()} for n, vals in row.items()}
(HIER / 'data' / ('e4_richtungen_%s.json' % '_'.join(sys.argv[1].split(',')))).write_text(json.dumps(out, indent=1), encoding='utf-8')
