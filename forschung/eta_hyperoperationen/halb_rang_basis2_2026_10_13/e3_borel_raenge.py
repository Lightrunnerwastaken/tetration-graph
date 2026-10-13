"""Etappe 3: Borel-Analyse der eps-Reihe je Rang (n=1,2): Gliedgroessen, Quotient Rang2/Rang1 je Ordnung,
Pole der Borel-Pade-Approximanten, Periodizitaet der spaeten Glieder."""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
import mpmath as mp
mp.mp.dps = 60
d = json.loads((HIER / 'data' / 'e3_raenge_J36.json').read_text())
J = d['J']
C = {n: {p: [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in d['raenge'][n][p]] for p in d['punkte']} for n in d['raenge']}
out = {}
for p in ['0.5', '0.25', '0.5+0.5j']:
    print('x =', p)
    for n in C:
        print('  Rang', n, 'log10|c_j|:', ' '.join('%.1f' % float(mp.log10(abs(c))) for c in C[n][p]))
    q = [C['2'][p][j] / C['1'][p][j] for j in range(J)]
    print('  c2/c1 je Ordnung:', ' '.join(mp.nstr(abs(v), 3) for v in q))
    out[p] = {'quot_rang2_rang1': [float(abs(v)) for v in q]}
# Periodizitaet der spaeten Glieder, Rang 2
for n in C:
    for j in [J - 6, J - 1]:
        r = C[n]['0.5'][j]
        print('Rang', n, 'Ordnung', j + 1, ' c(x)/c(0.5):', ' '.join('%s:%s' % (p, mp.nstr(C[n][p][j] / r, 4)) for p in d['punkte']))
# Borel-Pade-Pole
for n in C:
    for p in ['0.5', '0.5+1j']:
        c = [mp.mpc(complex(p)) + 1] + C[n][p]
        b = [c[j] / mp.factorial(j) for j in range(J + 1)]
        for L in [16, 17]:
            P_, Q = mp.pade(b, L, J - L)
            r = sorted(mp.polyroots(Q[::-1], maxsteps=400, extraprec=400), key=abs)[:6]
            print('Rang', n, 'x', p, '[%d/%d] Borel-Pole:' % (L, J - L), ' '.join(mp.nstr(z, 3) for z in r))
            out.setdefault('pole', {})['%s|%s|%d' % (n, p, L)] = [[float(mp.re(z)), float(mp.im(z))] for z in r]
(HIER / 'data' / 'e3_borel_raenge.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
