"""Auswertung: Reste nach Abzug der Implosionsvorhersage, Anpassung in eps.

r_unten(y) = slog_b(y) - slog_minus(y)
r_oben(y)  = slog_b(y) - slog_plus(y) - [Re(2 pi i / log lambda_b) + C_ch - C_T]
"""
import json
from pathlib import Path
import mpmath as mp

HIER = Path(__file__).resolve().parent
mp.mp.dps = 60
# (nenner, beste dps, Kontroll-dps)
WAHL = [(10, 30, None), (30, 30, None), (100, 30, None), (300, 30, None), (1000, 45, 30), (3000, 60, 45), (10000, 100, 80), (30000, 200, 150)]
CT = mp.mpf('2.7982481542313876624525831854690855')
CCH = mp.mpf('-1.2569611588725123189746884727585369')


def lade(n, d):
    return json.loads((HIER / 'cache' / f'phase_{d}_{n}.json').read_text(encoding='utf-8'))


def tabelle():
    zeilen = []
    for n, d, dk in WAHL:
        r = lade(n, d)
        eps = mp.mpf(r['eps'])
        R = mp.mpf(r['re_2pii_lnlam'])
        I = mp.mpf(r['im_2pii_lnlam'])
        lead = mp.mpf(r['lead'])
        z = {'n': n, 'eps': eps, 'R': R, 'I': I, 'lead': lead}
        for w in r['werte']:
            v = mp.mpf(w['diff'])
            if w['seite'] == 'oben':
                v = v - (R + CCH - CT)
            z[w['y']] = v
            if dk:
                v2 = mp.mpf([x for x in lade(n, dk)['werte'] if x['y'] == w['y']][0]['diff'])
                if w['seite'] == 'oben':
                    v2 = v2 - (R + CCH - CT)
                z['kontrolle_' + w['y']] = abs(v - v2)
        zeilen.append(z)
    return zeilen


def fit(xs, ys, exps):
    """kleinste Quadrate exakt (so viele Punkte wie Unbekannte: Interpolation der letzten Punkte)."""
    A = mp.matrix([[x**e for e in exps] for x in xs])
    b = mp.matrix(ys)
    return mp.lu_solve(A, b)


if __name__ == '__main__':
    Z = tabelle()
    print('n      eps          Im(2pi i/log lam)-pi/3   lead-R')
    for z in Z:
        print(z['n'], mp.nstr(z['eps'], 8), mp.nstr(z['I'] - mp.pi / 3, 8), mp.nstr(z['lead'] - z['R'], 8))
    for y in ['1/2', '2', 'exp(2)', '10', '2*exp(1)']:
        print('\n== y =', y, ' (r, r/eps, Kontrolle)')
        for z in Z:
            print(z['n'], mp.nstr(z[y], 15), mp.nstr(z[y] / z['eps'], 15), mp.nstr(z.get('kontrolle_' + y, mp.mpf(0)), 3))
        # Anpassung mit freier Konstante a0 und Potenzen eps^{1/2}
        for exps in ([0, 1, 2, 3, 4], [0, 0.5, 1, 1.5, 2, 2.5], [0, 1, 1.5, 2, 3]):
            sub = Z[-len(exps):]
            c = fit([z['eps'] for z in sub], [z[y] for z in sub], exps)
            print('  Potenzen', exps, '->', [mp.nstr(ci, 8) for ci in c])
