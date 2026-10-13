"""Etappe 4: seitliche Borel-Laplace-Summen der Rang-1-Reihe, Median, Stokes-Sprung, Vergleich mit Kneser.

S_pm(x) = int_0^inf e^{-t} B(t) dt laengs arg t = +-phi, B = Pade der Borel-Transformierten.
Median M = (S_+ + S_-)/2, Sprung D = S_+ - S_-.  Kneser - M wird auf theta-Moden projiziert.
Aufruf: python e4_theta.py e4_rang1_J100.json
"""
import sys, json
from pathlib import Path
HIER = Path(__file__).resolve().parent
import mpmath as mp
mp.mp.dps = 60


def lade_reihe(name):
    d = json.loads((HIER / 'data' / name).read_text())
    W = {}
    for k, row in d['werte'].items():
        z = complex(eval(k)) if not isinstance(k, complex) else k
        W[z] = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b, _ in row]
    return d['J'], W


def lade_kneser():
    K = {}
    for f in ['kneser_basis2.json', 'e4_kneser_raster.json']:
        for w in json.loads((HIER / 'data' / f).read_text())['werte']:
            K[complex(float(w['x'][0]), float(w['x'][1]))] = mp.mpc(mp.mpf(w['sexp'][0]), mp.mpf(w['sexp'][1]))
    return K


def borel_pade(c, L, M):
    b = [c[j] / mp.factorial(j) for j in range(L + M + 1)]
    p, q = mp.pade(b, L, M)
    return p, q


def laplace(p, q, phi):
    w = mp.expj(phi)
    f = lambda t: mp.exp(-t * w) * mp.polyval(p[::-1], t * w) / mp.polyval(q[::-1], t * w) * w
    return mp.quad(f, [0, 4, 8, 12, 20, 40, mp.inf])


def summen(x, terme, L, M, phi):
    c = [mp.mpc(x) + 1] + terme
    p, q = borel_pade(c, L, M)
    sp = laplace(p, q, phi); sm = laplace(p, q, -phi)
    return sp, sm


if __name__ == '__main__':
    name = sys.argv[1] if len(sys.argv) > 1 else 'e4_rang1_J100.json'
    J, W = lade_reihe(name)
    K = lade_kneser()
    phi = mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf('0.3')
    orders = [(n, n) for n in range(20, J // 2 + 1, 4)]
    res = {}
    for x in [0.25, 0.5, 0.75, -0.5, 1.5, 0.5 + 0.5j]:
        if x not in W:
            continue
        row = []
        for L, M in orders:
            try:
                sp, sm = summen(x, W[x], L, M, phi)
            except ZeroDivisionError:
                print('x=%s [%d/%d] Pade singulaer' % (x, L, M)); continue
            med = (sp + sm) / 2
            row.append((L, M, med - K[x], sp - sm))
            print('x=%-9s [%d/%d] Median-Kneser %s  Sprung S+-S- %s' % (x, L, M, mp.nstr(med - K[x], 6), mp.nstr(sp - sm, 6)), flush=True)
        res[str(x)] = [[L, M, [float(mp.re(a)), float(mp.im(a))], [float(mp.re(b)), float(mp.im(b))]] for L, M, a, b in row]
    (HIER / 'data' / ('e4_theta_' + name.replace('.json', '') + '_phi%s.json' % mp.nstr(phi, 3))).write_text(json.dumps(res, indent=1), encoding='utf-8')
