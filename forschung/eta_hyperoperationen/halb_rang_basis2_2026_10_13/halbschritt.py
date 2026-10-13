"""Rangfortsetzung n -> s der ganzzahligen Jets durch die Newtonreihe im Rang.

c_s = sum_p binom(s,p) Delta^p c(0)   (Auswahl: Rangfortsetzung polynomialen Wachstums, wie (A8)).
Gibt Partialsummen gegen die Abschneidung P aus.
"""
import sys, json, argparse
from pathlib import Path
import mpmath as mp
HIER = Path(__file__).resolve().parent


def newton_koeff(seq):
    """Vorwaertsdifferenzen Delta^p c(0), p=0..len-1."""
    d = list(seq); out = []
    while d:
        out.append(d[0]); d = [b - a for a, b in zip(d, d[1:])]
    return out


def newton_partial(seq, s, base=0):
    """Partialsummen der Newtonreihe zur Basis 'base' (Folge beginnt bei Rang base)."""
    co = newton_koeff(seq)
    tot = mp.mpc(0); b = mp.mpf(1); out = []
    t = mp.mpf(s) - base
    for p, c in enumerate(co):
        if p:
            b *= (t - p + 1) / p
        tot += b * c
        out.append(tot)
    return out


def lade(name):
    return json.loads((HIER / 'data' / f'raenge_{name}.json').read_text())


def folge(d, punkt, j, r=0):
    return [mp.mpc(mp.mpf(row[r][j][0]), mp.mpf(row[r][j][1])) for row in d['werte'][punkt]]


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('name'); ap.add_argument('--dps', type=int, default=300)
    ap.add_argument('--s', default='0.5')
    a = ap.parse_args(); mp.mp.dps = a.dps
    d = lade(a.name)
    for p in d['punkte']:
        print('x =', p)
        for j in range(d['J']):
            ps = newton_partial(folge(d, p, j), mp.mpf(a.s))
            N = len(ps) - 1
            idx = [N // 8, N // 4, N // 2, 3 * N // 4, N]
            print('  j=%d' % (j + 1), '  '.join('P%d:%s' % (i, mp.nstr(mp.re(ps[i]), 10)) for i in idx),
                  ' |dP|:', mp.nstr(abs(ps[N] - ps[N // 2]), 3))


def levin_wert(seq, s, variant='u'):
    """Levin-Transformation der Newton-Partialsummen (Rangfortsetzung bei s)."""
    ps = newton_partial(seq, s)
    while len(ps) > 1 and ps[0] == 0:
        ps = ps[1:]
    L = mp.levin(method='levin', variant=variant)
    v, e = L.update_psum(ps)
    return v, e
