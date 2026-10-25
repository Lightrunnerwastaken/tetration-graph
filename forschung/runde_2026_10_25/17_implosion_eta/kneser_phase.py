"""Kneser-slog nahe eta gegen die Implosionsvorhersage.

Aufruf: python kneser_phase.py DPS NENNER [NENNER ...]   (Basis b = eta + 1/NENNER)

Fuer jede Basis:
  A(y)  = slog_b(y) - slog_minus(y)                      (y < e, sollte -> 0)
  P(y)  = slog_b(y) - slog_plus(y)                       (y > e, sollte -> Phase P(b), unabhaengig von y)
  slog_b(y) fuer y>e: k-mal log_b bis unter 1, dann Motor-slog, plus k.
Ausgabe als JSON in cache/phase_<dps>_<nenner>.json
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import motor  # eigene Kopie, eigener Speicher
import eta_seite as es

Y_UNTEN = ['1/2', '2', '5/2']
Y_OBEN = ['exp(2)', '10', '2*exp(1)']


def lauf(dps, nenner):
    mp.mp.dps = dps + 15
    eta = mp.e ** (1 / mp.e)
    d = mp.mpf(1) / nenner
    b = eta + d
    lnb = mp.log(b)
    bgp = f"(exp(exp(-1))+1/{nenner})"
    ausdr = ['slog(1)', f'slog({bgp})']
    plan = []
    for ys in Y_UNTEN:
        y = mp.mpf(eval(ys.replace('exp', 'mp.exp'), {'mp': mp}))
        ausdr.append(f'slog({motor.gpzahl(y, dps + 10)})')
        plan.append(('unten', ys, y, 0))
    for ys in Y_OBEN:
        y0 = mp.mpf(eval(ys.replace('exp', 'mp.exp'), {'mp': mp}))
        y, k = y0, 0
        while y > 1:
            y = mp.log(y) / lnb
            k += 1
        ausdr.append(f'slog({motor.gpzahl(y, dps + 10)})')
        plan.append(('oben', ys, y0, k))
    t = time.time()
    w = motor.werte(bgp, ausdr, dps)
    zeit = time.time() - t
    eps = mp.e * lnb - 1
    lam = -mp.lambertw(-lnb)            # Multiplikator am Fixpunkt L (Im>0 Zweig)
    if mp.im(lam) < 0:
        lam = mp.conj(lam)
    res = {'nenner': nenner, 'dps': dps, 'eps': mp.nstr(eps, 30), 'zeit': zeit,
           'slog1': mp.nstr(w[0], 10), 'slogb_minus1': mp.nstr(w[1] - 1, 10),
           'lead': mp.nstr(mp.pi * mp.sqrt(2 / eps), 30),
           're_2pii_lnlam': mp.nstr(mp.re(2j * mp.pi / mp.log(lam)), 30),
           'im_2pii_lnlam': mp.nstr(mp.im(2j * mp.pi / mp.log(lam)), 30),
           'werte': []}
    for (seite, ys, y, k), v in zip(plan, w[2:]):
        s = mp.re(v) + k
        if seite == 'unten':
            ref = es.phi_unten(y, 3000) - es.phi_unten(mp.mpf(1), 3000)
        else:
            ref = es.phi_oben(y, 3000) - es.phi_oben(mp.e ** 2, 3000)
        res['werte'].append({'seite': seite, 'y': ys, 'k': k, 'slog_b': mp.nstr(s, dps), 'diff': mp.nstr(s - ref, dps),
                             'im': mp.nstr(mp.im(v), 5)})
    return res


if __name__ == '__main__':
    dps = int(sys.argv[1])
    for n in sys.argv[2:]:
        r = lauf(dps, int(n))
        (HIER / 'cache' / f'phase_{dps}_{n}.json').write_text(json.dumps(r, indent=1), encoding='utf-8')
        print(json.dumps(r, indent=1), flush=True)
    motor.schliessen()
