"""Randgesetz exakt: m(b) = lim_{y->oo} (slog_e(y) - slog_b(y)) gegen -Re(2 pi i/log lambda_b) + Konstante.

y = exp(exp(Y)):  slog_e(y) = 2 + slog_e(Y),  slog_b(y) = 2 + slog_b((Y - ln ln b)/ln b),  slog_plus(y) = 2 + slog_plus(e(Y+1)).
Aufruf: python mu_rand.py Y
"""
import json
import sys
from pathlib import Path
import mpmath as mp

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import motor
import eta_seite as es

WAHL = [(10, 30), (30, 30), (100, 30), (300, 30), (1000, 45), (3000, 60), (10000, 100), (30000, 200)]


def slog_mit_log(bgp, lnb, v, dps):
    k = 0
    while v > 1:
        v = mp.log(v) / lnb
        k += 1
    w = motor.werte(bgp, [f'slog({motor.gpzahl(v, dps + 10)})'], dps)[0]
    return mp.re(w) + k


if __name__ == '__main__':
    Y = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf(100)
    mp.mp.dps = 60
    CT = es.phi_unten(mp.mpf(1), 3000)
    CCH = es.phi_oben(mp.e ** 2, 3000)
    se = slog_mit_log('exp(1)', mp.mpf(1), Y, 45)
    mu_ech_a = se - (es.phi_oben(mp.e * (Y + 1), 3000) - CCH)
    v = mp.exp(mp.exp(mp.e - 1) - 1) - 1
    mu_ech_b = slog_mit_log('exp(1)', mp.mpf(1), v, 45) - 2
    print('mu_{e,cheta} (Y-Weg)    =', mp.nstr(mu_ech_a, 30))
    print('mu_{e,cheta} (exakt-Weg) =', mp.nstr(mu_ech_b, 30), ' Diff', mp.nstr(mu_ech_a - mu_ech_b, 3))
    K = mu_ech_b + CT - CCH
    print('Vorhersage-Konstante mu_{e,cheta} + C_T - C_ch =', mp.nstr(K, 30))
    aus = {'mu_e_cheta': mp.nstr(mu_ech_b, 40), 'K': mp.nstr(K, 40), 'zeilen': []}
    for n, d in WAHL:
        mp.mp.dps = d + 15
        eta = mp.e ** (1 / mp.e)
        b = eta + mp.mpf(1) / n
        lnb = mp.log(b)
        bgp = f"(exp(exp(-1))+1/{n})"
        sb = slog_mit_log(bgp, lnb, (Y - mp.log(lnb)) / lnb, d)
        m = se - sb
        eps = mp.e * lnb - 1
        lam = -mp.lambertw(-lnb)
        R = mp.re(2j * mp.pi / mp.log(lam))
        r = m + R - K
        print(n, 'm(b) =', mp.nstr(m, 20), '  m+R-K =', mp.nstr(r, 12), '  (m+R-K)/eps =', mp.nstr(r / eps, 10))
        aus['zeilen'].append({'n': n, 'm': mp.nstr(m, 30), 'rest': mp.nstr(r, 20), 'eps': mp.nstr(eps, 20)})
    (HIER / 'cache' / 'mu_rand.json').write_text(json.dumps(aus, indent=1), encoding='utf-8')
    motor.schliessen()
