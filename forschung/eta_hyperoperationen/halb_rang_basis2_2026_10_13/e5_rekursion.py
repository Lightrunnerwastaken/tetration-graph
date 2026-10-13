"""Etappe 5: Rangrekursion H_4.5(x+1) = H_3.5(H_4.5(x)) fuer die Rangfortsetzung der vollen Werte.
H_3.5 = Newton durch N=3..Nmax bei 3.5; H_4.5 = Newton durch N=4..Nmax bei 4.5.
Aufruf: python e5_rekursion.py Nmax dps"""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from e5_motor import motor
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
mp.mp.dps = dps + 10
def newton(seq, t):
    co = []; s = list(seq)
    while s:
        co.append(s[0]); s = [b - c for c, b in zip(s, s[1:])]
    tot = 0; b = mp.mpf(1)
    for p, c in enumerate(co):
        if p: b *= (t - p + 1) / p
        tot += b * c
    return tot
def gpz(z):
    z = mp.mpc(z)
    return '(%s)+(%s)*I' % (mp.nstr(z.real, dps + 5), mp.nstr(z.imag, dps + 5))
t0 = time.time(); gp = motor(dps, 1100); out = {}
try:
    gp.eval_batch('2', ['hinit(%d)' % Nmax])
    H = lambda z: gp.eval_batch('2', ['hval(%d,%s)' % (N, gpz(z)) for N in range(3, Nmax + 1)])
    for x in ['0.25', '0.5', '0.75', '0.25+0.25*I', '0.5+0.5*I']:
        z = mp.mpc(complex(x.replace('*I', 'j')))
        hx = H(z); hx1 = H(z + 1)
        h45x = newton(hx[1:], mp.mpf('0.5')); h45x1 = newton(hx1[1:], mp.mpf('0.5'))
        hw = H(h45x); h35w = newton(hw, mp.mpf('0.5'))
        # Kontrolle mit ganzzahligen Raengen: H_5(x+1) - H_4(H_5(x))
        kontrolle = hx1[2] - gp.eval_batch('2', ['hval(4,%s)' % gpz(hx[2])])[0]
        r = h45x1 - h35w
        out[x] = {'H45_x': str(h45x), 'H45_x1': str(h45x1), 'H35_w': str(h35w), 'residuum': [float(mp.re(r)), float(mp.im(r))], 'kontrolle_ganz': float(abs(kontrolle))}
        print('x=%-12s H4.5(x)=%s  Residuum %s  (ganzzahlig %.1e)  %.0fs' % (x, mp.nstr(h45x, 9), mp.nstr(r, 3), float(abs(kontrolle)), time.time() - t0), flush=True)
finally:
    gp.close()
    (HIER / 'data' / ('e5_rekursion_N%d.json' % Nmax)).write_text(json.dumps(out, indent=1), encoding='utf-8')
