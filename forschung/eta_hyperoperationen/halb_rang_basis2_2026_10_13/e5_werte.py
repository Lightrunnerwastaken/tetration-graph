"""Etappe 5: volle Rechenarten H_N(2,z), N=3..Nmax, an einem komplexen Raster (Motor der Vorrunde).
Aufruf: python e5_werte.py Nmax dps"""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from e5_motor import motor
Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 30
mp.mp.dps = dps + 10
XS = ['0.25', '0.5', '0.75', '1.25', '1.5', '1.75', '-0.5']
YS = (sys.argv[3].split(',') if len(sys.argv) > 3 else ['0', '0.25', '0.5', '1', '1.5', '2', '3', '4'])
t0 = time.time(); gp = motor(dps, 1100)
out = {'Nmax': Nmax, 'dps': dps, 'werte': {}}
try:
    v = gp.eval_batch('2', ['hinit(%d)' % Nmax, 'hrerror[%d]' % Nmax])
    out['init'] = [str(x) for x in v]; print('init', v, round(time.time() - t0, 1), flush=True)
    for x in XS:
        for y in YS:
            z = '%s+%s*I' % (x, y)
            ex = ['iferr(hval(%d,%s),E,10^99)' % (N, z) for N in range(3, Nmax + 1)]
            vals = gp.eval_batch('2', ex)
            out['werte'][x + '|' + y] = [[mp.nstr(mp.re(w), dps), mp.nstr(mp.im(w), dps)] for w in vals]
        print('x', x, round(time.time() - t0, 1), flush=True)
finally:
    gp.close()
    (HIER / 'data' / ('e5_werte_N%d.json' % Nmax)).write_text(json.dumps(out, indent=1), encoding='utf-8')
    print('fertig', round(time.time() - t0, 1))
