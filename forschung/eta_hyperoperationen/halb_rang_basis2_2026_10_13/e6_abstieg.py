"""Etappe 6: Halbrang durch Rangabstieg (Fatou-artig):  H_{3.5} = R^k(G),  R(S)(y)=S(S^{-1}(y)+1),
G = Newton-Interpolation der Raenge N0..9 bei N=3.5+k (moeglichst mittig).
R ist unabhaengig von der theta-Wahl der hohen Stufen: R(S o psi)=R(S) fuer psi(x+1)=psi(x)+1.
Werte der Raenge ueber den Motor der Vorrunde (e5_motor), nur reelle Hoehen.
Aufruf: python e6_abstieg.py N0 k"""
import sys, json, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
import mpmath as mp
from e5_motor import motor
mp.mp.dps = 35
FAELLE = [(int(a.split(':')[0]), int(a.split(':')[1])) for a in (sys.argv[1].split(',') if len(sys.argv) > 1 else ['5:4', '5:3', '4:3', '5:5'])]
XS = [mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), mp.mpf('1.5')]
t0 = time.time(); gp = motor(25, 1100)
cache = {}
def hv(N, x):
    key = (N, mp.nstr(x, 30))
    if key not in cache:
        cache[key] = mp.re(gp.eval_batch('2', ['hval(%d,%s)' % (N, mp.nstr(x, 30))])[0])
    return cache[key]
def gewichte(nodes, s):
    # Lagrange-Gewichte
    w = []
    for i, a in enumerate(nodes):
        p = mp.mpf(1)
        for j, b in enumerate(nodes):
            if j != i: p *= (s - b) / mp.mpf(a - b)
        w.append(p)
    return w
out = {}
try:
    gp.eval_batch('2', ['hinit(9)'])
    print('init', round(time.time() - t0), 's', flush=True)
    for N0, k in FAELLE:
        nodes = list(range(N0, 10)); s = mp.mpf(3.5) + k
        W = gewichte(nodes, s)
        G = lambda y: mp.fsum(w * hv(N, y) for w, N in zip(W, nodes))
        def Ginv(v, y0):
            y = y0
            for _ in range(40):
                g = G(y) - v
                d = (G(y + mp.mpf('1e-8')) - G(y - mp.mpf('1e-8'))) / mp.mpf('2e-8')
                y -= g / d
                if abs(g) < mp.mpf(10) ** -24: break
            return y
        def S(j, y):
            if j == 0: return G(y)
            return S(j - 1, Sinv(j - 1, y) + 1)
        def Sinv(j, y):
            if j == 0: return Ginv(y, y - 1)
            return S(j - 1, Sinv(j - 1, y) - 1)
        res = {}
        for x in XS:
            v = S(k, x)
            res[str(x)] = str(v)
        out['%d:%d' % (N0, k)] = {'s_start': str(s), 'knoten': nodes, 'H35': res}
        print('Knoten %d..9, Start bei N=%s, %d Abstiege: ' % (N0, mp.nstr(s, 3), k) + '  '.join('H3.5(%s)=%s' % (mp.nstr(x, 3), mp.nstr(mp.mpf(v), 10)) for x, v in res.items()), ' %.0fs' % (time.time() - t0), flush=True)
finally:
    gp.close()
    (HIER / 'data' / 'e6_abstieg.json').write_text(json.dumps(out, indent=1), encoding='utf-8')
