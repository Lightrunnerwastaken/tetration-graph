"""Halbrang H_{3.5} bei b=sqrt2: in zwei Familien (lin, geo) das t* mit Phi(g_t*)=3.5 suchen,
dann die Leitern T^n g_t* vergleichen.  python half2.py  (benutzt cache/ aus exp.py sqrt2 20 30 15 base half)"""
import sys, json
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
mp.mp.dps = 20
import phi as P

HERE = Path(__file__).resolve().parent
X, grad, R, b = 2, 30, 15, mp.sqrt(2)
tag = 'sqrt2_d20_g30_R15'
B = json.loads((P.CACHE / f'base_{tag}.json').read_text())
base = P.lnd(B['rows'], b, 3)


def chb(cs):
    c = P.Cheb.__new__(P.Cheb)
    c.a, c.b, c.n, c.c = mp.mpf(0), mp.mpf(X), len(cs), [mp.mpf(v) for v in cs]
    return c


H3, H4 = (lambda x: b ** x), chb(B['coeffs']['4'])


def interp_inv(base, target, guess, npt):
    rs = sorted(base)

    def L(t):
        near = sorted(rs, key=lambda s: abs(s - t))[:npt]
        tot = 0
        for i in near:
            w = 1
            for j in near:
                if j != i:
                    w *= (t - j) / (i - j)
            tot += w * base[i]
        return tot
    return mp.findroot(lambda t: L(t) - target, mp.mpf(guess))


def phi(rows, s0, npt):
    out = {}
    for r in rows:
        k = r['k']
        if s0 + k + 1 > max(base) - 2 or s0 + k < min(base) + 2:
            continue
        out[k] = interp_inv(base, mp.log(mp.mpf(r['p']) - b), s0 + k, npt) - k
    return out


D = json.loads(open(HERE / 'daten' / f'exp_{tag}_half.json').read())
erg = {}
for fam in 'lin', 'geo':
    ts = [mp.mpf(t) for t in ('0.3', '0.5', '0.7')]
    for npt in (4, 6):
        ph = []
        for t in ('0.3', '0.5', '0.7'):
            rows = json.loads((P.CACHE / f'half3_{fam}_{t}_{tag}.json').read_text())['rows']
            f = phi(rows, 3.5, npt)
            ph.append(f[max(f)] - mp.mpf(3.5))
        # quadratische Interpolation phi(t) und Nullstelle
        poly = lambda t: sum(ph[i] * mp.fprod((t - ts[j]) / (ts[i] - ts[j]) for j in range(3) if j != i) for i in range(3))
        ts_ = mp.findroot(poly, mp.mpf('0.62'))
        print(fam, 'npt', npt, 'Phi-3.5 bei t=.3,.5,.7:', [mp.nstr(v, 7) for v in ph], ' t* =', mp.nstr(ts_, 8), flush=True)
        erg[(fam, npt)] = ts_
    t6 = erg[(fam, 6)]
    if fam == 'lin':
        g = P.Cheb(lambda x: (1 - t6) * H3(x) + t6 * H4(x), 0, X, grad)
    else:
        g = P.Cheb(lambda x: H3(x) ** (1 - t6) * H4(x) ** t6, 0, X, grad)
    out = P.gecacht(f'half3_{fam}_star_{tag}', lambda: P.lauf(g, 11, X, grad))
    f4, f6 = phi(out['rows'], 3.5, 4), phi(out['rows'], 3.5, 6)
    print(' Phi(g*) (4 Pkt)', ' '.join(f'{k}:{mp.nstr(v, 8)}' for k, v in sorted(f4.items())))
    print(' Phi(g*) (6 Pkt)', ' '.join(f'{k}:{mp.nstr(v, 8)}' for k, v in sorted(f6.items())))
    print(' g*(0.5) =', mp.nstr(g(mp.mpf('0.5')), 12), ' g*(1.5) =', mp.nstr(g(mp.mpf('1.5')), 12),
          ' g*(0.25) =', mp.nstr(g(mp.mpf('0.25')), 12))
    erg[fam] = out['rows']

print('Vergleich der beiden Halbrang-Leitern (lin* gegen geo*): Rang 3.5+n')
for rl, rg in zip(erg['lin'], erg['geo']):
    pl, pg = mp.mpf(rl['p']), mp.mpf(rg['p'])
    vl, vg = mp.mpf(rl['v05']), mp.mpf(rg['v05'])
    print(f" n={rl['k']:2d}  p_lin={mp.nstr(pl, 14)}  p_lin-p_geo={mp.nstr(pl - pg, 4)}  "
          f"(p_lin-p_geo)/(p_lin-b)={mp.nstr((pl - pg) / (pl - b), 4)}  v05 diff={mp.nstr(vl - vg, 4)}")
