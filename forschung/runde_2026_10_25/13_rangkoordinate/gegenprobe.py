"""Gegenprobe beim Eintragen (25.10.2026): B2 (Gradient s1, c1 an H_4) und B7 (Halbrang-Leitern) frisch rechnen,
ohne den Zwischenspeicher cache/ (eigener Speicher gegenprobe_cache/). Benutzt phi.py unveraendert.
Aufruf: python gegenprobe.py  (b = sqrt2, dps 20, Grad 30, R 15 wie im Bericht)"""
import sys, time
from pathlib import Path
import mpmath as mp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import phi as P
P.CACHE = HERE / 'gegenprobe_cache'
P.CACHE.mkdir(exist_ok=True)
mp.mp.dps = 20
X, grad, R, b = 2, 30, 15, mp.sqrt(2)
t0 = time.time()

# Basisleiter frisch
f = lambda x: b ** x
H3 = f
rows, H = [], {}
for k in range(R - 3):
    g, p, lam = P.T(f, X, grad)
    rows.append(dict(k=k, p=str(p)))
    H[4 + k] = g
    f = g
base = P.lnd(rows, b, 3)
print('Basis', round(time.time() - t0, 1), 's; p_3..p_6 =', [mp.nstr(mp.mpf(r['p']), 12) for r in rows[:4]], flush=True)


def stoer(art, eps):
    e = mp.mpf(eps)
    if art == 's':
        psi = lambda x: x + e * mp.sin(2 * mp.pi * x)
    else:
        psi = lambda x: x + e * (1 - mp.cos(2 * mp.pi * x))
    return P.Cheb(lambda x: H[4](psi(x)), 0, X, grad)


for art in 's', 'c':
    fo = {}
    for sg in (+1, -1):
        out = P.lauf(stoer(art, sg * 0.005), R - 4 - 1, X, grad)
        fo[sg] = P.phi_folge(base, out, b, 4)
    ks = sorted(set(fo[1]) & set(fo[-1]))
    print(f'dPhi/d{art}1 an H4:', ' '.join(f'{k}:{mp.nstr((fo[1][k] - fo[-1][k]) / 0.01, 8)}' for k in ks),
          f'  ({round(time.time() - t0, 1)} s)', flush=True)


def interp_inv(target, guess, npt):
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


ps = {}
for fam, t in (('lin', '0.62762'), ('geo', '0.62300')):
    tt = mp.mpf(t)
    if fam == 'lin':
        g = P.Cheb(lambda x: (1 - tt) * H3(x) + tt * H[4](x), 0, X, grad)
    else:
        g = P.Cheb(lambda x: H3(x) ** (1 - tt) * H[4](x) ** tt, 0, X, grad)
    out = P.lauf(g, 11, X, grad)
    ph = {}
    for r in out:
        k = r['k']
        if 3.5 + k + 1 > max(base) - 2 or 3.5 + k < min(base) + 2:
            continue
        ph[k] = interp_inv(mp.log(mp.mpf(r['p']) - b), 3.5 + k, 6) - k
    ps[fam] = [mp.mpf(r['p']) for r in out]
    print(fam, 't*=', t, ' g*(0.5) =', mp.nstr(g(mp.mpf('0.5')), 10), ' g*(1.5) =', mp.nstr(g(mp.mpf('1.5')), 10),
          ' Phi (6 Pkt):', ' '.join(f'{k}:{mp.nstr(v, 8)}' for k, v in sorted(ph.items())), f'  ({round(time.time() - t0, 1)} s)', flush=True)
print('n  p_lin  p_lin-p_geo  relativ')
for n, (pl, pg) in enumerate(zip(ps['lin'], ps['geo'])):
    print(n, mp.nstr(pl, 10), mp.nstr(pl - pg, 3), mp.nstr((pl - pg) / (pl - b), 3))
print('gesamt', round(time.time() - t0, 1), 's')
