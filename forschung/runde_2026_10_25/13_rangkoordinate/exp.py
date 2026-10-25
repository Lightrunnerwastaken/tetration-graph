"""Experimente zur Rangkoordinate.  python exp.py <basis> <dps> <grad> <R> <experiment...>
basis: sqrt2 | 1.5    experiment: base | grad | add | rank | half | halfB
"""
import sys, json, time
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import phi as P

HERE = Path(__file__).resolve().parent
basis, dps, grad, R = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
exps = sys.argv[5:]
mp.mp.dps = dps
X = 2
b = mp.sqrt(2) if basis == 'sqrt2' else mp.mpf(basis)
tag = f'{basis}_d{dps}_g{grad}_R{R}'


def cheb_from_coeffs(cs):
    c = P.Cheb.__new__(P.Cheb)
    c.a, c.b, c.n, c.c = mp.mpf(0), mp.mpf(X), len(cs), [mp.mpf(v) for v in cs]
    return c


def seed():
    """Startglied und dessen Rang: sqrt2 -> H_3=b^x; 1.5 -> H_4 = Kneser (Werte aus cache/sexp_1.5_30.json)."""
    if basis == 'sqrt2':
        return (lambda x: b ** x), 3
    c = json.loads((HERE / 'cache' / f'sexp_{basis}_30.json').read_text())
    n = grad
    ts = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / n) for k in range(n)]
    xs = [X * (t + 1) / 2 for t in ts]
    with mp.workdps(30):
        keys = [mp.nstr(x, 38, min_fixed=-mp.inf, max_fixed=mp.inf) for x in xs]
    vals = [mp.mpf(c[k]) for k in keys]
    table = dict(zip([mp.nstr(x, 15) for x in xs], vals))
    return P.Cheb(lambda x: table[mp.nstr(x, 15)], 0, X, n), 4


# Basisleiter mit gespeicherten Stufen
cf = P.CACHE / f'base_{tag}.json'
f0, r0 = seed()
if cf.is_file():
    B = json.loads(cf.read_text())
else:
    t = time.time()
    rows, coeffs, f = [], {}, f0
    for k in range(R - r0):
        g, p, lam = P.T(f, X, grad)
        rows.append(dict(k=k, p=str(p), lam=str(lam), v025=str(f(mp.mpf('0.25'))), v05=str(f(mp.mpf('0.5'))),
                         v15=str(f(mp.mpf('1.5')))))
        coeffs[r0 + k + 1] = [str(v) for v in g.c]
        print(r0 + k, mp.nstr(p, 15), mp.nstr(lam, 10), flush=True)
        f = g
    B = dict(rows=rows, coeffs=coeffs, sec=round(time.time() - t, 1))
    cf.write_text(json.dumps(B, indent=0))
    print('base', B['sec'], 's')
H = {r0: f0}
for r, cs in B['coeffs'].items():
    H[int(r)] = cheb_from_coeffs(cs)
base = P.lnd(B['rows'], b, r0)
NST = R - r0 - 1


def modes(spec):
    """spec: Liste (art, n, eps) -> psi(x) = x + sum eps*phi."""
    def psi(x):
        s = x
        for art, n, e in spec:
            e = mp.mpf(e)
            s += e * (mp.sin(2 * mp.pi * n * x) if art == 's' else 1 - mp.cos(2 * mp.pi * n * x))
        return s
    return psi


def phi_von(name, f, s0, nsteps=None):
    nsteps = nsteps or (R - int(mp.ceil(s0)) - 1)
    out = P.gecacht(f'{name}_{tag}', lambda: P.lauf(f, nsteps, X, grad))
    fo = P.phi_folge(base, out['rows'], b, s0)
    return fo, out


def stoer(r, spec):
    f = H[r]
    psi = modes(spec)
    return P.Cheb(lambda x: f(psi(x)), 0, X, grad)


def zeige(label, fo):
    print(label.ljust(28), ' '.join(f'{k}:{mp.nstr(v, 9)}' for k, v in sorted(fo.items())), flush=True)


res = {}
for e in exps:
    if e == 'base':
        print('Basis', B['sec'], 's'); continue
    if e == 'grad':
        # Gradient von Phi nach Fourier-Moden der theta-Freiheit an H_4 (bzw. H_r0+1)
        rr = 4 if basis == 'sqrt2' else 5
        for art in 's', 'c':
            for n in 1, 2, 3:
                vals = {}
                for sg in (+1, -1):
                    eps = sg * 0.005
                    fo, _ = phi_von(f'r{rr}_{art}{n}_{eps}', stoer(rr, [(art, n, eps)]), rr)
                    zeige(f'H{rr} {art}{n} eps={eps}', fo)
                    vals[sg] = fo
                ks = sorted(set(vals[1]) & set(vals[-1]))
                g = {k: (vals[1][k] - vals[-1][k]) / 0.01 - 0 for k in ks}
                # Phi(stoer) - Phi(H_rr) = (s_n - n) - rr ; Gradient = Differenz / (2 eps)
                zeige(f'  dPhi/d{art}{n}', g)
                res[f'grad_{rr}_{art}{n}'] = {k: str(v) for k, v in g.items()}
    if e == 'grad1':
        # nur s1 und c1 an H_4 (fuer b=1.5: Kneser-Saat), eps = +-0.005
        for art in 's', 'c':
            vals = {}
            for sg in (+1, -1):
                eps = sg * 0.005
                fo, _ = phi_von(f'r4_{art}1_{eps}', stoer(4, [(art, 1, eps)]), 4)
                zeige(f'H4 {art}1 eps={eps}', fo)
                vals[sg] = fo
            ks = sorted(set(vals[1]) & set(vals[-1]))
            g = {k: (vals[1][k] - vals[-1][k]) / 0.01 for k in ks}
            zeige(f'  dPhi/d{art}1', g)
            res[f'grad1_4_{art}1'] = {k: str(v) for k, v in g.items()}
    if e == 'add':
        rr = 4 if basis == 'sqrt2' else 5
        for spec, lab in [([('s', 1, 0.02)], 's1'), ([('c', 1, 0.02)], 'c1'), ([('s', 1, 0.02), ('c', 1, 0.02)], 's1c1'),
                          ([('s', 2, 0.02)], 's2'), ([('s', 1, 0.02), ('s', 2, 0.02)], 's1s2'),
                          ([('s', 1, 0.1)], 's1_0.1'), ([('s', 1, -0.1)], 's1_-0.1')]:
            fo, _ = phi_von(f'r{rr}_{lab}', stoer(rr, spec), rr)
            zeige(f'H{rr} {lab}', fo)
            res[f'add_{lab}'] = {k: str(v - rr) for k, v in fo.items()}
    if e == 'rank':
        # gleiche Stoerung s1 eps=+-0.005 an verschiedenen Raengen
        for rr in range(r0 + 1 if basis == 'sqrt2' else r0, r0 + 6):
            vals = {}
            for sg in (+1, -1):
                fo, _ = phi_von(f'r{rr}_s1_{sg*0.005}', stoer(rr, [('s', 1, sg * 0.005)]), rr)
                vals[sg] = fo
            ks = sorted(set(vals[1]) & set(vals[-1]))
            g = {k: (vals[1][k] - vals[-1][k]) / 0.01 for k in ks}
            zeige(f'dPhi/ds1 an H{rr}', g)
            res[f'rank_{rr}'] = {k: str(v) for k, v in g.items()}
    if e in ('half', 'halfB'):
        # Halbrang zwischen H_a und H_{a+1}: Familien lin (1-t)H_a + t H_{a+1}, geo H_a^{1-t} H_{a+1}^t
        a = r0 if e == 'half' else r0 + 1
        Fa, Fb = H[a], H[a + 1]
        for fam in 'lin', 'geo':
            for t in ('0.3', '0.5', '0.7'):
                tt = mp.mpf(t)
                if fam == 'lin':
                    g = P.Cheb(lambda x: (1 - tt) * Fa(x) + tt * Fb(x), 0, X, grad)
                else:
                    g = P.Cheb(lambda x: Fa(x) ** (1 - tt) * Fb(x) ** tt, 0, X, grad)
                fo, out = phi_von(f'half{a}_{fam}_{t}', g, a + 0.5)
                zeige(f'{fam} t={t} (Phi-{a}.5)', {k: v - (a + 0.5) for k, v in fo.items()})
                res[f'half{a}_{fam}_{t}'] = dict(phi={k: str(v) for k, v in fo.items()},
                                                 g05=str(g(mp.mpf('0.5'))), g15=str(g(mp.mpf('1.5'))),
                                                 p=[r['p'] for r in out['rows']])
json.dump(res, open(HERE / 'daten' / f'exp_{tag}_{"_".join(exps)}.json', 'w'), indent=1)
