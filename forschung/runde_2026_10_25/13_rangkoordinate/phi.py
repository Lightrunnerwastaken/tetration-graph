"""Rangkoordinate Phi auf beschraenkten Rangleitern (Gebiet 13, Runde 25.10.2026).

Koenigs-Schritt T (wie runde_2026_10_18/07_raenge/koenigs_leiter.py): fuer f mit kleinstem
Fixpunkt p >= 1, lam = f'(p) < 1:  (T f)(x) = sigma^{-1}(lam^x sigma(1)),  T f(0)=1, Tf(x+1)=f(Tf(x)).

Phi(g) := lim_n [ s_n - n ],  s_n aus  ln(p(T^n g) - b) = L(s_n),  L = Interpolation (4 Punkte, in r)
der kanonischen Folge r -> ln(p_r - b), p_r = Fixpunkt von H_r (H_3 = b^x).  Normierung Phi(H_r) = r.

Aufruf: python phi.py <experiment>   (siehe unten).  Ergebnisse -> daten/<experiment>.json, Laeufe in cache/.
"""
import sys, json, os, time, hashlib
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
CACHE = HERE / 'cache'
CACHE.mkdir(exist_ok=True)


class Cheb:
    def __init__(self, f, a, b, n):
        self.a, self.b, self.n = mp.mpf(a), mp.mpf(b), n
        xs = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / n) for k in range(n)]
        fx = [f(self.a + (self.b - self.a) * (t + 1) / 2) for t in xs]
        self.c = []
        for j in range(n):
            s = mp.fsum(fx[k] * mp.cos(mp.pi * j * (k + mp.mpf(1) / 2) / n) for k in range(n))
            self.c.append(2 * s / n)
        self.c[0] /= 2

    def __call__(self, x):
        t = (2 * x - self.a - self.b) / (self.b - self.a)
        b1 = b2 = mp.mpf(0)
        for cj in reversed(self.c[1:]):
            b1, b2 = 2 * t * b1 - b2 + cj, b1
        return t * b1 - b2 + self.c[0]


def fixpunkt(f, tol):
    y = mp.mpf(1)
    for _ in range(4000):
        y2 = f(y)
        if abs(y2 - y) < mp.mpf(10) ** (-12):
            break
        y = y2
    for _ in range(60):
        g = f(y) - y
        d = mp.diff(f, y) - 1
        y2 = y - g / d
        if abs(y2 - y) < tol:
            return y2
        y = y2
    return y


def T(f, X, grad):
    """Ein Koenigs-Schritt. Rueckgabe (Tf als Cheb, p, lam) -- p, lam gehoeren zu f."""
    tol = mp.mpf(10) ** (-(mp.mp.dps - 5))
    p = fixpunkt(f, tol)
    lam = mp.diff(f, p)
    f2 = mp.diff(f, p, 2)
    f3 = mp.diff(f, p, 3)
    c2 = f2 / (2 * lam * (1 - lam))
    c3 = (f3 / 6 + 2 * c2 * lam * f2 / 2) / (lam - lam ** 3)
    eps = mp.mpf(10) ** (-(mp.mp.dps / 4.0))

    def sigma(yv):
        n = 0
        u = yv - p
        z = yv
        while abs(u) > eps:
            z = f(z)
            u = z - p
            n += 1
            if n > 100000:
                raise RuntimeError('keine Konvergenz')
        return (u + c2 * u ** 2 + c3 * u ** 3) / lam ** n

    s1 = sigma(mp.mpf(1))

    def Hn(x):
        w = lam ** x * s1
        if abs(w) < eps:
            u = w
            for _ in range(8):
                u = w - c2 * u ** 2 - c3 * u ** 3
            return p + u
        lo, hi = mp.mpf(1), p
        return mp.findroot(lambda yy: sigma(yy) - w, (lo + (hi - lo) * 0.3, lo + (hi - lo) * 0.7), solver='anderson')

    return Cheb(Hn, 0, X, grad), p, lam


def lauf(f0, n, X, grad, tag=None):
    """Leiter f0, T f0, ..., T^{n-1} f0: Liste von dicts (p, lam, Werte bei 0.25/0.5/1.5) je Glied."""
    rows = []
    f = f0
    for k in range(n):
        g, p, lam = T(f, X, grad)
        res = max(abs(g(mp.mpf(j) / 8 + 1) - f(g(mp.mpf(j) / 8))) for j in range(0, 9))
        rows.append(dict(k=k, p=str(p), lam=str(lam), v025=str(f(mp.mpf('0.25'))), v05=str(f(mp.mpf('0.5'))),
                         v15=str(f(mp.mpf('1.5'))), res_next=mp.nstr(res, 3)))
        f = g
    return rows


def gecacht(name, fn):
    cf = CACHE / f'{name}.json'
    if cf.is_file():
        return json.loads(cf.read_text())
    t = time.time()
    rows = fn()
    out = dict(name=name, dps=mp.mp.dps, sec=round(time.time() - t, 1), rows=rows)
    cf.write_text(json.dumps(out, indent=0))
    print(f'  [{name}] {out["sec"]} s', flush=True)
    return out


# ---------------------------------------------------------------- Phi-Auswertung
def lnd(rows, b, start):
    return {start + r['k']: mp.log(mp.mpf(r['p']) - b) for r in rows}


def interp_inv(base, target, guess):
    rs = sorted(base)

    def L(t):
        near = sorted(rs, key=lambda s: abs(s - t))[:4]
        tot = 0
        for i in near:
            w = 1
            for j in near:
                if j != i:
                    w *= (t - j) / (i - j)
            tot += w * base[i]
        return tot
    return mp.findroot(lambda t: L(t) - target, mp.mpf(guess))


def phi_folge(base, rows, b, s0):
    """s_n - n fuer die Leiter rows, die beim (nominellen) Rang s0 beginnt. Rueckgabe {n: Phi-Schaetzer}."""
    out = {}
    rmax = max(base)
    for r in rows:
        k = r['k']
        if s0 + k + 1 > rmax - 1 or s0 + k < min(base) + 1:
            continue
        v = mp.log(mp.mpf(r['p']) - b)
        try:
            s = interp_inv(base, v, s0 + k)
            out[k] = s - k
        except Exception:
            pass
    return out
