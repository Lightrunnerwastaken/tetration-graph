"""Etappe 4: Kneser-artige Superfunktion fuer f_eps(w) = w+1+eps(2^w-w-1) (Cauchy-Integral-Verfahren
nach Kouznetsovs Ansatz, eigene Umsetzung).

Gesucht F holomorph im Streifen, F(z+1)=f(F(z)), F(x+iy) -> L (oberer Fixpunkt) fuer y->+inf,
-> conj(L) fuer y->-inf, normiert F(0)=1.  Unbekannte: h(y)=F(iy) auf Gauss-Legendre-Knoten in [-A,A].
Cauchy-Formel auf dem Rechteck [-1,1]x[-A,A]:
  F(z) = 1/(2pi) int [ f(h(y))/(1+iy-z) - finv(h(y))/(-1+iy-z) ] dy
         + 1/(2pi i) [ L(log(-1+iA-z)-log(1+iA-z)) + Lb(log(1-iA-z)-log(-1-iA-z)) ].
Doppelte Genauigkeit (numpy).  Bei eps=1 Kontrolle gegen GP-Kneser.
"""
import sys, json, time, argparse
from pathlib import Path
import numpy as np
HIER = Path(__file__).resolve().parent
LN2 = np.log(2.0)


def fixpunkt(eps, start=0.8247 + 1.5674j):
    L = start
    for _ in range(200):
        g = L + 1 + eps * (2 ** L - L - 1) - L
        dg = eps * (2 ** L * LN2 - 1)
        L -= g / dg
    return L


def f(w, eps):
    return w + 1 + eps * (np.exp(LN2 * w) - w - 1)


def df(w, eps):
    return 1 + eps * (np.exp(LN2 * w) * LN2 - 1)


def finv(w, eps, u0):
    u = u0.copy()
    for _ in range(60):
        r = f(u, eps) - w
        u = u - r / df(u, eps)
        if np.max(np.abs(r)) < 1e-15:
            break
    return u


class Loeser:
    def __init__(self, eps, A=30.0, N=600, L=None):
        self.eps = eps; self.A = A
        self.L = L if L is not None else fixpunkt(eps)
        t, w = np.polynomial.legendre.leggauss(N)
        # Knoten verdichten bei kleinem |y|: y = A*sinh(c t)/sinh(c)
        c = 3.0
        self.y = A * np.sinh(c * t) / np.sinh(c)
        self.w = w * A * c * np.cosh(c * t) / np.sinh(c)

    def rand(self, z):
        L, Lb, A = self.L, np.conj(self.L), self.A
        return (L * (np.log(-1 + 1j * A - z) - np.log(1 + 1j * A - z)) +
                Lb * (np.log(1 - 1j * A - z) - np.log(-1 - 1j * A - z))) / (2j * np.pi)

    def auswerten(self, h, finvh, z):
        z = np.atleast_1d(np.asarray(z, dtype=complex))
        fr = f(h, self.eps)
        Y = self.y[None, :]
        I = (fr[None, :] / (1 + 1j * Y - z[:, None]) - finvh[None, :] / (-1 + 1j * Y - z[:, None])) @ self.w / (2 * np.pi)
        return I + self.rand(z)

    def loese(self, iters=4000, lam=0.3, tol=1e-15, verbose=False):
        L = self.L
        h = np.real(L) + 1j * np.imag(L) * np.tanh(self.y / 2) + (1 - np.real(L)) / np.cosh(self.y)
        u = h - 1.0
        # Fortsetzung fuer finv: von y=0 nach aussen
        u = self._finv_fortsetzung(h)
        z = 1j * self.y
        for k in range(iters):
            u = finv(h, self.eps, u)
            hn = self.auswerten(h, u, z)
            d = np.max(np.abs(hn - h))
            h = (1 - lam) * h + lam * hn
            if verbose and k % 200 == 0:
                print(k, d, flush=True)
            if d < tol:
                break
        self.h = h; self.u = finv(h, self.eps, u); self.iter = k; self.res = d
        return d

    def _finv_fortsetzung(self, h):
        order = np.argsort(np.abs(self.y))
        u = np.empty_like(h)
        prev_pos = h[order[0]] - 1; prev_neg = prev_pos
        for i in order:
            g = prev_pos if self.y[i] >= 0 else prev_neg
            ui = finv(np.array([h[i]]), self.eps, np.array([g]))[0]
            u[i] = ui
            if self.y[i] >= 0: prev_pos = ui
            else: prev_neg = ui
        return u

    def F(self, z):
        return self.auswerten(self.h, self.u, z)

    def normieren(self):
        c = 0.0
        for _ in range(50):
            v = self.F(c)[0] - 1
            dv = (self.F(c + 1e-5)[0] - self.F(c - 1e-5)[0]) / 2e-5
            c = c - np.real(v / dv)
            if abs(v) < 1e-15:
                break
        self.c = c
        return c

    def sexp(self, x):
        """Normierte Superfunktion fuer x in (-1-c, 1-c); groessere x ueber f."""
        x = np.atleast_1d(np.asarray(x, dtype=complex))
        out = []
        for xi in x:
            k = 0
            while np.real(xi) + self.c > 0.75:
                xi -= 1; k += 1
            v = self.F(xi + self.c)[0]
            for _ in range(k):
                v = f(v, self.eps)
            out.append(v)
        return np.array(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--eps', type=float, default=1.0)
    ap.add_argument('--A', type=float, default=30.0); ap.add_argument('--N', type=int, default=600)
    a = ap.parse_args()
    t0 = time.time()
    S = Loeser(a.eps, a.A, a.N)
    d = S.loese(verbose=True)
    c = S.normieren()
    print('eps', a.eps, 'L', S.L, 'Residuum', d, 'Iter', S.iter, 'c', c, round(time.time() - t0, 1), 's')
    if a.eps == 1.0:
        kn = json.loads((HIER / 'data' / 'e4_kneser_raster.json').read_text())
        err = 0
        for w in kn['werte']:
            x = complex(float(w['x'][0]), float(w['x'][1])); ref = complex(float(w['sexp'][0]), float(w['sexp'][1]))
            if abs(x.imag) > 0 or x.real > 1.5: continue
            err = max(err, abs(S.sexp(x)[0] - ref))
        print('max |F - GP-Kneser| auf [-0.5,1.5]:', err)
