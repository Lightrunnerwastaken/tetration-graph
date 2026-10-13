"""Etappe 6: Kneser-artige Superfunktion zu einer beliebigen unteren Funktion f (Cauchy-Verfahren).

Gesucht S mit S(x+1)=f(S(x)), S(0)=1, holomorph im Streifen, S(x+iy)->L (y->+inf), ->conj(L) (y->-inf),
L oberer Fixpunkt von f (f reell auf der reellen Achse).
Darstellung: F(w)=S(w-c) auf dem Rechteck [-1,1]x[-A,A] durch die Cauchy-Formel aus den Werten h(y)=F(iy);
Randwerte F(1+iy)=f(h), F(-1+iy)=f^{-1}(h).

Untere Funktionen implementieren: f(w), finv(w, u0) (vektorisiert, numpy complex).
"""
import numpy as np

LN2 = np.log(2.0)


class Exp2:
    name = 'H3'
    def f(self, w):
        return np.exp(LN2 * w)
    def finv(self, w, u0=None):
        return np.log(w) / LN2


class Superfunktion:
    def __init__(self, unten, L, A=30.0, N=600, kappa=3.0, name='S'):
        self.unten = unten; self.L = complex(L); self.A = A; self.name = name
        t, w = np.polynomial.legendre.leggauss(N)
        self.y = A * np.sinh(kappa * t) / np.sinh(kappa)
        self.w = w * A * kappa * np.cosh(kappa * t) / np.sinh(kappa)
        self.c = 0.0
        self.rmax = 0.75
        self.reparatur_alle = 0

    @classmethod
    def laden(cls, unten, d, name='S'):
        S = cls(unten, d['L'], A=d['A'], N=d['N'], name=name)
        S.h = d['h']; S.fr = d['fr']; S.u = d['u']; S.c = d['c']; S.res = d.get('res', 0)
        return S

    def speichern(self):
        return {'h': self.h, 'fr': self.fr, 'u': self.u, 'c': self.c, 'L': self.L, 'A': self.A, 'N': len(self.y), 'res': self.res}

    # ---------- Cauchy-Darstellung ----------
    def _rand(self, z):
        L, Lb, A = self.L, np.conj(self.L), self.A
        return (L * (np.log(-1 + 1j * A - z) - np.log(1 + 1j * A - z)) +
                Lb * (np.log(1 - 1j * A - z) - np.log(-1 - 1j * A - z))) / (2j * np.pi)

    def _F(self, z, h=None, fr=None, u=None):
        h = self.h if h is None else h
        fr = self.fr if fr is None else fr
        u = self.u if u is None else u
        z = np.atleast_1d(np.asarray(z, dtype=complex))
        Y = self.y[None, :]
        I = (fr[None, :] / (1 + 1j * Y - z[:, None]) - u[None, :] / (-1 + 1j * Y - z[:, None])) @ self.w / (2 * np.pi)
        return I + self._rand(z)

    def _finv_fortsetzung(self, h, start):
        # von grossem |y| (Naehe des Fixpunkts L der unteren Funktion) nach innen: Zweig der lokalen
        # Inversen am Fixpunkt (bei nicht injektiven unteren Funktionen wichtig)
        order = np.argsort(-np.abs(self.y))
        u = np.empty_like(h)
        prev = {1: self.L, -1: np.conj(self.L)}
        for i in order:
            s = 1 if self.y[i] >= 0 else -1
            u[i] = self.unten.finv(np.array([h[i]]), np.array([prev[s]]))[0]
            prev[s] = u[i]
        return u

    def _reparieren(self, h, u, faktor=8.0):
        """Zweigspruenge der linken Randwerte beseitigen: von grossem |y| nach innen; springt u staerker
        als faktor*(Sprung von h) + 0.05, wird u aus dem aeusseren Nachbarn neu per Newton bestimmt."""
        o_pos = [i for i in np.argsort(-self.y) if self.y[i] >= 0]
        o_neg = [i for i in np.argsort(self.y) if self.y[i] < 0]
        for ordn in (o_pos, o_neg):
            for a, b in zip(ordn, ordn[1:]):
                if abs(u[b] - u[a]) > faktor * abs(h[b] - h[a]) + 0.05:
                    u[b] = self.unten.finv(np.array([h[b]]), np.array([u[a]]))[0]
                    self.reparaturen = getattr(self, 'reparaturen', 0) + 1
        return u

    def loese(self, h0=None, iters=3000, lam=0.3, tol=1e-14, verbose=False):
        L = self.L
        h = h0 if h0 is not None else (np.real(L) + 1j * np.imag(L) * np.tanh(self.y / 2) + (1 - np.real(L)) / np.cosh(self.y))
        u = self._finv_fortsetzung(h, 0.0 + 0j)
        z = 1j * self.y
        d = np.inf
        for k in range(iters):
            fr = self.unten.f(h)
            u = self.unten.finv(h, u)
            if self.reparatur_alle and k % self.reparatur_alle == 0:
                u = self._reparieren(h, u)
            hn = self._F(z, h, fr, u)
            d = np.max(np.abs(hn - h))
            h = (1 - lam) * h + lam * hn
            if verbose and k % 200 == 0:
                print('   iter', k, '%.2e' % d, flush=True)
            if d < tol:
                break
        self.h = h; self.fr = self.unten.f(h); self.u = self.unten.finv(h, u)
        self.res = d; self.iter = k
        self._normieren()
        return d

    def _normieren(self):
        c = 0.0
        for _ in range(60):
            v = self._F(c)[0] - 1
            dv = (self._F(c + 1e-6)[0] - self._F(c - 1e-6)[0]) / 2e-6
            step = np.real(v / dv)
            c -= step
            if abs(step) < 1e-15:
                break
        self.c = c

    # ---------- Auswertung ueberall ----------
    def __call__(self, z):
        """S(z) fuer komplexe z (vektorisiert); Verschiebung in den Streifen |Re(z+c)|<=0.5."""
        z = np.atleast_1d(np.asarray(z, dtype=complex))
        re = np.real(z) + self.c
        k = np.where(re > self.rmax, np.ceil(re - self.rmax), np.where(re < -self.rmax, -np.ceil(-re - self.rmax), 0)).astype(int)
        w = z + self.c - k
        out = self._F(w)
        for kk in range(int(k.max(initial=0)), 0, -1):
            m = k >= kk
            out[m] = self.unten.f(out[m])
        # links: Inverse der unteren Funktion
        for kk in range(int(k.min(initial=0)), 0):
            m = k <= kk
            out[m] = self.unten.finv(out[m], out[m] - 1)
        # grosse |Im|: Rechteckgrenze beachten
        gross = np.abs(np.imag(z)) > 0.9 * self.A
        out[gross] = np.where(np.imag(z[gross]) > 0, self.L, np.conj(self.L))
        return out

    def ableitung(self, z, h=1e-6):
        return (self(z + h) - self(z - h)) / (2 * h)

    # ---------- als untere Funktion der naechsten Stufe ----------
    def f(self, w):
        return self(w)

    def finv(self, w, u0=None):
        w = np.atleast_1d(np.asarray(w, dtype=complex))
        u = np.array(w - 1 if u0 is None else u0, dtype=complex)
        for _ in range(100):
            r = self(u) - w
            d = self.ableitung(u)
            step = np.where(np.abs(d) > 1e-12, r / np.where(np.abs(d) > 1e-12, d, 1), 0.1 * r)
            big = np.abs(step) > 0.3
            step[big] = 0.3 * step[big] / np.abs(step[big])
            u = u - step
            if np.max(np.abs(r)) < 1e-14:
                break
        return u


def fixpunkt(funk, start, it=100):
    z = complex(start)
    for _ in range(it):
        v = funk(np.array([z]))[0] - z
        d = (funk(np.array([z + 1e-6]))[0] - funk(np.array([z - 1e-6]))[0]) / 2e-6 - 1
        z -= v / d
        if abs(v) < 1e-14:
            break
    return z
