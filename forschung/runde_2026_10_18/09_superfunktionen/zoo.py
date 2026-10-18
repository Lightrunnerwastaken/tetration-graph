"""Zoo der Kneser-Superfunktionen: allgemeiner Loeser fuer S(z+1) = f(S(z)).

Verfahren (Kouznetsov-Cauchy-Kontur, Gauss-Newton):
  S ist im Streifen |Re z - x0| < 1 holomorph, reell auf R, S -> L (oben), Lquer (unten).
  Cauchy ueber den unendlichen Streifen:
     S(z) = 1/(2 pi) Int_R [ f(S(x0+it)) / (1+it-zeta) - g(S(x0+it)) / (-1+it-zeta) ] dt,  zeta = z - x0,
  mit S(x0+1+it) = f(S(x0+it)), S(x0-1+it) = g(S(x0+it)) (g = Hauptast der Umkehrung, g(H) in H).
  Unbekannte: u_j = S(x0 + i t_j), t_j > 0 (Gauss-Legendre-Paneele); t < 0 ueber Spiegelung.
  Normierung: S(x0) = s0.

Nur numpy/scipy, doppelte Genauigkeit.  Zwischenspeicher: ./daten/*.npz (eigener Ordner).
"""
import numpy as np, os, json, math, cmath

HIER = os.path.dirname(os.path.abspath(__file__))
DATEN = os.path.join(HIER, 'daten')
os.makedirs(DATEN, exist_ok=True)


# ---------------------------------------------------------------- Abbildungen
class Abb:
    """f, f', g (Hauptast der Umkehrung auf H), Taylor-Koeffizienten von f bei L."""
    def __init__(self, name, f, fp, g, taylor, start=0.5 + 1.0j):
        self.name, self.f, self.fp, self.g, self.taylor = name, f, fp, g, taylor
        w = complex(start)
        for _ in range(400):
            w = complex(g(np.array([w]))[0])
        for _ in range(30):  # Newton auf f(w) = w
            fw = complex(f(np.array([w]))[0]); d = complex(fp(np.array([w]))[0])
            w = w - (fw - w) / (d - 1)
        self.L = w
        self.lam = complex(fp(np.array([w]))[0])
        self.loglam = cmath.log(self.lam)
        assert abs(self.lam) > 1 and self.loglam.imag > 0, (self.lam, self.loglam)
        self._schroeder(24)

    def gp(self, w):
        return 1.0 / self.fp(self.g(w))

    def _schroeder(self, K):
        """V(u) = Psi(u) - L = u + sum p_k u^k,  V(lam u) = f(L+V(u)) - L."""
        lam = self.lam
        fj = self.taylor(self.L, K)  # fj[j] = f^(j)(L)/j!,  j = 0..K
        p = np.zeros(K + 1, complex); p[1] = 1.0
        for k in range(2, K + 1):
            # [u^k] sum_{j>=2} fj V^j  mit V bis Ordnung k-1 bekannt
            V = p.copy(); V[k:] = 0
            acc = 0j; Vj = V.copy()
            for j in range(2, k + 1):
                Vj = np.convolve(Vj, V)[:K + 1]
                acc += fj[j] * Vj[k]
            p[k] = acc / (lam ** k - lam)
        self.p = p

    def V(self, s):
        r = np.zeros_like(s);
        for k in range(len(self.p) - 1, 0, -1):
            r = (r + self.p[k]) * s
        return r

    def Vp(self, s):
        r = np.zeros_like(s)
        for k in range(len(self.p) - 1, 0, -1):
            r = r * s + k * self.p[k]
        return r

    def _iter(self, w, rad=0.02, nmax=400):
        """g-Iteration je Punkt, bis |w - L| < rad. Liefert w_n, n (Feld), log Ableitung."""
        w = np.asarray(w, complex).copy(); n = np.zeros(w.shape, int); dl = np.zeros(w.shape, complex)
        akt = np.abs(w - self.L) >= rad
        for _ in range(nmax):
            if not akt.any(): break
            wa = w[akt]; gw = self.g(wa)
            dl[akt] -= np.log(self.fp(gw)); w[akt] = gw; n[akt] += 1
            akt = np.abs(w - self.L) >= rad
        v = w - self.L; s = v.copy()
        for _ in range(8):
            s = s - (self.V(s) - v) / self.Vp(s)
        return s, n, dl

    def logsigma(self, w, n=None):
        """log sigma(w) = n log lam + log s,  V(s) = g^n(w) - L  (Hauptwert log s), n je Punkt."""
        s, n, dl = self._iter(w)
        return n * self.loglam + np.log(s), n

    def h(self, w, n=None):
        ls, n = self.logsigma(w, n)
        return ls / self.loglam

    def hp(self, w, n=None):
        """h'(w) = sigma'/(log lam * sigma)."""
        s, n, dl = self._iter(w)
        return np.exp(dl) / self.Vp(s) / s / self.loglam


def taylor_poly(coef):
    """coef: Polynomkoeffizienten von f (aufsteigend)."""
    def t(L, K):
        c = np.array(coef, complex); out = np.zeros(K + 1, complex)
        # Taylor bei L: f(L+v)
        d = c.copy()
        for j in range(K + 1):
            if len(d) == 0: break
            out[j] = np.polyval(d[::-1], L) / math.factorial(j)
            d = np.array([k * d[k] for k in range(1, len(d))], complex)
        return out
    return t


def abb_quadrat(c, d=2):
    """f(x) = x^d + c (d gerade)."""
    coef = [c] + [0] * (d - 1) + [1]
    f = lambda w: w ** d + c
    fp = lambda w: d * w ** (d - 1)
    g = lambda w: (w - c + 0j) ** (1.0 / d)
    a = Abb(f'x^{d}+{c}', f, fp, g, taylor_poly(coef), start=0.5 + 1.0j)
    a.krit = dict(wc=0.0, v=float(c), d=d, a=1.0)
    return a


def abb_exp(lnb=1.0):
    """f(x) = b^x = exp(lnb x)."""
    f = lambda w: np.exp(lnb * w)
    fp = lambda w: lnb * np.exp(lnb * w)
    g = lambda w: np.log(w + 0j) / lnb
    def t(L, K):
        return np.array([cmath.exp(lnb * L) * lnb ** j / math.factorial(j) for j in range(K + 1)])
    a = Abb(f'exp({lnb}x)', f, fp, g, t, start=0.3 + 1.3j)
    a.krit = dict(wc=None, v=0.0, d=math.inf, a=None)
    a.lnb = lnb
    return a


def abb_cosh(c=0.0):
    """f(x) = cosh(x) + c."""
    f = lambda w: np.cosh(w) + c
    fp = lambda w: np.sinh(w)
    g = lambda w: np.arccosh(w - c + 0j)
    def t(L, K):
        return np.array([(cmath.cosh(L) + (c if j == 0 else 0) if j % 2 == 0 else cmath.sinh(L)) / math.factorial(j)
                         for j in range(K + 1)])
    a = Abb(f'cosh(x)+{c}', f, fp, g, t, start=0.5 + 1.0j)
    a.krit = dict(wc=0.0, v=1.0 + c, d=2, a=0.5)
    return a


# ---------------------------------------------------------------- Loeser
def gl_paneele(T, P, n):
    x, w = np.polynomial.legendre.leggauss(n)
    ts, ws = [], []
    a = 0.0
    while a < T - 1e-12:
        b = min(a + P, T)
        ts.append(a + (b - a) * (x + 1) / 2); ws.append(w * (b - a) / 2)
        a = b
    return np.concatenate(ts), np.concatenate(ws)


class Loesung:
    def __init__(self, abb, x0, s0, T=40.0, P=1.0, n=16, u=None, t=None):
        self.abb, self.x0, self.s0 = abb, x0, s0
        self.t, self.w = gl_paneele(T, P, n)
        self.T = T
        if u is None:
            L = abb.L
            u = L + (s0 - L) * np.exp(-0.7 * self.t)
        elif t is not None:  # Interpolation von altem Gitter
            u = np.interp(self.t, t, u.real) + 1j * np.interp(self.t, t, u.imag)
        self.u = np.asarray(u, complex)

    def _kern(self, zeta):
        t = self.t[None, :]; z = zeta[:, None]
        return (1 / (1 + 1j * t - z), 1 / (-1 + 1j * t - z), 1 / (1 - 1j * t - z), 1 / (-1 - 1j * t - z))

    def _tail(self, zeta):
        L, T = self.abb.L, self.T
        return (L * 1j * np.log((1 + 1j * T - zeta) / (-1 + 1j * T - zeta))
                - np.conj(L) * 1j * np.log((1 - 1j * T - zeta) / (-1 - 1j * T - zeta))) / (2 * np.pi)

    def S(self, z, u=None):
        """S(z) im Streifen |Re z - x0| < 1 (genau fuer |Re z - x0| <~ 0.7)."""
        u = self.u if u is None else u
        z = np.atleast_1d(np.asarray(z, complex))
        out = np.empty(z.shape, complex)
        flat = z.ravel(); res = np.empty(flat.shape, complex)
        R = self.abb.f(u); G = self.abb.g(u)
        for i0 in range(0, len(flat), 2000):
            zeta = flat[i0:i0 + 2000] - self.x0
            k1, k2, k3, k4 = self._kern(zeta)
            ww = self.w[None, :]
            res[i0:i0 + 2000] = ((ww * (k1 * R - k2 * G + k3 * np.conj(R) - k4 * np.conj(G))).sum(1) / (2 * np.pi)
                                 + self._tail(zeta))
        return res.reshape(z.shape)

    def residuum(self, u):
        z = np.concatenate([[self.x0 + 0j], self.x0 + 1j * self.t])
        C = self.S(z, u)
        return self.s0 - C[0], u - C[1:]

    def loese(self, maxit=40, tol=1e-15, info=True):
        abb = self.abb; N = len(self.t)
        zeta = np.concatenate([[0j], 1j * self.t])
        k1, k2, k3, k4 = self._kern(zeta)
        ww = self.w[None, :] / (2 * np.pi)
        alt = None
        for it in range(maxit):
            u = self.u
            r0, r = self.residuum(u)
            nr = math.sqrt(abs(r0) ** 2 + np.sum(abs(r) ** 2))
            if info:
                print(f'  it {it:2d}  |res| = {nr:.3e}', flush=True)
            if nr < tol:
                break
            fpu = abb.fp(u); gpu = abb.gp(u)
            A = ww * (k1 * fpu - k2 * gpu)
            B = ww * (k3 * np.conj(fpu) - k4 * np.conj(gpu))
            I = np.vstack([np.zeros((1, N)), np.eye(N)])
            Ja = I - A - B          # d r / d a   (komplex)
            Jb = 1j * (I - A + B)   # d r / d b
            J = np.vstack([np.hstack([Ja.real, Jb.real]), np.hstack([Ja.imag, Jb.imag])])
            rr = np.concatenate([[r0], r])
            rvec = np.concatenate([rr.real, rr.imag])
            # Zeile 0 (Normierung) staerker gewichten
            W = np.ones(len(rvec)); W[0] = 10.0
            dx = np.linalg.lstsq(J * W[:, None], -rvec * W, rcond=None)[0]
            du = dx[:N] + 1j * dx[N:]
            lam = 1.0
            for _ in range(20):
                r0n, rn = self.residuum(u + lam * du)
                nn = math.sqrt(abs(r0n) ** 2 + np.sum(abs(rn) ** 2))
                if nn < nr or nn < 1e-14:
                    break
                lam /= 2
            self.u = u + lam * du
            if alt is not None and nr > 0.5 * alt and nr < 1e-12:
                break
            alt = nr
        self.res = nr
        return nr

    def speichern(self, name):
        np.savez(os.path.join(DATEN, name + '.npz'), t=self.t, w=self.w, u=self.u, x0=self.x0, s0=self.s0, T=self.T)

    @staticmethod
    def laden(abb, name):
        d = np.load(os.path.join(DATEN, name + '.npz'))
        s = Loesung.__new__(Loesung)
        s.abb = abb; s.t = d['t']; s.w = d['w']; s.u = d['u']; s.x0 = float(d['x0']); s.s0 = float(d['s0']); s.T = float(d['T'])
        return s

    # ---- Fortsetzung ueber die Funktionalgleichung
    def Sweit(self, z):
        """S(z) fuer beliebiges Re z (Verschiebung in den Kernstreifen |Re-x0| <= 0.5)."""
        z = np.atleast_1d(np.asarray(z, complex))
        k = np.round(z.real - self.x0).astype(int)
        out = np.empty(z.shape, complex)
        for kk in np.unique(k):
            m = k == kk
            v = self.S(z[m] - kk)
            if kk > 0:
                for _ in range(kk): v = self.abb.f(v)
            else:
                for _ in range(-kk): v = self.abb.g(v)
            out[m] = v
        return out
