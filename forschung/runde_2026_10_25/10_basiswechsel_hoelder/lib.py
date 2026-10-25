"""Bausteine für den Basiswechsel auf Rang 4 und 5 (reelle Achse, Kneser-Tetration aus dem Motor).

Tet(basis)          Kneser-Tetration T_b über Tschebyschow-Interpolation auf [-1,0] (Motorwerte, cache/cheb_*.json),
                    T, T', slog, slog' auf der reellen Achse; Fixpunkt L in (-2,-1), Multiplikator lam = T'(L),
                    Schröder-Reihe sigma von T an L; Pentation P (regulär an L, P(0)=1) und spog = P^{-1}.
R(Tb, Td, u)        = slog_b(T_d(u)) exakt für beliebig grosse darstellbare u (Turmreste per Logarithmen-Abschälen).
phi4(Tb, Td, th)    Rang-4-Phase Phi_{b,d}(th) = lim slog_b(T_d(n+th)) - n - th.
H(Tb, Td, y)        Rang-5-Konjugator H = lim slog_b^n o T_d^n (y), T_b o H = H o T_d.
phi5(Tb, Td, th)    Rang-5-Phase Phi5(th) = lim spog_b(P_d(n+th)) - n - th = spog_b(H(P_d(th))) - th.
"""
import json
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
CDIR = HIER / "cache"
GROSS = mp.mpf(10) ** 7          # E_j mit log E > GROSS werden nicht mehr gebildet (Beitrag < 10^-10^6)


def _tag(basis):
    return basis.replace("(", "").replace(")", "").replace("/", "_")


class Tet:
    def __init__(self, basis, N=150, dps=80, K=48):
        d = json.loads((CDIR / f"cheb_{_tag(basis)}_{N}_{dps}.json").read_text(encoding="utf-8"))
        self.name = basis
        self.b = mp.mpf(mp.mpmathify(eval(basis.replace("exp(1)", "mp.e").replace("/", "*mp.mpf(1)/"))))
        self.lb = mp.log(self.b)
        f = [mp.mpf(v) for v in d["T"]]
        self.N = N
        # Tschebyschow-Koeffizienten (Knoten x_j = -1/2 + cos(pi(j+1/2)/N)/2, t = 2x+1)
        c = []
        for k in range(N):
            s = mp.fsum(f[j] * mp.cos(mp.pi * k * (j + mp.mpf(1) / 2) / N) for j in range(N))
            c.append(2 * s / N)
        c[0] /= 2
        self.c = c
        # Ableitungskoeffizienten (nach t)
        self.dc = self._cheb_abl(c)
        self.K = K
        self._fix()

    # ---------------------------------------------------------- Tschebyschow
    @staticmethod
    def _cheb_abl(c):
        n = len(c)
        d = [mp.mpf(0)] * (n + 1)
        for k in range(n - 1, 0, -1):
            d[k - 1] = d[k + 1] + 2 * k * c[k]
        d[0] /= 2
        return d[:n]

    @staticmethod
    def _clen(c, t):
        b1 = b2 = mp.mpf(0)
        for ck in reversed(c[1:]):
            b1, b2 = 2 * t * b1 - b2 + ck, b1
        return t * b1 - b2 + c[0]

    def T0(self, x):            # x in [-1, 0]
        return self._clen(self.c, 2 * x + 1)

    def dT0(self, x):
        return 2 * self._clen(self.dc, 2 * x + 1)

    # ---------------------------------------------------------- T, slog
    def T(self, x):
        m = int(mp.ceil(x))
        t = x - m
        y = self.T0(t)
        if m >= 0:
            for _ in range(m):
                y = mp.exp(self.lb * y)
        else:
            for _ in range(-m):
                y = mp.log(y) / self.lb
        return y

    def dT(self, x):
        m = int(mp.ceil(x))
        t = x - m
        y, dy = self.T0(t), self.dT0(t)
        if m >= 0:
            for _ in range(m):
                y = mp.exp(self.lb * y)
                dy = dy * self.lb * y
        else:
            for _ in range(-m):
                dy = dy / (y * self.lb)
                y = mp.log(y) / self.lb
        return dy

    def slog(self, y, mit_abl=False):
        n = 0
        fak = mp.mpf(1)
        while y > 1:
            fak /= (y * self.lb)
            y = mp.log(y) / self.lb
            n += 1
        while y <= 0:
            fak *= self.lb * mp.exp(self.lb * y)
            y = mp.exp(self.lb * y)
            n -= 1
        # T0(t) = y, t in (-1, 0]; Startwert linear
        t = y - 1
        for _ in range(200):
            dt = (self.T0(t) - y) / self.dT0(t)
            t -= dt
            if abs(dt) < mp.mpf(10) ** (-mp.mp.dps + 3):
                break
        if mit_abl:
            return n + t, fak / self.dT0(t)
        return n + t

    # ---------------------------------------------------------- Fixpunkt, Schröder, Pentation
    def _taylor_bei(self, x0, K):
        """Taylor-Koeffizienten von T0 (Tschebyschow-Polynom) um x0 in (-1,0)."""
        a = []
        c = self.c
        t0 = 2 * x0 + 1
        fak = mp.mpf(1)
        for k in range(K + 1):
            a.append(self._clen(c, t0) * fak)
            c = self._cheb_abl(c)
            fak = fak * 2 / (k + 1)
        return a

    def _fix(self):
        K = self.K
        lo, hi = mp.mpf(-2) + mp.mpf(10) ** -12, mp.mpf(-1)      # T(lo) - lo < 0 < T(hi) - hi
        for _ in range(60):
            mid = (lo + hi) / 2
            if self.T(mid) - mid < 0:
                lo = mid
            else:
                hi = mid
        L = (lo + hi) / 2
        for _ in range(100):
            g = self.T(L) - L
            dg = self.dT(L) - 1
            L -= g / dg
            if abs(g) < mp.mpf(10) ** (-mp.mp.dps + 5):
                break
        self.L = L
        self.lam = self.dT(L)
        # g(u) = T(L+u) - L = log_b(T0(L+1+u)) - L
        a = self._taylor_bei(L + 1, K)          # T0(L+1+u) = sum a_k u^k, a_0 = T0(L+1) = b^L
        # log einer Potenzreihe: log(a0) + log(1 + q)
        q = [mp.mpf(0)] + [a[k] / a[0] for k in range(1, K + 1)]
        lg = [mp.log(a[0])] + [mp.mpf(0)] * K
        # (log(1+q))' = q'/(1+q):  l_k = (k q_k - sum_{j=1}^{k-1} j l_j q_{k-j}) / k
        for k in range(1, K + 1):
            lg[k] = (k * q[k] - mp.fsum(j * lg[j] * q[k - j] for j in range(1, k))) / k
        g = [lg[k] / self.lb for k in range(K + 1)]
        g[0] -= L
        self.g = g                                # g[0] ~ 0, g[1] = lam
        # Schröder-Reihe sigma(L+u) = sum s_k u^k, sigma(T(z)) = lam sigma(z)
        gp = [None, g[:]]
        gp[1][0] = mp.mpf(0)
        for _ in range(2, K + 1):
            gp.append(self._mul(gp[-1], gp[1], K))
        s = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (K - 1)
        lam = self.lam
        for m_ in range(2, K + 1):
            s[m_] = mp.fsum(s[k] * gp[k][m_] for k in range(1, m_)) / (lam - lam ** m_)
        self.S = s
        self.dS = [k * s[k] for k in range(1, K + 1)]
        self.sig1 = self.sigma(mp.mpf(1))

    @staticmethod
    def _mul(a, b, K):
        c = [mp.mpf(0)] * (K + 1)
        for i in range(K + 1):
            if a[i] != 0:
                for j in range(K + 1 - i):
                    c[i + j] += a[i] * b[j]
        return c

    def sigma(self, y, mit_abl=False):
        n = 0
        fak = mp.mpf(1)
        tol = mp.mpf(10) ** (-mp.mp.dps)
        rho = self.L + 2                  # Konvergenzradius der Reihen an L (log-Singularität bei -2)
        while abs(y - self.L) > rho / 10 or (abs(y - self.L) / rho) ** (self.K + 1) > tol:
            if mit_abl:
                _, dsl = self.slog(y, True)
                fak *= dsl
            y = self.slog(y)
            n += 1
            if n > 400:
                raise RuntimeError("sigma: keine Konvergenz")
        u = y - self.L
        su = mp.polyval(self.S[::-1], u)
        if mit_abl:
            return self.lam ** n * su, self.lam ** n * fak * mp.polyval(self.dS[::-1], u)
        return self.lam ** n * su

    def spog(self, y, mit_abl=False):
        """Umkehrung der Pentation, spog(T(y)) = spog(y) + 1, spog(1) = 0."""
        if mit_abl:
            s, ds = self.sigma(y, True)
            return mp.log(s / self.sig1) / mp.log(self.lam), ds / (s * mp.log(self.lam))
        return mp.log(self.sigma(y) / self.sig1) / mp.log(self.lam)

    def P(self, x):
        """Pentation b^^^x (regulär an L, P(0) = 1)."""
        w = self.sig1 * mp.power(self.lam, x)
        n = 0
        rho = self.L + 2
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while abs(w) > rho / 10 or (abs(w) / rho) ** (self.K + 1) > tol:
            w /= self.lam
            n += 1
        # sigma(L+u) = w nach u auflösen
        u = w
        for _ in range(100):
            du = (mp.polyval(self.S[::-1], u) - w) / mp.polyval(self.dS[::-1], u)
            u -= du
            if abs(du) < mp.mpf(10) ** (-mp.mp.dps + 3):
                break
        y = self.L + u
        for _ in range(n):
            y = self.T(y)
        return y


# ------------------------------------------------------------------ Rang 4
def R(Tb, Td, u):
    """slog_b(T_d(u)) für reelles u (auch sehr gross), exakt bis auf Terme < 10^-10^6."""
    k = int(mp.floor(u))
    th = u - k
    j0 = 2                               # Abschälen endet bei E_{j0} (vermeidet log negativer Zahlen)
    if k <= 0:
        return Tb.slog(Td.T(u))
    v = Td.T(th)
    alpha = Td.lb / Tb.lb
    la = mp.log(alpha) / Tb.lb
    # E_j = exp_d^j(v), j = 1 .. k-1, solange darstellbar
    E = [v]
    while len(E) <= k - 1 and E[-1] * Td.lb < GROSS * mp.log(10):
        E.append(mp.exp(Td.lb * E[-1]))
    jmax = len(E) - 1
    j0 = min(j0, jmax, k - 1)
    a = mp.mpf(0) if jmax >= k - 1 else la
    for j in range(min(k - 1, jmax), j0, -1):
        a = la + mp.log1p(a / (alpha * E[j])) / Tb.lb
    # jetzt a + alpha*E_{j0} = log_b^{k-j0}(T_d(u))
    return (k - j0) + Tb.slog(a + alpha * E[j0])


def phi4(Tb, Td, th, N=12):
    return R(Tb, Td, N + th) - N - th


# ------------------------------------------------------------------ Rang 5
def H(Tb, Td, y, X=None):
    """H(y) = lim slog_b^n(T_d^n(y)); Abbruch, sobald T_d(U) > 10^X (Rest ~ 10^-X)."""
    if X is None:
        X = mp.mp.dps // 2
    grenze = mp.mpf(10) ** X
    A = Td.slog(grenze)                  # T_d(A) = 10^X
    U = [y]
    while U[-1] < A:
        U.append(Td.T(U[-1]))
        if len(U) > 60:
            raise RuntimeError("H: Kette zu lang")
    w = R(Tb, Td, U[-1])                 # = slog_b(T_d(U_last))  ~  slog_b(U_last+1 + c)
    for _ in range(len(U) - 1):
        w = Tb.slog(w)
    return w


def phi5(Tb, Td, th, X=None):
    y = Td.P(th)
    return Tb.spog(H(Tb, Td, y, X)) - th


def G5_endlich(Tb, Td, th, n):
    """spog_b(P_d(n+th)) - n - th direkt (nur für kleine n, solange P_d(n+th) darstellbar)."""
    return Tb.spog(Td.P(n + th)) - n - th
