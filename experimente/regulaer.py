"""Reguläre Iteration von exp_b für 1 < b < η über Koenigs-Reihen (mpmath).

Unterer Fixpunkt a (anziehend für f = b^x), oberer Fixpunkt A (anziehend für h = log_b).
alpha_unten: Abel-Funktion aus der Koenigs-Funktion an a, alpha_oben: aus der an A.
Beide erfüllen alpha(b^x) = alpha(x) + 1.
"""
import mpmath as mp


def fixpunkte(b):
    lb = mp.log(b)
    a = -mp.lambertw(-lb, 0) / lb
    A = -mp.lambertw(-lb, -1) / lb
    return mp.re(a), mp.re(A)


def schroeder_reihe(g, lam, N):
    """Koeffizienten s[0..N] von sigma mit sigma(g(w)) = lam*sigma(w), s[1] = 1. g[0] = 0, g[1] = lam."""
    potenzen = [None, list(g[: N + 1])]
    for j in range(2, N + 1):
        p, q = potenzen[-1], g
        r = [mp.mpf(0)] * (N + 1)
        for i in range(j - 1, N + 1):
            if p[i] == 0:
                continue
            for k in range(1, N + 1 - i):
                r[i + k] += p[i] * q[k]
        potenzen.append(r)
    s = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (N - 1)
    for n in range(2, N + 1):
        acc = mp.fsum(s[j] * potenzen[j][n] for j in range(1, n))
        s[n] = acc / (lam - lam ** n)
    return s


def horner(c, w):
    y = mp.mpf(0)
    for k in range(len(c) - 1, 0, -1):
        y = (y + c[k]) * w
    return y


class Regulaer:
    def __init__(self, b, N=90, tol=mp.mpf("2e-3")):
        self.b = mp.mpf(b)
        self.lb = mp.log(self.b)
        self.a, self.A = fixpunkte(self.b)
        self.lam_a = self.a * self.lb
        self.lam_A = self.A * self.lb
        self.N, self.tol = N, tol
        g = [mp.mpf(0)] + [self.a * self.lb ** k / mp.factorial(k) for k in range(1, N + 1)]
        self.s_a = schroeder_reihe(g, self.lam_a, N)
        h = [mp.mpf(0)] + [(-1) ** (k + 1) / (k * self.lb * self.A ** k) for k in range(1, N + 1)]
        self.s_A = schroeder_reihe(h, 1 / self.lam_A, N)

    def f(self, x):
        return mp.power(self.b, x)

    def h(self, x):
        return mp.log(x) / self.lb

    def alpha_unten(self, x):
        """Abel-Funktion vom Fixpunkt a; reell auf (-inf, a) und (a, A) (verschiedene Zweige)."""
        n = 0
        while abs(x - self.a) > self.tol:
            x = self.f(x)
            n += 1
            if n > 10000:
                raise RuntimeError("keine Konvergenz gegen a")
        s = horner(self.s_a, x - self.a)
        return mp.log(abs(s)) / mp.log(self.lam_a) - n

    def alpha_oben(self, x):
        """Abel-Funktion vom Fixpunkt A über Rückwärtsiteration log_b, für x in (a, A)."""
        n = 0
        while abs(x - self.A) > self.tol:
            x = self.h(x)
            n += 1
            if n > 10000:
                raise RuntimeError("keine Konvergenz gegen A")
        s = horner(self.s_A, x - self.A)
        return n + mp.log(abs(s)) / mp.log(self.lam_A)

    # --- Tetration (Zweig unter a), normiert auf T(0) = 1
    def slog(self, x):
        return self.alpha_unten(x) - self.alpha_unten(mp.mpf(1))

    def tet(self, t, x_lo=None, x_hi=None):
        """T(t) für reelles t > -2 durch Umkehren von slog (Bracketing auf (−∞, a))."""
        t = mp.mpf(t)
        k = int(mp.floor(t))
        frac = t - k
        # T(frac) liegt in [1, b) für frac in [0,1)
        x = mp.findroot(lambda y: self.slog(y) - frac, (mp.mpf(1), self.b), solver="illinois",
                       tol=mp.mpf(10) ** (-mp.mp.dps + 8))
        for _ in range(k):
            x = self.f(x)
        for _ in range(-k):
            x = self.h(x)
        return x
