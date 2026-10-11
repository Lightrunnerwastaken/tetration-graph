"""Hexation: Umkehrung der Pentation, Fixpunkte der Pentation, reguläre Hexation am anziehenden reellen Fixpunkt (b < η₅).

Baut auf stufen.pentation auf (Pentation regulär am negativen Fixpunkt L_T der Kneser-Tetration T_b).
"""
import mpmath as mp

from gemeinsam import kneser
from stufen import basis_text, pentation, poincare, psi_umkehr, psi_wert


def negativer_fixpunkt(b):
    """Fixpunkt von T_b in (−2, −1): Vorzeichenwechsel auf einem Raster, dann Sekanten."""
    xs = [mp.mpf("-1.97") + mp.mpf(k) / 50 for k in range(0, 36)]
    d = [mp.re(v) - x for v, x in zip(kneser("sexp", basis_text(b), xs), xs)]
    i = next(i for i in range(len(xs) - 1) if d[i] < 0 <= d[i + 1])
    x0, x1, f0, f1 = xs[i], xs[i + 1], d[i], d[i + 1]
    for _ in range(40):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
        f1 = mp.re(kneser("sexp", basis_text(b), [x1])[0]) - x1
        if abs(x1 - x0) < mp.mpf(10) ** -25:
            break
    return x1


class Pentation:
    """Pentation P_b mit Umkehrung P_b^{-1} (Hauptzweig: reell für reelle w > L_T)."""

    def __init__(self, b, M=8):
        self.b = b
        L0 = negativer_fixpunkt(b)
        self.L, self.lam, self.psiT, self.c, self.kette, self.P = pentation(b, L0, mp.mpf("0.08"), M=M)
        self.lnlam = mp.log(self.lam)

    def inv(self, ws, m=7):
        zs = list(ws)
        for _ in range(m):
            zs = kneser("slog", basis_text(self.b), zs)
        return [m + mp.log(psi_umkehr(self.L, self.psiT, z) / self.c) / self.lnlam for z in zs]

    def taylor(self, x0, r, N=64, K=30):
        zs = [x0 + r * mp.expjpi(2 * mp.mpf(j) / N) for j in range(N)]
        w = self.P(zs)
        return [mp.fsum(w[j] * mp.expjpi(-2 * m * mp.mpf(j) / N) for j in range(N)) / N / r ** m for m in range(K + 1)]

    def fixpunkt(self, z0, r=mp.mpf("0.15"), N=64, K=30):
        """Fixpunkt von P_b nahe z0 (reell oder komplex) und Taylorreihe dort (a[1] = Multiplikator)."""
        z = mp.mpc(z0)
        for _ in range(30):
            a = self.taylor(z, r, N, K)
            dz = (a[0] - z) / (a[1] - 1)
            z -= dz
            if abs(dz) < mp.mpf(10) ** -28:
                break
        return z, self.taylor(z, r, N, K)


class HexationRegulaer:
    """Reguläre Hexation H(x) = Ψ(μ^x σ(1)) am unteren, anziehenden reellen Fixpunkt p₁ der Pentation (nur für b < η₅)."""

    def __init__(self, b, KP=26):
        self.Pe = Pentation(b)
        P = self.Pe.P
        xs = [mp.mpf(1) + mp.mpf(k) / 20 for k in range(0, 41)]      # erster Vorzeichenwechsel von P(x) − x über 1 ist p₁
        d = [mp.re(v) - x for v, x in zip(P(xs), xs)]
        i = next(i for i in range(len(xs) - 1) if d[i] > 0 >= d[i + 1])
        p, a = self.Pe.fixpunkt((xs[i] + xs[i + 1]) / 2)
        self.p1, self.a = mp.re(p), [mp.re(c) for c in a]
        self.mu = self.a[1]
        self.psi = poincare(self.a, K=KP)
        self.rlok = min((mp.mpf(10) ** -32 / abs(self.psi[k])) ** (mp.mpf(1) / (k + 1)) for k in range(KP - 4, KP) if self.psi[k] != 0)
        y, n = mp.mpf(1), 0
        while abs(y - self.p1) > self.rlok / 2:
            y = mp.re(P([y])[0])
            n += 1
        self.sigma1 = psi_umkehr(self.p1, self.psi, y) / self.mu ** n

    def H(self, xs):
        """Werte H(x) für reelle x > s_b (Stapel); Ψ lokal, dann k-mal P^{-1}."""
        ks = []
        for x in xs:
            k = 0
            while abs(self.mu ** (x + k) * self.sigma1) > self.rlok:
                k += 1
            ks.append(k)
        ys = [psi_wert(self.p1, self.psi, self.mu ** (x + k) * self.sigma1) for x, k in zip(xs, ks)]
        for j in range(max(ks), 0, -1):
            idx = [i for i, k in enumerate(ks) if k >= j]
            neu = self.Pe.inv([ys[i] for i in idx])
            for i, v in zip(idx, neu):
                ys[i] = v
        return [mp.re(v) for v in ys]
