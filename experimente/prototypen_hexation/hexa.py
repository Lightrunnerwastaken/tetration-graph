"""Prototyp: Kneser-artige Hexation für b > η₅ — reguläre Superfunktion der Pentation am komplexen (anziehenden) Fixpunkt p."""
import sys
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/experimente")
import mpmath as mp
from gemeinsam import kneser
from stufen import pentation, poincare, psi_wert, psi_umkehr

mp.mp.dps = 40


class Hexa:
    def __init__(self, B, L0, p0, r=mp.mpf("0.15"), N=64, K=30, KP=25):
        self.B = B
        self.L, self.lam, self.psiT, self.c, _, self.P = pentation(B, L0, mp.mpf("0.08"))
        self.lnlam = mp.log(self.lam)
        p = mp.mpc(p0)
        for _ in range(20):                                       # Fixpunkt verfeinern über Taylorreihe
            a = self.taylorP(p, r, N, K)
            dp = (a[0] - p) / (a[1] - 1)
            p -= dp
            if abs(dp) < mp.mpf(10) ** -30:
                break
        self.p = p
        self.a = self.taylorP(p, r, N, K)
        self.mu = self.a[1]
        self.lnmu = mp.log(self.mu)
        self.psiP = poincare(self.a, K=KP)

    def taylorP(self, x0, r, N, K):
        zs = [x0 + r * mp.expjpi(2 * mp.mpf(j) / N) for j in range(N)]
        w = self.P(zs)
        return [mp.fsum(w[j] * mp.expjpi(-2 * m * mp.mpf(j) / N) for j in range(N)) / N / r ** m for m in range(K + 1)]

    def Pinv(self, ws, m=7, k=0):
        zs = list(ws)
        for _ in range(m):
            zs = kneser("slog", self.B, zs)
        return [m + (mp.log(psi_umkehr(self.L, self.psiT, z0) / self.c) + 2j * mp.pi * k) / self.lnlam for z0 in zs]

    def PsiP(self, ss, rlok=mp.mpf("0.05")):
        ns, ys = [], []
        for s in ss:
            n = 0
            while abs(self.mu ** n * s) > rlok:
                n += 1
            ns.append(n)
            ys.append(psi_wert(self.p, self.psiP, self.mu ** n * s))
        for k in range(max(ns), 0, -1):
            idx = [i for i, n in enumerate(ns) if n >= k]
            neu = self.Pinv([ys[i] for i in idx])
            for i, v in zip(idx, neu):
                ys[i] = v
        return ys, ns

    def Fp(self, zetas):
        return self.PsiP([mp.exp(z * self.lnmu) for z in zetas])
