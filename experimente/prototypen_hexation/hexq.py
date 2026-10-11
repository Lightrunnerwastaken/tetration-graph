"""Prototyp: Kneser-artige Hexation — reguläre Superfunktion der Pentation am abstossenden komplexen Fixpunkt q (nur Vorwärtsschritte)."""
import sys
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/experimente")
import mpmath as mp
from stufen import pentation, poincare, psi_wert
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/experimente/prototypen_hexation")
import cache_drossel  # noqa: F401

mp.mp.dps = 40


class HexQ:
    def __init__(self, B, L0, q0, r=mp.mpf("0.1"), N=64, K=30, KP=28):
        self.B = B
        self.L, self.lam, self.psiT, self.c, _, self.P = pentation(B, L0, mp.mpf("0.08"))
        q = mp.mpc(q0)
        for _ in range(20):
            a = self.taylorP(q, r, N, K)
            dq = (a[0] - q) / (a[1] - 1)
            q -= dq
            if abs(dq) < mp.mpf(10) ** -30:
                break
        self.q = q
        self.a = self.taylorP(q, r, N, K)
        self.mu = self.a[1]
        self.lnmu = mp.log(self.mu)
        self.psi = poincare(self.a, K=KP)
        # lokaler Radius: Glieder bis 1e-35
        self.rlok = min((mp.mpf(10) ** -35 / abs(self.psi[k])) ** (mp.mpf(1) / (k + 1)) for k in range(KP - 4, KP) if self.psi[k] != 0)

    def taylorP(self, x0, r, N, K):
        zs = [x0 + r * mp.expjpi(2 * mp.mpf(j) / N) for j in range(N)]
        w = self.P(zs)
        return [mp.fsum(w[j] * mp.expjpi(-2 * m * mp.mpf(j) / N) for j in range(N)) / N / r ** m for m in range(K + 1)]

    def Psi(self, ss):
        """Poincaré-Funktion: Ψ(s) = P^n(Ψ_lok(s/μ^n))."""
        ns, ys = [], []
        for s in ss:
            n = 0
            while abs(s / self.mu ** n) > self.rlok:
                n += 1
            ns.append(n)
            ys.append(psi_wert(self.q, self.psi, s / self.mu ** n))
        for k in range(1, max(ns) + 1):
            idx = [i for i, n in enumerate(ns) if n >= k]
            neu = self.P([ys[i] for i in idx])
            for i, v in zip(idx, neu):
                ys[i] = v
        return ys

    def F(self, zetas):
        return self.Psi([mp.exp(z * self.lnmu) for z in zetas])
