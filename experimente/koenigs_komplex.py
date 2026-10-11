"""Koenigs-Linearisierung für z ↦ exp(z·u) mit festem Logarithmus u = Log b (komplex), anziehender Fixpunkt.

Parametrisierung nach der Multiplikator-Rezeptur (Justus, Complex Tetration Dynamics, 2026):
  u = λ e^{−λ},  Fixpunkt L = e^{λ},  f'(L) = L·u = λ.
σ(z) = lim λ^{−n}(f^n(z) − L) über die Reihe am Fixpunkt; T(t) = σ^{-1}(λ^t σ(1)) durch Fortsetzung in t.
"""
import mpmath as mp

from regulaer import horner, schroeder_reihe


class KoenigsKomplex:
    def __init__(self, lam, N=40, tol="1e-3", log_lam=None):
        self.lam = mp.mpc(lam)
        self.u = self.lam * mp.exp(-self.lam)          # Log b
        self.L = mp.exp(self.lam)
        self.log_lam = mp.log(self.lam) if log_lam is None else mp.mpc(log_lam)
        self.tol = mp.mpf(tol)
        g = [mp.mpc(0)] + [self.L * self.u ** k / mp.factorial(k) for k in range(1, N + 1)]
        self.s = schroeder_reihe(g, self.lam, N)

    def f(self, z):
        return mp.exp(z * self.u)

    def sigma(self, z):
        n = 0
        while abs(z - self.L) > self.tol:
            z = self.f(z)
            n += 1
            if n > 20000:
                raise RuntimeError("keine Konvergenz gegen den Fixpunkt")
        return horner(self.s, z - self.L) / self.lam ** n

    def sigma_inv(self, w, start):
        """Lösung von σ(z) = w durch Newton, ausgehend von start."""
        z = mp.mpc(start)
        h = mp.mpf(10) ** (-mp.mp.dps // 2)
        for _ in range(60):
            fz = self.sigma(z) - w
            d = (self.sigma(z + h) - self.sigma(z - h)) / (2 * h)
            dz = fz / d
            z -= dz
            if abs(dz) < mp.mpf(10) ** (-mp.mp.dps + 8):
                break
        return z

    def tet(self, t, schritte=20):
        """T(t) mit T(0) = 1, stetig in s ∈ [0, t] fortgesetzt; λ^s = exp(s·log_lam)."""
        s1 = self.sigma(mp.mpc(1))
        z = mp.mpc(1)
        for j in range(1, schritte + 1):
            s = mp.mpf(t) * j / schritte
            z = self.sigma_inv(mp.exp(s * self.log_lam) * s1, z)
        return z
