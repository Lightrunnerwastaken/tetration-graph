"""Poincaré- und Schröder-Funktion an einem Fixpunkt L von f(z) = exp(u z) (anziehend oder abstossend)."""
import mpmath as mp

ZWEI_PI_I = 2j * mp.pi


class AmFixpunkt:
    """Poincaré-Funktion Ψ (Ψ(λw) = f(Ψ(w)), Ψ(0) = L, Ψ'(0) = 1) und Schröder-Funktion σ = Ψ^{-1} für f(z) = exp(u z)."""

    def __init__(self, u, L, N=70, log_lam=None):
        self.u, self.L = mp.mpc(u), mp.mpc(L)
        assert abs(mp.exp(self.u * self.L) - self.L) < mp.mpf(10) ** -30, "L ist kein Fixpunkt"
        self.lam = self.L * self.u
        self.log_lam = mp.log(self.lam) if log_lam is None else mp.mpc(log_lam)
        p, e = [mp.mpc(0), mp.mpc(1)], [mp.mpc(1), self.u]
        for k in range(2, N + 1):
            S = mp.fsum(j * self.u * p[j] * e[k - j] for j in range(1, k)) / k
            pk = self.L * S / (self.lam ** k - self.lam)
            p.append(pk)
            e.append(self.u * pk + S)
        self.p = p
        self.r0 = min((mp.mpf(10) ** -45 / abs(p[k])) ** (mp.mpf(1) / k) for k in range(N - 5, N + 1) if p[k] != 0)
        self.k_ast = int(mp.nint(((self.lam - mp.log(self.L)) / ZWEI_PI_I).real))   # Zweig von log mit g(L) = L

    def reihe(self, w):
        s = mp.mpc(0)
        for c in reversed(self.p[1:]):
            s = (s + c) * w
        return self.L + s

    def Psi(self, w):
        w, n = mp.mpc(w), 0
        anz = abs(self.lam) < 1                       # anziehend: Ψ(w) = g^n(Ψ(λ^n w)), abstossend: Ψ(w) = f^n(Ψ(w/λ^n))
        while abs(w) > self.r0:
            w = w * self.lam if anz else w / self.lam
            n += 1
        z = self.reihe(w)
        for _ in range(n):
            z = self.g(z) if anz else mp.exp(self.u * z)
        return z

    def g(self, z):
        return (mp.log(z) + ZWEI_PI_I * self.k_ast) / self.u

    def sigma(self, z, rho=None):
        rho = rho or self.r0 / 4
        z, n = mp.mpc(z), 0
        anz = abs(self.lam) < 1                       # anziehend: σ(z) = λ^{−n} σ(f^n z), abstossend: σ(z) = λ^n σ(g^n z)
        while abs(z - self.L) > rho:
            if z == 0 and not anz:
                raise ValueError("Rückwärtsbahn trifft 0")
            z = mp.exp(self.u * z) if anz else self.g(z)
            n += 1
            if n > 4000:
                raise ValueError("Bahn erreicht L nicht")
        w = z - self.L
        for _ in range(60):
            dw = (self.reihe(w) - z) / ((self.reihe(w + mp.mpf(10) ** -20) - self.reihe(w - mp.mpf(10) ** -20)) / (2 * mp.mpf(10) ** -20))
            w -= dw
            if abs(dw) < mp.mpf(10) ** (-mp.mp.dps + 6) * (1 + abs(w)):
                break
        return w * (self.lam ** (-n) if anz else self.lam ** n)
