"""Gebiet 20: Differenz der regulären Iterationen an den zwei reellen Fixpunkten a < A von b^x (1 < b < η).

Parametrisierung (wie Runde 18): a = (1+1/t)^t, A = (1+1/t)^(t+1), b = a^(1/a), λ_a = ln a, λ_A = ln A.
t -> 0: b -> 1;  t = 1: b = √2;  t -> ∞: b -> η.

Zwei periodische Differenzen (beide 1-periodisch):
  p(θ) = α_A(T_a(θ)) - θ     (Variable θ der unteren regulären Iteration, T_a(θ) = Ψ_a(λ_a^θ) in (a, A))
  q(s) = α_a(T_A(s)) - s     (Variable s der oberen regulären Iteration, T_A(s) = Ψ_A(-λ_A^s) in (a, A))
Fourier: c_k = (1/N) Σ_j f(j/N) e^{-2πi k j/N}.
"""
import mpmath as mp


def param(t):
    t = mp.mpf(t)
    l = mp.log1p(1 / t)
    la, lA = t * l, (t + 1) * l
    a, A = mp.exp(la), mp.exp(lA)
    return a, A, mp.exp(la / a), la, lA


def reihen(p, c, K):
    """Ψ_p(u) = p + Σ q_k u^k mit Ψ(λu) = exp(c Ψ(u)), λ = p c."""
    lam = p * c
    q = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (K - 1)
    for k in range(2, K + 1):
        ph = [mp.mpf(0)] + [c * q[j] for j in range(1, k)] + [mp.mpf(0)] * (K + 1 - k)
        e = [mp.mpf(1)] + [mp.mpf(0)] * k
        for m in range(1, k + 1):
            e[m] = mp.fsum(j * ph[j] * e[m - j] for j in range(1, m + 1)) / m
        q[k] = p * e[k] / (lam ** k - lam)
    return lam, q


class Paar:
    def __init__(self, t, K=40):
        self.t = mp.mpf(t)
        self.a, self.A, self.b, self.la, self.lA = param(t)
        self.c = mp.log(self.b)
        self.lam_a, self.qa = reihen(self.a, self.c, K)
        self.lam_A, self.qA = reihen(self.A, self.c, K)
        self.Lla, self.LlA = mp.log(self.lam_a), mp.log(self.lam_A)
        # Radius: Glieder klein genug
        self.eps = mp.mpf(10) ** (-(mp.mp.dps // K) - 1)
        self.dqa = [k * self.qa[k] for k in range(1, len(self.qa))]
        self.dqA = [k * self.qA[k] for k in range(1, len(self.qA))]

    # Grundbausteine
    def f(self, x):
        return mp.exp(self.c * x)

    def g(self, x):  # Umkehrung log_b
        return mp.log(x) / self.c

    @staticmethod
    def poly(q, u):
        return mp.polyval(q[::-1], u)

    def sigma(self, y, p, q, dq):
        u = y - p
        for _ in range(200):
            du = (p + self.poly(q, u) - y) / mp.polyval(dq[::-1], u)
            u -= du
            if abs(du) <= abs(u) * mp.mpf(10) ** (-mp.mp.dps + 3):
                break
        return u

    # untere reguläre Iteration
    def alpha_a(self, x):
        n = 0
        while abs(x - self.a) > self.eps:
            x = self.f(x)
            n += 1
        u = self.sigma(x, self.a, self.qa, self.dqa)
        return -n + mp.log(u) / self.Lla      # log: Hauptzweig (u<0 gibt iπ/ln λ_a)

    def T_a(self, th):
        n = 0
        u = mp.exp(th * self.Lla)
        while abs(u) > self.eps:
            u *= self.lam_a
            n += 1
        x = self.a + self.poly(self.qa, u)
        for _ in range(n):
            x = self.g(x)
        return x

    # obere reguläre Iteration
    def alpha_A(self, x):
        n = 0
        while abs(x - self.A) > self.eps:
            x = self.g(x)
            n += 1
        u = self.sigma(x, self.A, self.qA, self.dqA)
        return n + mp.log(-u) / self.LlA

    def T_A(self, s):
        n = 0
        u = -mp.exp(s * self.LlA)
        while abs(u) > self.eps:
            u /= self.lam_A
            n += 1
        x = self.A + self.poly(self.qA, u)
        for _ in range(n):
            x = self.f(x)
        return x

    def p(self, th):
        return self.alpha_A(self.T_a(th)) - th

    def q(self, s):
        return self.alpha_a(self.T_A(s)) - s

    def fourier(self, fun, N):
        vals = [fun(mp.mpf(j) / N) for j in range(N)]
        cs = []
        for k in range(N // 2):
            cs.append(mp.fsum(v * mp.expjpi(-2 * mp.mpf(k * j) / N) for j, v in enumerate(vals)) / N)
        return cs

    def S_a(self):
        return 2 * mp.pi ** 2 / abs(self.Lla)

    def S_A(self):
        return 2 * mp.pi ** 2 / abs(self.LlA)
