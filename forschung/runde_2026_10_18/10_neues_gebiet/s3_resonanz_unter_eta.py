"""s3: Zwei reelle Fixpunkte a < A von b^x (1 < b < η): reguläre Iteration an a gegen die an A.

Satz T3 (BERICHT): Ein gemeinsamer holomorpher Erzeuger an a und A ist nur möglich, wenn λ_a = λ_A^m (m ∈ Z),
mit λ_a = ln a, λ_A = ln A. Wegen ln(1+x) < x/sqrt(1+x) ist λ_a λ_A < 1 (m = -1 unmöglich);
für k >= 2 gibt es aber Resonanzbasen b_k mit λ_a λ_A^k = 1.
Hier: Abstand der beiden regulären Iterationen als 1-periodische Funktion
    P(s) = h_a(F_A(s)) - s,   h_a Abelfunktion an a,  F_A reguläre Superfunktion an A (von A nach a fallend),
und ihr erster Fourier-Koeffizient c_1 für viele Basen, insbesondere an den Resonanzen.
Gegenprobe: dps 60 gegen dps 90, N = 48 gegen N = 64 Stützstellen.
"""
import mpmath as mp


def param(t):
    """a = (1+1/t)^t, A = (1+1/t)^(t+1): a^A = A^a; Basis b = a^(1/a)."""
    l = mp.log1p(1 / t)
    la, lA = t * l, (t + 1) * l
    a, A = mp.exp(la), mp.exp(lA)
    return a, A, mp.exp(la / a), la, lA   # λ_a = ln a, λ_A = ln A


def reihen(p, c, K):
    """Ψ_p(u) = p + Σ q_k u^k mit Ψ(λu) = exp(c Ψ(u)), λ = p c."""
    lam = p * c
    q = [mp.mpf(0), mp.mpf(1)] + [mp.mpf(0)] * (K - 1)
    for k in range(2, K + 1):
        # [u^k] exp(c φ_{<k}(u)), φ = Σ_{j<k} q_j u^j
        ph = [mp.mpf(0)] + [c * q[j] for j in range(1, k)] + [mp.mpf(0)] * (K + 1 - k)
        e = [mp.mpf(1)] + [mp.mpf(0)] * k
        for m in range(1, k + 1):
            e[m] = mp.fsum(j * ph[j] * e[m - j] for j in range(1, m + 1)) / m
        q[k] = p * e[k] / (lam ** k - lam)
    return lam, q


class Paar:
    def __init__(self, t, K=40):
        self.a, self.A, self.b, self.la, self.lA = param(t)
        self.c = mp.log(self.b)
        self.lam_a, self.qa = reihen(self.a, self.c, K)
        self.lam_A, self.qA = reihen(self.A, self.c, K)
        self.eps = mp.mpf(10) ** (-(mp.mp.dps // K) - 1)

    def f(self, x):
        return mp.exp(self.c * x)

    def Psi(self, q, u):
        """Σ q_k u^k  (= Ψ_p(u) - p)"""
        return mp.polyval(q[::-1], u)

    def sigma_a(self, y):
        """σ_a(y), y nahe a: Newton für Ψ_a(u) = y."""
        u = y - self.a
        dq = [k * self.qa[k] for k in range(1, len(self.qa))]
        for _ in range(100):
            du = (self.a + self.Psi(self.qa, u) - y) / mp.polyval(dq[::-1], u)
            u -= du
            if abs(du) <= abs(u) * mp.mpf(10) ** (-mp.mp.dps + 2):
                break
        return u

    def h_a(self, x):
        n = 0
        while abs(x - self.a) > self.eps:
            x = self.f(x)
            n += 1
        return -n + mp.log(self.sigma_a(x)) / mp.log(self.lam_a)  # σ(x) = λ^{-n} σ_loc(f^n x); σ>0 für x>a

    def F_A(self, s):
        n = 0
        u = -mp.exp(s * mp.log(self.lam_A))
        while abs(u) > self.eps:
            u /= self.lam_A
            n += 1
        x = self.A + self.Psi(self.qA, u)
        for _ in range(n):
            x = self.f(x)
        return x

    def c1(self, N=48):
        vals = [self.h_a(self.F_A(mp.mpf(j) / N)) - mp.mpf(j) / N for j in range(N)]
        c0 = mp.fsum(vals) / N
        c1 = mp.fsum(v * mp.expjpi(-2 * mp.mpf(j) / N) for j, v in enumerate(vals)) / N
        return c0, c1


def resonanz(k):
    g = lambda t: (lambda a, A, b, la, lA: la * lA ** k - 1)(*param(t))
    # Vorzeichenwechsel suchen
    ts = [mp.mpf(10) ** (x / 4) for x in range(-40, 13)]
    for t0, t1 in zip(ts, ts[1:]):
        if g(t0) * g(t1) < 0:
            return mp.findroot(g, (t0, t1), solver="anderson")
    return None


if __name__ == "__main__":
    mp.mp.dps = 30
    print("Produkt λ_a λ_A (soll < 1):", [mp.nstr(param(mp.mpf(t))[3] * param(mp.mpf(t))[4], 8) for t in ["0.01", "0.5", "1", "2", "10", "100", "1e4"]])
    res = {}
    for k in [2, 3, 4]:
        res[k] = resonanz(k)
        a, A, b, la, lA = param(res[k])
        print(f"Resonanz k={k}: t={mp.nstr(res[k], 20)}  b={mp.nstr(b, 20)}  a={mp.nstr(a, 12)}  A={mp.nstr(A, 12)}")
    import sys
    tlist = [("sqrt2", mp.mpf(1))]  # t=1: a=2, A=4, b=sqrt2
    for k in [2, 3]:
        for d in ["-0.05", "0", "0.05"]:
            tlist.append((f"k={k}{'+' if d[0] != '-' else ''}{d}", res[k] * (1 + mp.mpf(d))))
    print("\n  Fall            b                    |c1| dps60/N48         |c1| dps90/N64     rel. Abw.")
    for name, t in tlist:
        out = []
        for dps, N in [(60, 48), (90, 64)]:
            mp.mp.dps = dps
            P = Paar(mp.mpf(t), K=40)
            c0, c1 = P.c1(N)
            out.append((P.b, c0, c1))
        mp.mp.dps = 30
        (b, _, c1a), (_, _, c1b) = out
        print(f"  {name:12s} {mp.nstr(b, 18):>20}  {mp.nstr(abs(c1a), 10):>18}  {mp.nstr(abs(c1b), 10):>18}  {mp.nstr(abs(c1a-c1b)/abs(c1b), 3)}")
        sys.stdout.flush()
