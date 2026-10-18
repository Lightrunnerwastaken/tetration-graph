"""Etappe 3: analytische Bausteine für Basis e (ohne Motor).

L      Fixpunkt e^L = L, Im L > 0 (Multiplikator λ = L, log λ = L)
sigma  Schröder-Funktion an L, σ(e^w) = L σ(w), σ'(L) = 1, Fortsetzung über Hauptzweig-Logarithmen
Psi    Umkehrung, Ψ(Lu) = e^{Ψ(u)}, ganz
F(s)   = Ψ(e^{Ls})  reguläre Superfunktion;  h = F^{-1}, e^{L h(w)} = σ(w)
theta  Knesers θ aus Fourier-Koeffizienten (etappe2, 44 Stellen, n ≤ 4096)
"""
import json
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E2 = HIER.parent.parent / "identitaeten_2026_10_14" / "etappe2"  # nur lesend


class Basis:
    def __init__(self, K=60):
        self.L = -mp.lambertw(-1)
        if mp.im(self.L) < 0:
            self.L = mp.conj(self.L)
        self.K = K
        L = self.L
        # Schröder-Reihe σ(L+u) = Σ s_k u^k  (wie etappe2 kern.Fix)
        g = [mp.mpc(0)] + [L / mp.factorial(k) for k in range(1, K + 1)]
        gp = [None, g]
        for _ in range(2, K + 1):
            gp.append(self._mul(gp[-1], g))
        s = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
        for m in range(2, K + 1):
            s[m] = mp.fsum(s[k] * gp[k][m] for k in range(1, m)) / (L - L ** m)
        self.S = s
        self.dS = [k * s[k] for k in range(1, K + 1)]
        # Ψ(u) = L + Σ p_k u^k mit Ψ(Lu) = exp(Ψ(u))
        p = [L, mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
        for k in range(2, K + 1):
            # [u^k] exp(Ψ(u)) = L * [u^k] exp(Σ_{j>=1} p_j u^j)
            q = [mp.mpc(0)] + p[1:k] + [mp.mpc(0)] * (K + 1 - k)
            ex = self._exp(q, k)
            p[k] = L * ex[k] / (L ** k - L)
        self.P = p
        self.dP = [j * p[j] for j in range(1, K + 1)]

    def _mul(self, a, b):
        K = self.K
        c = [mp.mpc(0)] * (K + 1)
        for i in range(K + 1):
            if a[i] != 0:
                for j in range(K + 1 - i):
                    c[i + j] += a[i] * b[j]
        return c

    @staticmethod
    def _exp(q, n):
        # exp einer Potenzreihe q (q[0] = 0) bis Ordnung n
        e = [mp.mpc(1)] + [mp.mpc(0)] * n
        for m in range(1, n + 1):
            e[m] = mp.fsum(j * q[j] * e[m - j] for j in range(1, m + 1)) / m
        return e

    # ---------- Ψ, F
    def Psi(self, u, ableitung=False):
        n = 0
        while abs(u) > mp.mpf("0.02"):
            u = u / self.L
            n += 1
        w = mp.polyval(self.P[::-1], u)
        d = mp.polyval(self.dP[::-1], u) if ableitung else None
        for _ in range(n):
            ew = mp.exp(w)
            if ableitung:
                d = d * ew / self.L
            w = ew
        return (w, d) if ableitung else w

    def F(self, s, ableitung=False):
        u = mp.exp(self.L * s)
        if not ableitung:
            return self.Psi(u)
        w, d = self.Psi(u, True)
        return w, d * self.L * u

    # ---------- σ, h (Hauptzweige), mit logarithmischer Ableitung
    def sigma(self, w, ableitung=False):
        """σ(w) = L^n σ(u_n), u_{k+1} = log u_k; dazu σ'(w)/σ(w)."""
        n = 0
        prod = mp.mpc(1)
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while abs(w - self.L) > mp.mpf("1e-3") or abs(w - self.L) ** (self.K + 1) > tol * abs(w - self.L):
            if ableitung:
                prod /= w
            w = mp.log(w)
            n += 1
            if n > 20000:
                raise RuntimeError("sigma: keine Konvergenz")
        u = w - self.L
        su = mp.polyval(self.S[::-1], u)
        val = self.L ** n * su
        if not ableitung:
            return val, n, su
        dsu = mp.polyval(self.dS[::-1], u)
        return val, n, su, prod * dsu / su

    def h(self, w):
        _, n, su = self.sigma(w)
        return n + mp.log(su) / self.L

    def h_mit_ableitung(self, w):
        _, n, su, dlog = self.sigma(w, True)
        return n + mp.log(su) / self.L, dlog / self.L


def lade_theta(datei=E2 / "daten" / "a_e_80_0.003_16384.json", nmax=None):
    d = json.loads(Path(datei).read_text(encoding="utf-8"))
    a = [mp.mpc(mp.mpf(x[0]), mp.mpf(x[1])) for x in d["a"]]
    if nmax:
        a = a[: nmax + 1]
    return a


class ThetaReihe:
    def __init__(self, a, B):
        self.a = a
        self.B = B

    def theta(self, z, ableitung=0):
        q = mp.expjpi(2 * z)
        # Horner für Σ a_n q^n und Ableitungen nach z
        n = len(self.a) - 1
        th = mp.mpc(0); d1 = mp.mpc(0); d2 = mp.mpc(0)
        qn = mp.mpc(1)
        tpi = 2j * mp.pi
        res0 = mp.fsum([])
        out0 = []; out1 = []; out2 = []
        for k, ak in enumerate(self.a):
            t = ak * qn
            out0.append(t)
            if ableitung >= 1:
                out1.append(t * tpi * k)
            if ableitung >= 2:
                out2.append(t * (tpi * k) ** 2)
            qn *= q
            if abs(qn) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        r = [mp.fsum(out0)]
        if ableitung >= 1:
            r.append(mp.fsum(out1))
        if ableitung >= 2:
            r.append(mp.fsum(out2))
        return r if ableitung else r[0]

    def sexp(self, z, ableitung=False):
        """sexp(z) = F(z + θ(z)) für Im z > 0 (ausreichend gross für die Reihenkonvergenz)."""
        if not ableitung:
            return self.B.F(z + self.theta(z))
        th, th1 = self.theta(z, 1)
        w, dF = self.B.F(z + th, True)
        return w, dF * (1 + th1)
