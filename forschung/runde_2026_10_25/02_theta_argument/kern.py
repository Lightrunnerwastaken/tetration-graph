"""Etappe 2: gemeinsamer Kern (Motor-Zugang, Fixpunkt, Schröder-Reihe, θ aus sexp-Werten).

Motor: fatou.gp-Fork wie in ../motor.py, aber mit eigenem Zwischenspeicher etappe2/cache/<basis>_<dps>.json
(der Etappe-1-Speicher bleibt unangetastet).
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
CDIR = HIER / "cache"
_motoren = {}


def _cfile(basis, dps):
    tag = "e" if basis == "exp(1)" else basis.replace("(", "").replace(")", "")
    return CDIR / f"{tag}_{dps}.json"


def motor(basis, dps):
    key = (basis, dps)
    if key not in _motoren:
        if str(GP_SRC) not in sys.path:
            sys.path.insert(0, str(GP_SRC))
        from fatou_backend import FatouGP
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=900.0, init_timeout=1100.0)
        t = time.time()
        s = m.session(basis)
        s.eval_batch(["1"])
        print(f"[motor] {basis} dps={dps} bereit nach {time.time()-t:.1f}s", flush=True)
        _motoren[key] = s
    return _motoren[key]


def gpzahl(x, digits):
    x = mp.mpc(x)
    f = lambda v: mp.nstr(v, digits, min_fixed=-mp.inf, max_fixed=mp.inf)
    if mp.im(x) == 0:
        return f"({f(mp.re(x))})"
    return f"({f(mp.re(x))}+({f(mp.im(x))})*I)"


def werte(basis, ausdruecke, dps, chunk=200):
    CDIR.mkdir(exist_ok=True)
    cf = _cfile(basis, dps)
    c = json.loads(cf.read_text(encoding="utf-8")) if cf.is_file() else {}
    fehlt = [a for a in dict.fromkeys(ausdruecke) if a not in c]
    for i in range(0, len(fehlt), chunk):
        teil = fehlt[i:i + chunk]
        res = motor(basis, dps).eval_batch(teil)
        for a, w in zip(teil, res):
            c[a] = [mp.nstr(mp.re(w), dps + 10), mp.nstr(mp.im(w), dps + 10)]
        cf.write_text(json.dumps(c, indent=0), encoding="utf-8")
    with mp.workdps(dps + 20):
        return [mp.mpc(mp.mpf(c[a][0]), mp.mpf(c[a][1])) for a in ausdruecke]


def sexp(basis, zs, dps):
    return werte(basis, [f"sexp({gpzahl(z, dps + 10)})" for z in zs], dps)


class Fix:
    """Fixpunkt L von b^w (Im L > 0) und Schröder-Funktion σ dort, σ'(L) = 1."""

    def __init__(self, b, K=40):
        self.b = mp.mpf(b) if not isinstance(b, mp.mpf) else b
        lb = mp.log(self.b)
        self.lb = lb
        # b^L = L  <=>  L = -W(-ln b)/ln b ; Zweig mit Im L > 0
        Lc = -mp.lambertw(-lb) / lb
        self.L = mp.conj(Lc) if mp.im(Lc) < 0 else Lc
        self.lam = self.L * lb                      # Multiplikator  (b^w)' = ln b · b^L = λ
        self.K = K
        self.S = self._reihe(K)

    def _reihe(self, K):
        L, lam = self.L, self.lam
        # g(u) = b^{L+u} − L = L(e^{lb u} − 1)
        g = [mp.mpc(0)] + [L * self.lb ** k / mp.factorial(k) for k in range(1, K + 1)]

        def mul(a, b_):
            c = [mp.mpc(0)] * (K + 1)
            for i in range(K + 1):
                if a[i] != 0:
                    for j in range(K + 1 - i):
                        c[i + j] += a[i] * b_[j]
            return c
        gp = [None, g]
        for _ in range(2, K + 1):
            gp.append(mul(gp[-1], g))
        s = [mp.mpc(0), mp.mpc(1)] + [mp.mpc(0)] * (K - 1)
        for m in range(2, K + 1):
            rest = mp.fsum(s[k] * gp[k][m] for k in range(1, m))
            s[m] = rest / (lam - lam ** m)
        return s

    def logb(self, w):
        return mp.log(w) / self.lb

    def h(self, w, schritte=None):
        """Abelfunktion der regulären Superfunktion: F(h(w)) = w, F(s) = σ⁻¹(λ^s); Hauptzweige der Logarithmen."""
        n = 0
        while abs(w - self.L) > mp.mpf("1e-3") or n < 5:
            w = self.logb(w)
            n += 1
            if n > 5000:
                raise RuntimeError("keine Konvergenz")
        u = w - self.L
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while abs(u) ** (self.K + 1) > tol * abs(u):
            u = self.logb(self.L + u) - self.L
            n += 1
        su = mp.polyval(self.S[::-1], u)
        return n + mp.log(su) / mp.log(self.lam)
