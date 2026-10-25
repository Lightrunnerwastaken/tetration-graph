"""Gebiet 05: Knesers Welle für komplexe Basen (Paulsen–Cowgill).

Motor: fatou.gp-Fork (nur lesend), eigener Zwischenspeicher cache/<tag>_<dps>.json in diesem Ordner, state_cache=False.
Basis b komplex, ln b = Hauptzweig (wie im Motor: kc = log(log(B)) + 1).

Fix(b, Lstart): Fixpunkt L von b^w (Newton ab Lstart), Multiplikator λ = L ln b, Schröder-Reihe σ (σ'(L)=1),
Zweig k des Umkehrzweigs log_b^{(k)} w = (log w + 2πik)/ln b mit log_b^{(k)}(L) = L.
h(w) = Abelfunktion der regulären Superfunktion F(s) = σ⁻¹(λ^s) an L, über Iteration von log_b^{(k)}.
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


def _tag(basis):
    return "".join(ch if ch.isalnum() or ch in ".-" else "_" for ch in basis)


def motor(basis, dps):
    key = (basis, dps)
    if key not in _motoren:
        if str(GP_SRC) not in sys.path:
            sys.path.insert(0, str(GP_SRC))
        from fatou_backend import FatouGP
        m = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=280.0, init_timeout=280.0, state_cache=False)
        t = time.time()
        s = m.session(basis)
        s.eval_batch(["1"])
        print(f"[motor] {basis} dps={dps} bereit nach {time.time()-t:.1f}s", flush=True)
        _motoren[key] = (m, s)
    return _motoren[key][1]


def schliessen():
    for m, s in _motoren.values():
        try:
            m.close()
        except Exception:
            pass
    _motoren.clear()


def gpzahl(x, digits):
    x = mp.mpc(x)
    f = lambda v: mp.nstr(v, digits, min_fixed=-mp.inf, max_fixed=mp.inf)
    return f"({f(mp.re(x))}+({f(mp.im(x))})*I)"


def werte(basis, ausdruecke, dps, chunk=128):
    CDIR.mkdir(exist_ok=True)
    cf = CDIR / f"{_tag(basis)}_{dps}.json"
    c = json.loads(cf.read_text(encoding="utf-8")) if cf.is_file() else {}
    fehlt = [a for a in dict.fromkeys(ausdruecke) if a not in c]
    for i in range(0, len(fehlt), chunk):
        teil = fehlt[i:i + chunk]
        res = motor(basis, dps).eval_batch(teil)
        for a, w in zip(teil, res):
            c[a] = [mp.nstr(mp.re(w), dps + 10), mp.nstr(mp.im(w), dps + 10)]
        for _ in range(20):
            try:
                tmp = cf.with_suffix(".tmp")
                tmp.write_text(json.dumps(c, indent=0), encoding="utf-8")
                tmp.replace(cf)
                break
            except OSError:
                time.sleep(0.5)
    with mp.workdps(dps + 20):
        return [mp.mpc(mp.mpf(c[a][0]), mp.mpf(c[a][1])) for a in ausdruecke]


def sexp(basis, zs, dps):
    return werte(basis, [f"sexp({gpzahl(z, dps + 10)})" for z in zs], dps)


class Fix:
    def __init__(self, b, Lstart, K=40):
        self.b = mp.mpc(b)
        lb = mp.log(self.b)
        self.lb = lb
        L = mp.mpc(Lstart)
        for _ in range(200):
            f = mp.exp(L * lb) - L
            d = lb * mp.exp(L * lb) - 1
            L1 = L - f / d
            if abs(L1 - L) < mp.mpf(10) ** (-mp.mp.dps + 5):
                L = L1
                break
            L = L1
        self.L = L
        self.lam = L * lb
        self.k = int(mp.nint(mp.im(self.lam - mp.log(L)) / (2 * mp.pi)))
        self.K = K
        self.S = self._reihe(K)
        self.per = 2j * mp.pi / mp.log(self.lam)

    def _reihe(self, K):
        L, lam = self.L, self.lam
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
        return (mp.log(w) + 2j * mp.pi * self.k) / self.lb

    def h(self, w, nmax=20000):
        n = 0
        while abs(w - self.L) > mp.mpf("1e-3") or n < 5:
            w = self.logb(w)
            n += 1
            if n > nmax:
                raise RuntimeError("keine Konvergenz")
        u = w - self.L
        tol = mp.mpf(10) ** (-mp.mp.dps)
        while abs(u) ** (self.K + 1) > tol * abs(u):
            u = self.logb(self.L + u) - self.L
            n += 1
            if n > nmax:
                raise RuntimeError("keine Konvergenz (Feinschritt)")
        su = mp.polyval(self.S[::-1], u)
        return n + mp.log(su) / mp.log(self.lam)


def fft(v):
    n = len(v)
    if n == 1:
        return v[:]
    ev, od = fft(v[0::2]), fft(v[1::2])
    out = [None] * n
    for k in range(n // 2):
        tw = mp.expjpi(-2 * mp.mpf(k) / n) * od[k]
        out[k] = ev[k] + tw
        out[k + n // 2] = ev[k] - tw
    return out


def entwirren(vals, per):
    for j in range(1, len(vals)):
        m = min(range(-6, 7), key=lambda m: abs(vals[j] + m * per - vals[j - 1]))
        vals[j] += m * per
    return vals
