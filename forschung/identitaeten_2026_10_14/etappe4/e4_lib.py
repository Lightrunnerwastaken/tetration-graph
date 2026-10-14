"""Etappe 4: gestörte Tetrationen T(x) = sexp(x + p(x)) auf der reellen Achse (Basis e).

θ_T(x) = h(T x) − x  (stetiger Zweig auf (0,1)),  1 + θ_T'(x) = T'(x) h'(T x).
Gebraucht wird nur Kneser auf der reellen Achse (Taylor-Daten aus etappe3) und h, σ (lib3.Basis).

Störung p: 1-periodisch, reell, p(0) = 0, 1 + p' > 0. Formate:
  ("trig", [(k, a_k, b_k), ...])  p = Σ a_k sin 2πkx + b_k (cos 2πkx − 1)
  ("moeb", a)                       Möbius-Kreisdiffeo: φ(x) = x + p(x) mit e^{2πiφ} = e^{iα}(q−a)/(1−āq), φ(0)=0
  ("lin", None)                     lineare Tetration (Hooshmand): T(x) = e^x auf [0,1]  (kein p nötig)
"""
import json
import sys
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
E3 = HIER.parent / "etappe3"
sys.path.insert(0, str(E3))
from lib3 import Basis  # noqa: E402

_T = json.loads((E3 / "daten" / "taylor_130_256.json").read_text(encoding="utf-8"))
_ZEN_STR = _T["zentren"]
_ZEN = None
B = None
L = None
PER = None


def init(dps):
    global _ZEN, B, L, PER
    mp.mp.dps = dps
    _ZEN = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in _ZEN_STR.items()}
    B = Basis(70)
    L = B.L
    PER = 2j * mp.pi / L


def sexp_reihe(u):
    """(sexp(u) − 1, sexp'(u)) für u in [−0.55, 0.55] (Taylor um −0.4…0.4)."""
    c = min(_ZEN, key=lambda c: abs(u - c))
    co = _ZEN[c]
    v = u - c
    if c == 0:
        s1 = mp.polyval(co[:0:-1], v) * v
    else:
        s1 = mp.polyval(co[::-1], v) - 1
    d = mp.polyval([k * co[k] for k in range(len(co) - 1, 0, -1)], v)
    return s1, d


def knoten(k, tmax=5):
    h = mp.mpf(2) ** (-k)
    n = int(tmax / h)
    pts = []
    for j in range(-n, n + 1):
        t = j * h
        a = mp.pi * mp.sinh(t)
        dl = 1 / (1 + mp.exp(-a))
        dr = 1 / (1 + mp.exp(a))
        pts.append((dl, dr, h * mp.pi * mp.cosh(t) * dl * dr))
    return pts


_KCACHE = {}


def knoten_c(k):
    key = (k, mp.mp.dps)
    if key not in _KCACHE:
        _KCACHE[key] = knoten(k)
    return _KCACHE[key]


def entwirren(werte, per):
    out = [werte[0]]
    for v in werte[1:]:
        r = (out[-1] - v) / per
        mr = int(mp.nint(mp.re(r)))
        cand = [mr + m for m in range(-2, 3)] + list(range(-3, 4))
        best = min(cand, key=lambda m: abs(v + m * per - out[-1]))
        out.append(v + best * per)
    return out


# ---------------- Störungen
def p_wert(stoer, x):
    """(p(x), p'(x)) für x in (−1, 1)."""
    art, par = stoer
    if art == "trig":
        p = mp.mpf(0); dp = mp.mpf(0)
        for k, a, b in par:
            w = 2 * mp.pi * k
            s, c = mp.sin(w * x), mp.cos(w * x)
            p += a * s + b * (c - 1)
            dp += w * (a * c - b * s)
        return p, dp
    if art == "moeb":
        a = mp.mpc(par)
        q = mp.expjpi(2 * x)
        ab = mp.conj(a)
        # φ(x) − x = −(1/π)[arg(1 − ā q) − arg(1 − ā)]
        p = -(mp.arg(1 - ab * q) - mp.arg(1 - ab)) / mp.pi
        dp = (1 - abs(a) ** 2) / abs(1 - ab * q) ** 2 - 1
        return p, dp
    raise ValueError(art)


def auswerten(stoer, k=5):
    """Knoten X (in (0,1)), Gewichte W, θ_T (stetiger Zweig), lg = log(1+θ_T') (stetiger Zweig)."""
    pts = knoten_c(k)
    th, lg = [], []
    for dl, dr, w in pts:
        if stoer[0] == "lin":
            # T(x) = e^x auf (0,1): θ = h(e^x) − x = 1 + h(x) − x,  1+θ' = h'(x)
            if dl <= mp.mpf("0.5"):
                l1 = mp.log(dl)
                hv, hd = B.h_mit_ableitung(l1)          # h(x) = 1 + h(log x), h'(x) = h'(log x)/x
                hx, hdx = 1 + hv, hd / dl
            else:
                l1 = mp.log1p(-dr)
                hv, hd = B.h_mit_ableitung(l1)
                hx, hdx = 1 + hv, hd / (1 - dr)
            th.append(1 + hx - dl)
            lg.append(mp.log(hdx))
            continue
        xs = dl if dl <= mp.mpf("0.5") else -dr
        p, dp = p_wert(stoer, xs)
        u = xs + p
        s1, ds = sexp_reihe(u)
        l1 = mp.log1p(s1)
        hv, hd = B.h_mit_ableitung(l1)
        th.append(1 + hv - xs)
        lg.append(mp.log(ds * (1 + dp) * hd / (1 + s1)))
    mid = len(pts) // 2
    th = entwirren(th[mid::-1], PER)[::-1][:-1] + entwirren(th[mid:], PER)
    lg = entwirren(lg[mid::-1], 2j * mp.pi)[::-1][:-1] + entwirren(lg[mid:], 2j * mp.pi)
    X = [p[0] for p in pts]
    W = [p[2] for p in pts]
    return X, W, th, lg


def kenngroessen(stoer, k=5, moden=6, momente=6):
    X, W, th, lg = auswerten(stoer, k)
    M = lambda vals: mp.fsum(w * v for w, v in zip(W, vals))  # noqa: E731
    mt = M(th)
    r = {"Mtheta": mt}
    r["a"] = M([t ** 2 for t in th]) - mt ** 2                                    # M[θ²] − M[θ]²
    r["b"] = M([mp.exp(L * t) for t in th]) * M([mp.exp(-L * t) for t in th]) - 1  # σ-Produkt − 1
    r["c"] = M([mp.exp(-v) for v in lg]) - 1                                       # ∫ slog'²/h' − 1
    r["d"] = M(lg)                                                                 # M[log(T' h'(T))]
    eks = [mp.expjpi(2 * x) for x in X]
    for m in range(1, moden + 1):
        r[f"neg{m}"] = M([t * e ** m for t, e in zip(th, eks)])   # Koeffizient von e^{−2πimx}
    for j in range(2, momente + 1):
        r[f"V{j}"] = M([(t - mt) ** j for t in th])               # zentrale Momente (Kneser: 0)
    return r
