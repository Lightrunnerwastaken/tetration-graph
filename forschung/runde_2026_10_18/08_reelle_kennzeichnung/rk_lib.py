"""Runde 18.10.2026, Agent 08: Taylorreihen der Kneser-Tetration (Basis e) und gestörter Tetrationen an reellen Stellen.

Quelle: Taylor-Daten von sexp_K um c = -0.4, -0.2, 0, 0.2, 0.4 (128 Koeffizienten, 130 Stellen,
forschung/identitaeten_2026_10_14/etappe3/daten/taylor_130_256.json, nur gelesen). Kein Motor, kein Zwischenspeicher.

Reihen sind Listen [a_0, a_1, ..., a_{N-1}] (Koeffizienten von u^k), abgeschnitten bei N.
"""
import json
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).resolve().parent
QUELLE = HIER.parent.parent / "identitaeten_2026_10_14" / "etappe3" / "daten" / "taylor_130_256.json"
_ROH = json.loads(QUELLE.read_text(encoding="utf-8"))["zentren"]
_Z = None


def init(dps=60):
    global _Z
    mp.mp.dps = dps
    _Z = {mp.mpf(k): [mp.mpf(x) for x in v] for k, v in _ROH.items()}


# ---------------- Reihenarithmetik
def mul(a, b, N):
    c = [mp.mpf(0)] * N
    for i, ai in enumerate(a[:N]):
        if ai == 0:
            continue
        for j in range(min(len(b), N - i)):
            c[i + j] += ai * b[j]
    return c


def compose(c, v, N):
    """Σ c_k v^k mit v[0] = 0 (Horner)."""
    assert v[0] == 0
    r = [mp.mpf(0)] * N
    for ck in reversed(c[:N]):
        r = mul(r, v, N)
        r[0] += ck
    return r


def exp_s(a, N):
    """exp einer Reihe: e = exp(a), e' = a' e."""
    e = [mp.mpf(0)] * N
    e[0] = mp.exp(a[0])
    da = [k * a[k] for k in range(1, N)] + [mp.mpf(0)]
    for m in range(1, N):
        e[m] = mp.fsum(da[j - 1] * e[m - j] for j in range(1, m + 1)) / m
    return e


def log_s(a, N):
    """log einer Reihe mit a[0] > 0: l' = a'/a."""
    l = [mp.mpf(0)] * N
    l[0] = mp.log(a[0])
    for m in range(1, N):
        # m a_m = Σ_{j=1}^{m} j l_j a_{m-j}
        s = m * a[m] - mp.fsum(j * l[j] * a[m - j] for j in range(1, m))
        l[m] = s / (m * a[0])
    return l


def shift(c, h, N):
    """Koeffizienten um c0+h aus Koeffizienten um c0 (Taylor-Verschiebung, Horner-Schema)."""
    b = list(c)
    M = len(b)
    for i in range(min(N, M)):
        for j in range(M - 2, i - 1, -1):
            b[j] += h * b[j + 1]
    return b[:N]


# ---------------- Kneser
def kneser(x, N, erlaube_exp=False):
    """Taylor-Koeffizienten von sexp_K um die reelle Stelle x (x > -2), N Stück."""
    x = mp.mpf(x)
    if x > mp.mpf("0.5"):
        if not erlaube_exp:
            raise ValueError("kneser: x > 0.5 nur mit erlaube_exp (hohe Ordnungen unzuverlässig, siehe BERICHT)")
        return exp_s(kneser(x - 1, N, erlaube_exp), N)
    if x < mp.mpf("-0.5"):
        return log_s(kneser(x + 1, N, erlaube_exp), N)
    c0 = min(_Z, key=lambda c: abs(x - c))
    return shift(_Z[c0], x - c0, N)


# ---------------- Störungen φ(x) = x + p(x)
def phi_reihe(stoer, x, N):
    """Taylor-Koeffizienten von φ(x+u) = x + u + p(x+u)."""
    art, par = stoer
    x = mp.mpf(x)
    out = [mp.mpf(0)] * N
    out[0] = x
    if N > 1:
        out[1] = mp.mpf(1)
    if art == "kneser":
        return out
    if art == "trig":
        # p = Σ a sin(2πk x) + b (cos(2πk x) - 1)
        for k, a, b in par:
            w = 2 * mp.pi * k
            s, c = mp.sin(w * x), mp.cos(w * x)
            for j in range(N):
                # d^j/du^j sin(w(x+u)) = w^j sin(wx + jπ/2)
                fj = w ** j / mp.factorial(j)
                sj = s if j % 4 == 0 else (c if j % 4 == 1 else (-s if j % 4 == 2 else -c))
                cj = c if j % 4 == 0 else (-s if j % 4 == 1 else (-c if j % 4 == 2 else s))
                out[j] += fj * (a * sj + b * cj)
            out[0] -= b
        return out
    if art == "moeb":
        # φ − x = −(1/π)[arg(1 − ā q) − arg(1 − ā)],  q = e^{2πi(x+u)}
        # arg(1 − ā q) = Im log(1 − ā q),  log(1 − ā q) = −Σ_{n≥1} (ā q)^n / n  (|a| < 1)
        a = mp.mpc(par)
        ab = mp.conj(a)
        q0 = mp.expjpi(2 * x)
        # Reihe von log(1 - ā q0 e^{2πiu}) in u:  -Σ_n (ā q0)^n/n · e^{2πinu}, Koeffizient von u^j: -Σ_n (ā q0)^n/n (2πin)^j/j!
        la = mp.log(abs(a))
        def lt(n):  # log der Termgrösse für die höchste Ordnung
            return n * la + (N - 1) * mp.log(2 * mp.pi * n) - mp.loggamma(N) - mp.log(n)
        peak = max(1, int((N - 1) / max(-la, mp.mpf('1e-3'))))
        top = max(lt(n) for n in range(1, peak + 2))
        nmax = peak + 1
        while lt(nmax) > top - (mp.mp.dps + 10) * mp.log(10) and nmax < 200000:
            nmax += 1 + nmax // 20
        base = [(ab * q0) ** n / n for n in range(1, nmax + 1)]
        for j in range(N):
            s = -mp.fsum(base[n - 1] * (2j * mp.pi * n) ** j for n in range(1, nmax + 1)) / mp.factorial(j)
            out[j] += -mp.im(s) / mp.pi
        out[0] += mp.im(mp.log(1 - ab)) / mp.pi
        return out
    raise ValueError(art)


def tet(stoer, x, N):
    """Taylor-Koeffizienten von T = sexp_K ∘ φ um x."""
    ph = phi_reihe(stoer, x, N)
    K = kneser(ph[0], N)
    v = [mp.mpf(0)] + ph[1:]
    return compose(K, v, N)
