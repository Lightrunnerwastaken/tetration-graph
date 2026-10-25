"""Gebiet 21: schnelle (Gleitkomma-)Auswertung von Knesers sexp zur Basis e, ohne Motor.

Quellen (nur lesen):
  Taylor um 0, 128 Koeffizienten, 130 Stellen: identitaeten_2026_10_14/etappe3/daten/taylor_130_256.json
  theta-Fourier-Koeffizienten a_n, n <= 16384, 44+ Stellen: identitaeten_2026_10_14/etappe2/daten/a_e_80_0.003_16384.json
  Psi (Umkehrung der Schroeder-Funktion an L): Reihe aus lib3.Basis (K=60), hier in complex128.

sexp(z):
  |z| <= R_T            : Taylor um 0
  |Im z| >= Y_TH        : F(z + theta(z)),  F(s) = Psi(e^{L s})  (untere Halbebene per Spiegelung)
  sonst                 : sexp(z) = exp(sexp(z-1)) (Re z > 0)  bzw.  Log(sexp(z+1)) (Re z < 0)
"""
import json
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

HIER = Path(__file__).resolve().parent
ID = HIER.parent.parent / "identitaeten_2026_10_14"
sys.path.insert(0, str(ID / "etappe3"))

R_T = 1.6
Y_TH = 0.2

_d = json.loads((ID / "etappe3" / "daten" / "taylor_130_256.json").read_text(encoding="utf-8"))
TAYLOR = np.array([float(mp.mpf(x)) for x in _d["zentren"]["0.0"]])  # c_0 .. c_127
_t = json.loads((ID / "etappe2" / "daten" / "a_e_80_0.003_16384.json").read_text(encoding="utf-8"))
THETA = np.array([complex(float(mp.mpf(a)), float(mp.mpf(b))) for a, b in _t["a"]])

mp.mp.dps = 30
import lib3  # noqa: E402

_B = lib3.Basis(K=60)
L = complex(_B.L)
PSI = np.array([complex(p) for p in _B.P])  # Psi(u) = sum PSI[k] u^k, |u| klein
SIG = np.array([complex(s) for s in _B.S])  # sigma(L+u) = sum SIG[k] u^k
P_PER = 2j * np.pi / L  # Periode von F


def psi(u):
    u = np.asarray(u, dtype=complex).copy()
    n = np.zeros(u.shape, dtype=int)
    while True:
        m = np.abs(u) > 0.02
        if not m.any():
            break
        u[m] /= L
        n[m] += 1
    w = np.polyval(PSI[::-1], u)
    for k in range(1, n.max() + 1 if n.size else 1):
        m = n >= k
        w[m] = np.exp(w[m])
    return w


def theta(z, nmax=None):
    z = np.asarray(z, dtype=complex)
    q = np.exp(2j * np.pi * z)
    a = THETA if nmax is None else THETA[: nmax + 1]
    # Horner, mit Abbruch: |q|^n < 1e-18
    qmax = np.abs(q).max()
    if qmax > 0:
        nn = int(min(len(a) - 1, max(10, np.ceil(-41.0 / np.log(qmax))))) if qmax < 1 else len(a) - 1
    else:
        nn = 10
    return np.polyval(a[: nn + 1][::-1], q)


def tau(z):
    return np.asarray(z, dtype=complex) + theta(z)


def F(s):
    return psi(np.exp(L * np.asarray(s, dtype=complex)))


def sexp(z, tiefe=0):
    z = np.asarray(z, dtype=complex)
    flach = z.ndim == 0
    z = np.atleast_1d(z).copy()
    out = np.empty(z.shape, dtype=complex)
    a = np.abs(z) <= R_T
    out[a] = np.polyval(TAYLOR[::-1], z[a])
    b = (~a) & (z.imag >= Y_TH)
    if b.any():
        out[b] = F(tau(z[b]))
    c = (~a) & (z.imag <= -Y_TH)
    if c.any():
        out[c] = np.conj(F(tau(np.conj(z[c]))))
    r = (~a) & (np.abs(z.imag) < Y_TH)
    if r.any():
        if tiefe > 60:
            raise RuntimeError("sexp: Rekursion zu tief")
        zr = z[r]
        rechts = zr.real > 0
        tmp = np.empty(zr.shape, dtype=complex)
        if rechts.any():
            tmp[rechts] = np.exp(sexp(zr[rechts] - 1, tiefe + 1))
        if (~rechts).any():
            tmp[~rechts] = np.log(sexp(zr[~rechts] + 1, tiefe + 1))
        out[r] = tmp
    return out[0] if flach else out


def dsexp(z, h=1e-5):
    z = np.asarray(z, dtype=complex)
    return (sexp(z + h) - sexp(z - h)) / (2 * h)


# ---------- Geometrie: Doppelpunkte geschlossener / offener Polygone
def doppelpunkte(w, geschlossen=True, max_treffer=5):
    """Liefert Paare (i, j) sich schneidender, nicht benachbarter Strecken des Polygons w."""
    w = np.asarray(w, dtype=complex)
    if geschlossen:
        a = w
        b = np.roll(w, -1)
    else:
        a = w[:-1]
        b = w[1:]
    n = len(a)
    ax, ay, bx, by = a.real, a.imag, b.real, b.imag
    treffer = []
    for i in range(n):
        j = np.arange(i + 2, n)
        if geschlossen and i == 0:
            j = j[j != n - 1]
        if j.size == 0:
            continue
        # Bounding-Box-Vorfilter
        m = (np.maximum(ax[j], bx[j]) >= min(ax[i], bx[i])) & (np.minimum(ax[j], bx[j]) <= max(ax[i], bx[i])) & \
            (np.maximum(ay[j], by[j]) >= min(ay[i], by[i])) & (np.minimum(ay[j], by[j]) <= max(ay[i], by[i]))
        j = j[m]
        if j.size == 0:
            continue
        d1 = _kreuz(a[i], b[i], a[j])
        d2 = _kreuz(a[i], b[i], b[j])
        d3 = _kreuz(a[j], b[j], a[i])
        d4 = _kreuz(a[j], b[j], b[i])
        hit = (d1 * d2 < 0) & (d3 * d4 < 0)
        for jj in j[hit]:
            treffer.append((i, int(jj)))
            if len(treffer) >= max_treffer:
                return treffer
    return treffer


def _kreuz(p, q, r):
    return ((q - p).conjugate() * (r - p)).imag


def schnitt_zweier(w1, w2, max_treffer=3):
    """Schnittpunkte zweier offener Polygone."""
    a1, b1 = w1[:-1], w1[1:]
    treffer = []
    a2, b2 = w2[:-1], w2[1:]
    for i in range(len(a1)):
        m = (np.maximum(a2.real, b2.real) >= min(a1[i].real, b1[i].real)) & \
            (np.minimum(a2.real, b2.real) <= max(a1[i].real, b1[i].real)) & \
            (np.maximum(a2.imag, b2.imag) >= min(a1[i].imag, b1[i].imag)) & \
            (np.minimum(a2.imag, b2.imag) <= max(a1[i].imag, b1[i].imag))
        j = np.nonzero(m)[0]
        if j.size == 0:
            continue
        d1 = _kreuz(a1[i], b1[i], a2[j]); d2 = _kreuz(a1[i], b1[i], b2[j])
        d3 = _kreuz(a2[j], b2[j], a1[i]); d4 = _kreuz(a2[j], b2[j], b1[i])
        hit = (d1 * d2 < 0) & (d3 * d4 < 0)
        for jj in j[hit]:
            treffer.append((i, int(jj)))
            if len(treffer) >= max_treffer:
                return treffer
    return treffer
