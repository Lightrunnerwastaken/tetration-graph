"""Gemeinsame Bausteine (nur Lesen der e6-Daten aus eta_hyperoperationen/halb_rang_basis2_2026_10_13, doppelte Genauigkeit, Basis 2).
H4 = Kneser-sexp_2 (Cauchy-Loeser e6), H5 = Kneser-artige Pentation (e6), L = anziehender Fixpunkt L4 von H4,
Koenigs-Reihe s (sigma_lok) an L, h(w) = log sigma(w)/log mu (Abel-Funktion), theta5(z) = h(H5 z) - z.
klasse(w): Becken-Klassifikation unter Vorwaertsiteration von H4."""
import sys, warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
HIER = Path(__file__).resolve().parent
E6 = HIER.parents[1] / "eta_hyperoperationen" / "halb_rang_basis2_2026_10_13"
sys.path.insert(0, str(E6))
import e6_lade
from e6_superfunktion import fixpunkt
R = e6_lade.lade_raenge(5); H4, H5 = R[4], R[5]
L = fixpunkt(H4, 1.47 + 1.66j)
L3 = R[4].L          # Fixpunkt von 2^z, Grenzwert von sexp fuer Im -> +inf
M_, r_ = 128, 0.4
_ph = np.exp(2j * np.pi * np.arange(M_) / M_)
_c = np.fft.fft(H4(L + r_ * _ph)) / M_ / r_ ** np.arange(M_)
K = 30
_f = _c[:K + 1].copy(); _f[0] = 0
mu = _f[1]
_fp = [None, _f.copy()]
for k in range(2, K + 1): _fp.append(np.convolve(_fp[-1], _f)[:K + 1])
s = np.zeros(K + 1, complex); s[1] = 1
for m in range(2, K + 1):
    s[m] = sum(s[k] * _fp[k][m] for k in range(1, m)) / (mu - mu ** m)
logmu = np.log(mu)
per = 2j * np.pi / logmu


def klasse(w, nmax=400, eps=1e-6):
    """1: Becken von L4, 2: Becken von conj(L4), 3: Re -> +inf / nan (entkommt), 0: unentschieden. Dazu Schrittzahl."""
    w = np.array(w, complex).ravel().copy()
    kl = np.zeros(len(w), int); n = np.zeros(len(w), int)
    for k in range(nmax):
        m = kl == 0
        if not m.any(): break
        w[m] = H4(w[m]); n[m] += 1
        kl[m & (np.abs(w - L) < eps)] = 1
        kl[m & (np.abs(w - np.conj(L)) < eps)] = 2
        kl[m & (~np.isfinite(w) | (np.real(w) > 60))] = 3
    return kl, n


def h(w, rad=0.05, nmax=400):
    """Abel-Funktion an L4 (Koenigs): h(w) = log sigma(sexp^n w)/log mu - n; ok = Bahn erreicht die Scheibe."""
    w = np.array(w, complex).ravel().copy(); n = np.zeros(len(w), int); ok = np.ones(len(w), bool)
    for _ in range(nmax):
        m = ok & (np.abs(w - L) > rad)
        if not m.any(): break
        bad = m & (~np.isfinite(w) | (np.real(w) > 60))
        ok[bad] = False
        m2 = m & ~bad
        w[m2] = H4(w[m2]); n[m2] += 1
    ok &= np.abs(w - L) <= rad
    sig = np.polyval(s[::-1], w - L)
    return np.log(sig) / logmu - n, ok


def theta_linie(y, N):
    """theta5 auf der Linie x + iy, x = j/N, stetig gemacht (Periode per), dazu ok-Maske."""
    x = np.arange(N) / N
    th, ok = h(H5(x + 1j * y)); th = th - (x + 1j * y)
    ks = np.arange(-5, 6)
    for j in range(1, N):
        th[j] += per * ks[np.argmin(np.abs(th[j] + ks * per - th[j - 1]))]
    return th, ok
