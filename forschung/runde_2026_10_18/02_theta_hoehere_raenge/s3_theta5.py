"""Rang 5: theta5 der Kneser-artigen Pentation (Cauchy-Loeser aus eta_hyperoperationen/.../e6, doppelte Genauigkeit, Basis 2)
bezogen auf die regulaere Loesung am ANZIEHENDEN Fixpunkt L4 von sexp (Multiplikator mu, |mu| = 0.71):
sigma(w) = mu^-n sigma_lok(sexp^n w) (Koenigs), h = log sigma / log mu, theta5(z) = h(H5(z)) - z.
FFT auf mehreren Hoehen y; gleiche a_n auf allen Hoehen <=> theta5 dazwischen holomorph.
Aufruf: python s3_theta5.py"""
import sys, json, time, warnings
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
# Taylor von H4 um L (Cauchy, Radius 0.4)
M, r = 128, 0.4
ph = np.exp(2j * np.pi * np.arange(M) / M)
c = np.fft.fft(H4(L + r * ph)) / M / r ** np.arange(M)
K = 30
f = c[:K + 1].copy(); f[0] = 0
mu = f[1]
# Schroeder-Reihe
fp = [None, f.copy()]
for k in range(2, K + 1): fp.append(np.convolve(fp[-1], f)[:K + 1])
s = np.zeros(K + 1, complex); s[1] = 1
for m in range(2, K + 1):
    s[m] = sum(s[k] * fp[k][m] for k in range(1, m)) / (mu - mu ** m)
logmu = np.log(mu)
print("L4 =", L, " mu =", mu, " |s_k| r^k bei r=0.05:", abs(s[K]) * 0.05 ** K)
def h(w, rad=0.05, nmax=300):
    w = np.array(w, complex); n = np.zeros(len(w), int); ok = np.ones(len(w), bool)
    for _ in range(nmax):
        m = ok & (np.abs(w - L) > rad)
        if not m.any(): break
        m2 = m & (np.abs(w) < 40)
        ok[m & ~m2] = False
        w[m2] = H4(w[m2]); n[m2] += 1
    ok &= np.abs(w - L) <= rad
    sig = np.polyval(s[::-1], w - L)
    return np.log(sig) / logmu - n, ok
per = 2j * np.pi / logmu
erg = {}
for y in [1.5, 1.0, 0.7, 0.5, 0.35]:
    N = 128
    x = np.arange(N) / N
    th, ok = h(H5(x + 1j * y)); th -= x + 1j * y
    if not ok.all():
        print(y, "nicht konvergiert:", int((~ok).sum())); continue
    for j in range(1, N):
        m = np.round(((th[j - 1] - th[j]) / per).real) if False else 0
        ks = np.arange(-5, 6); th[j] += per * ks[np.argmin(np.abs(th[j] + ks * per - th[j - 1]))]
    F = np.fft.fft(th) / N
    a = F[:N // 2] * np.exp(2 * np.pi * np.arange(N // 2) * y)
    aneg = F[-1:-9:-1] * np.exp(-2 * np.pi * np.arange(1, 9) * y)
    erg[y] = a
    print(f"y={y}: a0={a[0]:.10f}  |a_-1..-3|={np.abs(aneg[:3])}")
    print("    |a_n| n=1..10:", " ".join("%.3e" % v for v in np.abs(a[1:11])))
ys = sorted(erg)
print("\nVergleich der Hoehen (relativer Unterschied zu y=%.2f):" % ys[-1])
for yy in ys[:-1]:
    print("  y=%.2f:" % yy, " ".join("%.1e" % v for v in np.abs(erg[yy][1:11] - erg[ys[-1]][1:11]) / np.abs(erg[ys[-1]][1:11])))
a = erg[ys[0]]
print("\nlog|a_n| Schritte (exponentiell <=> konstant): ", " ".join("%.3f" % v for v in np.diff(np.log(np.abs(a[1:12])))))
print("arg a_n:", " ".join("%.3f" % v for v in np.angle(a[1:12])))
json.dump({str(k): [[float(v.real), float(v.imag)] for v in erg[k][:20]] for k in erg}, open(HIER / "daten" / "theta5_e6.json", "w"), indent=0)
