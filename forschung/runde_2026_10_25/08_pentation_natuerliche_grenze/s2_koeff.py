"""Schritt 2: Wachstum der theta5-Koeffizienten gegen den obersten Randpunkt (Cauchy-Hadamard), und Fortsetzung UNTER den Konvergenzkreis.
(a) a_n aus zwei Linien y1, y2 (FFT, N Punkte): rho_n = ln|a_n|/(2 pi n), Wurzeltest, Phasenschritt arg(a_n/a_{n-1}) -> -2 pi x_tip.
(b) Senkrechte Linie x = X von 0.45 bis Y herab: alle Punkte im Becken (Klasse 1)? Dann Mittelwerttest theta5(z1) = Mittel auf Kreis um z1
    (holomorph unterhalb y_max) und Fourier-Reihe dort divergent.
Aufruf: python s2_koeff.py y1 y2 N X Y"""
import sys, json, time
import numpy as np
from lib import H5, klasse, h, theta_linie, HIER, per
t0 = time.time()
y1, y2, N = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
X, Y = float(sys.argv[4]), float(sys.argv[5])
A = {}
for y in (y1, y2):
    th, ok = theta_linie(y, N)
    print("y=%.3f: alle Punkte konvergiert: %s" % (y, ok.all()))
    A[y] = np.fft.fft(th)[:N // 2] / N * np.exp(2 * np.pi * np.arange(N // 2) * y)
    neg = np.fft.fft(th)[-1:-6:-1] / N
    print("   |negative Moden 1..5| (roh) =", " ".join("%.1e" % v for v in np.abs(neg)))
a, b = A[y1], A[y2]
print(" n    |a_n|(y1)     rel.Diff(y2)  ln|a_n|/(2pi n)  (ln|a_n|-ln|a_{n-8}|)/(16pi)   -arg(a_n/a_{n-1})/2pi mod 1")
for n in list(range(2, 21, 2)) + list(range(24, N // 4, 8)):
    rd = abs(a[n] - b[n]) / abs(a[n])
    loc = (np.log(abs(a[n])) - np.log(abs(a[n - 8]))) / (16 * np.pi) if n > 8 else np.nan
    ph = (-np.angle(a[n] / a[n - 1]) / (2 * np.pi)) % 1
    print("%3d  %.4e   %.1e      %.4f          %.4f                        %.4f" % (n, abs(a[n]), rd, np.log(abs(a[n])) / (2 * np.pi * n), loc, ph))
json.dump({str(y): [[float(v.real), float(v.imag)] for v in A[y]] for y in A}, open(HIER / "cache" / ("s2_a_%g_%g_%d.json" % (y1, y2, N)), "w"))
# (b)
ys = np.linspace(0.45, Y - 0.04, 200)
kl, n = klasse(H5(X + 1j * ys), nmax=300)
print("\nsenkrechte Linie x=%.3f, y von 0.45 bis %.3f: Klassen %s" % (X, Y - 0.04, np.bincount(kl, minlength=4).tolist()))
z1 = X + 1j * Y
for rho in (0.02, 0.03):
    M = 64
    zc = z1 + rho * np.exp(2j * np.pi * np.arange(M) / M)
    v, ok = h(H5(np.concatenate([[z1], zc])))
    v = v - np.concatenate([[z1], zc])
    vc = v[1:].copy(); ks = np.arange(-5, 6)
    for j in range(1, M): vc[j] += per * ks[np.argmin(np.abs(vc[j] + ks * per - vc[j - 1]))]
    v0 = v[0] + per * ks[np.argmin(np.abs(v[0] + ks * per - vc[0]))]
    print("Kreis rho=%.2f um z1=%s: alle ok %s, |theta5(z1) - Kreismittel| = %.2e, Umlauf %.2e, theta5(z1) = %s" %
          (rho, z1, ok.all(), abs(v0 - vc.mean()), abs(vc[0] + (vc[-1] - vc[-2]) - vc[0] - (vc[0] - vc[-1]) + (vc[-1] - vc[0])) if False else abs(vc[-1] - vc[0]), v0))
q = np.exp(2j * np.pi * z1)
terme = np.abs(a[:N // 4] * q ** np.arange(N // 4))
print("Fourier-Terme |a_n q^n| bei z1 (n=10,20,40,60):", " ".join("%.2e" % terme[k] for k in (10, 20, 40, 60) if k < len(terme)))
print("Zeit %.1f s" % (time.time() - t0))
