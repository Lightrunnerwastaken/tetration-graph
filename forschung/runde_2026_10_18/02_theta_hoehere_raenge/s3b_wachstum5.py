"""Rang 5: Wachstum der a_n von theta5 (Koenigs an L4) bis n = 24 aus den Linien y = 0.30 und 0.35 (N = 256):
|a_n| ~ e^{2 pi n y_c}, y_c = Hoehe des obersten Punktes der Beckengrenze (s4_becken5: ~0.26-0.27).
Aufruf: python s3b_wachstum5.py"""
import runpy, sys, json, numpy as np
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.argv = ["x"]
g = runpy.run_path(str(HIER / "s3_theta5.py"), run_name="lib")   # laedt H4, H5, h (und rechnet die Grundtabelle)
h, H5 = g["h"], g["H5"]; per = g["per"]
res = {}
for y in [0.30, 0.35, 0.45]:
    N = 256; x = np.arange(N) / N
    th, ok = h(H5(x + 1j * y)); th -= x + 1j * y
    if not ok.all(): print(y, "nicht konvergiert:", int((~ok).sum())); continue
    for j in range(1, N):
        ks = np.arange(-5, 6); th[j] += per * ks[np.argmin(np.abs(th[j] + ks * per - th[j - 1]))]
    a = np.fft.fft(th)[:30] / N * np.exp(2 * np.pi * np.arange(30) * y)
    res[y] = a
b = res[0.35]
print(" n   |a_n|(y=.35)   rel.Unterschied zu y=.30   log|a_n|/(2 pi n)   (log|a_n|-log|a_n-4|)/(8 pi)")
for n in range(4, 25, 2):
    print("%3d  %.4e   %.1e   %.4f   %.4f" % (n, abs(b[n]), abs(b[n] - res[0.30][n]) / abs(b[n]), np.log(abs(b[n])) / (2 * np.pi * n),
          (np.log(abs(b[n])) - np.log(abs(b[n - 4]))) / (8 * np.pi)))
json.dump({str(k): [[float(v.real), float(v.imag)] for v in res[k]] for k in res}, open(HIER / "daten" / "theta5_wachstum.json", "w"), indent=0)
