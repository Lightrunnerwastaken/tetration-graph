"""Frage 2, Fortsetzung: c_{-n} von psi bis n = 14 aus den Linien y = 1.5 und 1.6 (N = 128) und Vergleich mit dem
Verzweigungspunkt z_L3 (K5(z_L3) = L3, Fixpunkt von 2^z): Erwartung c_{-n} ~ C n^beta e^{2 pi i n z_L3}.
Aufruf: python s5b_psi5_hoch.py"""
import sys, json, warnings
import numpy as np
from pathlib import Path
warnings.filterwarnings("ignore")
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[1] / "eta_hyperoperationen" / "halb_rang_basis2_2026_10_13"))
import e6_lade, motor
R = e6_lade.lade_raenge(5); K5 = R[5]
per = complex(motor.roh("2", ["pentx(0.5)", "pperiod"], 30)[1])
L3 = 0.8246785461420744 + 1.5674321238496478j
z = -0.6181507 + 1.6902j
for _ in range(30):
    v = K5(np.array([z]))[0] - L3; dv = (K5(np.array([z + 1e-6]))[0] - K5(np.array([z - 1e-6]))[0]) / 2e-6; z -= v / dv
zL3 = z
print("z_L3 =", zL3, " |K5(z_L3)-L3| =", abs(K5(np.array([zL3]))[0] - L3))
N = 128; res = {}
for y in [1.5, 1.6]:
    x = np.arange(N) / N
    w = K5(x + 1j * y)
    psi = np.array([complex(t) for t in motor.roh("2", [f"ipentx({motor.gpzahl(complex(wi), 20)})" for wi in w], 30)]) - (x + 1j * y)
    for j in range(1, N):
        ks = np.arange(-3, 4); psi[j] += per * ks[np.argmin(np.abs(psi[j] + ks * per - psi[j - 1]))]
    F = np.fft.fft(psi) / N
    n = np.arange(1, 17)
    res[y] = F[-n] * np.exp(-2 * np.pi * n * y)
a, b = res[1.5], res[1.6]
print(" n   |c_-n|        rel.Unterschied 1.5/1.6   |c_-n| e^{2 pi n d}*n   Phase(c_-n e^{-2 pi i n x0})   ")
d, x0 = zL3.imag, zL3.real % 1
zeilen = []
for k in range(14):
    nn = k + 1
    g = b[k] * np.exp(-2j * np.pi * nn * zL3)          # c_{-n} e^{-2 pi i n z_L3}
    print("%2d  %.6e   %.1e    %.6f   %+.5f" % (nn, abs(b[k]), abs(a[k] - b[k]) / abs(b[k]), abs(g) * nn, np.angle(g)))
    zeilen.append([nn, abs(b[k]), abs(a[k] - b[k]) / abs(b[k]), abs(g) * nn, float(np.angle(g))])
r = b[1:14] / b[0:13]
print("Verhaeltnisse -> e^{2 pi i z_L3} =", np.exp(2j * np.pi * zL3))
print("   ", " ".join("%.4e%+.4ej" % (t.real, t.imag) for t in r[-4:]))
json.dump({"z_L3": [zL3.real, zL3.imag], "zeilen": zeilen}, open(HIER / "daten" / "psi5_hoch.json", "w"), indent=1)
