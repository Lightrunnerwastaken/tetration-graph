"""Frage 2: Welle zwischen regulaerer Pentation P (Motor, am reellen Fixpunkt p) und Kneser-artiger Pentation K5
(Cauchy-Loeser e6, doppelt, Fixpunkt L4).  psi(z) = P^{-1}(K5(z)) - z ist 1-periodisch, reell auf R.
Fourier-Koeffizienten c_n (n in Z) aus Linien y = 0, 0.8, 1.2, 1.5, 1.7: c_{-n} e^{2 pi n y} und c_n e^{-2 pi n y}.
Aufruf: python s5_psi5.py"""
import sys, json, time, warnings
import numpy as np, mpmath as mp
from pathlib import Path
warnings.filterwarnings("ignore")
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[1] / "eta_hyperoperationen" / "halb_rang_basis2_2026_10_13"))
import e6_lade
import motor
R = e6_lade.lade_raenge(5); K5 = R[5]
mp.mp.dps = 30
per = complex(motor.roh("2", ["pentx(0.5)", "pperiod"], 30)[1])
print("Periode der regulaeren Pentation:", per)
erg = {}
N = 64
for y in [0.0, 0.8, 1.2, 1.5, 1.7, 1.8]:
    x = np.arange(N) / N
    w = K5(x + 1j * y)
    v = motor.roh("2", [f"iferr(ipentx({motor.gpzahl(complex(wi), 20)}),E,1000)" for wi in w], 30)
    psi = np.array([complex(t) for t in v]) - (x + 1j * y)
    if np.any(np.abs(psi) > 500):
        print(y, "ipent scheitert an", int(np.sum(np.abs(psi) > 500)), "Punkten"); continue
    for j in range(1, N):
        ks = np.arange(-3, 4); psi[j] += per * ks[np.argmin(np.abs(psi[j] + ks * per - psi[j - 1]))]
    F = np.fft.fft(psi) / N
    n = np.arange(1, 6)
    cpos = F[1:6] * np.exp(2 * np.pi * n * y)        # c_n
    cneg = F[-1:-6:-1] * np.exp(-2 * np.pi * n * y)   # c_{-n}
    erg[y] = {"c0": F[0] - 0, "cpos": cpos, "cneg": cneg, "sprung": abs(psi[0] + per * 0 - psi[-1])}
    print(f"y={y}: c0={F[0]:.4e}  |c_-n| n=1..5: " + " ".join("%.4e" % abs(t) for t in cneg) +
          "   |c_n|: " + " ".join("%.2e" % abs(t) for t in cpos), flush=True)
out = {str(k): {"c0": [v["c0"].real, v["c0"].imag], "cneg": [[t.real, t.imag] for t in v["cneg"]], "cpos": [[t.real, t.imag] for t in v["cpos"]]} for k, v in erg.items()}
json.dump(out, open(HIER / "daten" / "psi5.json", "w"), indent=1)
