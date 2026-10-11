"""Strang A, Runde 2: Fehler E_N = P_N − f_K auf dem Rand |z| = |L| (P_N = Π_N f_K).

Koeffizienten von E_N: (b_k − a_k) für k ≤ N, −a_k für N < k ≤ K (K = max(4N, 4000), a_k ab 701 aus der Darboux-Reihe);
Werte per FFT (Differenzen sind klein ⇒ double genügt) plus feines Gitter um arg L.
Ausgaben je N: max|E| auf dem Rand und Ort (θ − arg L)·N, gewichtete Maxima max|E|·(|z−L||z−L̄|)^{−a}·… (siehe unten),
max|P_N| auf dem Rand in der Nähe von L (|θ − arg L| ≤ 0.3) und global.
"""
import json
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

from lin_kern import DATEN, kneser_koeff, lade  # noqa: E402

mp.mp.dps = 60
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
argL = float(mp.arg(L))
Lc = complex(L)
Ldat = lade()
NMAX = max(Ldat)
K = max(4 * NMAX, 4000)
t0 = time.time()
a = kneser_koeff(K)
aL = np.array([float(x * RL ** k) for k, x in enumerate(a)])
print(f"Koeffizienten bis {K} ({time.time() - t0:.0f} s)", flush=True)
pfad = DATEN / "lin_kreisfehler.json"
aus = json.loads(pfad.read_text(encoding="utf-8")) if pfad.is_file() else {}
for N in sorted(Ldat):
    if str(N) in aus:
        continue
    b = [mp.mpf(x) for x in Ldat[N]["b"]]
    KN = max(4 * N, 4000)
    e = np.zeros(KN + 1)
    e[:N + 1] = [float((b[k] - a[k]) * RL ** k) for k in range(N + 1)]
    e[N + 1:] = -aL[N + 1:KN + 1]
    Mg = 1 << int(np.ceil(np.log2(4 * (KN + 1))))
    th1 = 2 * np.pi * np.arange(Mg) / Mg
    E1 = np.fft.ifft(np.concatenate([e, np.zeros(Mg - KN - 1)])) * Mg
    th2 = argL + np.linspace(-40.0 / N, 40.0 / N, 4001)
    E2 = np.polyval(e[::-1], np.exp(1j * th2))
    th = np.concatenate([th1, th2])
    E = np.concatenate([E1, E2])
    z = float(RL) * np.exp(1j * th)
    w = np.abs(z - Lc) * np.abs(z - np.conj(Lc))
    i = int(np.argmax(np.abs(E)))
    tt = (th[i] + np.pi) % (2 * np.pi) - np.pi
    res = {"maxE": float(abs(E[i])), "ort_N": float((abs(tt) - argL) * N)}
    # P_N selbst: |P| = |f + E|; f auf dem Rand aus den Kneser-Koeffizienten (FFT)
    Fv = np.fft.ifft(np.concatenate([aL[:KN + 1], np.zeros(Mg - KN - 1)])) * Mg
    Pv = Fv + E1
    nahe = np.abs(((th1 + np.pi) % (2 * np.pi) - np.pi) - argL) <= 0.3
    res["maxP_global"] = float(np.max(np.abs(Pv)))
    res["maxP_nahL"] = float(np.max(np.abs(Pv[nahe])))
    for a_ in (1, 2, 3, 3.4):
        y = np.abs(E) * w ** (-a_ / 10)          # nur zur Lage; Hauptgrösse unten
    # gewichteter Fehler: max |E|·|z − L|^{−s} für s = 0.5, 1, 2 (wie stark E an L verschwindet)
    for s in (1, 2, 3):
        y = np.abs(E2) / np.abs(float(RL) * np.exp(1j * th2) - Lc) ** (-s)
        j = int(np.argmax(y))
        res[f"maxE_mal_abst^{s}"] = [float(y[j]), float((th2[j] - argL) * N)]
    aus[str(N)] = res
    pfad.write_text(json.dumps(aus, indent=0), encoding="utf-8")
    print(f"N = {N}: max_Rand|P_N − f_K| = {res['maxE']:.4g} bei (θ − arg L)·N = {res['ort_N']:.3g};  max|P_N| global "
          f"{res['maxP_global']:.6f}, nahe L {res['maxP_nahL']:.6f}  ({time.time() - t0:.0f} s)", flush=True)
