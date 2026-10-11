"""Strang C, Runde 4, c26: Gegenprobe (b) bis N = 1000 aus Strang As Daten — sieht δ1 die Eckspitze?

Exakt (BERICHT §12.6): mit A_f = P_N J_∞(f) P_N ist δ1^f := A_f⁻¹B_N t_N = Π_N f_K − q_N, also
α_N^f = |Π_N f_K − q_N|_{ℓ¹_ρ} (c21: = α_N mit A = J_N(q_N) auf alle gedruckten Stellen). Π_N f_K liegt bei Strang A
(gemeinsam/lin_PiN_fK_bis_1000.json, 60 Stellen, Präzision 0.72N+200 Bit), die a_k in a_moden/daten/kneser_koeff_lang.json
(dieselben a_k, die A benutzt hat). Kein Motor, keine Matrix: nur Differenzen und FFT.
Ausgabe: α_N (s = 0 und s = 4), Fit log α = A + b log N + cN, Profil |δ1(ρe^{iθ})| (Lage des Maximums, Werte neben L,
bei −ρ), Anteil der Koeffizientenmasse in den Bändern k/N.
"""
import json
import os
import sys

import mpmath as mp
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
mp.mp.dps = 80
L = mp.mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)
ARGL = float(mp.arg(L))

PI = json.load(open(os.path.join(G, "lin_PiN_fK_bis_1000.json"), encoding="utf-8"))
AK = json.load(open(os.path.join(HIER, "..", "a_moden", "daten", "kneser_koeff_lang.json"), encoding="utf-8"))
NS = sorted(int(k) for k in PI if k != "beschreibung")

M = 1 << 15                     # FFT-Länge: Winkelauflösung 0.011°
erg = []
for N in NS:
    b = PI[str(N)]["b"]
    d = [(mp.mpf(b[k]) - mp.mpf(AK[k])) * RHO ** k for k in range(N + 1)]    # δ1_k ρ^k
    dl = np.array([float(x) for x in d])
    al0 = float(np.sum(np.abs(dl)))
    al4 = float(np.sum(np.abs(dl) * (1.0 + np.arange(N + 1)) ** 4))
    werte = np.fft.ifft(np.concatenate([dl, np.zeros(M - N - 1)])) * M   # δ1(ρ e^{iθ_m}), θ_m = 2πm/M
    betr = np.abs(werte)
    m_max = int(np.argmax(betr[: M // 2 + 1]))
    th_max = 360.0 * m_max / M

    def bei(th_grad):
        m = int(round(th_grad / 360.0 * M)) % M
        return float(betr[m])

    gL = ARGL * 180 / np.pi
    k = np.arange(N + 1) / N
    baender = [(0, .5), (.5, .8), (.8, .9), (.9, .95), (.95, 1.0001)]
    anteile = [float(np.sum(np.abs(dl)[(k >= lo) & (k < hi)])) / al0 for lo, hi in baender]
    erg.append((N, al0, al4, th_max, float(betr[m_max]), bei(gL - 1), bei(gL + 1), bei(gL + 5), bei(180.0), bei(0.0), anteile))

print("N | α_N (s=0) | α_N (s=4) | Max |δ1| auf |x|=ρ bei θ | Wert | θ=argL−1° | argL+1° | argL+5° | −ρ | +ρ | Bänder k/N "
      "[0,.5) [.5,.8) [.8,.9) [.9,.95) [.95,1]")
for N, a0, a4, th, mx, m1, p1, p5, w180, w0, an in erg:
    print(f"{N} | {a0:.4e} | {a4:.4e} | {th:.3f}° | {mx:.3e} | {m1:.2e} | {p1:.2e} | {p5:.2e} | {w180:.2e} | {w0:.2e} | "
          + " ".join(f"{x:.3f}" for x in an))

# lokale Exponenten und Fit mit Exponentialterm
Ns = np.array([e[0] for e in erg], dtype=float)
for name, idx in (("α_N (s=0)", 1), ("α_N (s=4)", 2), ("|δ1(−ρ)|", 8), ("|δ1| bei argL+5°", 7)):
    y = np.log(np.array([e[idx] for e in erg]))
    lok = [(Ns[i], (y[i + 1] - y[i]) / (np.log(Ns[i + 1]) - np.log(Ns[i]))) for i in range(len(Ns) - 1) if Ns[i + 1] - Ns[i] >= 25]
    print(f"{name}: lokale Exponenten " + " ".join(f"{int(n)}:{v:.2f}" for n, v in lok))
    for nmin in (300, 500):
        sel = Ns >= nmin
        X = np.vstack([np.ones(sel.sum()), np.log(Ns[sel]), Ns[sel]]).T
        koef, res, *_ = np.linalg.lstsq(X, y[sel], rcond=None)
        r = y[sel] - X @ koef
        s2 = float(r @ r) / max(1, sel.sum() - 3)
        cov = s2 * np.linalg.inv(X.T @ X)
        print(f"   Fit N ≥ {nmin}: log = {koef[0]:.3f} + ({koef[1]:.3f}) log N + ({koef[2]:.2e}) N   ± c: {np.sqrt(cov[2, 2]):.1e}"
              f"   (Eckrate D: 0.0120)")
