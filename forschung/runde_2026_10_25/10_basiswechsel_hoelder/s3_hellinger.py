"""S3: Hellinger-Affinität der Ziffernverteilung je Ebene bei Auflösung 1/K (Satz B, Beweis über das Martingal Z_n).

Unter Lebesgue ist der Bruchteil t_n = {U_n} (in K Fächer eingeteilt) fast gleichverteilt: q_i = 1/K.
Unter dem Mass mu = dH ist er nach nu = d(t + P4(t)) verteilt: p_i = 1/K + P4((i+1)/K) - P4(i/K).
A_K = sum_i sqrt(p_i q_i) < 1 genau dann, wenn P4 nicht 1/K-periodisch ist; dann E_lambda sqrt(Z_n) <= prod (A_K + delta_j) -> 0.
Gegenprobe: Fourier-Daten mit M = 64 und M = 128 (Runde 18).
Aufruf: python s3_hellinger.py
"""
import json
from pathlib import Path

import numpy as np

ALT = Path(__file__).resolve().parents[2] / "runde_2026_10_18" / "01_basiswechsel_raenge" / "daten"


def lade(name):
    dat = json.loads((ALT / name).read_text(encoding="utf-8"))
    return np.array([complex(float(re), float(im)) for re, im in dat["f4"]][1:-1])


def P4(ak, t):
    k = np.arange(1, len(ak) + 1)
    return 2 * np.real(np.exp(2j * np.pi * np.outer(np.atleast_1d(t), k)) @ ak)


for name in ["e1_exp1_2_M128.json", "e1_exp1_2_M64.json", "e1_exp1_3_M64.json"]:
    ak = lade(name)
    zeile = []
    for K in [2, 3, 4, 8, 16, 64]:
        g = P4(ak, np.arange(K + 1) / K)
        p = 1 / K + np.diff(g)
        A = np.sum(np.sqrt(p / K))
        zeile.append(f"K={K}: 1-A_K={1 - A:.4e}")
    G = 20000
    t = np.arange(G) / G
    k = np.arange(1, len(ak) + 1)
    q = 2 * np.real(np.exp(2j * np.pi * np.outer(t, k)) @ (ak * 2j * np.pi * k))
    chi = -np.mean(np.log1p(q))
    print(f"{name}: " + "  ".join(zeile) + f"   (K->inf: 1-int sqrt(1+P4') = {1 - np.mean(np.sqrt(1 + q)):.6e}; chi = {chi:.8e})")
