import json, sys
import mpmath as mp
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_80.json"))]
for N in (20, 40, 60, 80, 100, 120, 140, 160, 180, 194):
    aN = P[N][N]
    print(N, "a_N^(N) =", mp.nstr(aN, 5), " |a_N|^(-1/N) =", mp.nstr(abs(aN) ** (-mp.mpf(1) / N), 6),
          " |a_N/K_N| =", mp.nstr(abs(aN / K[N]), 4) if N <= 80 else "")
# Nullstellen von p_N: Lage
for N in (60, 120, 194):
    c = [float(x) for x in reversed(P[N])]
    r = np.roots(c)
    print(N, "Nullstellen |z| min/median/max:", np.round(np.min(abs(r)), 3), np.round(np.median(abs(r)), 3), np.round(np.max(abs(r)), 3))
    np.save(f"b05_nullstellen_{N}.npy", r)
