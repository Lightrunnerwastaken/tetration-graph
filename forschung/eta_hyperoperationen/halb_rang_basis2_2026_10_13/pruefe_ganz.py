"""Pruefung der Rangschritte an ganzen Hoehen: c_n(k+1) aus c_n(k) und den Ableitungen von c_{n-1} bei k+1.

S_n(k+1) = S_{n-1}(S_n(k)),  S_n(k) = k+1+e_n(k);  Taylor um k+1 mit e_n(k) als Reihe in eps.
"""
import sys
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
from halbschritt import lade
import mpmath as mp
mp.mp.dps = 80
d = lade(sys.argv[1]); J = d['J']; K = d['K']


def mul(a, b):
    c = [mp.mpf(0)] * (J + 1)
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(J + 1 - i):
            c[i + j] += x * b[j]
    return c


worst = 0
for n in range(1, d['N'] + 1):
    for k in range(1, K):
        G = d['ganz']
        e = [mp.mpf(0)] + [mp.mpf(G[str(k)][n][0][j]) for j in range(J)]  # e_n(k) als eps-Reihe
        # S_{n-1}(k+1+e) = k+2+e + sum_j eps^j sum_r c_{n-1}^{(j),(r)}(k+1) e^r/r!
        pred = [mp.mpf(0)] * (J + 1)
        for i in range(J + 1): pred[i] += e[i]
        powe = [mp.mpf(1)] + [mp.mpf(0)] * J
        for r in range(J):
            for j in range(1, J + 1):
                c = mp.mpf(G[str(k + 1)][n - 1][r][j - 1]) / mp.factorial(r)
                for i in range(J + 1 - j):
                    pred[i + j] += c * powe[i]
            powe = mul(powe, e)
        ist = [mp.mpf(0)] + [mp.mpf(G[str(k + 1)][n][0][j]) for j in range(J)]
        rel = max(abs(pred[j] - ist[j]) / (1 + abs(ist[j])) for j in range(1, J + 1))
        worst = max(worst, rel)
print('groesste relative Abweichung (alle n, k=2..K):', mp.nstr(worst, 3))
