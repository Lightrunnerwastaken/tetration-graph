"""z2: Schwelle β* für Türme mit Basen ~ U(0, β): ab wann konvergiert der Zufallsturm?
Lücke |T_N − T_{N+1}| für N = 4000 und 16000 (M = 1000 Folgen), χ mit M = 20000 (N = 2000).
Aufruf: python z2_schwelle.py
"""
import sys
import numpy as np
from z1_zufallsturm import lauf

if __name__ == "__main__":
    print(" β       P(L>1e-8) N=4000  N=16000   mittl. Lücke N=16000    χ (M=20000, Saat 3)   χ (Saat 4)")
    for beta in np.arange(0.020, 0.0501, 0.0025):
        p = []
        for N, seed in [(4000, 5), (16000, 6)]:
            l, _ = lauf(1e-300, beta, N, 1000, seed)
            p.append((np.mean(l > 1e-8), np.mean(l)))
        _, c1 = lauf(1e-300, beta, 2000, 20000, 3)
        _, c2 = lauf(1e-300, beta, 2000, 20000, 4)
        print(f" {beta:.4f}   {p[0][0]:6.3f}   {p[1][0]:6.3f}     {p[1][1]:10.3e}         {c1:+.4f}            {c2:+.4f}")
        sys.stdout.flush()
