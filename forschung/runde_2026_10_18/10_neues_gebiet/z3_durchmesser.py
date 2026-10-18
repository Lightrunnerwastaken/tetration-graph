"""z3: Durchmesser D_N = |F_N(1) − F_N(0)| des Bildes von [0,1] unter F_N = f_1∘…∘f_N (Basen ~ U(0, β)).
D_N fällt monoton in N (geschachtelte Intervalle); der Turm konvergiert genau dann, wenn D_N -> 0.
Gleiche Basen-Folge für alle N (Präfix). M Folgen; Ausgabe: Anteil D_N > 1e-8 und Median log10 D_N.
Gegenprobe: zwei Saaten. Aufruf: python z3_durchmesser.py
"""
import sys
import numpy as np

def durchmesser(beta, Ns, M, seed):
    rng = np.random.default_rng(seed)
    la = np.log(rng.uniform(1e-300, beta, size=(M, max(Ns))))
    out = []
    for N in Ns:
        x = np.ones(M); y = np.zeros(M)
        for k in range(N - 1, -1, -1):
            x = np.exp(la[:, k] * x); y = np.exp(la[:, k] * y)
        D = np.abs(x - y)
        out.append((N, np.mean(D > 1e-8), np.median(np.log10(np.maximum(D, 1e-300)))))
    return out

if __name__ == "__main__":
    Ns = [1000, 4000, 16000, 64000, 256000]
    for beta in [0.015, 0.02, 0.025, 0.03, 0.035]:
        for seed in [11, 12]:
            r = durchmesser(beta, Ns, 200, seed)
            print(f"β={beta:.3f} Saat {seed}: " + "  ".join(f"N={N}: {p:.3f}/{m:7.1f}" for N, p, m in r))
            sys.stdout.flush()
