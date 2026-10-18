"""z1: Zufällige Türme a_1^(a_2^(a_3^...)) mit unabhängig gleichverteilten Basen in (0, 1).

Für jede Basen-Folge: T_N = f_1∘…∘f_N(1) und T_{N+1} (f_a(x) = a^x). Lücke = |T_N − T_{N+1}|.
Konvergiert der Turm, geht die Lücke gegen 0; schwingt er (wie unter e^{-e}), bleibt sie gross.
Dazu der Ljapunow-Exponent der Zufalls-Fixlösung  χ = E[ log( T · |ln a_1| ) ]  (T = Turm, a_1 erste Basis);
deterministisch ist χ = log|λ| und χ = 0 genau bei e^{-e}.
Familie: a ~ Gleichverteilung auf (α, β).
Gegenprobe: N = 3000 gegen N = 6000, zwei Zufallssaaten.
"""
import sys
import numpy as np


def lauf(alpha, beta, N, M, seed):
    rng = np.random.default_rng(seed)
    a = rng.uniform(alpha, beta, size=(M, N + 1))
    la = np.log(a)
    # rückwärts: T_N = f_1(...f_N(1)), T_{N+1} = f_1(...f_{N+1}(1))
    x = np.ones(M)            # Kette für T_N  (Start oben bei Index N)
    y = np.exp(la[:, N])      # f_{N+1}(1)
    for k in range(N - 1, -1, -1):
        x = np.exp(la[:, k] * x)
        y = np.exp(la[:, k] * y)
    # Achtung: x nutzt Basen 0..N-1 (N Stück), y Basen 0..N (N+1 Stück)
    luecke = np.abs(x - y)
    # χ: Turm ab Index 1 (Basen 1..N) als S, dann T = a_0^S
    S = np.ones(M)
    for k in range(N, 0, -1):
        S = np.exp(la[:, k] * S)
    T = np.exp(la[:, 0] * S)
    chi = np.mean(np.log(T * np.abs(la[:, 0])))
    return luecke, chi


if __name__ == "__main__":
    print("Familie U(0, β):  Anteil Lücke > 1e-8,  mittlere Lücke,  χ   (N=3000/6000, M=400, Saaten 1/2)")
    for beta in [0.02, 0.05, 0.066, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0]:
        r = []
        for N, seed in [(3000, 1), (6000, 2)]:
            l, chi = lauf(1e-300, beta, N, 400, seed)
            r.append((np.mean(l > 1e-8), np.mean(l), chi))
        print(f"  β={beta:6.3f}:  " + "   |   ".join(f"{p:5.3f} {m:9.3e} χ={c:+.4f}" for p, m, c in r))
        sys.stdout.flush()
