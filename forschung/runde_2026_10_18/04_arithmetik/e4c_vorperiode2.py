"""E4c: Beginn der Periodizitaet von EF(n) mod 2^k (Periode 2): N(2^k) gegen k."""
from e4_expfakultaet import EF_mod
for k in (1, 2, 3, 4, 5, 18, 19, 20, 21, 22, 40):
    m = 2 ** k
    vals = [EF_mod(n, m) for n in range(1, 30)]
    N = 1
    while not all(vals[i] == vals[i + 2] for i in range(N - 1, len(vals) - 2)):
        N += 1
    print(k, N)
