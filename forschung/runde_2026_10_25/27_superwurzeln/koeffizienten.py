# Taylor-Koeffizienten von eps_n und omega_n = eps_n^{-1} um 0, exakt (python-flint fmpq_series).
# eps_0 = x, eps_{n+1} = x*exp(eps_n).  omega_n = Reversion von eps_n.
# Ausgabe: a_n(k) = k! [y^k] omega_n (ganze Zahlen), Stabilisierung gegen omega_inf = y e^{-y},
# Wurzeltest-Schaetzung des Konvergenzradius R_n.
import sys, math
import flint
from flint import fmpq_series, fmpq, fmpz

K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 6

flint.ctx.cap = K + 1
x = fmpq_series([0, 1], prec=K + 1)

def eps(n):
    e = x
    for _ in range(n):
        e = x * e.exp()
    return e

def coeffs(s, K):
    c = s.coeffs()
    c = list(c) + [fmpq(0)] * (K + 1 - len(c))
    return c[:K + 1]

res = {}
for n in range(1, NMAX + 1):
    en = eps(n)
    wn = en.reversion()
    c = coeffs(wn, K)
    a = [c[k] * math.factorial(k) for k in range(K + 1)]
    assert all(v.q == 1 for v in a), "nicht ganz"
    a = [int(v.p) for v in a]
    res[n] = a
    # Ableitung: omega_inf hat k![y^k] = k*(-1)^(k-1)
    ainf = [0] + [k * (-1) ** (k - 1) for k in range(1, K + 1)]
    firstdiff = next((k for k in range(K + 1) if a[k] != ainf[k]), None)
    # Wurzeltest |c_k|^{-1/k} fuer die letzten k
    est = []
    for k in range(K - 4, K + 1):
        ck = abs(float(c[k])) if c[k] != 0 else 0.0
        est.append(ck ** (-1.0 / k) if ck > 0 else float('inf'))
    print(f"n={n}: a_n(1..12) = {a[1:13]}")
    print(f"      erste Abweichung von omega_inf bei k={firstdiff}; Wurzeltest R ~ {[round(e,5) for e in est]}")

# a_n(k) fuer n=2 als Folge (OEIS-Suche)
print("omega_2 k![y^k]:", res[2][1:16])
if NMAX >= 3:
    print("omega_3 k![y^k]:", res[3][1:16])
