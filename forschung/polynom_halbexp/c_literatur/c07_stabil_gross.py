"""Strang C, c07: Stabilität für grosse N und für das lineare Modell.

(1) |J_N(p_N)^{-1}| in l^1 mit w_k = |L|^k (Urbild = Bild), N bis 194 (arb, mittlere Werte).
(2) Lineares Modell (Paulsen 2016, Julia-Gleichung lambda(e^x) = e^x lambda(x), c0 = 1 normiert):
    Unbekannte c_1..c_N, Zeilen i = 1..N:  [x^i] (lambda∘exp - e^x lambda) = 0 mit lambda = 1 + sum c_j x^j.
    Matrix M[i][j] = [x^i] (e^{x})^j - [x^{i-j}] e^x = j^i/i! - 1/(i-j)!  (j ≤ i), rechte Seite -([x^i]1∘exp - e^x·1) = -(0 - 1/i!) ...
    Hier nur die Norm von M_N^{-1} in l^1_{|L|} (Stabilität), und c_0-Normierung wie bei Paulsen: lambda(0) = c0 frei,
    wir setzen den konstanten Term 1 und lösen die übrigen (entspricht Paulsens System bis auf Skalierung).
    Ebenso Walkers Abel-System G(e^x) = G(x) + 1: M[i][j] = j^i/i! - delta_ij, rechte Seite e_0.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, factorial
from flint import arb_mat, arb, ctx

mp.dps = 160
ctx.prec = 700
HIER = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HIER, "..", "gemeinsam", "polynome_bis_194.json"), encoding="utf-8"))
L_abs = mpf("1.3745570107437075")


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def jacobi(a):
    n = len(a) - 1
    pot = [[mpf(1)] + [mpf(0)] * n]
    for _ in range(1, n + 1):
        pot.append(mul(pot[-1], a, n))
    dp = [k * a[k] for k in range(1, n + 1)]
    dpp = [mpf(0)] * (n + 1)
    for j, c in enumerate(dp):
        for k in range(n + 1):
            dpp[k] += c * pot[j][k]
    return [[pot[j][k] + (dpp[k - j] if j <= k else mpf(0)) for j in range(n + 1)] for k in range(n + 1)]


def invnorm(M, w):
    n = len(M)
    A = arb_mat([[arb(mp.nstr(x, 150)) for x in z] for z in M])
    Ai = A.inv()
    best = 0.0
    for j in range(n):
        s = sum(abs(float(Ai[i, j].mid())) * float(w[i]) for i in range(n)) / float(w[j])
        best = max(best, s)
    return best


print("(1) Polynom-Wurzel, |J_N^-1| in l1_|L|:")
for N in (60, 80, 100, 120, 150, 194):
    a = [mpf(x) for x in P[str(N)]]
    w = [L_abs ** k for k in range(N + 1)]
    print(f"   N={N:3d}: {invnorm(jacobi(a), w):.3f}")
    sys.stdout.flush()

print("(2) Lineare Modelle, |M_N^-1| in l1_|L| (Indizes j = 1..N, Gewicht |L|^j):")
for N in (20, 40, 60, 80, 120):
    # Julia: lambda = sum_{j>=0} c_j x^j, lambda(e^x) - e^x lambda(x); Spalte j: (e^x)^j - e^x x^j
    Mj = [[mpf(j) ** i / factorial(i) - (mpf(1) / factorial(i - j) if j <= i else 0) for j in range(1, N + 1)]
          for i in range(1, N + 1)]
    # Abel (Walker): G(e^x) - G(x); Spalte j: (e^x)^j - x^j, Zeilen i = 0..N-1
    Ma = [[mpf(j) ** i / factorial(i) - (1 if i == j else 0) for j in range(1, N + 1)] for i in range(0, N)]
    wj = [L_abs ** j for j in range(1, N + 1)]
    wa_out = [L_abs ** i for i in range(0, N)]
    # Norm Urbild-Gewicht |L|^j (j=1..N), Bild-Gewicht |L|^i (Zeilenindex)
    def invn(M, w_in, w_out):
        n = len(M)
        A = arb_mat([[arb(mp.nstr(x, 150)) for x in z] for z in M])
        Ai = A.inv()
        return max(sum(abs(float(Ai[r, c].mid())) * float(w_in[r]) for r in range(n)) / float(w_out[c]) for c in range(n))
    print(f"   N={N:3d}: Julia {invn(Mj, wj, [L_abs ** i for i in range(1, N + 1)]):.3e}   Abel {invn(Ma, wj, wa_out):.3e}")
    sys.stdout.flush()
