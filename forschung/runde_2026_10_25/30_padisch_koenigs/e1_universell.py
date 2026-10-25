"""E1: kappa = alpha * sigma_lam(-lam) = K(log a) fuer a ≡ 1 -- reell UND p-adisch dieselbe Reihe.
(a) Ganzzahligkeit m! K_m (bis M)
(b) reell: K(ln a) gegen direkte Grenzwertbildung (a^^n - alpha)/lam^n, mpmath
(c) p-adisch: K(log_p a) gegen Turm-Grenzwert (unabhaengig, ueber exakte Tuerme mod p^N)
Aufruf: python e1_universell.py [M]
"""
import sys, time
from math import factorial
from fractions import Fraction as Fr
from kr_lib import *
import ar_lib
import mpmath as mp

t0 = time.time()
M = int(sys.argv[1]) if len(sys.argv) > 1 else 40
K = K_series(M)
ints = [K[m] * factorial(m) for m in range(M + 1)]
print("(a) m!K_m ganzzahlig bis M =", M, ":", all(x.denominator == 1 for x in ints))
print("    m!K_m, m=0..20:", [int(x) for x in ints[:21]])
print("    Zeit", round(time.time() - t0, 1), "s")

# (b) reell
mp.mp.dps = 90
for a in ("1.02", "1.05", "1.1", "1.15"):
    a = mp.mpf(a); L = mp.log(a)
    al = -mp.lambertw(-L) / L
    lam = L * al
    x = mp.mpf(1)  # a^^0
    n = 0
    vals = []
    while n < 400:
        x = a ** x; n += 1
        vals.append((x - al) / lam ** n)
        if abs(lam) ** n < mp.mpf(10) ** (-42):
            break
    kd = vals[-1]
    ks = mp.fsum(mp.mpf(K[m].numerator) / K[m].denominator * L ** m for m in range(M + 1))
    print(f"(b) a={mp.nstr(a,4)}  kappa direkt={mp.nstr(kd,25)}  K(ln a)={mp.nstr(ks,25)}  Abw={mp.nstr(abs(kd-ks),3)}  (e*L)^M={mp.nstr((mp.e*L)**M,3)}")

# (c) p-adisch
def kappa_direct(a, p, N):
    """kappa_p(a) fuer a ≡ 1 mod p (p=2: mod 4) direkt: (x_n - alpha)/lam^n, Rueckgabe (R, E) als R/p^E mod p^Nk."""
    q = p ** N
    al = ar_lib.tower_mod(a, 10 * N, q)
    L = ar_lib.padic_log(a, p, N)
    lam = al * L % q
    V = vp(lam, p)
    lu = lam // p ** V
    best = None
    # x_n - alpha hat Bewertung V(n+1); kappa_n = (x_n-alpha)/lam^n hat Bewertung V
    # Genauigkeit von kappa_n: min(N - nV, V(n+2)) -> n ~ (N-2V)/(2V)
    n = max(1, (N - 2 * V) // (2 * V))
    x = ar_lib.tower_mod(a, n, q)
    d = (x - al) % q
    v = vp(d, p)
    prec = min(N - n * V, V * (n + 2)) - 2
    ku = (d // p ** v) * pow(pow(lu, n, p ** prec), -1, p ** prec) % p ** prec
    kv = v - n * V
    return ku, kv, prec


print("(c) p-adisch: K(log_p a) gegen direkten Grenzwert")
bad = 0; cnt = 0
for p in (2, 3, 5, 7, 11, 13, 101):
    for a in range(2, 200):
        if (p > 2 and a % p != 1) or (p == 2 and a % 4 != 1):
            continue
        N = 90
        ku, kv, prec = kappa_direct(a, p, N)
        L = ar_lib.padic_log(a, p, N + 40)
        R, E = eval_padic(K, L, p, prec + kv + 5)
        # Wert = R/p^E; Bewertung und Einheitsteil vergleichen
        vR = vp(R, p) - E
        uR = R // p ** vp(R, p)
        P = min(prec, 30)
        ok = (vR == kv) and (uR - ku) % p ** P == 0
        cnt += 1
        if not ok:
            bad += 1
            print("   ABW", p, a, vR, kv, vp(uR - ku, p))
print(f"    Faelle {cnt}, Abweichungen {bad} (Vergleich auf min(prec,30) Stellen des Einheitsteils, M={M})")
print("Zeit", round(time.time() - t0, 1), "s")
