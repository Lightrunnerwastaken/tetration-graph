"""E6: (a) Frischer Block = Koenigs-Gesetz: F_j(H) ≡ unit(kappa) * u^H (mod p^j), u = unit(lambda)
       (b) vorhergesagter frischer Block der Pentation: unit(kappa) * u^(alpha^(a)), alpha^ = profiniter Turmgrenzwert
       (c) profiniter Turmgrenzwert: alpha(0)=0, alpha(1)=1, alpha(-1)=-1 (a ≡ x mod Lambda(m), a gross)
"""
from ar_lib import *
from math import gcd

def unit_ord(u, p, j):
    o, x = 1, u % p ** j
    while x != 1:
        x = x * u % p ** j; o += 1
    return o

print("(a)+(b)")
for p, a, j in ((5, 2, 4), (5, 3, 4), (3, 2, 5), (7, 3, 3), (11, 2, 2), (2, 3, 6), (2, 5, 6)):
    N = 160
    q = p ** N
    seq = tower_run(a, p, N, 6, 6 + N + 10)
    al = seq[-1]
    lam = al * padic_log(a, p, N) % q
    V = vp(lam, p); u = lam // p ** V
    # kappa aus spaeter Hoehe
    H = 40
    d = (seq[H - 6] - al) % q; v = vp(d, p)
    ku = (d // p ** v) * pow(pow(u, H, p ** (N - v)), -1, p ** (N - v)) % p ** (N - v)
    bad = 0; tested = 0
    for H in range(10, 120):
        d = (seq[H - 6] - al) % q; v = vp(d, p)
        if v + j > N: break
        F = (d // p ** v) % p ** j
        if H * V > j + 5:  # erst ab hier ist der Fehler kleiner als p^j
            tested += 1
            if F != ku * pow(u, H, p ** j) % p ** j: bad += 1
    o = unit_ord(u, p, j)
    ahat = tower_mod(a, 200, o)  # alpha^ mod ord(u)
    beta = ku * pow(u, ahat, p ** j) % p ** j
    print(f"  p={p} a={a}: F_j(H)=unit(kappa)u^H mod p^{j}: {tested} Hoehen, Abweichungen {bad};"
          f" Pentationsblock beta_{j} = {beta} (ord u = {o}, alpha^ mod ord = {ahat})")

print("(c) profiniter Turmgrenzwert")
for m in (10 ** 6, 2 ** 20 * 3 ** 5, 7 ** 5 * 11 ** 3):
    L = 1; mm = m
    while mm > 1:
        L = L * mm // gcd(L, mm); mm = carmichael(mm)
    for x in (0, 1, -1, 2):
        a = x % L + 7 * L
        print(f"  m={m}: alpha({x}) mod m = {tower_mod(a, 60, m)}  (also {tower_mod(a, 60, m) - m if tower_mod(a, 60, m) > m // 2 else tower_mod(a, 60, m)})")
