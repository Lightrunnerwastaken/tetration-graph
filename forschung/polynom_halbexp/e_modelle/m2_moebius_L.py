"""Modell 2: elliptische Möbius-Abbildung mit Fixpunkten A ± iB und Drehwinkel c:
   g(x) = A + B·tan(arctan((x − A)/B) + c),  f_M(x) = A + B·tan(arctan((x − A)/B) + c/2).
Spezialfall „Möbius-exp“: A + iB = L (Fixpunkt von exp), c so, dass g(0) = 1.
Aufruf: python m2_moebius_L.py NMAX A B c   (c = 'exp' wählt g(0) = 1)."""
import json, sys, time
from pathlib import Path
import mpmath as mp
from flint import arb
from kern import zweig, horner

mp.mp.dps = 80
HIER = Path(__file__).parent
NMAX = int(sys.argv[1])
A, B = mp.mpf(sys.argv[2]), mp.mpf(sys.argv[3])
c = mp.atan((1 - A) / B) + mp.atan(A / B) if sys.argv[4] == "exp" else mp.mpf(sys.argv[4])
BITS = int(sys.argv[5]) if len(sys.argv) > 5 else 900
g = lambda x: A + B * mp.tan(mp.atan((x - A) / B) + c)
fM = lambda x: A + B * mp.tan(mp.atan((x - A) / B) + c / 2)
gk_mp = mp.taylor(g, 0, NMAX)
fk = mp.taylor(fM, 0, NMAX)
gk = lambda N: [arb(mp.nstr(x, 75)) for x in gk_mp[:N + 1]]
pol = A + B * mp.cot(c / 2)
print(f"A={mp.nstr(A, 6)} B={mp.nstr(B, 6)} c={mp.nstr(c, 8)}: g(0)={mp.nstr(g(0), 8)} g'(0)={mp.nstr(gk_mp[1], 8)} f_M(0)={mp.nstr(fM(0), 8)} "
      f"Pol von f_M bei {mp.nstr(pol, 6)}, |Fixpunkt|={mp.nstr(mp.sqrt(A * A + B * B), 6)}")
t0 = time.time()
try:
    loes = zweig(gk, NMAX, bits=BITS)
except RuntimeError as e:
    print("ABBRUCH", e)
    sys.exit()
d = {N: loes[N][0] - fk[0] for N in loes}
e1 = {N: horner(loes[N], mp.mpf(1)) - fM(mp.mpf(1)) for N in loes}
for N in (1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100):
    if N in loes:
        print(f"N={N:3d} a0-f(0)={mp.nstr(d[N], 4):>11}  p(1)-f(1)={mp.nstr(e1[N], 4):>11}")
hull = lambda dd, N: max(abs(dd[n]) for n in range(N - 5, N + 1))
for N1, N2 in ((10, 30), (30, 50), (50, 70), (70, 90)):
    if N2 in loes:
        print(f"Rate (Einhüllende a0) {N1}->{N2}: {mp.nstr((hull(d, N2) / hull(d, N1)) ** (mp.mpf(1) / (N2 - N1)), 5)}")
print(f"Zeit {time.time() - t0:.0f}s")
