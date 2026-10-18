"""A12: M = max |g| auf |z| = 2.5 (g = sexp - log(z+2), holomorph in |z|<3), für |d_n| <= M 2.5^{-n}.
Halbstellen-Gitter (kein Punkt genau auf dem Schnitt; dort hängt der Motor), N = 400 und 300, dps 40 (bei dps 30 hängt der Motor bei Im z ~ 2.5).
"""
import mpmath as mp
from motor import sexp
mp.mp.dps = 40
r = mp.mpf("2.5")
for N in (400, 300):
    pts = [r * mp.expjpi((j + mp.mpf(0.5)) / N) for j in range(N)]
    v = sexp(pts, 40)
    g = [abs(a - mp.log(z + 2)) for a, z in zip(v, pts)]
    i = max(range(len(g)), key=lambda k: g[k])
    print(N, "max|g| =", mp.nstr(g[i], 12), "bei arg/pi =", mp.nstr((i + 0.5) / N, 5), " |g| nahe arg=pi:", mp.nstr(g[-1], 6), flush=True)
print("g(2.5) reell:", mp.nstr(sexp([mp.mpf(2.5)], 40)[0] - mp.log(4.5), 12))
M = 180
n0 = next(n for n in range(1, 1000) if M * n * (mp.mpf(2) / r) ** n < mp.mpf("0.5"))
print("M =", M, ": |d_n| < 1/(2 n 2^n) sicher für n >=", n0)
