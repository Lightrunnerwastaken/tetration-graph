"""Modell 1b: elliptische Möbius-Abbildung bei grösserem Drehwinkel c — wo bricht der Zweig ab und warum?"""
import sys, time
import mpmath as mp
from kern import zweig, g_moebius_ell, horner

mp.mp.dps = 60
NMAX = int(sys.argv[1])
c = mp.mpf(sys.argv[2])
tau = mp.tan(c / 2)
fM = lambda z: (z + tau) / (1 - tau * z)
fk = [tau] + [(1 + tau ** 2) * tau ** (k - 1) for k in range(1, NMAX + 1)]
print(f"c = {c}, tau = {mp.nstr(tau, 6)}, Pol 1/tau = {mp.nstr(1 / tau, 6)}")


def zeige(N, a):
    d = horner(a, mp.mpf("0.5")) - fM(mp.mpf("0.5"))
    dk = [mp.nstr(a[k] - fk[k], 3) for k in range(min(4, N + 1))]
    print(f"N={N:3d} a0={mp.nstr(a[0], 10)} p(.5)-f(.5)={mp.nstr(d, 3)} da_0..3={dk} a_N={mp.nstr(a[N], 4)}", flush=True)


try:
    zweig(g_moebius_ell(c), NMAX, bits=600, ausgabe=zeige, still=False)
except RuntimeError as e:
    print("ABBRUCH:", e)
