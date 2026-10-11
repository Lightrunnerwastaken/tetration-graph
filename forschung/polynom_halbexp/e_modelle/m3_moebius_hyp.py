"""Modell 3: hyperbolische Möbius-Abbildung g(x) = tanh(artanh x + c) (reelle Fixpunkte ±1), Wurzel f_M = tanh(artanh x + c/2).
Zweig von N = 1 aus und Fehler zu f_M; Prony-Mode der a_0-Fehler."""
import sys, time
import mpmath as mp
from kern import zweig, g_moebius_hyp, horner

mp.mp.dps = 100
NMAX = int(sys.argv[1])
CS = [mp.mpf(x) for x in sys.argv[2].split(",")]
for c in CS:
    t0 = time.time()
    T = mp.tanh(c / 2)
    fM = lambda z: (z + T) / (1 + T * z)
    try:
        loes = zweig(g_moebius_hyp(c), NMAX, bits=1000)
    except RuntimeError as e:
        print("c =", c, e)
        continue
    d = {N: loes[N][0] - T for N in loes}
    e5 = {N: horner(loes[N], mp.mpf("0.5")) - fM(mp.mpf("0.5")) for N in loes}
    em = {N: horner(loes[N], mp.mpf("-0.9")) - fM(mp.mpf("-0.9")) for N in loes}
    zs = {}
    for N in range(10, NMAX - 3):
        A = mp.matrix([[d[N + 1], -d[N]], [d[N + 2], -d[N + 1]]])
        b = mp.matrix([d[N + 2], d[N + 3]])
        try:
            s, p = mp.lu_solve(A, b)
            zs[N] = (s + mp.sqrt(s * s - 4 * p)) / 2
        except ZeroDivisionError:
            pass
    print(f"\n== c={mp.nstr(c, 4)} T=tanh(c/2)={mp.nstr(T, 6)} tanh(c)={mp.nstr(mp.tanh(c), 6)} Pol von f bei {mp.nstr(-1 / T, 5)}  ({time.time() - t0:.0f}s)")
    for N in (1, 2, 3, 5, 10, 20, 30, 40, 60, 80, 100):
        if N in loes:
            zz = zs.get(N - 4)
            zt = f"  z(N)={mp.nstr(abs(zz), 5)}∠{mp.nstr(mp.arg(zz) * 180 / mp.pi, 5)}" if zz is not None else ""
            print(f"N={N:3d} a0-T={mp.nstr(d[N], 4):>11} p(.5)-f(.5)={mp.nstr(e5[N], 4):>11} p(-.9)-f(-.9)={mp.nstr(em[N], 4):>11}{zt}")
