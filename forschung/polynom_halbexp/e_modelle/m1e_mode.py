"""Modell 1e: genaue Bestimmung der dominanten Fehlermode z (e_N ≈ Re(A N^β z^N)) für die elliptische Möbius-Abbildung.
Gleitende Prony-Anpassung (2 Moden, konjugiert) an d_N = a_0^(N) − τ, dann Extrapolation in 1/N."""
import json, sys, time
from pathlib import Path
import mpmath as mp
from kern import zweig, g_moebius_ell, horner

mp.mp.dps = 80
HIER = Path(__file__).parent
NMAX = int(sys.argv[1])
CS = [mp.mpf(x) for x in sys.argv[2].split(",")]
BITS = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
aus = {}
for c in CS:
    t0 = time.time()
    tau = mp.tan(c / 2)
    loes = zweig(g_moebius_ell(c), NMAX, bits=BITS)
    d = {N: loes[N][0] - tau for N in loes}
    zs = {}
    for N in range(10, NMAX - 3):
        # e_{n+2} = s e_{n+1} - p e_n  für n = N, N+1  ->  z, zbar Wurzeln von t^2 - s t + p
        A = mp.matrix([[d[N + 1], -d[N]], [d[N + 2], -d[N + 1]]])
        b = mp.matrix([d[N + 2], d[N + 3]])
        try:
            s, p = mp.lu_solve(A, b)
        except ZeroDivisionError:
            continue
        disc = s * s - 4 * p
        z = (s + mp.sqrt(disc)) / 2
        zs[N] = z
    # Richardson in 1/N: z(N) ≈ z + u/N + v/N^2
    Ns = sorted(zs)
    n1, n2, n3 = Ns[-1], Ns[-1] - 10, Ns[-1] - 20
    M = mp.matrix([[1, 1 / mp.mpf(n), 1 / mp.mpf(n) ** 2] for n in (n1, n2, n3)])
    r = mp.lu_solve(M, mp.matrix([zs[n1], zs[n2], zs[n3]]))
    zinf = r[0]
    print(f"c={mp.nstr(c, 4)} tau={mp.nstr(tau, 8)}  z(N) bei N={n3},{n2},{n1}: " +
          ", ".join(f"{mp.nstr(abs(zs[n]), 6)}∠{mp.nstr(mp.arg(zs[n]) * 180 / mp.pi, 6)}" for n in (n3, n2, n1)) +
          f"  extrapoliert |z|={mp.nstr(abs(zinf), 6)} arg={mp.nstr(mp.arg(zinf) * 180 / mp.pi, 6)}°  ({time.time() - t0:.0f}s)", flush=True)
    aus[mp.nstr(c, 6)] = {"tau": mp.nstr(tau, 12), "z_N": {n: [mp.nstr(mp.re(zs[n]), 10), mp.nstr(mp.im(zs[n]), 10)] for n in Ns},
                          "z_inf": [mp.nstr(mp.re(zinf), 8), mp.nstr(mp.im(zinf), 8)], "d_N": {N: mp.nstr(d[N], 6) for N in d}}
json.dump(aus, open(HIER / "daten" / f"m1e_mode_{sys.argv[2].replace(',', '_')}.json", "w"), indent=0)
