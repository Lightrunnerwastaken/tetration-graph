"""Abstand Boettcher <-> Kneser fuer x^d + c: Fourier-Koeffizienten von p (S_B(x) = S_K(x-1+p(x)))."""
import numpy as np, json, os, sys, math, time
from zoo import *
from s7_boettcher import logphi, green, SB, AK
erg = []
ALLE = [(2, 0.3), (2, 0.5), (2, 1.0), (2, 2.0), (2, 4.0), (4, 1.0), (4, 2.0), (6, 1.0), (6, 3.0)]
teil = int(sys.argv[1]); faelle = ALLE[:5] if teil == 1 else (ALLE[5:] if teil == 2 else [(2, 0.3)])
for d, c in faelle:
    t0 = time.time()
    a = abb_quadrat(c, d)
    s0 = c ** d + c
    for (T, P, n) in ([(110, 1.0, 16), (110, 1.0, 12)] if teil == 3 else ([(40, 1.0, 16), (50, 0.75, 16)] if c == 1.0 else [(40, 1.0, 16), (40, 1.0, 12)])):
        s = Loesung(a, 0.0, s0, T=T, P=P, n=n); res = s.loese(maxit=60, info=False)
        Gc = green(c, d); M = 256; x = np.arange(M) / M
        wK = s.Sweit(x - 1 + 0j); wB = SB(x, wK, c, d, Gc)
        zK = AK(s, wB.real + 0j, x - 1 + 0j); p = zK.real + 1 - x
        ph = np.fft.fft(p) / M
        if T == 40: pA = ph
    dAB = abs(ph[1] - pA[1])
    z = dict(d=d, c=c, L=[a.L.real, a.L.imag], res=res, p1=abs(ph[1]), p2=abs(ph[2]), p3=abs(ph[3]), dAB=dAB,
             ln_p1=math.log(abs(ph[1])), ziel=-math.pi ** 2 / math.log(d), maxp=float(abs(p).max()))
    erg.append(z)
    print(f'd={d} c={c}: |p^(1)|={z["p1"]:.4e} (A-B {dAB:.0e}) |p^(2)|={z["p2"]:.3e} |p^(3)|={z["p3"]:.2e}  ln|p^(1)|={z["ln_p1"]:.3f}  -pi^2/ln d={z["ziel"]:.3f}  ln|p2/p1| = {math.log(z["p2"]/z["p1"]):.3f}  ({time.time()-t0:.0f}s)', flush=True)
json.dump(erg, open(os.path.join(DATEN, f'boettcher_scan_{teil}.json'), 'w'), indent=1)
