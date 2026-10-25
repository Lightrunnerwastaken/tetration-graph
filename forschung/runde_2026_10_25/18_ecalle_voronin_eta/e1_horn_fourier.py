# Fourier-Koeffizienten der Hornabbildungen von e^u-1 auf mehreren Hoehen Y.
import mpmath as mp, sys, json
from fatou import horn
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
M = int(sys.argv[2]) if len(sys.argv) > 2 else 32
Ys = [float(y) for y in sys.argv[3:]] or [1.0, 1.5]
c_plus = -1j * mp.pi / 3
out = {}
for Y in Ys:
    vals = []
    for j in range(M):
        W = mp.mpc(mp.mpf(j) / M, Y)
        cp = c_plus if Y > 0 else -c_plus
        vals.append(horn(W, K=36, delta=0.04, R=50) - W - cp)
    res = []
    sgn = 1 if Y > 0 else -1
    for n in range(0, M // 2):
        s = mp.fsum(vals[j] * mp.expjpi(-2 * sgn * n * mp.mpf(j) / M) for j in range(M)) / M
        an = s * mp.exp(2 * mp.pi * n * abs(Y))
        res.append(an)
        print(f"Y={Y} n={n:2d}  a_n = {mp.nstr(an, 18)}   |a_n e^-2pinY| = {mp.nstr(abs(s), 3)}")
    out[str(Y)] = [[str(mp.re(a)), str(mp.im(a))] for a in res]
json.dump(out, open(f"cache/e1_{mp.mp.dps}_{M}.json", "w"), indent=1)
