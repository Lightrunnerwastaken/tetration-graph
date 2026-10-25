# a_1..a_4 der oberen Hornabbildung hochgenau, zwei Hoehen + zwei Abbruch-Einstellungen
import mpmath as mp, sys, json
from fatou import horn
mp.mp.dps = int(sys.argv[1]); K = int(sys.argv[2]); dl = float(sys.argv[3]); M = int(sys.argv[4]); Y = mp.mpf(sys.argv[5])
cp = -1j * mp.pi / 3
vals = [horn(mp.mpc(mp.mpf(j) / M, Y), K=K, delta=dl, R=60) - mp.mpc(mp.mpf(j) / M, Y) - cp for j in range(M)]
res = {}
for n in range(0, 6):
    s = mp.fsum(vals[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M
    an = s * mp.exp(2 * mp.pi * n * Y)
    res[n] = an
    print(n, mp.nstr(an, 45))
json.dump({n: [str(mp.re(a)), str(mp.im(a))] for n, a in res.items()}, open(f"cache/a_{mp.mp.dps}_{K}_{dl}_{M}_{sys.argv[5]}.json", "w"))
