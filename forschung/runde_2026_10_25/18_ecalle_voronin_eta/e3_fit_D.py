# Bestimme D aus c_k, c_{k+1} und vergleiche mit a_1 e^{-pi^2/3} / (2 pi i) (Vorzeichen/Phase offen)
import mpmath as mp, json, sys
mp.mp.dps = 60
c = [None] + [mp.mpf(x) for x in json.load(open(sys.argv[1]))]
K = len(c) - 1
a1 = mp.mpc('10.7488070931472942', '63.2462981521074514')
def G(k): return mp.power(1j, -k) * mp.gamma(k - 2j * mp.pi / 3) * mp.power(4 * mp.pi, -k + 2j * mp.pi / 3)
Ds = {}
for k in list(range(20, K, 10)) + [K - 1]:
    g0, g1 = G(k), G(k + 1) if k + 1 <= K else None
    if k + 1 > K: continue
    A = mp.matrix([[2 * mp.re(g0), -2 * mp.im(g0)], [2 * mp.re(g1), -2 * mp.im(g1)]])
    x = mp.lu_solve(A, mp.matrix([c[k], c[k + 1]]))
    D = mp.mpc(x[0], x[1]); Ds[k] = D
    print(k, mp.nstr(D, 15), ' |D| =', mp.nstr(abs(D), 12))
ref = a1 * mp.exp(-mp.pi**2 / 3) / (2 * mp.pi)
print('a1 e^{-pi^2/3}/(2pi) =', mp.nstr(ref, 15), '|.|=', mp.nstr(abs(ref), 12))
k = max(Ds); print('D/ref =', mp.nstr(Ds[k] / ref, 12))
# Richardson in 1/k
ks = sorted(Ds)[-4:]
