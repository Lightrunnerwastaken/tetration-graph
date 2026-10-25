# c_k = 2 Re[ sum_{m=0}^{M} d_m G_{k-m} ],  G_j = i^{-j} Gamma(j - 2 pi i/3) (4 pi)^{-j + 2 pi i/3}
import mpmath as mp, json, sys
mp.mp.dps = 80
c = [None] + [mp.mpf(x) for x in json.load(open(sys.argv[1]))]
K = len(c) - 1
a1 = mp.mpc('10.7488070931472942', '63.2462981521074514')
def G(j): return mp.power(1j, -j) * mp.gamma(j - 2j * mp.pi / 3) * mp.power(4 * mp.pi, -j + 2j * mp.pi / 3)
ref = -1j * a1 * mp.exp(-mp.pi**2 / 3) / (2 * mp.pi)
for M in [3, 6, 9, 12, 15]:
    for kend in [K, K - 10]:
        ks = list(range(kend - 2 * (M + 1) + 1, kend + 1))
        rows = []
        for k in ks:
            row = []
            for m in range(M + 1):
                g = G(k - m); row += [2 * mp.re(g), -2 * mp.im(g)]
            rows.append(row)
        x = mp.lu_solve(mp.matrix(rows), mp.matrix([c[k] for k in ks]))
        d0 = mp.mpc(x[0], x[1]); d1 = mp.mpc(x[2], x[3])
        print(M, kend, 'd0 =', mp.nstr(d0, 20), ' d0/ref-1 =', mp.nstr(d0 / ref - 1, 3), ' d1/d0 =', mp.nstr(d1 / d0, 10))
