import mpmath as mp, json, sys
mp.mp.dps = 60
c = [None] + [mp.mpf(x) for x in json.load(open(sys.argv[1]))]
K = len(c) - 1
a1 = mp.mpc('10.7488070931472942', '63.2462981521074514')
def G(k): return mp.power(1j, -k) * mp.gamma(k - 2j * mp.pi / 3) * mp.power(4 * mp.pi, -k + 2j * mp.pi / 3)
def Dk(k):
    g0, g1 = G(k), G(k + 1)
    A = mp.matrix([[2 * mp.re(g0), -2 * mp.im(g0)], [2 * mp.re(g1), -2 * mp.im(g1)]])
    x = mp.lu_solve(A, mp.matrix([c[k], c[k + 1]]))
    return mp.mpc(x[0], x[1])
ref = -1j * a1 * mp.exp(-mp.pi**2 / 3) / (2 * mp.pi)
for order in [2, 4, 6, 8, 10]:
    for step in [4, 6]:
        ks = [K - 1 - step * j for j in range(order + 1)]
        # Fit D_k = sum_{m=0}^{order} e_m / k^m (Interpolation)
        A = mp.matrix([[mp.mpf(k) ** (-m) for m in range(order + 1)] for k in ks])
        br = mp.lu_solve(A, mp.matrix([mp.re(Dk(k)) for k in ks]))
        bi = mp.lu_solve(A, mp.matrix([mp.im(Dk(k)) for k in ks]))
        Dinf = mp.mpc(br[0], bi[0])
        print(order, step, mp.nstr(Dinf, 16), ' D/ref-1 =', mp.nstr(Dinf / ref - 1, 3), ' 1/k-Koeff:', mp.nstr(mp.mpc(br[1], bi[1]) / ref, 8))
