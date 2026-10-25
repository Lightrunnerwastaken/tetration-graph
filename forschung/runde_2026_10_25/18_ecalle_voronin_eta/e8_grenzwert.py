# Grenzwert s -> 0 der Entfaltungs-Koeffizienten: |t1| e^{2pi^2/s} -> |aI_1| = |a1| e^{-2pi^2/3}? t2/t1^2 -> 2 pi i - a2/a1^2 ?
import mpmath as mp, sys, json
from e6_entfaltung import setup, alpha_a, alpha_r, psi_a
mp.mp.dps = 200; M = 12
ss = [mp.mpf(x) for x in sys.argv[1:]]
absa1 = mp.mpf('64.15318451854875769074032856744156932425')
aI1 = absa1 * mp.exp(-2 * mp.pi**2 / 3)
L2 = mp.mpc('0.0253204930997575512295974960086', '4.70896624194763')
rows = []
for s in ss:
    P = setup(s); lam, mu, lnb, a, r = P
    w0 = alpha_a((a + r) / 2, P)
    th = [alpha_r(psi_a(w0 + mp.mpf(j) / M, P), P) - (w0 + mp.mpf(j) / M) for j in range(M)]
    t = {n: mp.fsum(th[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M * mp.expjpi(-2 * n * w0) for n in (1, 2)}
    x = abs(t[1]) * mp.exp(2 * mp.pi**2 / s) / aI1
    y = t[2] / t[1]**2
    rows.append((s, x, y, mp.log(mu)))
    print(mp.nstr(s, 4), ' |t1|e^{2pi^2/s}/|aI1| =', mp.nstr(x, 15), '  t2/t1^2 - L2 =', mp.nstr(y - L2, 6), '  1/ln(mu)-1/s =', mp.nstr(1 / mp.log(mu) - 1 / s, 10), flush=True)
# Polynom-Extrapolation in s (Grad len-1 und len-2)
for deg in [len(rows) - 1, len(rows) - 2, len(rows) - 3]:
    use = rows[-(deg + 1):]
    A = mp.matrix([[r[0] ** m for m in range(deg + 1)] for r in use])
    cx = mp.lu_solve(A, mp.matrix([r[1] for r in use]))
    cyr = mp.lu_solve(A, mp.matrix([mp.re(r[2]) for r in use])); cyi = mp.lu_solve(A, mp.matrix([mp.im(r[2]) for r in use]))
    cz = mp.lu_solve(A, mp.matrix([1 / r[3] - 1 / r[0] for r in use]))
    print('Grad', deg, ': x(0) =', mp.nstr(cx[0], 12), ' dx/ds =', mp.nstr(cx[1], 6), '  t2/t1^2(0) =', mp.nstr(mp.mpc(cyr[0], cyi[0]), 12), '  (1/s_r-1/s)(0) =', mp.nstr(cz[0], 12))
