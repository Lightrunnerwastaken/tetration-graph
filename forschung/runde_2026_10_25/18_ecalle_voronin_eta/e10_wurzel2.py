# Beispiel b = sqrt(2): Fixpunkte 2 und 4. Welle theta zwischen den regulaeren Abelfunktionen an 2 und 4.
import mpmath as mp
from e6_entfaltung import setup, alpha_a, alpha_r, psi_a
mp.mp.dps = 80; M = 12
s = -mp.log(mp.log(2)); P = setup(s); lam, mu, lnb, a, r = P
print('b =', mp.nstr(mp.exp(lnb), 15), ' a =', mp.nstr(a, 15), ' r =', mp.nstr(r, 15))
w0 = alpha_a(mp.mpf(3), P)
th = [alpha_r(psi_a(w0 + mp.mpf(j) / M, P), P) - (w0 + mp.mpf(j) / M) for j in range(M)]
t = {n: mp.fsum(th[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M * mp.expjpi(-2 * n * w0) for n in (0, 1, 2)}
for n in t: print('t', n, mp.nstr(t[n], 15), ' |.|=', mp.nstr(abs(t[n]), 10))
pred = mp.mpf('64.15318451854875769') * mp.exp(-2 * mp.pi**2 / s - 2 * mp.pi**2 / 3)
print('Vorhersage |a1| e^{-2pi^2/s-2pi^2/3} =', mp.nstr(pred, 10), '  Verhaeltnis =', mp.nstr(abs(t[1]) / pred, 10))
print('max|theta - t0| ~ 2|t1| =', mp.nstr(2 * abs(t[1]), 6))
