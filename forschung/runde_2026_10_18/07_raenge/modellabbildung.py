"""Plateau-Modellabbildung der Rangleiter (beschraenkter Bereich):
    mu_r = lam_r^(a+delta_{r+1}),  delta_{r+1} = delta_r (1-mu_r),  lam_{r+1} = delta_r * ln(1/lam_r) * mu_r
in Log-Variablen (D = -ln delta, L = -ln lam), iteriert bis r = 10^7.
Start: letzte Zeile einer Koenigs-Leiter-Datei. Gibt 1/mu-Zuwachs, r*mu, delta*r^k*ln r aus.
python modellabbildung.py daten/koenigs_....json
"""
import sys, json, math
import mpmath as mp

d = json.load(open(sys.argv[1]))
b = float(eval(d['b'], {'sqrt': math.sqrt, 'e': math.e, 'exp': math.exp}))
a = b - 1
k = (2 - b) / (b - 1)
last = d['rows'][-1]
r = last['r']
D = -math.log(float(mp.mpf(last['p'])) - b)
L = -math.log(float(mp.mpf(last['lam'])))
A = a * k ** k
print(f'b={b} k={k:.6f} 1/k={1/k:.6f} A=a k^k={A:.6f}; Start r={r}')
checkpoints = {10 ** j for j in range(2, 8)} | {20, 30, 50, 200, 500, 2000, 5000}
mu_prev = None
while r < 10 ** 7:
    # mu = lam^(a+delta'), delta' haengt schwach von mu ab: zwei Fixpunktschritte
    delta = math.exp(-D)
    mu = math.exp(-a * L)
    for _ in range(3):
        dn = delta * (1 - mu)
        mu = math.exp(-(a + dn) * L)
    Dn = D - math.log1p(-mu)
    Ln = D - math.log(L) - math.log(mu)   # lam' = delta L mu
    if r + 1 in checkpoints:
        inc = 1 / mu - 1 / mu_prev
        rr = r + 1
        print(f'r={rr:>8}  1/mu-Zuwachs={inc:.6f}  r*mu={rr*mu:.5f}  delta*r^k*ln r={math.exp(-Dn)*rr**k*math.log(rr):.5f}'
              f'  L/ln r={Ln/math.log(rr):.5f} (Ziel 1/a={1/a:.4f})')
    mu_prev = mu
    D, L = Dn, Ln
    r += 1
