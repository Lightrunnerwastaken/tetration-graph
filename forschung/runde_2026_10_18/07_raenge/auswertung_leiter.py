"""Auswertung der Koenigs-Leiter: Plateau-Modell und Asymptotik.
python auswertung_leiter.py daten/koenigs_....json
"""
import sys, json
import mpmath as mp
mp.mp.dps = 30

d = json.load(open(sys.argv[1]))
b = mp.mpf(eval(d['b'], {'sqrt': mp.sqrt, 'e': mp.e, 'exp': mp.exp}))
a = b - 1
k = (2 - b) / (b - 1)
rows = d['rows']
p = {r['r']: mp.mpf(r['p']) for r in rows}
lam = {r['r']: mp.mpf(r['lam']) for r in rows}
prof = {r['r'] + 1: {x: mp.mpf(r[key]) for x, key in
                     [(0.25, 'next_quarter'), (0.5, 'next_half'), (0.75, 'next_3q'), (1.5, 'next_1p5')]}
        for r in rows}
dl = {r: p[r] - b for r in p}
print(f'b={mp.nstr(b,10)} a={mp.nstr(a,6)} k=(2-b)/(b-1)={mp.nstr(k,6)}')
print(' r   delta_r        lam_r        mu_r=1-d_{r+1}/d_r  1/mu_{r+1}-1/mu_r  (Ziel 1/k=%s)  Modellfehler d, lam' % mp.nstr(1 / k, 5))
rs = sorted(p)
mu = {}
for r in rs[:-1]:
    mu[r] = 1 - dl[r + 1] / dl[r]
for r in rs[:-1]:
    L = -mp.log(lam[r])
    # Plateau-Modell: ein Schritt
    m = lam[r] ** (a + dl[r + 1])
    dpred = dl[r] * (1 - m)
    lpred = dl[r] * L * m
    inc = (1 / mu[r + 1] - 1 / mu[r]) if (r + 1) in mu else None
    print(f'{r:2d} {mp.nstr(dl[r],10):>14} {mp.nstr(lam[r],10):>14} {mp.nstr(mu[r],8):>12} '
          f'{mp.nstr(inc,6) if inc is not None else "":>12}     '
          f'{mp.nstr(dpred/dl[r+1]-1,3):>10} {mp.nstr(lpred/lam[r+1]-1,3):>10}')
print('\nProfil H_{r+1}(x) gegen Grenzprofil min(x+1,b) und Plateau-Formel')
print(' r+1  x   H        min(x+1,b)   Plateau-Formel 1-(-ln(1-rho(1-lam_r^x)))/L_{r-1}')
for r in rs[1:-1]:
    rho = dl[r] / dl[r - 1]
    Lm = -mp.log(lam[r - 1])
    line = []
    for x in (0.25, 0.5, 0.75):
        g = prof[r + 1][x]
        form = 1 - mp.log(1 - rho * (1 - lam[r] ** x)) / Lm
        line.append(f'x={x}: {mp.nstr(g,8)} / {mp.nstr(min(x+1,b),6)} / {mp.nstr(form,8)}')
    print(r + 1, ' | '.join(line))
