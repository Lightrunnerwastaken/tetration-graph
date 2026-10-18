"""Rang 6: Fourier-Koeffizienten von theta6 fuer die Kneser-Hexation des Motors (hexx, Basis 2),
bezogen auf die regulaere Loesung am abstossenden komplexen Fixpunkt kp der Pentation (derselbe, den der Motor nutzt).
theta6(z) = h(hex z) - z,  h(w) = n + log(sigma(ipent^n w))/log mu.   a_n = e^{2 pi n y} (1/N) sum theta(x_j+iy) e^{-2 pi i n x_j}.
Aufruf: python s1_theta6.py DPS Y N   -> daten/theta6_<dps>_<y>_<N>.json"""
import json, sys, time
from pathlib import Path
import mpmath as mp
import motor
from lib6 import Rang6
DPS, Y, N = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
t0 = time.time()
R = Rang6("2", DPS, start=mp.mpc(-1.99294843593232, 1.44512849663921), rho=mp.mpf("0.2"), rad=mp.mpf("0.04"))
mp.mp.dps = DPS + 10
y = mp.mpf(Y)
per = 2j * mp.pi / R.logmu
xs = [mp.mpf(j) / N for j in range(N)]
zs = [x + 1j * y for x in xs]
ws = motor.werte("2", [f"iferr(hexx({motor.gpzahl(z, DPS + 10)}),E,0)" for z in zs], DPS)
fehl = [j for j, w in enumerate(ws) if w == 0]
if fehl: sys.exit("hexx scheitert an %d Punkten, z.B. x=%s" % (len(fehl), [float(xs[j]) for j in fehl[:5]]))
print(f"[{time.time()-t0:.0f}s] hex fertig", flush=True)
hs, ns = R.h6(ws)
th = [h - z for h, z in zip(hs, zs)]
spr = 0
for j in range(1, N):
    m = min(range(-5, 6), key=lambda m: abs(th[j] + m * per - th[j - 1]))
    if m: spr += 1
    th[j] += m * per
schritt = max(abs(th[j] - th[j - 1]) for j in range(1, N))
rand = min(abs(th[0] + m * per - th[-1]) for m in range(-5, 6))
print(f"[{time.time()-t0:.0f}s] theta fertig: Zweigkorrekturen {spr}, max Schritt {mp.nstr(schritt,3)}, Randschritt {mp.nstr(rand,3)}", flush=True)

def fft(a):
    n = len(a)
    if n == 1: return a[:]
    ev, od = fft(a[0::2]), fft(a[1::2]); out = [None] * n
    for k in range(n // 2):
        tw = mp.expjpi(-2 * mp.mpf(k) / n) * od[k]; out[k] = ev[k] + tw; out[k + n // 2] = ev[k] - tw
    return out
F = fft(th)
a = [F[n] / N * mp.exp(2 * mp.pi * n * y) for n in range(N // 2)]
am = [F[N - n] / N * mp.exp(-2 * mp.pi * n * y) for n in range(1, 9)]   # negative Moden a_{-n}
out = {"dps": DPS, "y": Y, "N": N, "kp": [mp.nstr(R.L.real, DPS), mp.nstr(R.L.imag, DPS)],
       "mu": [mp.nstr(R.mu.real, DPS), mp.nstr(R.mu.imag, DPS)],
       "a": [[mp.nstr(v.real, DPS), mp.nstr(v.imag, DPS)] for v in a],
       "a_neg": [[mp.nstr(v.real, 8), mp.nstr(v.imag, 8)] for v in am], "iter_n": [min(ns), max(ns)]}
fn = Path(__file__).resolve().parent / "daten" / f"theta6_{DPS}_{Y}_{N}.json"
fn.write_text(json.dumps(out, indent=0))
print("a0 =", mp.nstr(a[0], 20))
print("negative Moden |a_-n| (n=1..8):", [mp.nstr(abs(v), 3) for v in am])
x0 = mp.mpf("-2.765022050525603451110465872431272359") % 1
for n in [1, 2, 3, 4, 8, 16, 32, 64, 128, 256]:
    if n < N // 2:
        ph = mp.arg(a[n] / a[n - 1]) if n > 1 else 0
        print(f"n={n:4d} |a_n|={mp.nstr(abs(a[n]),10):>14}  arg a_n={mp.nstr(mp.arg(a[n]),8):>11}  arg(a_n/a_n-1)={mp.nstr(ph,8):>11}  (-2pi x*={mp.nstr(-2*mp.pi*x0 + 2*mp.pi*0,8)})")
print(f"[{time.time()-t0:.0f}s] fertig")
