"""Rang 8 (Basis 2): dieselbe Pruefung wie s1_theta6 fuer die Kneser-Oktation des Motors (hyper(8,.)) am abstossenden
komplexen Fixpunkt kp8 von H7.  Vorhersage: theta8 singulaer genau bei x8* + Z, x8* = ihyper(8, p6), p6 = reeller
Fixpunkt von H6 (Wert, den H7 fuer Re z -> -inf annimmt);  arg(a_n/a_{n-1}) -> -2 pi x8*.
Aufruf: python s6_theta8.py DPS Y N"""
import sys, json, time
from pathlib import Path
import mpmath as mp
import motor
from lib_r import RangR
DPS, Y, N = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
t0 = time.time()
r0 = motor.roh("2", ["hyper(8,0.5)", "kp", "klam", "hyp[7]", "hylam[7]", "hyperdigits(8)"], DPS)
print("H8(0.5) =", r0[0], " kp8 =", r0[1], " klam =", r0[2], " p6 =", r0[3], " lam6 =", r0[4], " Stellen", r0[5], f"{time.time()-t0:.0f}s", flush=True)
p6 = r0[3]
x8 = motor.roh("2", [f"ihyper(8,{motor.gpzahl(p6, DPS)})"], DPS)[0]
print("x8* = ihyper(8, p6) =", x8, flush=True)
R = RangR("2", DPS, r0[1], 8, rho=mp.mpf(sys.argv[4]) if len(sys.argv) > 4 else mp.mpf("0.2"), rad=mp.mpf("0.02"))
mp.mp.dps = DPS + 10
print("kp8 =", R.L, " mu =", R.mu, abs(R.mu), flush=True)
y = mp.mpf(Y); per = 2j * mp.pi / R.logmu
xs = [mp.mpf(j) / N for j in range(N)]; zs = [x + 1j * y for x in xs]
ws = motor.werte("2", [f"iferr(hyper(8,{motor.gpzahl(z, DPS + 10)}),E,0)" for z in zs], DPS)
if any(w == 0 for w in ws): sys.exit("hyper(8) scheitert bei x = %s" % [float(xs[j]) for j, w in enumerate(ws) if w == 0][:8])
hs, ns = R.h6(ws)
th = [h - z for h, z in zip(hs, zs)]
for j in range(1, N):
    m = min(range(-5, 6), key=lambda m: abs(th[j] + m * per - th[j - 1])); th[j] += m * per
print("max Schritt", mp.nstr(max(abs(th[j] - th[j - 1]) for j in range(1, N)), 3), " Randschritt", mp.nstr(abs(th[0] - th[-1]), 3), flush=True)
a = [mp.fsum(th[j] * mp.expjpi(-2 * mp.mpf(j * n) / N) for j in range(N)) / N * mp.exp(2 * mp.pi * n * y) for n in range(-4, N // 2)]
neg, a = a[:4], a[4:]
print("negative Moden:", [mp.nstr(abs(v), 3) for v in neg])
x0 = mp.re(x8) % 1
zeilen = []
for n in range(1, min(N // 2, 25)):
    ph = mp.arg(a[n] / a[n - 1]) if n > 1 else mp.mpf(0)
    print(f"n={n:3d} |a_n|={mp.nstr(abs(a[n]),8):>12}  arg(a_n/a_n-1)={mp.nstr(ph,7):>10}   -2 pi x8* = {mp.nstr(-2*mp.pi*x0 + (2*mp.pi if -2*mp.pi*x0 < -mp.pi else 0),7)}")
    zeilen.append([n, float(abs(a[n])), float(ph)])
json.dump({"x8": str(x8), "kp8": str(R.L), "mu": str(R.mu), "a0": str(a[0]), "zeilen": zeilen}, open(Path(__file__).resolve().parent / "daten" / f"theta8_{DPS}_{Y}_{N}.json", "w"), indent=1)
print(f"[{time.time()-t0:.0f}s] fertig")
