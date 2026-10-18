"""Frage 3: Basisabhaengigkeit der Kneser-Welle.  a_n(b) fuer mehrere Basen (Motor, dps 30, Linie y=0.1, N=256)
und A_b = sexp'(0); Vergleich mit dem fuehrenden Term T_b(n) = -h_b'((c_b - ln n)/ln b)/(n ln b).
Aufruf: python d2_basen.py 1.5 1.7 2 exp(1) 4 10 100"""
import json, sys, time
from pathlib import Path
import mpmath as mp
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parents[1] / "identitaeten_2026_10_14" / "etappe2"))
from kern import Fix
import motor
DPS, Y, N, NMAX = 30, mp.mpf("0.1"), 256, 64
mp.mp.dps = DPS + 15

def fft(a):
    n = len(a)
    if n == 1: return a[:]
    ev, od = fft(a[0::2]), fft(a[1::2]); out = [None] * n
    for k in range(n // 2):
        tw = mp.expjpi(-2 * mp.mpf(k) / n) * od[k]; out[k] = ev[k] + tw; out[k + n // 2] = ev[k] - tw
    return out

def hprime(fx, v):
    acc = mp.mpc(1)
    while abs(v) > 30:
        acc /= v * fx.lb; v = mp.log(v) / fx.lb
    return acc * mp.diff(fx.h, v)

res = {}
for B in sys.argv[1:]:
    t0 = time.time()
    b = mp.e if B == "exp(1)" else mp.mpf(B)
    fx = Fix(b, K=40); per = 2j * mp.pi / mp.log(fx.lam)
    ts = [mp.mpf(j) / N - mp.mpf(1) / 2 for j in range(N)]
    zs = [t + 1j * Y for t in ts]
    K = 32; rho = mp.mpf("0.1"); pts = [rho * mp.expjpi(2 * mp.mpf(j) / K) for j in range(K)]
    vals = motor.werte(B, [f"sexp({motor.gpzahl(z, DPS + 10)})" for z in zs + pts], DPS)
    ws, cv = vals[:N], vals[N:]
    A = mp.re(mp.fsum(v * mp.expjpi(-2 * mp.mpf(j) / K) for j, v in enumerate(cv)) / K / rho)
    th = [fx.h(w) - z for w, z in zip(ws, zs)]
    for j in range(1, N):
        m = min(range(-5, 6), key=lambda m: abs(th[j] + m * per - th[j - 1])); th[j] += m * per
    rand = min(abs(th[0] + m * per - th[-1] - 0) for m in range(-5, 6))
    F = fft(th)
    a = [(-1) ** n * F[n] / N * mp.exp(2 * mp.pi * n * Y) for n in range(NMAX + 1)]
    c = mp.log(A / (2 * mp.pi * fx.lb)) + 1j * mp.pi / 2
    zeile = {"L": [float(fx.L.real), float(fx.L.imag)], "A": mp.nstr(A, 20), "a0": [float(a[0].real), float(a[0].imag)], "a": {}}
    print(f"Basis {B}: L={mp.nstr(fx.L,8)}, A={mp.nstr(A,12)}, a0={mp.nstr(a[0],10)}, Randsprung {mp.nstr(rand,2)}, {time.time()-t0:.0f}s")
    for n in [1, 2, 4, 8, 16, 32, 64]:
        Tn = -hprime(fx, (c - mp.log(n)) / fx.lb) / fx.lb / n
        zeile["a"][n] = [float(abs(a[n])), float(mp.arg(a[n])), float(abs(Tn)), float(mp.arg(Tn))]
        print(f"   n={n:3d} |a_n|={mp.nstr(abs(a[n]),8):>12} arg={mp.nstr(mp.arg(a[n]),6):>9}  fuehrend |T|={mp.nstr(abs(Tn),8):>12} arg={mp.nstr(mp.arg(Tn),6)}")
    res[B] = zeile
    fn = HIER / "daten" / "d2_basen.json"
    alt = json.loads(fn.read_text()) if fn.exists() else {}
    alt.update(res); fn.write_text(json.dumps(alt, indent=1))
