"""Fourier-Koeffizienten a_n von θ bis n = nmax, Abtastung auf Im z = y mit N Punkten (N Zweierpotenz).
a_n = e^{2πny}·(1/N)Σ_j θ(t_j+iy) e^{−2πi n t_j}.  Fehler durch Überlagerung ≈ |a_{n+N}| e^{−2πNy}·…, Rundung ≈ 10^{−dps}·e^{2πny}.
Aufruf: python p2_koeff.py <basis: e|2> <dps> <y> <N> <nmax>
Ausgabe: daten/a_<basis>_<dps>_<y>_<N>.json
"""
import json
import sys
import time
from pathlib import Path

import mpmath as mp
from kern import Fix, sexp

B, DPS, Y, N, NMAX = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
mp.mp.dps = DPS + 20
basis = "exp(1)" if B == "e" else B
fx = Fix(mp.e if B == "e" else mp.mpf(B), K=40)
per = 2j * mp.pi / mp.log(fx.lam)
y = mp.mpf(Y)
t0 = time.time()


def fft(a):
    n = len(a)
    if n == 1:
        return a[:]
    ev, od = fft(a[0::2]), fft(a[1::2])
    out = [None] * n
    for k in range(n // 2):
        tw = mp.expjpi(-2 * mp.mpf(k) / n) * od[k]
        out[k] = ev[k] + tw
        out[k + n // 2] = ev[k] - tw
    return out


ts = [mp.mpf(j) / N - mp.mpf(1) / 2 for j in range(N)]
zs = [t + 1j * y for t in ts]
ws = sexp(basis, zs, DPS)
print(f"[{time.time()-t0:.0f}s] sexp fertig", flush=True)
th = [fx.h(w) - z for w, z in zip(ws, zs)]
for j in range(1, N):
    m = min(range(-5, 6), key=lambda m: abs(th[j] + m * per - th[j - 1]))
    th[j] += m * per
schritt = max(abs(th[j] - th[j - 1]) for j in range(1, N))
rand = min(abs(th[0] + m * per - th[-1]) for m in range(-5, 6))
print(f"[{time.time()-t0:.0f}s] theta fertig, max Schritt {mp.nstr(schritt,3)}, Randsprung {mp.nstr(rand,3)}", flush=True)
# t_j = j/N − ½  ⇒  e^{−2πin t_j} = (−1)^n e^{−2πi nj/N}
F = fft(th)
a = []
for n in range(NMAX + 1):
    a.append((-1) ** n * F[n] / N * mp.exp(2 * mp.pi * n * y))
out = {"basis": B, "dps": DPS, "y": Y, "N": N, "L": mp.nstr(fx.L, 40),
       "a": [[mp.nstr(mp.re(v), DPS), mp.nstr(mp.im(v), DPS)] for v in a]}
Path("daten").mkdir(exist_ok=True)
fn = Path("daten") / f"a_{B}_{DPS}_{Y}_{N}.json"
fn.write_text(json.dumps(out, indent=0), encoding="utf-8")
print(f"[{time.time()-t0:.0f}s] geschrieben {fn}; a_1 = {mp.nstr(a[1], 25)}", flush=True)
for n in (1, 2, 64, 256, 512, 1024, 2048):
    if n <= NMAX:
        print(f"  a_{n} = {mp.nstr(a[n], 20)}")
