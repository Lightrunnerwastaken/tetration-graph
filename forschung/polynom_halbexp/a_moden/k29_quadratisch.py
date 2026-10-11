"""Strang A, Runde 2: nichtlineare Korrektur Δ_N = p_N − Π_N f_K (p_N aus k28, Newton von Π_N f_K aus) für N bis 1000.

Grössen: max_{|z|=|L|} |Δ_N| (FFT + feines Gitter um arg L) und Ort, |Δ_N(L)|, max_{|z|=1}|Δ_N|, Verhältnis zu |E_N(L)|²
(E_N = Π_N f_K − f_K), Modenamplitude von p_N gegen die von Π_N f_K, p_N(1) − T(½).
"""
import json
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

from lin_kern import DATEN, lade  # noqa: E402

mp.mp.dps = 60
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
argL = float(mp.arg(L))
P = lade()
pn = json.loads((DATEN / "pN_gross.json").read_text(encoding="utf-8"))
A = json.loads((DATEN / "lin_auswertung.json").read_text(encoding="utf-8"))
ld = json.loads((DATEN / "linie48.json").read_text(encoding="utf-8"))
ws = [mp.mpf(x) for x in ld["w"]]
zs = [mp.mpf(x) for x in ld["z"]]
Tp = [mp.mpf(x) for x in ld["Tp"]]
T12 = mp.mpf("1.6463542337511945809719240315921145182")
cx = lambda p: mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))  # noqa: E731


def hv(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


print(" N    max_Rand|Δ|  Ort(θ−argL)·N   |Δ(L)|      max_{|z|=1}|Δ|  |Δ(L)|/|E(L)|²  rel. Δc/c      p_N(1)−T(½)   Newton-Korr.")
for N in sorted(int(k) for k in pn):
    p = [mp.mpf(x) for x in pn[str(N)]["p"]]
    b = [mp.mpf(x) for x in P[N]["b"]]
    d = [x - y for x, y in zip(p, b)]
    e = np.array([float(x * RL ** k) for k, x in enumerate(d)])
    Mg = 1 << int(np.ceil(np.log2(16 * (N + 1))))
    v1 = np.fft.ifft(e, Mg) * Mg
    th2 = argL + np.linspace(-40.0 / N, 40.0 / N, 4001)
    v2 = np.polyval(e[::-1], np.exp(1j * th2))
    th = np.concatenate([2 * np.pi * np.arange(Mg) / Mg, th2])
    vv = np.concatenate([v1, v2])
    i = int(np.argmax(np.abs(vv)))
    tt = (th[i] + np.pi) % (2 * np.pi) - np.pi
    e1 = np.array([float(x) for x in d])
    u1 = np.fft.ifft(e1, Mg) * Mg
    dL = abs(hv(d, L))
    EL = abs(cx(A[str(N)]["PL"]))
    th_p = [(hv(p, z) - hv(b, z)) / t for z, t in zip(zs, Tp)]
    dc = mp.fsum(t * mp.expjpi(2 * w) for t, w in zip(th_p, ws)) / len(ws)
    c = cx(A[str(N)]["c"])
    print(f"{N:4d}  {abs(vv[i]):10.3e}  {(abs(tt) - argL) * N:10.3f}  {mp.nstr(dL, 4):>10}  {np.max(np.abs(u1)):12.3e}  "
          f"{mp.nstr(dL / EL ** 2, 4):>12}  {mp.nstr(abs(dc) / abs(c), 3):>10}  {mp.nstr(mp.fsum(p) - T12, 5):>12}  {pn[str(N)]['korr']}")
