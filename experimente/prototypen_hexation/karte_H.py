"""Phasenbild der Kneser-artigen Hexation H(z) = F_q(x* + z + θ(z)) (θ aus hex_lauf_2.json, K = 64) über dem Streifen:
 Vermutung: Singularitäten bei s_b + 2πi/ln μ_q + m (m ≥ 0), also etwa 2.94 über der reellen Achse."""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hexq import HexQ
mp.mp.dps = 30
t0 = time.time()
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
K = int(sys.argv[1]) if len(sys.argv) > 1 else 64
alt = json.load(open(HIER + "/prototypen_hexation/hex_lauf_2.json"))
g = [mp.mpf(v) for v in alt[str(K)]["g"]]
def dft(f):
    M = len(f)
    return [mp.fsum(f[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M for n in range(M)]
gh = dft(g)
def g_at(x):
    return mp.re(gh[0] + sum(2 * mp.re(gh[n] * mp.expjpi(2 * n * x)) for n in range(1, K // 2)) + mp.re(gh[K // 2] * mp.cos(mp.pi * K * x)))
xs = [mp.mpf(j) / K for j in range(K)]
u = [mp.mpf(0)] * K
for it in range(300):
    v = [g_at(x + uu) for x, uu in zip(xs, u)]
    vh = dft(v)
    a = [1j * vh[0]] + [2j * vh[n] for n in range(1, K // 2)]
    un = [mp.re(mp.fsum(a[n] * mp.expjpi(2 * n * x) for n in range(len(a)))) for x in xs]
    if max(abs(p - r) for p, r in zip(u, un)) < mp.mpf(10) ** -22:
        break
    u = un
theta = lambda z: mp.fsum(a[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(a)))
XS = mp.mpf(alt[str(K)]["x_stern"])
periode = 2j * mp.pi / HQ.lnmu
print("Periode 2πi/ln μ =", mp.nstr(periode, 10), " θ̂₀ =", mp.nstr(a[0], 8), flush=True)
# Kontrolle: H(0) = 1, H(1) = 2
print("H(0), H(1):", [mp.nstr(h, 10) for h in HQ.F([XS + z + theta(XS + z) for z in (mp.mpf(0), mp.mpf(1))])], flush=True)
# s_b: H(t0) = L_T auf (−3, −2), s_b = t0 − 1
def H(zs):
    return HQ.F([XS + z + theta(XS + z) for z in zs])
ts = [mp.mpf("-3") + mp.mpf(j) / 100 for j in range(101)]
hs = [mp.re(h) - HQ.L for h in H(ts)]
j = next(j for j in range(100) if hs[j] < 0 <= hs[j + 1])
a0, a1, f0, f1 = ts[j], ts[j + 1], hs[j], hs[j + 1]
for _ in range(8):
    t = a1 - f1 * (a1 - a0) / (f1 - f0)
    a0, f0, a1, f1 = a1, f1, t, mp.re(H([t])[0]) - HQ.L
s_b = a1 - 1
print("s_b =", mp.nstr(s_b, 10), " vorhergesagte Singularitäten:", [mp.nstr(s_b + periode + m, 6) for m in range(4)], f"[{time.time() - t0:.0f}s]", flush=True)
gx = np.arange(-1.5, 2.5001, 0.04)
gy = np.arange(1.6, 4.4001, 0.04)
pkt = [mp.mpc(x, y) for y in gy for x in gx]
Z = np.array([complex(h) for h in H(pkt)]).reshape(len(gy), len(gx))
np.save(HIER + f"/prototypen_hexation/karte_H_{K}.npy", Z)
print(f"H auf {len(pkt)} Punkten [{time.time() - t0:.0f}s]", flush=True)
fig, ax = plt.subplots(1, 2, figsize=(16, 6.5))
for a_, f, titel, cm in ((ax[0], np.angle(Z), "arg H(z)", "hsv"), (ax[1], np.log10(np.abs(Z) + 1e-300).clip(-1, 3), "log10 |H(z)|", "viridis")):
    im = a_.imshow(f, origin="lower", extent=[gx[0], gx[-1], gy[0], gy[-1]], cmap=cm, aspect="equal")
    for m in range(4):
        zz = complex(s_b + periode + m)
        a_.plot([zz.real], [zz.imag], "kx", ms=10, mew=2)
    a_.set_title(titel + f"  (K = {K}; x = vorhergesagt s_b + 2πi/ln μ + m)")
    plt.colorbar(im, ax=a_, shrink=0.7)
plt.tight_layout()
plt.savefig(HIER + f"/prototypen_hexation/karte_H_{K}.png", dpi=90)
print("Bild gespeichert", flush=True)
