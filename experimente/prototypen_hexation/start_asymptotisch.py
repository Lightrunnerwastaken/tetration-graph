"""Startwert ohne F_q-Kette oberhalb y = 1.5: grobe Lösung (Hauptzweig, dort stetig) für y ≤ 1.5, darüber q + C·μ^{i(y−1.5)} (glatt angesetzt)."""
import sys, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
from hexq import HexQ
mp.mp.dps = 40
A, dy = float(sys.argv[1]), float(sys.argv[2])
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
alt = json.load(open(HIER + "/prototypen_hexation/hex_lauf_2.json"))
g = [mp.mpf(v) for v in alt["32"]["g"]]
K = 32
def dft(f):
    M = len(f)
    return [mp.fsum(f[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M for n in range(M)]
gh = dft(g)
def g_at(x):
    return mp.re(gh[0] + sum(2 * mp.re(gh[n] * mp.expjpi(2 * n * x)) for n in range(1, K // 2)) + mp.re(gh[K // 2] * mp.cos(mp.pi * K * x)))
xs = [mp.mpf(j) / K for j in range(K)]
u = [mp.mpf(0)] * K
for it in range(200):
    v = [g_at(x + uu) for x, uu in zip(xs, u)]
    vh = dft(v)
    a = [1j * vh[0]] + [2j * vh[n] for n in range(1, K // 2)]
    un = [mp.re(mp.fsum(a[n] * mp.expjpi(2 * n * x) for n in range(len(a)))) for x in xs]
    if max(abs(p - r) for p, r in zip(u, un)) < mp.mpf(10) ** -22:
        break
    u = un
theta = lambda z: mp.fsum(a[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(a)))
XS = mp.mpf(alt["32"]["x_stern"])
y0 = 1.5
ys_u = [mp.mpf(k) * mp.mpf(dy) for k in range(int(y0 / dy) + 1)]
Hu = HQ.F([XS + 1j * y + theta(XS + 1j * y) for y in ys_u])
q = complex(HQ.q)
lnmu = complex(HQ.lnmu)
C = complex(Hu[-1]) - q
ys_o = np.arange(y0 + dy, A + 1e-9, dy)
Ho = q + C * np.exp(-1.068 * (ys_o - y0))          # ohne Drehung: Weg direkt nach q (Spirale umläuft einen Verzweigungspunkt)
y_all = [float(y) for y in ys_u] + list(ys_o)
H_all = [complex(h) for h in Hu] + list(Ho)
print("Übergang bei y = 1.5:", complex(Hu[-1]), " Ansatz oben bis", H_all[-1], " q =", q)
json.dump({"A": A, "dy": dy, "y": y_all, "H": [[h.real, h.imag] for h in H_all]}, open(HIER + f"/prototypen_hexation/start_verfolgt_{int(A)}.json", "w"))
print("gespeichert")
