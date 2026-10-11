"""Kneser-artige Hexation zur Basis 2 mit Kouznetsovs Cauchy-Integral-Methode.

 H holomorph im Streifen −1 ≤ Re z ≤ 1, H(z) → q für Im z → +∞, H(z̄) = conj H(z).
 Auf der imaginären Achse: H(iy) = (1/2π)∫ [P(H(iy'))/(1 + iy' − iy) − P^{-1}(H(iy'))/(−1 + iy' − iy)] dy' + Deckel/Boden mit H ≈ q, q̄.
 Iteration mit Dämpfung; Startwert aus der groben θ-Lösung (K = 32). Aufruf: hex_kouz.py A N iter
"""
import sys, time, json
sys.path.insert(0, r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad")
sys.path.insert(0, r"D:/Users/janis/Documents/Matura/tetration-graph/experimente")
import numpy as np
import mpmath as mp
from hexq import HexQ
from hexation import Pentation
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

A, N, ITER = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
Pe = Pentation("2")
q = HQ.q
log("q", mp.nstr(q, 20))

# Gauss–Legendre auf [−A, A]; symmetrisch, wir rechnen nur mit y > 0 (N gerade)
xg, wg = np.polynomial.legendre.leggauss(N)
y = A * xg
wy = A * wg
pos = np.where(y > 0)[0]
yp = y[pos]

# Startwert: grobe θ-Lösung (K = 32) aus hex_lauf_2.json neu aufbauen (Theodorsen, ohne neue F-Auswertungen)
alt = json.load(open(r"C:/Users/janis/AppData/Local/Temp/claude/D--Users-janis-Documents-Matura/16e6767d-2c78-46bf-8ec0-6069646a66c1/scratchpad/hex_lauf_2.json"))
g = [mp.mpf(v) for v in alt["32"]["g"]]
K = 32
def dft(f):
    M = len(f)
    return [mp.fsum(f[j] * mp.expjpi(-2 * n * mp.mpf(j) / M) for j in range(M)) / M for n in range(M)]
gh = dft(g)
def g_at(x):
    s = gh[0]
    for n in range(1, K // 2):
        s += 2 * mp.re(gh[n] * mp.expjpi(2 * n * x))
    s += mp.re(gh[K // 2] * mp.cos(mp.pi * K * x))
    return mp.re(s)
xs = [mp.mpf(j) / K for j in range(K)]
u = [mp.mpf(0)] * K
for it in range(200):
    v = [g_at(x + uu) for x, uu in zip(xs, u)]
    vh = dft(v)
    a = [1j * vh[0]] + [2j * vh[n] for n in range(1, K // 2)]
    un = [mp.re(mp.fsum(a[n] * mp.expjpi(2 * n * x) for n in range(len(a)))) for x in xs]
    if max(abs(p - r) for p, r in zip(u, un)) < mp.mpf(10) ** -22:
        u = un
        break
    u = un
def theta(z):
    return mp.fsum(a[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(a)))
XS = mp.mpf(alt["32"]["x_stern"])
Hy = HQ.F([XS + 1j * mp.mpf(float(yy)) + theta(XS + 1j * mp.mpf(float(yy))) for yy in yp])
log("Startwert H(iy) bei y = ", [f"{yy:.2f}" for yy in yp[::max(1, len(yp) // 6)]], [mp.nstr(h, 6) for h in Hy[::max(1, len(yp) // 6)]])

qc = complex(q)
# Deckel- und Bodenbeitrag für z = i y_k (H ≈ q oben, q̄ unten)
def deckel(z):
    return qc / (2j * np.pi) * (np.log(-1 + 1j * A - z) - np.log(1 + 1j * A - z)) + np.conj(qc) / (2j * np.pi) * (np.log(1 - 1j * A - z) - np.log(-1 - 1j * A - z))

ziel = 1j * yp
ziel_all = 1j * y
D = deckel(ziel)
omega = 0.5
for it in range(ITER):
    Pf = Pe.P(Hy)
    Pinv = Pe.inv(Hy)
    hy = np.array([complex(h) for h in Hy])
    pf = np.array([complex(v) for v in Pf])
    pi_ = np.array([complex(v) for v in Pinv])
    # Werte auf allen Knoten (unten: konjugiert)
    full_P = np.empty(N, complex); full_I = np.empty(N, complex)
    full_P[pos] = pf; full_I[pos] = pi_
    neg = np.where(y < 0)[0]
    # Zuordnung y_neg = −y_pos
    for j in neg:
        k = pos[np.argmin(np.abs(y[pos] + y[j]))]
        full_P[j] = np.conj(full_P[k]); full_I[j] = np.conj(full_I[k])
    rechts = (wy[None, :] * full_P[None, :] / (1 + 1j * y[None, :] - ziel[:, None])).sum(axis=1) / (2 * np.pi)
    links = (wy[None, :] * full_I[None, :] / (-1 + 1j * y[None, :] - ziel[:, None])).sum(axis=1) / (2 * np.pi)
    neu = rechts - links + D
    d = np.max(np.abs(neu - hy))
    hy2 = (1 - omega) * hy + omega * neu
    Hy = [mp.mpc(complex(h)) for h in hy2]
    log(f"Iteration {it}: max |H_neu − H| = {d:.3e}, H(i·{yp[0]:.3f}) = {hy2[0]:.10f}")
    if d < 1e-11:
        break

# H im Inneren (reelle Achse) per Cauchy-Formel aus dem letzten Randwertsatz
def H_innen(zs):
    out = []
    for z in zs:
        rechts = (wy * full_P / (1 + 1j * y - z)).sum() / (2 * np.pi)
        links = (wy * full_I / (-1 + 1j * y - z)).sum() / (2 * np.pi)
        out.append(rechts - links + deckel(np.array([z]))[0])
    return out
xx = np.linspace(-0.9, 0.9, 19)
hv = H_innen(xx)
log("H auf der reellen Achse:", [(f"{x:.2f}", f"{h.real:.10f}", f"{h.imag:.1e}") for x, h in zip(xx[::3], hv[::3])])
# Normierung H(x0) = 1
from scipy.optimize import brentq
f = lambda x: H_innen([x])[0].real - 1
x0 = brentq(f, -0.9, 0.9, xtol=1e-14)
h05 = H_innen([x0 + 0.5])[0]
h15 = complex(Pe.P([mp.mpc(h05)])[0])
log(f"x0 = {x0:.12f}; 2↑↑↑↑0.5 = {h05.real:.12f} (Im {h05.imag:.1e}); 2↑↑↑↑1.5 = {h15.real:.12f} (Im {h15.imag:.1e})  [Levenstein ≈ 2.6729]")
hm05 = H_innen([x0 - 0.5])[0]
log(f"2↑↑↑↑(−0.5) = {hm05.real:.12f}, H(x0 − 1) = {H_innen([x0 - 0.999999])[0].real:.8f} (soll ≈ 0)")
json.dump({"A": A, "N": N, "x0": x0, "H05": [h05.real, h05.imag], "H15": [h15.real, h15.imag], "Hm05": [hm05.real, hm05.imag],
           "Hy": [[complex(h).real, complex(h).imag] for h in Hy], "yp": list(map(float, yp))},
          open(sys.argv[0].replace(".py", f"_{int(A)}_{N}.json"), "w"), indent=1)
log("fertig")
