"""Startwert für die Hexation: grobe θ-Lösung (K = 32), H(iy) = F_q(x* + iy + θ(x* + iy)) — jetzt mit Blattverfolgung entlang y.

 F_q(ζ) = P^n(Ψ_lok(μ^ζ / μ^n)) mit festem n; jede Pentation der Kette wird entlang y fortgesetzt (Start y = 0: Hauptzweig).
 Ergebnis: H(iy) auf y ∈ [0, A] (fein), dazu P(H(iy)) verfolgt (Kontrolle: oben P(H) ≈ H ≈ q?). Aufruf: start_verfolgt.py A dy
"""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import mpmath as mp
from blatt import BlattPentation
from hexq import HexQ
from stufen import psi_wert
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

A, dy = mp.mpf(sys.argv[1]), mp.mpf(sys.argv[2])
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
BP = BlattPentation("2")
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

from blatt import verfolge


def kette(param):
    zet = [XS + 1j * y + theta(XS + 1j * y) for y in param]
    sl = [mp.exp(z * HQ.lnmu) for z in zet]
    w = [psi_wert(HQ.q, HQ.psi, s / HQ.mu ** n) for s in sl]
    sp = [mp.mpf(0)] * len(param)
    for k in range(n):
        w, s = BP.P_pfad(w)
        sp = [max(a, b) for a, b in zip(sp, s)]
    return w, sp


ys0 = [dy * k for k in range(int(A / dy) + 1)]
ss0 = [mp.exp((XS + 1j * y + theta(XS + 1j * y)) * HQ.lnmu) for y in ys0]
n = max(0, max(int(mp.ceil(mp.log(abs(s) / HQ.rlok) / mp.log(abs(HQ.mu)))) for s in ss0))
log(f"{len(ys0)} Wegpunkte am Anfang, Kettenlänge n = {n}")
ys, H, smax = verfolge(ys0, kette, log=log)
zetas = [XS + 1j * y + theta(XS + 1j * y) for y in ys]
log(f"nach Verfeinerung {len(ys)} Punkte, grösster Winkelsprung {mp.nstr(smax, 3)}")
haupt = HQ.F(zetas[::10])
abw = [ys[10 * i] for i in range(len(haupt)) if abs(haupt[i] - H[10 * i]) > 1e-15]
log("Hauptzweig weicht ab bei y ∈", (mp.nstr(abw[0], 4), mp.nstr(abw[-1], 4)) if abw else "—")
log("H(iy) verfolgt:", [(mp.nstr(ys[i], 3), mp.nstr(H[i], 8)) for i in range(0, len(ys), max(1, len(ys) // 12))])
PH, _ = BP.P_pfad(H)
log("P(H(iy)) verfolgt am oberen Ende:", mp.nstr(PH[-1], 12), " H(iA) =", mp.nstr(H[-1], 12), " q =", mp.nstr(HQ.q, 12))
log("|P(H) − H| oben (sollte klein sein, wenn q der richtige Fixpunkt ist):", mp.nstr(abs(PH[-1] - H[-1]), 3), " grösster Sprung:", mp.nstr(smax, 3))
json.dump({"A": float(A), "dy": float(dy), "y": [float(y) for y in ys], "H": [[float(mp.re(h)), float(mp.im(h))] for h in H],
           "PH": [[float(mp.re(h)), float(mp.im(h))] for h in PH], "x_stern": float(XS)},
          open(HIER + f"/prototypen_hexation/start_verfolgt_{int(A)}.json", "w"))
log("fertig")
