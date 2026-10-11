"""Ist die Kneser-artige Hexation auf dem Streifen S = {a ≤ Re z ≤ a + 1} (links der ersten Verzweigung bei Re −0.741) holomorph und injektiv?
 H(z) = F_q(x* + z + θ(x* + z)) mit dem groben θ (hex_lauf_2.json, K = 128), F_q mit Blattverfolgung auf senkrechten Wegen von oben.
 Ränder γ = H(a + iy) und P(γ) = H(a + 1 + iy), y ≥ 0 (unten konjugiert). Prüft: Ränder einfach und disjunkt, w_b im Inneren.
 Aufruf: streifen.py a
"""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from blatt import BlattPentation, verfolge
from hexq import HexQ
from stufen import psi_wert
mp.mp.dps = 30
t0 = time.time()
a = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf("-1.8")
HQ = HexQ("2", mp.mpf("-1.7439"), mp.mpc("-1.99294843593232", "1.44512849663921"))
BP = BlattPentation("2")
K = 128
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
    koef = [1j * vh[0]] + [2j * vh[n] for n in range(1, K // 2)]
    un = [mp.re(mp.fsum(koef[n] * mp.expjpi(2 * n * x) for n in range(len(koef)))) for x in xs]
    if max(abs(p - r) for p, r in zip(u, un)) < mp.mpf(10) ** -22:
        break
    u = un
theta = lambda z: mp.fsum(koef[n] * mp.exp(2j * mp.pi * n * z) for n in range(len(koef)))
XS = mp.mpf(alt[str(K)]["x_stern"])


def H_weg(zs):
    zetas = [XS + z + theta(XS + z) for z in zs]
    n = max(0, max(int(mp.ceil(mp.re((mp.log(abs(mp.exp(c * HQ.lnmu))) - mp.log(HQ.rlok)) / mp.re(HQ.lnmu)))) for c in zetas)) + 1
    w = [psi_wert(HQ.q, HQ.psi, mp.exp((c - n) * HQ.lnmu)) for c in zetas]
    sp = [mp.mpf(0)] * len(zs)
    for _ in range(n):
        w, s_ = BP.P_pfad(w)
        sp = [max(x, y) for x, y in zip(sp, s_)]
    return w, sp


raender = {}
for kante in (a, a + 1):
    ys = [mp.mpf(8) * (mp.mpf("0.0025") / 8) ** (mp.mpf(k) / 300) for k in range(301)]     # 8 … 0.0025 geometrisch
    p, w, sm = verfolge([kante + 1j * y for y in ys], H_weg, runden=8)
    raender[float(kante)] = (np.array([complex(z) for z in p]), np.array([complex(x) for x in w]))
    print(f"[{time.time() - t0:5.0f}s] Rand Re z = {mp.nstr(kante, 4)}: {len(p)} Punkte, Sprung {mp.nstr(sm, 3)}, "
          f"oben H = {mp.nstr(w[0], 8)} (q = {mp.nstr(HQ.q, 8)}), unten H = {mp.nstr(w[-1], 10)}", flush=True)
links, rechts = raender[float(a)][1], raender[float(a + 1)][1]
print("  Kontrolle P(H(a + iy)) = H(a + 1 + iy) oben:", mp.nstr(abs(BP.Pe.P([mp.mpc(links[0])])[0] - rechts[0]), 3))


def schnitte(c1, c2, selbst=False):
    n = 0
    for i in range(len(c1) - 1):
        for j in range(i + 2 if selbst else 0, len(c2) - 1):
            p1, p2, p3, p4 = c1[i], c1[i + 1], c2[j], c2[j + 1]
            d = (p2 - p1).real * (p4 - p3).imag - (p2 - p1).imag * (p4 - p3).real
            if d == 0:
                continue
            t = ((p3 - p1).real * (p4 - p3).imag - (p3 - p1).imag * (p4 - p3).real) / d
            s = ((p3 - p1).real * (p2 - p1).imag - (p3 - p1).imag * (p2 - p1).real) / d
            if 1e-9 < t < 1 - 1e-9 and 1e-9 < s < 1 - 1e-9:
                n += 1
    return n


rand = np.concatenate([links[::-1][:-1], np.conj(links), np.conj(rechts[::-1]), rechts])   # geschlossen über q̄ und q (Enden nahe q)
print(f"  Selbstschnitte γ: {schnitte(links, links, True)}, P(γ): {schnitte(rechts, rechts, True)}, γ mit P(γ): {schnitte(links, rechts)}")
wb = complex(-1.9786, 1.91185)
pfad = rand
wind = np.sum(np.angle(np.roll(pfad, -1)[:-1] - wb) - np.angle(pfad[:-1] - wb))
umlauf = np.sum(np.angle((np.roll(pfad, -1)[:-1] - wb) / (pfad[:-1] - wb))) / (2 * np.pi)
print(f"  Umlaufzahl des Randes um w_b ≈ {umlauf:.3f} (≠ 0: w_b liegt im Bild H(S))")
fig, ax = plt.subplots(figsize=(7, 7))
for c, f, name in ((links, "#3d6fb6", f"γ = H({float(a):.1f} + iy)"), (rechts, "#b5523b", f"P(γ) = H({float(a) + 1:.1f} + iy)")):
    ax.plot(c.real, c.imag, color=f, label=name); ax.plot(c.real, -c.imag, color=f, ls=":")
ax.plot([complex(HQ.q).real], [complex(HQ.q).imag], "ko"); ax.plot([wb.real], [wb.imag], "k^")
ax.set_aspect("equal"); ax.legend(); ax.set_title("Bild des Streifens S unter H (Rand)")
plt.savefig(HIER + "/prototypen_hexation/streifen.png", dpi=90)
print(f"fertig [{time.time() - t0:.0f}s]")
