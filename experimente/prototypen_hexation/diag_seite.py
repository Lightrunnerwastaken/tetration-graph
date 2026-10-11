"""Auf welcher Seite der Singularität w_s (Linie Im w = π/ln λ) muss der Weg H₁ = P(H₀) vorbei, damit P(H₁) oben gegen q geht?
 H₁ als Polygonzug (Schritt 0.01 in y), bei der Kreuzung ein Umweg links bzw. rechts von w_s."""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
from blatt import BlattPentation, verfolge
from stufen import psi_wert
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

BP = BlattPentation("2")
Pe = BP.Pe
q = mp.mpc("-1.992948435932322332933652", "1.445128496639205463533289")
st = json.load(open(HIER + "/prototypen_hexation/start_verfolgt_12.json"))
sy = np.array(st["y"]); sH = np.array([complex(a, b) for a, b in st["H"]])
def h0(y):
    return complex(np.interp(y, sy, sH.real), np.interp(y, sy, sH.imag))
ys = [mp.mpf(k) / 100 for k in range(0, 1201)]
H1, s = BP.P_pfad([mp.mpc(h0(float(t))) for t in ys])
log("H₁ berechnet, Sprung", mp.nstr(max(s), 3), " oben", mp.nstr(H1[-1], 10))
im_linie = mp.pi / mp.re(Pe.lnlam)
# Kreuzung von H₁ mit der Linie
i0 = next(i for i in range(len(H1) - 1) if (mp.im(H1[i]) - im_linie) * (mp.im(H1[i + 1]) - im_linie) <= 0)
a, b = H1[i0], H1[i0 + 1]
t = (im_linie - mp.im(a)) / (mp.im(b) - mp.im(a))
kreuz = a + t * (b - a)
log(f"Kreuzung bei y ≈ {mp.nstr(ys[i0], 6)}: w = {mp.nstr(kreuz, 12)} (Linie Im = {mp.nstr(im_linie, 12)})")


# Singularität auf der Linie: Stufe j der reellen Kette = −2; Suche über Re w
def kette_min(x):
    w = mp.mpc(x, im_linie)
    v = psi_wert(Pe.L, Pe.psiT, Pe.c * Pe.lam ** (w - BP.M))
    werte = [v]
    for _ in range(BP.M):
        if mp.re(v) <= -2:
            break
        v = BP.T_haupt([mp.re(v)])[0]
        werte.append(v)
    return werte
for x in (mp.mpf("-0.3"), mp.mpf("-0.1"), mp.mpf(0), mp.re(kreuz), mp.mpf("0.1"), mp.mpf("0.3")):
    log(f"  Kette auf der Linie bei Re w = {mp.nstr(x, 8)}:", [mp.nstr(mp.re(v), 8) for v in kette_min(x)])


def P_polygon(punkte, name):
    pts = [mp.mpc(p) for p in punkte]
    def weg(tt):
        out = []
        for t_ in tt:
            j = min(int(mp.floor(t_)), len(pts) - 2)
            out.append(pts[j] + (t_ - j) * (pts[j + 1] - pts[j]))
        return BP.P_pfad(out)
    p, w, sm = verfolge([mp.mpf(j) for j in range(len(pts))], weg, runden=14)
    log(f"{name}: {len(p)} Punkte, Sprung {mp.nstr(sm, 3)}, oben P(H₁) = {mp.nstr(w[-1], 12)}, Abstand zu q {mp.nstr(abs(w[-1] - q), 3)}")
    return w


vor, nach = H1[:i0 - 4], H1[i0 + 6:]
for name, dx in (("links vorbei (Re w_s − 0.05)", mp.mpf("-0.05")), ("rechts vorbei (Re w_s + 0.05)", mp.mpf("0.05"))):
    umweg = [mp.mpc(mp.re(kreuz) + dx, im_linie)]
    P_polygon(vor + umweg + nach, name)
