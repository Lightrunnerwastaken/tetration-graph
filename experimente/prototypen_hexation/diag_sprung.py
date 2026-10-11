"""Wo genau springt P(P(H₀)) bei y ≈ 2.8453? Weg 0 → 2.86 (fein um die Stelle), Kette der zweiten Pentation Stufe für Stufe."""
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
st = json.load(open(HIER + "/prototypen_hexation/start_verfolgt_12.json"))
sy = np.array(st["y"]); sH = np.array([complex(a, b) for a, b in st["H"]])
def h0(y):
    return complex(np.interp(y, sy, sH.real), np.interp(y, sy, sH.imag))
ys = [mp.mpf(k) / 20 for k in range(0, 57)] + [mp.mpf("2.8") + mp.mpf(k) / 4000 for k in range(1, 241)]
ys = sorted(set(ys))
w1, s1 = BP.P_pfad([mp.mpc(h0(float(t))) for t in ys])
log("Stufe 1 max Sprung", mp.nstr(max(s1), 3))
Pe = BP.Pe
werte = [psi_wert(Pe.L, Pe.psiT, Pe.c * Pe.lam ** (w - BP.M)) for w in w1]
j = min(range(len(ys)), key=lambda i: abs(ys[i] - mp.mpf("2.8453")))
for stufe in range(BP.M):
    werte, s = BP.T_pfad(werte)
    i = max(range(len(s)), key=lambda t: s[t])
    log(f"T-Stufe {stufe + 1}: max Sprung {mp.nstr(s[i], 3)} bei y = {mp.nstr(ys[i], 8)}; Werte dort",
        [mp.nstr(werte[t], 6) for t in range(max(0, i - 2), min(len(ys), i + 3))])
    log("      um y = 2.8453:", [mp.nstr(werte[t], 6) for t in range(j - 2, j + 3)])
log("P(H₀) um die Stelle:", [(mp.nstr(ys[t], 7), mp.nstr(w1[t], 10)) for t in range(j - 3, j + 4)])
