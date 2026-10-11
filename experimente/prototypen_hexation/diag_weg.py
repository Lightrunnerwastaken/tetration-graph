"""Wo löst sich der Weg y ↦ H₀(y) (Startwert) nicht auf? Stufen H₀ → P(H₀) → P(P(H₀)) einzeln, mit Protokoll."""
import sys, time, json
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import numpy as np
import mpmath as mp
from blatt import BlattPentation, verfolge
mp.mp.dps = 40
t0 = time.time()
def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)

BP = BlattPentation("2")
st = json.load(open(HIER + "/prototypen_hexation/start_verfolgt_12.json"))
sy = np.array(st["y"]); sH = np.array([complex(a, b) for a, b in st["H"]])
def h0(y):
    return complex(np.interp(y, sy, sH.real), np.interp(y, sy, sH.imag))
ys = [mp.mpf(k) / 20 for k in range(0, 241)]
log("Stufe 1: P(H₀)")
p1, w1, s1 = verfolge(ys, lambda tt: BP.P_pfad([mp.mpc(h0(float(t))) for t in tt]), log=log)
log("Stufe 2: P(P(H₀))")
def zwei(tt):
    a, sa = BP.P_pfad([mp.mpc(h0(float(t))) for t in tt])
    b, sb = BP.P_pfad(a)
    return b, [max(x, y) for x, y in zip(sa, sb)]
p2, w2, s2 = verfolge(ys, zwei, log=log)
for name, p, w in (("P(H₀)", p1, w1), ("P(P(H₀))", p2, w2)):
    log(name, [(mp.nstr(p[i], 3), mp.nstr(w[i], 6)) for i in range(0, len(p), max(1, len(p) // 16))])
