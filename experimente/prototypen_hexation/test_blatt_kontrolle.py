"""Kontrolle von Test (3): P entlang 1 → 1 + 3i auf drei Arten — Vorgabe, eine Stufe mehr Verschiebung, ohne analytischen Logarithmus."""
import sys, time
HIER = r"D:/Users/janis/Documents/Matura/tetration-graph/experimente"
sys.path.insert(0, HIER + "/prototypen_hexation")
sys.path.insert(0, HIER)
import mpmath as mp
from blatt import BlattPentation, verfolge
mp.mp.dps = 40
t0 = time.time()
senk = [mp.mpf(1) + 1j * mp.mpf(k) / 100 for k in range(0, 301)]
for name, XA, kx in (("Vorgabe", mp.mpf("1.5"), 0), ("k + 1", mp.mpf("1.5"), 1), ("ohne analyt. Log", mp.mpf(100), 0), ("ohne, k + 1", mp.mpf(100), 1)):
    BP = BlattPentation("2")
    BP.XA, BP.k_extra = XA, kx
    p, v, s = verfolge(senk, BP.P_pfad, runden=16)
    werte = {mp.nstr(mp.im(p[i]), 3): mp.nstr(v[i], 10) for i in range(len(p)) if mp.im(p[i]) * 2 % 1 == 0}
    print(f"{name}: {len(p)} Punkte, Sprung {mp.nstr(s, 3)}, Werte {werte} [{time.time() - t0:.0f}s]", flush=True)
