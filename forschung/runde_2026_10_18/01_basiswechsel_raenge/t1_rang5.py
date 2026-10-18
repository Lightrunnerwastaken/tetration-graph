"""Rang 5: Konvergenz von spog_b(P_d(n+th)) - n - th, Phase Phi5, Vergleich mit Rang 4."""
import sys
import time
import mpmath as mp
import lib

mp.mp.dps = 150
bb, dd = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("exp(1)", "2")
t = time.time()
Tb, Td = lib.Tet(bb), lib.Tet(dd)
print(f"Aufbau {time.time()-t:.1f}s")
th = mp.mpf("0.3")
print("endliche n:", [mp.nstr(lib.G5_endlich(Tb, Td, th, n), 25) for n in (0, 1)])
for X in (40, 60, 75):
    print(f"Phi5(0.3), X={X}:", mp.nstr(lib.phi5(Tb, Td, th, X), 40))
t = time.time()
M = 32
v = [lib.phi5(Tb, Td, mp.mpf(j) / M) for j in range(M)]
print(f"Phi5 Mittel (M={M})", mp.nstr(mp.fsum(v) / M, 40), f"  {time.time()-t:.1f}s")
v2 = [lib.phi5(Tb, Td, mp.mpf(j) / 16) for j in range(16)]
print(f"Phi5 Mittel (M=16)", mp.nstr(mp.fsum(v2) / 16, 40))
print("max-min Rang 5", mp.nstr(max(v) - min(v), 12))
w = [lib.phi4(Tb, Td, mp.mpf(j) / M) for j in range(M)]
print("Mittel Rang 4", mp.nstr(mp.fsum(w) / M, 30), " max-min Rang 4", mp.nstr(max(w) - min(w), 12))
for j in range(0, M, 4):
    print(f"  th={j/M:.4f}  Phi5-mu5={mp.nstr(v[j]-mp.fsum(v)/M, 12):>20}   Phi4-mu4={mp.nstr(w[j]-mp.fsum(w)/M, 12):>20}")
