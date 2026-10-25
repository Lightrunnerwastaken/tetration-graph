"""E1: Satz R (Ordnungs-Rekursion) gegen direkte Rechnung; Schranken l-1 <= t <= l (+1 nur a=2)."""
import sys, time
from collections import Counter
from tlib import t_direct, t_rec, ell

A = int(sys.argv[1]) if len(sys.argv) > 1 else 60
M = int(sys.argv[2]) if len(sys.argv) > 2 else 600
t0 = time.time()
bad = 0; n = 0; over = Counter(); excess = []
for m in range(2, M + 1):
    l = ell(m)
    for a in range(1, A + 1):
        td = t_direct(a, m); tr = t_rec(a, m)
        n += 1
        if td != tr:
            bad += 1
            if bad < 10:
                print("ABWEICHUNG", a, m, td, tr)
        over[td - l] += 1
        if td > l:
            excess.append((a, m))
        if a >= 20 and t_rec(a, m, large=True) != td:
            print("large-Abweichung", a, m)
print(f"{n} Paare (a<= {A}, m<= {M}), Abweichungen Rekursion/direkt: {bad}")
print("Verteilung t - l:", sorted(over.items()))
print("a mit t > l:", sorted(set(a for a, _ in excess)), " Anzahl", len(excess), "Beispiele", excess[:12])
print(f"{time.time()-t0:.1f}s")
