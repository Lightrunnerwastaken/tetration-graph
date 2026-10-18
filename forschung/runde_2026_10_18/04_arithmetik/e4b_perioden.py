"""E4b: Struktur der kleinsten Periode P(m) von n -> EF(n) mod m."""
from e4_expfakultaet import EF_mod, Lam
from ar_lib import carmichael, fac
from math import gcd
import functools
def lcm(a,b): return a*b//gcd(a,b)
@functools.lru_cache(None)
def P(m):
    if m == 1: return 1
    L = Lam(m)
    N0 = 12
    seq = [EF_mod(n, m) for n in range(N0, N0 + 2 * L + 1)]
    for d in sorted(x for x in range(1, L + 1) if L % x == 0):
        if all(seq[i] == seq[i + d] for i in range(len(seq) - d)):
            return d
rows=[]
for m in range(2, 260):
    if Lam(m) > 40000: continue
    pm = P(m); pl = P(carmichael(m))
    rows.append((m, pm, pl, lcm(m, pl)))
bad=[r for r in rows if r[1]!=r[3]]
print("P(m) == kgV(m, P(lambda m)):", len(rows)-len(bad), "von", len(rows))
print("Abweichungen (m, P(m), P(lambda m), kgV):", bad[:40])

def tilde(x):
    while x % 4 == 0:
        x //= 2
    return x
def Lam_t(m):
    L = 1
    while m > 1:
        L = lcm(L, tilde(m)); m = carmichael(m)
    return L
ok = sum(1 for r in rows if r[1] == Lam_t(r[0]))
print("P(m) == Lambda~(m) (Kette mit gekapptem 2-Anteil):", ok, "von", len(rows))
# groesserer Bereich, nur Formeltest fuer Lambda~ <= 30000 (unabhaengig von Lambda)
cnt = bad2 = 0
for m in range(260, 1500):
    Lt = Lam_t(m)
    if Lt > 30000: continue
    seq = [EF_mod(n, m) for n in range(12, 12 + 2 * Lt + 1)]
    per = min(d for d in range(1, Lt + 1) if Lt % d == 0 and all(seq[i] == seq[i + d] for i in range(len(seq) - d)))
    # Lt selbst muss Periode sein (Teil des Beweises), und minimal
    cnt += 1
    if per != Lt: bad2 += 1; print("  Abweichung", m, per, Lt)
print("m in [260,1500), Lambda~<=30000:", cnt, "Faelle, Abweichungen", bad2)
# Vorperiode ab n = 1
from e4_expfakultaet import EF_small
mx = []
for m in range(2, 400):
    Lt = Lam_t(m)
    if Lt > 30000: continue
    vals = [EF_mod(n, m) for n in range(1, 40 + 2 * Lt)]
    N = 1
    while not all(vals[i] == vals[i + Lt] for i in range(N - 1, len(vals) - Lt)):
        N += 1
    mx.append((N, m))
from collections import Counter
print("Beginn der Periodizitaet N(m) (n ab 1), Haeufigkeit:", sorted(Counter(n for n, _ in mx).items()), " Beispiele N>=5:", [x for x in mx if x[0] >= 5][:12])
