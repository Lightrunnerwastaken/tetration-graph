"""V2: kleinste Moduln a(n) mit Vorperiode N = n (ueber delta-Rekursion, Primzahlen <= PMAX und Primzahlpotenzen),
Gegenprobe der kleinen Werte mit dem direkten Turm-Algorithmus; Statistik N(p) gegen N(p-1)+1; 2er-Potenzen bis k = 40.
Laufzeit ~1-3 min fuer PMAX = 10^6."""
import sys, time
from collections import Counter
from ef_lib import Delta, N_formel, N_brute, factor, Lam_t, EF_mod, HUGE, INF

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
NMAX = 26
t0 = time.time()
D = Delta(PMAX, NMAX)
print("Delta fertig: %d Primzahlen, n <= %d, %.1f s" % (len(D.primes), NMAX, time.time() - t0))
# Kontrolle: in der letzten Zeile (n = NMAX) darf kein delta < 2 auftreten (sonst NMAX zu klein)
small_last = [p for p in D.primes if D.d[(p, NMAX)] < 3]
print("Primzahlen mit delta_p(NMAX) < 3 (muss leer sein):", small_last[:5])

Np = {p: D.bad_max(p, 1) + 1 for p in D.primes}
cnt = Counter(Np.values())
print("Haeufigkeit N(p), p <= %d:" % PMAX, sorted(cnt.items()))
# kleinste Primzahl und kleinste Primzahlpotenz je N
first_p = {}
for p in D.primes:
    first_p.setdefault(Np[p], p)
first_pp = {}
for p in D.primes:
    k = 1
    while p ** k <= PMAX:
        n = D.bad_max(p, k) + 1
        first_pp.setdefault(n, p ** k)
        first_pp[n] = min(first_pp[n], p ** k)
        k += 1
print("kleinste Primzahl je N:", sorted(first_p.items()))
print("kleinste Primzahlpotenz (= kleinster Modul) je N:", sorted(first_pp.items()))
# Gegenprobe direkt fuer die kleinsten Moduln bis 10^5 (Periode Lam_t kann gross sein; Turm-Algorithmus braucht nur die Kettenlaenge)
for n, m in sorted(first_pp.items()):
    if m <= 200000:
        print("  direkt N(%d) = %d (Formel %d)" % (m, N_brute(m, nmax=NMAX + 2), n))
# N(p) gegen N(p-1)+1
le = sum(1 for p in D.primes if p > 2 and Np[p] <= N_formel(p - 1, D) + 1)
eq = sum(1 for p in D.primes if p > 2 and Np[p] == N_formel(p - 1, D) + 1)
print("p ungerade: N(p) <= N(p-1)+1 in %d von %d; Gleichheit in %d" % (le, len(D.primes) - 1, eq))
# Cunningham-Ketten p -> 2p+1
cc = [(p, 2 * p + 1) for p in D.primes if 2 * p + 1 <= PMAX and Np.get(2 * p + 1)]
ok = sum(1 for p, q in cc if Np[q] >= Np[p] + 1 or Np[p] == 1)
print("q = 2p+1 prim: N(q) >= N(p)+1 (oder N(p)=1) in %d von %d" % (ok, len(cc)))
# Wachstum: max N(p) gegen log2 p
import math
worst = max((Np[p] - math.log2(p + 1), p) for p in D.primes)
print("max (N(p) - log2(p+1)) = %.3f bei p = %d" % worst)
# 2er-Potenzen direkt bis k = 40 gegen geschlossene Formel N(2^k) = min{n : delta_2(n) >= k}
d2 = {2: 1, 3: 3, 4: 18, 5: 20, 6: HUGE}
def N2(k):
    if k <= 1:
        return 1
    return min(n for n in d2 if d2[n] >= k)
bad = [(k, N_brute(2 ** k, nmax=20, P=2), N2(k)) for k in range(1, 41) if N_brute(2 ** k, nmax=20, P=2) != N2(k)]
print("2^k, k <= 40, direkt gegen Formel: Abweichungen", bad)
print("N(2^k) fuer k=1..25:", [N_brute(2 ** k, nmax=20, P=2) for k in range(1, 26)])
# 3er- und 5er-Potenzen
for p in (3, 5, 7, 11):
    row = []
    for k in range(1, 9):
        m = p ** k
        row.append((k, N_formel(m, D), N_brute(m, nmax=20) if Lam_t(m) < 10 ** 12 else None))
    print("p=%d: (k, Formel, direkt)" % p, row)
print("Zeit %.1f s" % (time.time() - t0))
