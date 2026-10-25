"""V1: Vorperiode N(m) direkt (Turm-Algorithmus) gegen die delta-Formel, alle m <= MMAX.
Ausgabe: Abweichungen, Haeufigkeit von N, Beispiele. Laufzeit ~1-2 min."""
import sys, time
from collections import Counter
from ef_lib import N_brute, N_formel, Delta, Lam_t, phi, carmichael, EF_mod

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
t0 = time.time()
D = Delta(MMAX, 40)
bad, cnt = [], Counter()
ex = {}
for m in range(2, MMAX + 1):
    nb = N_brute(m, nmax=30)
    nf = N_formel(m, D)
    cnt[nb] += 1
    ex.setdefault(nb, []).append(m)
    if nb != nf:
        bad.append((m, nb, nf))
print("m <= %d: Abweichungen Formel/direkt: %d %s" % (MMAX, len(bad), bad[:10]))
print("Haeufigkeit N(m):", sorted(cnt.items()))
for k in sorted(ex):
    print("  N=%d: erste m %s" % (k, ex[k][:12]))
# Gegenprobe mit phi statt lambda (anderer Reduktionsweg, andere Periodenvielfache)
gp = 0
for m in list(range(2, 400)) + [47, 94, 283, 563, 1129]:
    if m > MMAX:
        continue
    P = Lam_t(m)
    nb2 = 1 + max([0] + [n for n in range(1, 30) if EF_mod(n, m, phi) != EF_mod(n + P, m, phi)])
    if nb2 != N_formel(m, D):
        gp += 1
        print("  phi-Abweichung", m, nb2)
print("phi-Gegenprobe Abweichungen:", gp)
print("Zeit %.1f s" % (time.time() - t0))
