"""Strang A, Runde 2: Π_N f_K für eine Liste von N berechnen und in daten/lin_PN.json speichern (fortsetzbar, Zeitlimit).
Aufruf: python k21_lauf.py N1,N2,... [Zeitlimit_s]. Präzision 0.72·N + 200 Bit (geprüft: +150 Bit ändert 50 Stellen nicht).
"""
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

from lin_kern import kneser_koeff, lade, pi_N, speichere  # noqa: E402

mp.mp.dps = 70
Ns = [int(x) for x in sys.argv[1].split(",")]
limit = float(sys.argv[2]) if len(sys.argv) > 2 else 1080
fertig = lade()
a = kneser_koeff(max(Ns) + 1)
t0 = time.time()
for N in Ns:
    if N in fertig:
        continue
    if time.time() - t0 > limit:
        print("Zeitlimit erreicht", flush=True)
        break
    prec = int(0.72 * N) + 200
    t1 = time.time()
    b = pi_N(a, N, prec)
    speichere(N, prec, b)
    print(f"N = {N}: {time.time() - t1:.1f} s (prec {prec})", flush=True)
print(f"fertig ({time.time() - t0:.0f} s)")
