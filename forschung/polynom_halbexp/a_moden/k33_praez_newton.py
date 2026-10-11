"""Strang A: Präzisionstest des Newton-Schritts bei grossem N: ein Schritt von Π_N f_K (60-stellig gespeichert) aus mit prec Bit."""
import sys, time, json
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from lin_kern import DATEN, lade, newton_N
N, prec = int(sys.argv[1]), int(sys.argv[2])
t0 = time.time()
p, korr = newton_N(lade()[N]["b"], N, prec, schritte=1)
st = json.loads((DATEN / "pN_gross_schritte.json").read_text(encoding="utf-8"))
mp.mp.dps = 60
L = -mp.lambertw(-1); RL = abs(L)
ref = [mp.mpf(x) for x in st[str(N)]["iterate"][0]]
d = max(abs(mp.mpf(x) - y) * RL ** k for k, (x, y) in enumerate(zip(p, ref)))
print(f"N = {N}, prec {prec}: Korrektur {mp.nstr(korr[0], 4)}; Abstand zum Iterat bei 0.72N+200 Bit: {mp.nstr(d, 3)} ({time.time() - t0:.0f} s)")
