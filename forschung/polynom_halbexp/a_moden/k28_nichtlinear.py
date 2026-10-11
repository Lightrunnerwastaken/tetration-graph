"""Strang A, Runde 2: die echten Lösungen p_N bis N = 1000 per Newton von Π_N f_K aus; nichtlineare Korrektur p_N − Π_N f_K.
Aufruf: python k28_nichtlinear.py N1,N2,... [Zeitlimit_s]. Ergebnis: daten/pN_gross.json (Koeffizienten, Korrekturen)."""
import json, sys, time
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from lin_kern import DATEN, lade, newton_N
mp.mp.dps = 60
Ns = [int(x) for x in sys.argv[1].split(",")]
limit = float(sys.argv[2]) if len(sys.argv) > 2 else 1080
P = lade()
pfad = DATEN / "pN_gross.json"
aus = json.loads(pfad.read_text(encoding="utf-8")) if pfad.is_file() else {}
t0 = time.time()
for N in Ns:
    if str(N) in aus or N not in P:
        continue
    if time.time() - t0 > limit:
        print("Zeitlimit erreicht", flush=True); break
    t1 = time.time()
    p, korr = newton_N(P[N]["b"], N, int(0.72 * N) + 200, schritte=2)
    aus[str(N)] = {"p": p, "korr": [mp.nstr(x, 4) for x in korr]}
    pfad.write_text(json.dumps(aus), encoding="utf-8")
    print(f"N = {N}: Newton-Korrekturen (max |d_k||L|^k) = {[mp.nstr(x, 3) for x in korr]}  ({time.time() - t1:.0f} s)", flush=True)
print("fertig")
