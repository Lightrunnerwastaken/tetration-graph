"""Strang A: Newton direkt von Knesers q_N aus (4 Schritte) für N = 300, 400, 500; Vergleich mit p_N (Newton von Π_N f_K aus)."""
import json, sys, time
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp
from lin_kern import DATEN, kneser_koeff, newton_N
mp.mp.dps = 60
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
pn = json.loads((DATEN / "pN_gross.json").read_text(encoding="utf-8"))
a = kneser_koeff(1001)
for N in [int(x) for x in sys.argv[1].split(",")]:
    t0 = time.time()
    p, korr = newton_N([mp.nstr(x, 55) for x in a[:N + 1]], N, int(0.72 * N) + 200, schritte=4)
    ref = [mp.mpf(x) for x in pn[str(N)]["p"]]
    d = max(abs(mp.mpf(x) - y) * RL ** k for k, (x, y) in enumerate(zip(p, ref)))
    print(f"N = {N}: Korrekturen von q_N aus {[mp.nstr(x, 3) for x in korr]}; Abstand zum p_N aus Π_N f_K: {mp.nstr(d, 3)} ({time.time() - t0:.0f} s)", flush=True)
