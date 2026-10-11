"""Strang A, Test: Zeit von slog bei dps 80 an Punkten w = f_K(z), |z| = 1.34, direkt und über log(w)."""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 90
d = json.loads((Path(__file__).resolve().parent / "daten" / "kreis_r1.34_M4096.json").read_text(encoding="utf-8"))
f = [mp.mpc(*x) for x in d["f"]]
A = [mp.mpc(*x) for x in d["A"]]
idx = list(range(0, 2049, 128))
gp = FatouGP(dps=80, fatou_gp="fork")
s = gp.session("exp(1)")
for i in idx:
    t0 = time.time()
    a = s.slog_batch([f[i]])[0]
    t1 = time.time()
    b = s.slog_batch([mp.log(f[i])])[0] + 1
    t2 = time.time()
    print(f"j={i:5d} w={mp.nstr(f[i], 6):>28} direkt {t1-t0:6.2f}s  log {t2-t1:6.2f}s  |A(w)−A(z)−½| = {mp.nstr(abs(a - A[i] - 0.5), 3)}"
          f"  |direkt−log| = {mp.nstr(abs(a - b), 3)}", flush=True)
gp.close()
