"""Strang A, Test: welche dps-80-slog-Aufrufe sind langsam? (z auf |z| = 1.34, e^z, f(z)), je Punkt gemessen."""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 90
r = mp.mpf("1.34")
M = 8192
gp = FatouGP(dps=80, fatou_gp="fork", eval_timeout=60)
s = gp.session("exp(1)")
for j in list(range(0, 4097, 256)):
    z = r * mp.expjpi(2 * mp.mpf(j) / M)
    out = []
    for name, p in (("z", z), ("e^z", mp.exp(z))):
        t0 = time.time()
        try:
            s.slog_batch([p])
            out.append(f"{name}: {time.time() - t0:6.2f}s")
        except Exception as e:
            out.append(f"{name}: FEHLER {str(e).splitlines()[-1][:60]}")
    print(f"j={j:5d} z={mp.nstr(z, 5):>24}  " + "  ".join(out), flush=True)
gp.close()
