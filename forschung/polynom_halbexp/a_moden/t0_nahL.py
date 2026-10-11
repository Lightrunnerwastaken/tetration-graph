"""Strang A, Test: dps-80-slog nahe arg L (Punkte j ≈ 1743 von M = 8192) — welche Punkte hängen?"""
import sys
import time

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 90
M = 8192
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 80
r = mp.mpf(sys.argv[1])
for j in range(1664, 1792, 8):
    z = r * mp.expjpi(2 * mp.mpf(j) / M)
    gp = FatouGP(dps=dps, fatou_gp="fork", eval_timeout=15)
    s = gp.session("exp(1)")
    t0 = time.time()
    try:
        a = s.slog_batch([z])[0]
        out = f"{time.time() - t0:6.2f}s A = {mp.nstr(a, 8)}"
    except Exception as e:
        out = f"FEHLER/Zeit {time.time() - t0:.0f}s {str(e).splitlines()[-1][:50]}"
    gp.close()
    print(f"r={r} j={j} |z−L|={mp.nstr(abs(z - (-mp.lambertw(-1)).conjugate()), 3)}: {out}", flush=True)
