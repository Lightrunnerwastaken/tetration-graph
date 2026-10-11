"""Strang A, Test: im Paket j = 1664..1791 (r = 1.34, M = 8192) jeden dps-80-slog-Aufruf einzeln mit Zeitlimit."""
import sys
import time

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 90
M = 8192
r = mp.mpf("1.34")
js = list(range(1664, 1792))
zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in js]
g57 = FatouGP(dps=57, fatou_gp="fork")
s57 = g57.session("exp(1)")
A57 = s57.slog_batch(zs)
w0 = s57.sexp_batch([a + mp.mpf("0.5") for a in A57])
g57.close()
gp = None
langsam = []
for j, z, w in zip(js, zs, w0):
    for name, p in (("z", z), ("w0", w), ("e^z", mp.exp(z))):
        if gp is None:
            gp = FatouGP(dps=80, fatou_gp="fork", eval_timeout=10)
            s = gp.session("exp(1)")
        t0 = time.time()
        try:
            s.slog_batch([p])
            dt = time.time() - t0
            if dt > 0.5:
                langsam.append((j, name, round(dt, 2)))
        except Exception as e:
            langsam.append((j, name, "FEHLER " + str(e).splitlines()[-1][:40]))
            try:
                gp.close()
            except Exception:
                pass
            gp = None
print("langsam/Fehler:", langsam)
