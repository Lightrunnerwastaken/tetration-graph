"""Strang A, Test: Motor direkt (ohne Zwischenspeicher) bei höherem dps."""
import sys
import time

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 100
dps = int(sys.argv[1])
cache = sys.argv[2] == "1" if len(sys.argv) > 2 else True
t0 = time.time()
gp = FatouGP(dps=dps, fatou_gp="fork", state_cache=cache)
s = gp.session("exp(1)")
print("init", time.time() - t0)
pts = [mp.mpf("0.5"), mp.mpc("0.5", "0.1"), mp.mpc("0.5943548392272753163375129499607193", "-0.8953258429"),
       mp.mpc("0.59", "-0.9"), mp.mpc("0.59", "0.9"), mp.mpc("0.6904", "0.41"), mp.mpc("-1.3", "1.2"), mp.mpc("-1.3", "2.2"), mp.mpc("-1.0", "4.0")]
for p in pts:
    try:
        t0 = time.time()
        v = s.sexp_batch([p])[0]
        print(mp.nstr(p, 8), mp.nstr(v, 30), f"{time.time()-t0:.2f}s")
    except Exception as e:
        print(mp.nstr(p, 8), "FEHLER", str(e).splitlines()[-1])
gp.close()
