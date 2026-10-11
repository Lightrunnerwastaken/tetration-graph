"""Strang A, Test: nur slog bei hohem dps: A(e^z) − A(z) − 1 und Vergleich mit dps 57."""
import sys
import time

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 130
r = mp.mpf(sys.argv[1])
n = int(sys.argv[2])
dpsl = [int(x) for x in sys.argv[3].split(",")]
zs = [r * mp.expjpi((mp.mpf(j) + mp.mpf(1) / 3) / n) for j in range(n)]
res = {}
for d in dpsl:
    gp = FatouGP(dps=d, fatou_gp="fork")
    s = gp.session("exp(1)")
    t0 = time.time()
    A = s.slog_batch(zs)
    t1 = time.time()
    Ae = s.slog_batch([mp.exp(z) for z in zs])
    gp.close()
    res[d] = A
    fe = max(abs(b - a - 1) for a, b in zip(A, Ae))
    print(f"dps {d}: slog {t1-t0:.1f}s/{n}; max|A(e^z)-A(z)-1| = {mp.nstr(fe,3)}")
for d1, d2 in zip(dpsl[:-1], dpsl[1:]):
    print(f"dps {d1} gegen {d2}: max|ΔA| = {mp.nstr(max(abs(a - b) for a, b in zip(res[d1], res[d2])), 3)}")
