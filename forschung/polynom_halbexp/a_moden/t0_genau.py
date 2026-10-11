"""Strang A, Test: Genauigkeit von slog/sexp des Motors bei dps 45/50/55 auf |z| = r (obere Hälfte).

Prüfungen: (1) Vergleich zweier dps, (2) slog(e^z) − slog(z) − 1, (3) sexp(slog z) − z, (4) sexp(A+1) − e^z.
"""
import sys
import time

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

mp.mp.dps = 110
r = mp.mpf(sys.argv[1])
n = int(sys.argv[2])
dpsl = [int(x) for x in sys.argv[3].split(",")]
zs = [r * mp.expjpi((mp.mpf(j) + mp.mpf(1) / 3) / n) for j in range(n)]  # obere Hälfte, Winkel in (0, π)
res = {}
for d in dpsl:
    gp = FatouGP(dps=d, fatou_gp="fork")
    s = gp.session("exp(1)")
    t0 = time.time()
    A = s.slog_batch(zs)
    t1 = time.time()
    Ae = s.slog_batch([mp.exp(z) for z in zs])
    t2 = time.time()
    f = s.sexp_batch([a + mp.mpf("0.5") for a in A])
    t3 = time.time()
    rund = s.sexp_batch(A)
    f1 = s.sexp_batch([a + 1 for a in A])
    gp.close()
    res[d] = (A, f)
    fe = max(abs(b - a - 1) for a, b in zip(A, Ae))
    rr = max(abs(w - z) for w, z in zip(rund, zs))
    ee = max(abs(w - mp.exp(z)) for w, z in zip(f1, zs))
    print(f"dps {d}: slog {t1-t0:.1f}s/{n}, sexp {t3-t2:.1f}s/{n}; max|A(e^z)-A(z)-1| = {mp.nstr(fe,3)}, "
          f"max|T(A(z))-z| = {mp.nstr(rr,3)}, max|T(A+1)-e^z| = {mp.nstr(ee,3)}")
for d1, d2 in zip(dpsl[:-1], dpsl[1:]):
    dA = max(abs(a - b) for a, b in zip(res[d1][0], res[d2][0]))
    df = max(abs(a - b) for a, b in zip(res[d1][1], res[d2][1]))
    print(f"dps {d1} gegen {d2}: max|ΔA| = {mp.nstr(dA,3)}, max|Δf| = {mp.nstr(df,3)}")
