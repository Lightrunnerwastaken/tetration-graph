"""Strang A, Test: Geschwindigkeit und Genauigkeit des Kneser-Motors bei komplexen Punkten auf |z| = 1.25."""
import sys
import time

sys.path.insert(0, r"D:\Users\janis\Documents\Matura\tetration-graph\experimente")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

from gemeinsam import kneser  # noqa: E402

mp.mp.dps = 70
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
r = mp.mpf("1.25")
zs = [r * mp.expjpi(2 * mp.mpf(j) / n + mp.mpf(1) / 1000) for j in range(n)]
t0 = time.time()
A = kneser("slog", "exp(1)", zs, dps=dps)
t1 = time.time()
f = kneser("sexp", "exp(1)", [a + mp.mpf("0.5") for a in A], dps=dps)
t2 = time.time()
# Rundreise: sexp(slog(z)) = z
back = kneser("sexp", "exp(1)", A, dps=dps)
t3 = time.time()
print(f"dps={dps} n={n}: slog {t1-t0:.1f}s, sexp {t2-t1:.1f}s, sexp(rund) {t3-t2:.1f}s")
for z, a, w, b in zip(zs, A, f, back):
    print(mp.nstr(z, 8), mp.nstr(a, 12), mp.nstr(w, 12), "Rundreise", mp.nstr(abs(b - z), 3))
