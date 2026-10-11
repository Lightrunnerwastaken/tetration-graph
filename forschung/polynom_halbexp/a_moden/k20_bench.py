"""Strang A, Runde 2: Zeit- und Genauigkeitstest für Π_N f_K (zwei Präzisionen, Rauschen in den Eingangskoeffizienten,
Vergleich mit p_N bei N = 194 bzw. 232)."""
import json
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402

from lin_kern import GEM, DATEN, kneser_koeff, pi_N  # noqa: E402

mp.mp.dps = 70
Ns = [int(x) for x in sys.argv[1].split(",")]
a = kneser_koeff(max(Ns) + 1)
P194 = [mp.mpf(x) for x in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8"))["194"]]
P232 = [mp.mpf(x) for x in json.loads((GEM / "polynome_195_bis_232.json").read_text(encoding="utf-8"))["232"]]


def wert(b, z):
    s = mp.mpc(0)
    for x in reversed(b):
        s = s * z + mp.mpf(x.mid().str(50, radius=False))
    return s


for N in Ns:
    p1 = int(0.72 * N) + 200
    t0 = time.time()
    b1 = pi_N(a, N, p1)
    t1 = time.time()
    b2 = pi_N(a, N, p1 + 150)
    t2 = time.time()
    # Rauschen in a: 4e-48·1.34^{−k} (Grösse des Datenfehlers)
    import random
    random.seed(1)
    rs = [mp.mpf(4e-48) * mp.mpf("1.34") ** (-k) * (2 * random.random() - 1) for k in range(N + 2)]
    b3 = pi_N(a, N, p1, rauschen=rs)
    zs = [mp.mpf("0.5"), mp.mpf(1), mp.mpc(0, 1), mp.mpc("-0.9", "0.5")]
    d12 = max(abs(wert(b1, z) - wert(b2, z)) for z in zs)
    d13 = max(abs(wert(b1, z) - wert(b3, z)) for z in zs)
    f1 = mp.fsum(a[:1001]) if N < 1000 else None
    extra = ""
    if N in (194, 232):
        P = P194 if N == 194 else P232
        extra = f"; |Π_N f(1) − p_N(1)| = {mp.nstr(abs(wert(b1, 1) - mp.fsum(P)), 3)}, Π_N f(1) − T(½) = {mp.nstr(mp.re(wert(b1, 1)) - f1, 4)}, p_N(1) − T(½) = {mp.nstr(mp.fsum(P) - f1, 4)}"
    print(f"N = {N}: prec {p1}: {t1 - t0:.1f} s, prec {p1 + 150}: {t2 - t1:.1f} s; Präzisionsabstand {mp.nstr(d12, 3)}, "
          f"Rauschwirkung {mp.nstr(d13, 3)}{extra}", flush=True)
