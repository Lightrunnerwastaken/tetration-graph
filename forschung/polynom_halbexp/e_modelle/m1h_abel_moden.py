"""Modell 1h: Ist der Fehler p_N − f_M der Möbius-Drehung in Abel-Koordinaten EIN Kernmodus (wie bei exp, Strang A/B)?
θ_N(w) := α(p_N(α⁻¹(w))) − w − ½, α(x) = arctan(x)/c. Fourier-Koeffizienten von θ_N auf w ∈ [−½, ½] (Höhen v = 0, 0.1).
Kernmoden: ungerade k (θ(w+½) = −θ(w)); gerade k (inkl. 0) sind keine Kernrichtungen."""
import sys
import mpmath as mp
from kern import zweig, g_moebius_ell, horner

mp.mp.dps = 60
c = mp.mpf(sys.argv[1])
NS = [int(x) for x in sys.argv[2].split(",")]
loes = zweig(g_moebius_ell(c), max(NS), bits=900)
M = 64
for N in NS:
    a = loes[N]
    for v in (mp.mpf(0), mp.mpf("0.1")):
        th = []
        for j in range(M):
            w = mp.mpf(j) / M - mp.mpf("0.5") + 1j * v
            x = mp.tan(c * w)
            y = horner(a, mp.mpc(x))
            th.append(mp.atan(y) / c - w - mp.mpf("0.5"))
        koef = {}
        for k in range(-4, 5):
            koef[k] = mp.fsum(th[j] * mp.expjpi(-2 * k * (mp.mpf(j) / M - mp.mpf("0.5"))) for j in range(M)) / M
        print(f"c={mp.nstr(c, 3)} N={N} v={mp.nstr(v, 2)}: " + "  ".join(f"k={k}: {mp.nstr(abs(koef[k]), 3)}" for k in range(-4, 5)))
