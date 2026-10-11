"""Modell 6: parabolischer Fixpunkt im Entwicklungspunkt, g(x) = e^x − 1.
Satz: der Zweig von p_1 = x aus ist p_N = T_N h (h = formale reguläre halbe Iterierte, a_0 = 0, Gleichungssystem dreieckig).
h divergiert (Baker 1958; Écalle: Gevrey-1) ⇒ p_N(x) divergiert für jedes x ≠ 0 — Gegenbeispiel zur Konvergenz."""
import sys
import mpmath as mp
from kern import zweig, g_expm1, horner

mp.mp.dps = 60
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
loes = zweig(g_expm1(), NMAX, bits=500)
a = loes[NMAX]
# formale Reihe dreieckig: [x^n] h∘h = 1/n!, h_1 = 1, Koeffizient von h_n in [x^n] h∘h ist 2 (n ≥ 2)
h = [mp.mpf(0), mp.mpf(1)]
for n in range(2, NMAX + 1):
    h.append(mp.mpf(0))
    # [x^n] Σ_m h_m h(x)^m mit h_n = 0
    pw = [mp.mpf(1)] + [mp.mpf(0)] * n
    s = mp.mpf(0)
    for m in range(1, n + 1):
        pw = [mp.fsum(pw[i] * h[j - i] for i in range(j + 1)) for j in range(n + 1)]
        s += h[m] * pw[n]
    h[n] = (1 / mp.factorial(n) - s) / 2
print("max_N max_k |a_k^(N) − h_k| =", mp.nstr(max(abs(loes[N][k] - h[k]) for N in loes for k in range(N + 1)), 3))
print("a_0^(N) = 0 für alle N:", all(abs(loes[N][0]) < mp.mpf(10) ** -100 for N in loes))
for n in (10, 20, 30, 40, 50, 60):
    if n <= NMAX:
        q = abs(h[n]) / (mp.factorial(n - 1) / (2 * mp.pi) ** n)
        print(f"n={n}: h_n = {mp.nstr(h[n], 6)},  |h_n| / ((n−1)!/(2π)^n) = {mp.nstr(q, 5)}")
for x in ("0.5", "1"):
    print(f"p_N({x}) für N = 10, 20, 30, 40, 50, 60:", [mp.nstr(horner(loes[N], mp.mpf(x)), 8) for N in (10, 20, 30, 40, 50, 60) if N <= NMAX])
