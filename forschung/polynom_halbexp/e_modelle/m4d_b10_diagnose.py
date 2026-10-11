"""Modell 4d: b = 10 — warum scheitert der Newton-Schritt bei N = 99? Newton-Verlauf, |J⁻¹| und kleinster Singulärwert-Ersatz."""
import sys, time
import mpmath as mp
from flint import arb, ctx
from kern import zweig, g_exp_basis, F_J

mp.mp.dps = 40
B = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf(10)
N0 = int(sys.argv[2]) if len(sys.argv) > 2 else 98
BITS = int(sys.argv[3]) if len(sys.argv) > 3 else 900
t0 = time.time()
loes = zweig(g_exp_basis(B), N0, bits=BITS)
print("Zweig bis", N0, "in", round(time.time() - t0), "s", flush=True)
ctx.prec = BITS
g = [arb(x) for x in g_exp_basis(B)(N0 + 6)]
for N in range(N0 - 2, N0 + 5):
    if N <= N0:
        continue
    a = [arb(mp.nstr(x, 250)) for x in loes[N - 1]] + [arb(0)]
    verlauf = []
    for it in range(25):
        F, J = F_J(a, g[:N + 1], N)
        d = J.solve(F)
        nd = max(abs(mp.mpf(d[i, 0].mid().str(20, radius=False))) for i in range(N + 1))
        Ji = J.inv()
        nJi = max(sum(abs(mp.mpf(Ji[i, j].mid().str(15, radius=False))) for i in range(N + 1)) for j in range(N + 1))
        verlauf.append((mp.nstr(nd, 3), mp.nstr(nJi, 3)))
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
        if nd < mp.mpf(10) ** (-200):
            break
    print(f"N={N}: Newton-Verlauf (|Schritt|, |J⁻¹|_1): {verlauf[:12]} … ({len(verlauf)} Schritte)", flush=True)
    loes[N] = [mp.mpf(x.mid().str(250, radius=False)) for x in a]
