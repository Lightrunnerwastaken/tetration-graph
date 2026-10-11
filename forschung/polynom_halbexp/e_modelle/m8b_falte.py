"""Modell 8b: elliptische Drehung c = 1.055 (> π/3): wo und warum bricht die Fortsetzung des Zweigs ab?
Newton vom aufgefüllten Vorgänger bei N0+1, …: Schrittgrössen, |J⁻¹|; dazu kleinste Singulärwerte der Jacobi-Matrix (mpmath)."""
import sys, time
import mpmath as mp
from flint import arb, ctx
from kern import zweig, g_moebius_ell, F_J

mp.mp.dps = 40
c = mp.mpf(sys.argv[1])
N0 = int(sys.argv[2])
BITS = int(sys.argv[3]) if len(sys.argv) > 3 else 1200
t0 = time.time()
loes = zweig(g_moebius_ell(c), N0, bits=BITS)
print("Zweig bis", N0, "in", round(time.time() - t0), "s; a0 − τ =", mp.nstr(loes[N0][0] - mp.tan(c / 2), 5), flush=True)
ctx.prec = BITS
g = [arb(x) for x in g_moebius_ell(c)(N0 + 6)]
for N in range(N0 + 1, N0 + 4):
    a = [arb(mp.nstr(x, 300)) for x in loes[N - 1]] + [arb(0)]
    verlauf = []
    ok = False
    for it in range(30):
        F, J = F_J(a, g[:N + 1], N)
        d = J.solve(F)
        nd = max(abs(mp.mpf(d[i, 0].mid().str(20, radius=False))) for i in range(N + 1))
        verlauf.append(mp.nstr(nd, 2))
        a = [arb((a[i] - d[i, 0]).mid()) for i in range(N + 1)]
        if nd < mp.mpf(10) ** (-150):
            ok = True
            break
    # kleinste Singulärwerte (gewichtet ρ = 0.5) an der Endstelle
    _, J = F_J(a, g[:N + 1], N)
    W = mp.matrix(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            W[k, j] = mp.mpf(J[k, j].mid().str(30, radius=False)) * mp.mpf("0.5") ** k / mp.mpf("0.5") ** j
    S = sorted(mp.svd_r(W, compute_uv=False))
    a0f = mp.mpf(a[0].mid().str(30, radius=False)) - mp.tan(c / 2)
    print(f"N={N}: Newton {'konvergiert' if ok else 'NICHT konvergiert'}; a0 − τ = {mp.nstr(a0f, 4)}; Schritte {verlauf[:8]}; kleinste SW {[mp.nstr(s, 3) for s in S[:3]]}", flush=True)
    if not ok:
        break
    loes[N] = [mp.mpf(x.mid().str(300, radius=False)) for x in a]
