"""Strang C, Runde 2, c10: Grösse und Profil des Residuums F_N(q_N) = -[t_N ∘ q_N]_{≤N}.

Druckt |F|_{ρ,s} = Σ_j |F_j| ρ^j (1+j)^s für s = 0, 4 (ρ = |L|), das Profil |F_j| ρ^j an einigen j, und zum Vergleich
α_N = |J(q_N)^{-1}F|_{ρ,s}. Prüft ausserdem die Identität F_N(q_N) = -[t∘q]_{≤N} (t = Kneser-Schwanz ab N+1, bis k ≤ 999).
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf
from flint import arb_mat, arb, ctx

mp.dps = 80
ctx.prec = 600
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
RHO = mpf("1.3745570107437075")


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def potenzen(a, n, m):
    pot = [[mpf(1)] + [mpf(0)] * n]
    for _ in range(1, m + 1):
        pot.append(mul(pot[-1], a, n))
    return pot


for N in (20, 40, 80, 120, 160):
    q = KN[: N + 1]
    pot = potenzen(q, N, N)
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    F = [qq[k] - mpf(1) / mp.factorial(k) for k in range(N + 1)]
    # -[t∘q]_{≤N} mit t-Koeffizienten k = N+1 .. KMAX
    KMAX = min(999, 6 * N + 200)
    tq = [mpf(0)] * (N + 1)
    cur = pot[N]
    for k in range(N + 1, KMAX + 1):
        cur = mul(cur, q, N)
        fk = KN[k]
        for j in range(N + 1):
            tq[j] += fk * cur[j]
    diff = max(abs(F[j] + tq[j]) * RHO ** j for j in range(N + 1))
    # J und alpha
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    J = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 70)) for j in range(N + 1)] for k in range(N + 1)])
    sF = J.solve(arb_mat([[arb(mp.nstr(x, 70))] for x in F]))
    zeile = []
    for s in (0, 4):
        w = [RHO ** j * (1 + j) ** s for j in range(N + 1)]
        nF = sum(abs(F[j]) * w[j] for j in range(N + 1))
        al = sum(abs(mpf(sF[j, 0].mid().str(60, radius=False))) * w[j] for j in range(N + 1))
        zeile.append(f"s={s}: |F|={float(nF):.3e} alpha={float(al):.3e}")
    prof = [float(abs(F[j]) * RHO ** j) for j in range(N + 1)]
    jmax = max(range(N + 1), key=lambda j: prof[j])
    print(f"N={N:3d}: {' | '.join(zeile)} | Identitaet max|F+[t∘q]|ρ^j = {float(diff):.1e} | "
          f"Profil |F_j|ρ^j: j=0:{prof[0]:.1e} j=N/4:{prof[N//4]:.1e} j=N/2:{prof[N//2]:.1e} j=3N/4:{prof[3*N//4]:.1e} "
          f"j=N:{prof[N]:.1e} max bei j={jmax}: {prof[jmax]:.1e}")
    sys.stdout.flush()
