"""Strang C, Runde 2, c11: Antwort-Operator für die Konsistenz (b).

Exakt (Lemma D1a): F_N(q_N) = -B_N t_N,  B_N = ([q_N^k]_{≤N})_{k>N},  t_N = Kneser-Schwanz (f_k)_{k>N}.
Daher α_N = |A^{-1}B_N t_N|_s ≤ R_N |t_N|_s mit
   R_N := sup_{k>N} |A^{-1} b_k|_s / (ρ^k (1+k)^s),   A = J_N(q_N),  |v|_s = Σ_j |v_j| ρ^j (1+j)^s.
Druckt R_N, das Spaltenprofil r_k = |A^{-1}b_k|_s/(ρ^k(1+k)^s) für k = N+1 … KMAX, |t_N|_s und α_N/|t_N|_s.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf
from flint import arb_mat, arb, ctx

mp.dps = 90
ctx.prec = 640
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


def mid(x):
    return mpf(x.mid().str(80, radius=False))


for N in (20, 40, 80, 120, 160):
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    dpp = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]
    A = arb_mat([[arb(mp.nstr(pot[j][k] + (dpp[k - j] if j <= k else 0), 85)) for j in range(N + 1)] for k in range(N + 1)])
    Ai = A.inv()
    KMAX = min(999, 5 * N + 100)
    cols = []
    cur = pot[N]
    for k in range(N + 1, KMAX + 1):
        cur = mul(cur, q, N)
        cols.append(cur[:])
    Bm = arb_mat([[arb(mp.nstr(cols[c][j], 85)) for c in range(len(cols))] for j in range(N + 1)])
    R = Ai * Bm
    tvec = arb_mat([[arb(mp.nstr(KN[k], 85))] for k in range(N + 1, KMAX + 1)])
    d = R * tvec  # = A^{-1} B_N t_N  (= -A^{-1}F)
    for s in (0, 4):
        w = [RHO ** j * (1 + j) ** s for j in range(N + 1)]
        prof = []
        for c in range(len(cols)):
            k = N + 1 + c
            prof.append(sum(abs(mid(R[j, c])) * w[j] for j in range(N + 1)) / (RHO ** k * (1 + k) ** s))
        RN = max(prof)
        kmax = N + 1 + prof.index(RN)
        tn = sum(abs(KN[k]) * RHO ** k * (1 + k) ** s for k in range(N + 1, 1000))
        al = sum(abs(mid(d[j, 0])) * w[j] for j in range(N + 1))
        ausw = [N + 1, N + 2, N + 5, int(1.2 * N), int(1.5 * N), 2 * N, 3 * N]
        pr = " ".join(f"k={k}:{float(prof[k - N - 1]):.2e}" for k in ausw if k <= KMAX)
        print(f"N={N:3d} s={s}: R_N={float(RN):.3e} (bei k={kmax}) |t_N|_s={float(tn):.3e} alpha={float(al):.3e} "
              f"alpha/|t|={float(al / tn):.3e} | Profil r_k: {pr}")
    sys.stdout.flush()
