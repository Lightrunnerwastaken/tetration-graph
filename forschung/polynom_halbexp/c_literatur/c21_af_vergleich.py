"""Strang C, Runde 3, c21: J_N(q_N) gegen die abgeschnittene exakte Linearisierung A_f.

A δ   = [δ∘q + (q'∘q) δ]_{≤N}          (Jacobi-Matrix am Kneser-Polynom q = q_N)
A_f δ = [δ∘q + (f'∘f) δ]_{≤N},  f'∘f = e^x / f'(x)   (P_N J_∞(f) P_N; Unterschied A_f − A = Multiplikation mit [t'∘q]_{≤N})
Druckt: |[q'∘q]_{≤N}|_0 und |[f'∘f]_{≤N}|_0 (Multiplikatornormen), R_N(L) und R^f_N(L), α_N = |A^{-1}B_N t| und α^f_N = |A_f^{-1}B_N t|.
"""
import json, sys, os
sys.stdout.reconfigure(encoding="utf-8")
from mpmath import mp, mpf, mpc
from flint import arb_mat, arb, ctx

mp.dps = 130
ctx.prec = 950
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
KN = [mpf(x) for x in json.load(open(os.path.join(G, "kneser_halbexp_koeff_1000.json"), encoding="utf-8"))["a"]]
L = mpc("0.31813150520476413", "1.3372357014306895")
RHO = abs(L)


def mul(a, b, n):
    c = [mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def cmul(a, b, n):
    c = [mpc(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        for j, bj in enumerate(b[: n + 1 - i]):
            c[i + j] += ai * bj
    return c


def inv_reihe(c, n):
    r = [type(c[0])(0)] * (n + 1)
    r[0] = 1 / c[0]
    for k in range(1, n + 1):
        r[k] = -sum(c[j] * r[k - j] for j in range(1, k + 1)) / c[0]
    return r


def norm0(v):
    return sum(abs(v[j]) * RHO ** j for j in range(len(v)))


for N in [int(x) for x in (sys.argv[1].split(",") if len(sys.argv) > 1 else "40,80,120,160,200")]:
    q = KN[: N + 1]
    pot = [[mpf(1)] + [mpf(0)] * N]
    for _ in range(N + 1):
        pot.append(mul(pot[-1], q, N))
    dq = [k * q[k] for k in range(1, N + 1)]
    qpq = [sum(dq[j] * pot[j][k] for j in range(len(dq))) for k in range(N + 1)]          # [q'∘q]_{≤N}
    fprime = [(k + 1) * KN[k + 1] for k in range(N + 1)]                                    # f' bis Grad N
    expx = [mpf(1) / mp.factorial(k) for k in range(N + 1)]
    fpf = mul(expx, inv_reihe(fprime, N), N)                                                # [e^x/f']_{≤N} = [f'∘f]_{≤N}
    def matrix_mit(mult):
        return arb_mat([[arb(mp.nstr(pot[j][k] + (mult[k - j] if j <= k else 0), 125)) for j in range(N + 1)] for k in range(N + 1)])
    A, Af = matrix_mit(qpq), matrix_mit(fpf)
    # rechte Seiten: B_N t (= -F(q)) und Pol bei L
    qq = [sum(q[j] * pot[j][k] for j in range(N + 1)) for k in range(N + 1)]
    Bt = [-(qq[k] - mpf(1) / mp.factorial(k)) for k in range(N + 1)]
    nenner = [L - q[0]] + [-q[k] for k in range(1, N + 1)]
    vL = cmul([mpc(x) for x in pot[N + 1]], inv_reihe(nenner, N), N)
    rhs = arb_mat([[arb(mp.nstr(Bt[j], 125)), arb(mp.nstr(vL[j].real, 125)), arb(mp.nstr(vL[j].imag, 125))] for j in range(N + 1)])
    out = []
    for name, M in (("A", A), ("A_f", Af)):
        s = M.solve(rhs)
        al = norm0([mpf(s[j, 0].mid().str(120, radius=False)) for j in range(N + 1)])
        rl = norm0([mpc(mpf(s[j, 1].mid().str(120, radius=False)), mpf(s[j, 2].mid().str(120, radius=False))) for j in range(N + 1)]) / abs(L) ** (N + 1)
        out.append(f"{name}: α={mp.nstr(al, 4)} R(L)={mp.nstr(rl, 5)}")
    print(f"N={N}: |[q'∘q]|_0={mp.nstr(norm0(qpq), 4)} |[f'∘f]|_0={mp.nstr(norm0(fpf), 4)} | " + " | ".join(out))
    sys.stdout.flush()
