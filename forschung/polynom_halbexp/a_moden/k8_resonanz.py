"""Strang A: Woher kommt der Faktor N in c_N = C·N^{1−2γ}? Test „fast-singuläre Richtung mit σ ~ 1/N in |L|-Gewichten“.

J = J_N(f_K) an q_N, S = diag(|L|^k) (Urbild = Bild, l²). Weil |J| ~ 1.63^N, aber |J⁻¹| mässig (Strang C), wird J⁻¹ in hoher
Genauigkeit gebildet und die SVD von W⁻¹ = S J⁻¹ S⁻¹ betrachtet: grösste Singulärwerte s_i = 1/σ_i(W), Richtung u_i (im Urbild).
Vergleich von u_1 mit d_N = p_N − q_N und mit den abgeschnittenen δ_{±1} (Re/Im), alle gewichtet mit |L|^k.
"""
import json
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb_mat, arb, arb_mat, arb_poly, ctx  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
mp.mp.dps = 50
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
RL = abs(L)
q = json.loads((GEM / "kneser_halbexp_koeff_1000.json").read_text(encoding="utf-8"))["a"]
P = {int(k): v for k, v in json.loads((GEM / "polynome_bis_194.json").read_text(encoding="utf-8")).items()}
D = json.loads((GEM / "stoerrichtungen_delta.json").read_text(encoding="utf-8"))["delta"]["1"]


def jac(a, N):
    f = arb_poly(a)
    Pw = [arb_poly([1])]
    for _ in range(N):
        Pw.append((Pw[-1] * f).truncate(N + 1))
    fpf = arb_poly([0])
    for m in range(1, N + 1):
        fpf += m * a[m] * Pw[m - 1]
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            J[k, j] = Pw[j][k] + (fpf[k - j] if j <= k else 0)
    return J


def fl(x):
    return float(x.mid())


def ortho(vs):
    out = []
    for v in vs:
        for u in out:
            p = sum(a * b for a, b in zip(u, v))
            v = [b - p * a for a, b in zip(u, v)]
        n = sum(a * a for a in v) ** 0.5
        out.append([a / n for a in v])
    return out


print(" N   s1=1/σ_min  s2        s3        s1/N    cos∠(u1, d_N)  cos∠(u1, span δ±1)  cos∠(d_N, span δ±1)  |J⁻¹|_l¹")
t0 = time.time()
for N in (20, 30, 45, 60, 80, 100, 120, 150, 190):
    ctx.prec = 300 + 6 * N
    a = [arb(x) for x in q[:N + 1]]
    J = jac(a, N)
    Ji = J.inv()
    w = [arb(mp.nstr(RL ** k, 60)) for k in range(N + 1)]
    Wi = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            Wi[k, j] = Ji[k, j] * w[k] / w[j]
    l1 = max(sum(abs(fl(Wi[k, j])) for k in range(N + 1)) for j in range(N + 1))
    ctx.prec = 200
    Wf = arb_mat([[arb(Wi[k, j].mid()) for j in range(N + 1)] for k in range(N + 1)])
    E, R = acb_mat(Wf.transpose() * Wf).eig(right=True, algorithm="approx")
    ev = sorted(range(N + 1), key=lambda i: -fl(E[i].real))
    s = [fl(E[i].real) ** 0.5 for i in ev[:3]]
    v1 = arb_mat([[R[j, ev[0]].real] for j in range(N + 1)])
    u1 = Wf * v1
    u1 = [fl(u1[k, 0]) for k in range(N + 1)]
    nu = sum(x * x for x in u1) ** 0.5
    u1 = [x / nu for x in u1]
    dw = [float(mp.mpf(P[N][k]) - mp.mpf(q[k])) * float(RL) ** k for k in range(N + 1)]
    nd = sum(x * x for x in dw) ** 0.5
    dw = [x / nd for x in dw]
    B = ortho([[float(D["re"][k]) * float(RL) ** k for k in range(N + 1)], [float(D["im"][k]) * float(RL) ** k for k in range(N + 1)]])
    cud = abs(sum(x * y for x, y in zip(u1, dw)))
    cub = sum(sum(x * y for x, y in zip(u1, b)) ** 2 for b in B) ** 0.5
    cdb = sum(sum(x * y for x, y in zip(dw, b)) ** 2 for b in B) ** 0.5
    print(f"{N:3d}  {s[0]:9.4f}  {s[1]:8.4f}  {s[2]:8.4f}  {s[0] / N:7.4f}  {cud:12.4f}  {cub:16.4f}  {cdb:18.4f}  {l1:9.3f}"
          f"   ({time.time() - t0:.0f} s)", flush=True)
