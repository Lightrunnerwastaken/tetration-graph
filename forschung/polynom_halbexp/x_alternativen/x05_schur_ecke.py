"""X05 — Reduktion der Eckmoden auf eine d×d-Matrix (kubisches Modell), über rechteckige Abschnitte.

R := T_{N+d} J_∞|_{P_N} = [J_N; B]  ((N+1+d)×(N+1), ℓ²_ρ-skaliert), R⁺ = (RᵀR)⁻¹Rᵀ = [R⁺_top, R⁺_bot] (Spalten N+1 | d).
Exakt (bewiesen, elementar): R⁺R = I ⇒ J_N⁻¹ = R⁺_top + R⁺_bot S⁻¹ B R⁺_top,  S := I_d − B R⁺_bot  (d×d),
und σ_{d+1}(J_N⁻¹) ≤ ‖R⁺‖ = 1/σ_min(R).
Ausgabe je N: σ_min(R) für d = 1…8; für d = dfix: Singulärwerte von S, Normen ‖R⁺_bot‖, ‖B R⁺_top‖, Kontrolle σ₁(J_N⁻¹) über die Formel.
Aufruf: python x05_schur_ecke.py N[,N2,…] [dfix]
"""
import sys
import time

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
NS = [int(x) for x in sys.argv[1].split(",")]
DFIX = int(sys.argv[2]) if len(sys.argv) > 2 else 6
DLIST = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else list(range(0, 9))
DMAX = max(max(DLIST), DFIX)


def kubik():
    Lr = arb("0.31813150520476413531265425158766451720351761387139986692237")
    Li = arb("1.33723570143068940890116214319371061253950213846051241687")
    R2 = Lr * Lr + Li * Li
    rho = R2.sqrt()
    mu_r = ((rho + Lr) / 2).sqrt()
    mu_i = ((rho - Lr) / 2).sqrt()
    rL_re = mu_i / (2 * Li)
    rL_im = -(mu_r - 1) / (2 * Li)
    beta = rL_im / Li
    alpha = rL_re - beta * Lr
    c = [R2 * alpha, 1 + R2 * beta - 2 * Lr * alpha, alpha - 2 * Lr * beta, beta]
    return arb_poly(c), arb_poly([c[1], 2 * c[2], 3 * c[3]]), rho


def to_np(Mt):
    return np.array([[float(Mt[i, j].mid()) for j in range(Mt.ncols())] for i in range(Mt.nrows())])


for N in NS:
    t0 = time.time()
    ctx.prec = int(2.6 * N) + 500
    f, fp, rho = kubik()
    m = fp(f)
    M = N + DMAX
    A = arb_mat(M + 1, N + 1)
    rk = [arb(1)]
    for k in range(1, M + 1):
        rk.append(rk[-1] * rho)
    pot = arb_poly([1])
    for j in range(N + 1):
        col = pot + m * arb_poly([0] * j + [1])
        inv_rj = 1 / rk[j]
        for k in range(min(col.degree(), M) + 1):
            A[k, j] = col[k] * rk[k] * inv_rj
        pot = (pot * f).truncate(M + 1)
    # quadratisch
    JN = arb_mat(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1):
            JN[i, j] = A[i, j]
    G = JN.transpose() * JN
    out = []
    Gd = G
    dl = sorted(set(DLIST + [DFIX]))
    for d in range(0, DMAX + 1):
        if d > 0:
            row = arb_mat(1, N + 1)
            for j in range(N + 1):
                row[0, j] = A[N + d, j]
            Gd = Gd + row.transpose() * row
        if d not in dl:
            continue
        Gi = Gd.inv()
        ew = np.sort(np.linalg.eigvalsh((to_np(Gi) + to_np(Gi).T) / 2))[::-1]
        out.append((d, 1 / np.sqrt(ew[:3])))
        if d == DFIX:
            Gi_fix = Gi
    print(f"N={N}: σ_min(R_d): " + " | ".join(f"{d}:{o[0]:.4f},{o[1]:.4f}" for d, o in out)
          + f"  [{time.time() - t0:.0f}s]", flush=True)
    # Schur-Reduktion für d = DFIX
    d = DFIX
    R = arb_mat(N + 1 + d, N + 1)
    for i in range(N + 1 + d):
        for j in range(N + 1):
            R[i, j] = A[i, j]
    Rp = Gi_fix * R.transpose()                          # (N+1)×(N+1+d)
    Rtop = arb_mat(N + 1, N + 1)
    Rbot = arb_mat(N + 1, d)
    for i in range(N + 1):
        for j in range(N + 1):
            Rtop[i, j] = Rp[i, j]
        for j in range(d):
            Rbot[i, j] = Rp[i, N + 1 + j]
    B = arb_mat(d, N + 1)
    for i in range(d):
        for j in range(N + 1):
            B[i, j] = A[N + 1 + i, j]
    I_d = arb_mat(d, d)
    for i in range(d):
        I_d[i, i] = 1
    S = I_d - B * Rbot
    mp.mp.dps = 30
    Sm = mp.matrix([[mp.mpf(S[i, j].mid().str(25, radius=False)) for j in range(d)] for i in range(d)])
    sv = sorted([abs(x) for x in mp.svd_r(Sm, compute_uv=False)])
    BRt = to_np(B * Rtop)
    nb = np.linalg.norm(to_np(Rbot), 2)
    nbr = np.linalg.norm(BRt, 2)
    # Kontrolle: J_N⁻¹ über die Formel, grösster Singulärwert
    try:
        Sinv = S.inv()
    except ZeroDivisionError:
        print('   S in arb nicht invertierbar (Präzision) — Formelkontrolle übersprungen', flush=True)
        continue
    corr = Rbot * (Sinv * (B * Rtop))
    Jinv = Rtop + corr
    sJ = np.linalg.svd(to_np(Jinv), compute_uv=False)
    print(f"   d={d}: σ(S) = " + ", ".join(mp.nstr(s, 5) for s in sv) + f";  ‖R⁺_bot‖ = {nb:.4g}, ‖B R⁺_top‖ = {nbr:.4g};"
          f"  σ(J_N⁻¹) per Formel: {sJ[0]:.5g}, {sJ[1]:.5g}, {sJ[2]:.5g}, {sJ[3]:.5g};  ‖R⁺_top‖ = {np.linalg.norm(to_np(Rtop), 2):.5g}"
          f"  [{time.time() - t0:.0f}s]", flush=True)
