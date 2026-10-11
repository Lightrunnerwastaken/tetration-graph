"""O22 — Rigorose (arb) untere Schranken für c_N := inf_{δ∈P_N} ‖J_∞δ‖_{H²(D)}/‖δ‖_{H²(D)} bei KNESER über den Rechteck-Abschnitt
R_N = T_{3N} J_∞|_{P_N} (σ_min(R_N) ≤ c_N, da Abschneiden die Norm nur verkleinert).
Zertifikat mit Vorkonditionierung: Q = genäherte orthonormale Eigenbasis von G = RᵀR (mpmath, 90 Stellen), M = Qᵀ(G − c²I)Q in arb;
M strikt diagonaldominant mit positiver Diagonale (Gershgorin) ⇒ M ≻ 0 ⇒ G − c²I ≻ 0 (Kongruenz) ⇒ σ_min(R) > c.
Eingabe: Kneser-Koeffizienten aus gemeinsam/kneser_halbexp_koeff_1000.json als arb-Bälle a_k ± 10·fehler_abs_k (Fehler = Strang-A-Schätzung,
nicht selbst zertifiziert), m = e^x / f' als Reihe. Zertifikat: Cholesky von R^T R − c² I in Ball-Arithmetik gelingt ⇒ σ_min(R) > c.
Aufruf: python o22_kneser_rechteck.py N1,N2,...  [c_kandidaten]
"""
import json
import os
import sys
import time

from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HIER, "..", "gemeinsam", "kneser_halbexp_koeff_1000.json")))
Ns = [int(s) for s in sys.argv[1].split(",")]
cs = [float(s) for s in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.70, 0.72, 0.73, 0.7315]
ctx.prec = 500
ctx.cap = 4000
L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"
rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()


def cholesky_ok(G, c2):
    """Cholesky-Zerlegung von G − c2·I in arb; True, wenn alle Pivots sicher > 0."""
    n = G.nrows()
    A = [[G[i, j] - (c2 if i == j else 0) for j in range(n)] for i in range(n)]
    Lm = [[arb(0)] * n for _ in range(n)]
    for j in range(n):
        s = A[j][j]
        for k in range(j):
            s -= Lm[j][k] * Lm[j][k]
        if not (s > 0):
            return False, j
        Lm[j][j] = s.sqrt()
        for i in range(j + 1, n):
            t = A[i][j]
            for k in range(j):
                t -= Lm[i][k] * Lm[j][k]
            Lm[i][j] = t / Lm[j][j]
    return True, n


for N in Ns:
    t0 = time.time()
    K = 3 * N
    a = [arb(d["a"][k], 10 * float(d["fehler_abs"][k])) for k in range(K + 2)]
    f = arb_poly(a[:K + 1])
    fp = [a[k + 1] * (k + 1) for k in range(K + 1)]
    ex = [arb(1)]
    for k in range(1, K + 1):
        ex.append(ex[-1] / k)
    mk = [arb(0)] * (K + 1)
    for k in range(K + 1):
        s = ex[k]
        for j in range(k):
            s -= mk[j] * fp[k - j]
        mk[k] = s / fp[0]
    pw = [arb_poly([1])]
    for _ in range(N):
        pw.append((pw[-1] * f).truncate(K + 1))
    rp = [rho ** k for k in range(K + 1)]
    R = arb_mat(K + 1, N + 1)
    for j in range(N + 1):
        for k in range(K + 1):
            v = pw[j][k] if k <= pw[j].degree() else arb(0)
            if k >= j:
                v = v + mk[k - j]
            R[k, j] = v * rp[k] / rp[j]
    G = R.transpose() * R
    rad = max(float(G[i, j].rad()) for i in range(N + 1) for j in range(N + 1))
    import mpmath as mp
    mp.mp.dps = 90
    n = N + 1
    Gm = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Gm[i, j] = mp.mpf(G[i, j].mid().str(85, radius=False))
    ev, Qm = mp.eigsy(Gm)
    lam_min = min(ev)
    Q = arb_mat(n, n, [arb(mp.nstr(Qm[i, j], 85)) for i in range(n) for j in range(n)])
    QtGQ = Q.transpose() * G * Q
    QtQ = Q.transpose() * Q
    out = [f"σ_min(mid) = {mp.nstr(mp.sqrt(lam_min), 8)}"]
    for c in cs:
        c2 = arb(c) ** 2
        ok = True
        for i in range(n):
            diag = QtGQ[i, i] - c2 * QtQ[i, i]
            off = arb(0)
            for j in range(n):
                if j != i:
                    off += abs(QtGQ[i, j] - c2 * QtQ[i, j])
            if not (diag - off > 0):
                ok = False
                break
        out.append(f"c = {c}: {'σ_min > c ZERTIFIZIERT' if ok else 'nicht zertifiziert'}")
    print(f"Kneser N = {N} ({time.time()-t0:.0f}s, max Radius in G {rad:.1e}): " + "; ".join(out), flush=True)
