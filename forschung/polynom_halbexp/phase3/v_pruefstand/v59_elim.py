"""Strang V — v59 (aus v58/v48): adjungierte Köpfe ℓ_n als VEKTOREN für das gewichtete CF-Problem (Auftrag Leitung 13:00) und
primale Köpfe am Steg (Kopf×Kopf-Abgleich mit Ks (KK), KERNPROBLEM §4).

J̃_N = D J_N D⁻¹ (x-Rahmen), LU ohne Pivot. Mitgeführt:
  y = L⁻¹[ẽ_L, ẽ_t(Steg)]          ⇒ Λ_n(L), Λ_n(t) für alle n (wie v58),
  R = [ε_L; ε_z(Steg)]·U⁻¹ (Zeilen)  ⇒ δ_n(L), δ_n(z) für alle n (primal am Steg, z = t = −ρ(1+S)),
  Spalten von L gespeichert            ⇒ ℓ_n = e_nᵀL_n⁻¹ als Vektor (Rückwärts-Substitution), für n ∈ V59_NL.
Ausgabe v59_ell_<MODELL>_N<N>.json (ℓ_n exakt als Mantisse/Exponent, Kopf×Kopf-Daten lp, la); Auswertung auf der Linse in v59_eval.py.
Aufruf: python v59_elim.py kubisch|kneser N [P]
"""
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
from flint import arb, acb, arb_mat, arb_poly, acb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parent.parent / "gemeinsam"
MODELL, N = sys.argv[1], int(sys.argv[2])
P = int(sys.argv[3]) if len(sys.argv) > 3 else int(1.3 * N) + 400
NL = [int(x) for x in os.environ.get("V59_NL", "200,400,600,800,1000,1200,1400,1500,1600,1800,2000").split(",") if int(x) <= N]
LINSEN = os.environ.get("V59_S", "0.0166,0.022,0.03").split(",")
MB = int(os.environ.get("V59_MB", "16384"))
SS = (0.0, 0.0126, 0.0137, 0.0147, 0.0166, 0.022, 0.026, 0.03)
L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"



def modell(N):
    L = acb(arb(L_RE), arb(L_IM))
    if MODELL == "kubisch":
        mu = L.sqrt()
        rL = (mu - 1) / (L - L.conjugate())
        beta = rL.imag / L.imag
        alpha = rL.real - beta * L.real
        R2 = L.real ** 2 + L.imag ** 2
        c = [R2 * alpha, 1 + R2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]
        f = arb_poly(c)
        m = f.derivative()(f)
        fk = [f[k] if k <= f.degree() else arb(0) for k in range(N + 1)]
        mk = [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
        return fk, mk
    z = json.loads((GEM / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))["a"]
    g = json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]
    a = []
    for k in range(N + 2):
        x = z[k] if k < len(z) and k <= 640 else g[k]
        a.append(arb(x[0] if isinstance(x, list) else x))
    fk = a[:N + 1]
    fp = [a[k + 1] * (k + 1) for k in range(N + 1)]
    ex = [arb(1)]
    for k in range(1, N + 1):
        ex.append(ex[-1] / k)
    mk = [arb(0)] * (N + 1)
    for k in range(N + 1):
        s_ = ex[k]
        for j in range(k):
            s_ -= mk[j] * fp[k - j]
        mk[k] = s_ / fp[0]
    return fk, mk


def lnabs(x):
    return float(abs(x).log().mid())


def main():
    t0 = time.time()
    ctx.prec = P
    fk, mk = modell(N)
    f = arb_poly(fk)
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    L = acb(arb(L_RE), arb(L_IM))
    rp = [rho ** k for k in range(N + 1)]
    A = arb_mat(N + 1, N + 1)
    col = arb_poly([1])
    for j in range(N + 1):
        for k in range(N + 1):
            ck = col[k] if k <= col.degree() else arb(0)
            mm = mk[k - j] if k >= j else arb(0)
            A[k, j] = ((ck + mm) * rp[k] / rp[j]).mid()
        col = (col * f).truncate(N + 1)
    print(f"[{MODELL}, N = {N}, P = {P} Bit, NL = {NL}, Linsen = {LINSEN}, MB = {MB}] J̃ aufgebaut ({time.time() - t0:.0f}s)", flush=True)
    q = acb(rho) / L
    e_re, e_im = [], []
    qk = acb(1)
    for k in range(N + 1):
        e_re.append(qk.real.mid())
        e_im.append(qk.imag.mid())
        qk = qk * q
    nb = 2 + len(SS)
    b = arb_mat(N + 1, nb)
    R = arb_mat(nb, N + 1)
    for k in range(N + 1):
        b[k, 0], b[k, 1] = e_re[k], e_im[k]
        R[0, k], R[1, k] = e_re[k], -e_im[k]          # ε_L = conj ẽ_L
    for a_, S_ in enumerate(SS):
        qq = arb(-1) / (1 + arb(S_))                    # ẽ_t, t = −ρ(1+S): (ρ/t)^k
        qz = -(1 + arb(S_))                             # ε_z, z = −ρ(1+S): (z/ρ)^k
        q1, q2 = arb(1), arb(1)
        for k in range(N + 1):
            b[k, 2 + a_] = q1.mid()
            R[2 + a_, k] = q2.mid()
            q1, q2 = q1 * qq, q2 * qz
    off = 0
    m = N + 1
    yL, wL, ys, ws = [], [], [], []
    Lcol = []
    for k in range(N + 1):
        i = k - off
        p = A[i, i]
        yL.append(acb(b[i, 0], b[i, 1]))
        ys.append([b[i, 2 + a_] for a_ in range(len(SS))])
        wcol = arb_mat(nb, 1)
        for r in range(nb):
            wcol[r, 0] = (R[r, i] / p).mid()
        wL.append(acb(wcol[0, 0], wcol[1, 0]))
        ws.append([wcol[2 + a_, 0] for a_ in range(len(SS))])
        cvec = arb_mat(m, 1)
        rvec = arb_mat(1, m)
        cl = []
        for ii in range(i + 1, m):
            v = (A[ii, i] / p).mid()
            cvec[ii, 0] = v
            cl.append(v)
        Lcol.append(cl)                                  # L[k+1.., k]
        for jj in range(i + 1, m):
            rvec[0, jj] = A[i, jj]
        A = (A - cvec * rvec).mid()
        bro = arb_mat(1, nb)
        for jj in range(nb):
            bro[0, jj] = b[i, jj]
        b = (b - cvec * bro).mid()
        R = (R - wcol * rvec).mid()
        if (k + 1) % 250 == 0 and k < N:
            noff = k + 1
            mn = N + 1 - noff
            d = noff - off
            A2 = arb_mat(mn, mn)
            b2 = arb_mat(mn, nb)
            R2 = arb_mat(nb, mn)
            for ii in range(mn):
                for jj in range(mn):
                    A2[ii, jj] = A[ii + d, jj + d]
                for jj in range(nb):
                    b2[ii, jj] = b[ii + d, jj]
                    R2[jj, ii] = R[jj, ii + d]
            A, b, R, off, m = A2, b2, R2, noff, mn
            print(f"  Schritt {k + 1}/{N + 1} ({time.time() - t0:.0f}s)", flush=True)
    del A, b, R
    Pn = [wL[n] * yL[n] for n in range(N + 1)]
    for n in (1000, 1500, 2000):
        if n <= N:
            print(f"  |P_{n}| = {float(abs(Pn[n]).mid()):.5f}", flush=True)
    # Kopf×Kopf am Steg: lnR_prim(n) = ln|δ_n(z)/δ_n(L)|, lnR_adj(n) = ln|Λ_n(t)/Λ_n(L)|, z = t = −ρ(1+S)
    lp = [[lnabs(ws[n][a_]) - lnabs(wL[n]) for a_ in range(len(SS))] for n in range(N + 1)]
    la = [[lnabs(ys[n][a_]) - lnabs(yL[n]) for a_ in range(len(SS))] for n in range(N + 1)]
    print("Kopf×Kopf am Steg (z = t = −ρ(1+S)), S = " + str(SS), flush=True)
    for n in [300, 600, 1000, 1500, 2000]:
        if n <= N:
            print(f"  n = {n}: ln|δ_n(z)/δ_n(L)| = " + " ".join(f"{v:+.3f}" for v in lp[n]) + " | ln|Λ_n(t)/Λ_n(L)| = " +
                  " ".join(f"{v:+.3f}" for v in la[n]) + " | Summe = " + " ".join(f"{lp[n][a_] + la[n][a_]:+.3f}" for a_ in range(len(SS))), flush=True)
    for n1, n2 in [(300, 600), (600, 1000), (1000, 1500), (1500, 2000), (1600, 2000)]:
        if n2 <= N:
            rp_ = [(lp[n2][a_] - lp[n1][a_]) / (n2 - n1) for a_ in range(len(SS))]
            ra_ = [(la[n2][a_] - la[n1][a_]) / (n2 - n1) for a_ in range(len(SS))]
            print(f"  Raten {n1}→{n2}: primal " + " ".join(f"{v:+.5f}" for v in rp_) + " | adjungiert " + " ".join(f"{v:+.5f}" for v in ra_) +
                  " | Summe " + " ".join(f"{rp_[a_] + ra_[a_]:+.5f}" for a_ in range(len(SS))), flush=True)
    for n in (300, 600, 1000):
        if n <= N:
            print(f"  mittlere Raten (ln R)/n bei n = {n} (Ks h-Definition in k52): adjungiert " + " ".join(f"{v / n:+.5f}" for v in la[n]) +
                  " | primal " + " ".join(f"{v / n:+.5f}" for v in lp[n]), flush=True)
    # ℓ_n als Vektoren
    ells = {}
    for n in NL:
        ell = [arb(0)] * (n + 1)
        ell[n] = arb(1)
        for k in range(n - 1, -1, -1):
            cl = Lcol[k]
            s_ = arb(0)
            for i_ in range(k + 1, n + 1):
                s_ += ell[i_] * cl[i_ - k - 1]
            ell[k] = (-s_).mid()
        LamL = sum((ell[k] * acb(e_re[k], e_im[k]) for k in range(n + 1)), acb(0))
        rel = float(abs(LamL / yL[n] - 1).mid())
        mx = max(float(abs(x).mid()) for x in ell)
        ells[n] = (ell, LamL, mx)
        print(f"  ℓ_{n}: Λ_n(L) aus Vektor / aus Elimination − 1 = {rel:.1e}; max|ℓ_k| = {mx:.3e}, |Λ_n(L)| = {float(abs(LamL).mid()):.3e} ({time.time() - t0:.0f}s)", flush=True)
    del Lcol
    ausg = dict(modell=MODELL, N=N, P=P, SS=list(SS), lp=lp, la=la,
                ell={str(n): [[str(x.mid().man_exp()[0]), int(x.mid().man_exp()[1])] for x in ells[n][0]] for n in NL},
                LamL={str(n): [[str(v.mid().man_exp()[0]), int(v.mid().man_exp()[1])] for v in (ells[n][1].real, ells[n][1].imag)] for n in NL},
                maxell={str(n): ells[n][2] for n in NL})
    json.dump(ausg, open(HIER / f"v59_ell_{MODELL}_N{N}.json", "w", encoding="utf-8"))
    print(f"fertig ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
