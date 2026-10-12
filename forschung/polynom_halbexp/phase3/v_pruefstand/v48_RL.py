"""Strang V — (R_L) und Eckprodukt für alle n ≤ N aus EINER Vorwärts-Elimination (Auftrag Leitung 09:30, Ks §24).

J̃_N = D J_N D⁻¹ (D = diag ρ^k), J_N[k][j] = [x^k] f^j + [x^{k−j}] m, m = f'∘f. LU ohne Pivot, J̃_n = L_n U_n für alle führenden Blöcke.
  y := L⁻¹ẽ_L (Vorwärts-Substitution mitgeführt), w := ε_Lᵀ U⁻¹ (Rekursion r ← r − w_k·(Zeile k von U)), ẽ_L = ((ρ/L)^k), ε_L = ((L/ρ)^k).
Dann (Satz K5/Ks §24.1): primaler Kopf c_n = U_n⁻¹e_n, adjungierter ℓ_n = e_nᵀL_n⁻¹ (⟨ℓ_n, J̃c_n⟩ = 1) und
  δ_n(L) = ε_Lᵀc_n = w_n,  Λ_n(L) = ℓ_nᵀẽ_L = y_n,  P_n = w_n y_n,  K_n(L,L) = ε_Lᵀ J̃_n⁻¹ ẽ_L = Σ_{k≤n} w_k y_k.
Optional (V48_NORM = "n1,n2,…"): ‖c_n‖_{ℓ²} per eigener Lösung J̃_n c = e_n ⇒ B_n = √n‖c_n‖|Λ_n(L)|, Q_n = ‖c_n‖/(|δ_n(L)||d_n|^{1/2}).
Arithmetik: Mittelpunkt mit P Bit (Standard 1.3N + 400), Rang-1-Schritte als C-Matrixoperationen, alle 250 Schritte Block verkleinert.
Aufruf: python v48_RL.py kubisch|kneser N [P]
Ausgabe: v48_RL_<MODELL>_N<N>.json (n, K_n, P_n, δ_n(L), Λ_n(L)) und Kurzprotokoll.
"""
import cmath
import json
import math
import os
import sys
import time
from pathlib import Path

from flint import arb, acb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parent.parent / "gemeinsam"
MODELL, N = sys.argv[1], int(sys.argv[2])
P = int(sys.argv[3]) if len(sys.argv) > 3 else int(1.3 * N) + 400
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
    print(f"[{MODELL}, N = {N}, P = {P} Bit] J̃ aufgebaut ({time.time() - t0:.0f}s)", flush=True)
    # ẽ_L = ((ρ/L)^k), ε_L = conj(ẽ_L)
    q = acb(rho) / L
    e_re, e_im = [], []
    qk = acb(1)
    for k in range(N + 1):
        e_re.append(qk.real.mid())
        e_im.append(qk.imag.mid())
        qk = qk * q
    b = arb_mat(N + 1, 2)
    for k in range(N + 1):
        b[k, 0], b[k, 1] = e_re[k], e_im[k]
    r_re = list(e_re)                       # ε_L = conj ẽ_L: r = ε_Lᵀ
    r_im = [-x for x in e_im]
    off = 0
    m = N + 1
    w_l, y_l = [], []
    for k in range(N + 1):
        i = k - off
        p = A[i, i]
        yk = acb(b[i, 0], b[i, 1])
        wk = acb(r_re[k], r_im[k]) / p
        w_l.append(wk)
        y_l.append(yk)
        # r ← r − w_k·(Zeile k von U), Spalten > k
        for jj in range(i + 1, m):
            u = A[i, jj]
            r_re[off + jj] = (r_re[off + jj] - wk.real * u).mid()
            r_im[off + jj] = (r_im[off + jj] - wk.imag * u).mid()
        # Rang-1-Schritt auf A und b
        cvec = arb_mat(m, 1)
        rvec = arb_mat(1, m)
        for ii in range(i + 1, m):
            cvec[ii, 0] = (A[ii, i] / p).mid()
        for jj in range(i + 1, m):
            rvec[0, jj] = A[i, jj]
        A = (A - cvec * rvec).mid()
        bro = arb_mat(1, 2)
        bro[0, 0], bro[0, 1] = b[i, 0], b[i, 1]
        b = (b - cvec * bro).mid()
        if (k + 1) % 250 == 0 and k < N:
            noff = k + 1
            mn = N + 1 - noff
            A2 = arb_mat(mn, mn)
            b2 = arb_mat(mn, 2)
            d = noff - off
            for ii in range(mn):
                for jj in range(mn):
                    A2[ii, jj] = A[ii + d, jj + d]
                b2[ii, 0], b2[ii, 1] = b[ii + d, 0], b[ii + d, 1]
            A, b, off, m = A2, b2, noff, mn
            print(f"  Schritt {k + 1}/{N + 1} ({time.time() - t0:.0f}s)", flush=True)
    Lc = complex(float(L_RE), float(L_IM))
    rhof = abs(Lc)
    a_exp = 1 + 1j * Lc / (4 * math.pi)
    Ks, Ps, out = acb(0), [], []
    for n in range(N + 1):
        Pn = w_l[n] * y_l[n]
        Ks = Ks + Pn
        Pc = complex(float(Pn.real.mid()), float(Pn.imag.mid()))
        Kc = complex(float(Ks.real.mid()), float(Ks.imag.mid()))
        dc = complex(float(w_l[n].real.mid()), float(w_l[n].imag.mid()))
        lc = complex(float(y_l[n].real.mid()), float(y_l[n].imag.mid()))
        out.append([n, [Kc.real, Kc.imag], [Pc.real, Pc.imag], [dc.real, dc.imag], [lc.real, lc.imag]])
    zeilen = []
    for n in [50, 100, 200, 300, 400, 600, 800, 1000, 1200, 1500, 1800, 2000]:
        if n <= N:
            _, K_, P_, d_, l_ = out[n]
            Kc, Pc = complex(*K_), complex(*P_)
            dn = rhof * n ** (-a_exp)
            Cn = Pc * n * dn / (a_exp * Lc)
            zeilen.append(f"n = {n}: K_n = {Kc.real:.4f}{Kc.imag:+.4f}i, |K_n|/n^Re a = {abs(Kc) / n ** a_exp.real:.5f}, P_n = {Pc.real:.5f}{Pc.imag:+.5f}i, "
                          f"P_n·n·d_n/(aL) = {Cn.real:.5f}{Cn.imag:+.5f}i")
    print("\n".join(zeilen), flush=True)
    json.dump(dict(modell=MODELL, N=N, P=P, zeilen=out), open(HIER / f"v48_RL_{MODELL}_N{N}.json", "w", encoding="utf-8"))
    print(f"fertig ({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
