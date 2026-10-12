"""k55 — Vortest (b): Lemma ρ für die zusammengesetzte Parametrix der Eck-Pol-Antwort (BERICHT §37). Kubisches Modell.
u_n = J_n⁻¹T_nk_L (k_L = 1/(L − x)) ist exakt Π_n s mit s = J_∞⁻¹k_L = c_s/(x − L) + s_reg, c_s = −1/(μ + 1/μ).
Parametrizes (ρ-skalierte Koeffizientenvektoren in ℂ^{n+1}):
  P0: ũ = 0                         (Kontrolle: ρ ≡ 1)
  P1: ũ = T_n s = T_n[s_reg] + T_n[c_s/(x − L)]       (nur äussere Lösung)
  P2: ũ = T_n[s_reg] + T_n[(C/d_n)·V((x − L)/d_n)]    (äussere Lösung + inneres Profil statt Pol)
s_reg ≈ Π_{2n} s_reg = J_{2n}⁻¹T_{2n}g, g = J s_reg = 1/(L − x) − c_s[1/(f − L) + m/(x − L)] (regulär in L).
Lemma ρ: ρ_{m,n} = ⟨ℓ_m, r̃⟩_m/⟨ℓ_m, b⟩_m = (L⁻¹r̃)_m/(L⁻¹b)_m (J̃_n = LU ohne Pivot), r̃ = b − J̃_nũ, b = ρ-Koeff. von 1/(L − x).
Ausgabe: max_m|ρ_{m,n}| und Lage, Profil |ρ_{m,n}| über m/n, Fehler in L: e_n = ε_Lᵀ(J̃⁻¹b − ũ) relativ zu u_n(L).
Aufruf: python k55_lemma_rho.py 200 [P2-Parameter: Rc-Faktor]
"""
import math
import os
import sys
import time

import mpmath as mp

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, acb_mat, ctx  # noqa: E402
from k_kern import J_tilde, L_RE, L_IM, kubik_koeff, modell  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)


def acb2mp(z, nd):
    return mp.mpc(mp.mpf(z.real.mid().str(nd, radius=False)), mp.mpf(z.imag.mid().str(nd, radius=False)))


def main():
    n = int(sys.argv[1])
    t0 = time.time()
    bits = int(1.3 * 2 * n) + 400
    ctx.prec = bits
    mp.mp.dps = int(bits / 3.32) + 10
    nd = mp.mp.dps
    L = mp.mpc(L_RE, L_IM); Lb = mp.conj(L); rho = abs(L)
    al = L / (4 * mp.pi); a = 1 + 1j * al; mu = mp.sqrt(L)
    cs = -1 / (mu + 1 / mu)
    bP = 0.5 + 1j * al
    lamP = mp.mpc("0.603424", "1.165020"); lamPt = -lamP * rho ** 2 / L ** 2
    Ck = -lamPt * mp.gamma(0.5 + 1j * al) / ((mu + 1 / mu) * mp.gamma(0.5 + 2j * al))
    dn = rho * mp.power(n, -a)
    # ---- g = J s_reg, Taylor-Koeffizienten bis N = 2n (exakt über Rekursionen, mpmath)
    N = 2 * n
    ck = [mp.mpc(str(x.mid())) if False else mp.mpf(x.mid().str(nd, radius=False)) for x in kubik_koeff()]
    alpha = (ck[2] + 2 * mp.re(L) * ck[3])          # aus kubik_koeff: c2 = α − 2Re(L)β, c3 = β
    beta = ck[3]
    # Q(x) = 1 + (x − L̄)(α + βx) = (1 − L̄α) + (α − L̄β)x + βx²
    Q = [1 - Lb * alpha, alpha - Lb * beta, beta]
    invQ = [mp.mpc(0)] * (N + 1)
    invQ[0] = 1 / Q[0]
    for k in range(1, N + 1):
        s_ = Q[1] * invQ[k - 1] + (Q[2] * invQ[k - 2] if k >= 2 else 0)
        invQ[k] = -s_ / Q[0]
    invxL = [-1 / L ** (k + 1) for k in range(N + 1)]      # 1/(x − L)
    B = [mp.fsum(invxL[j] * invQ[k - j] for j in range(k + 1)) for k in range(N + 1)]   # 1/(f − L) = 1/((x−L)Q)
    fk, mk = modell("kubisch", 12)
    m = [mp.mpf(x.mid().str(nd, radius=False)) for x in mk[:7]]
    while len(m) < 7:
        m.append(mp.mpf(0))
    mL = mp.polyval(m[::-1], L)
    # (m(x) − m(L))/(x − L): synthetische Division
    qd = [mp.mpc(0)] * 6
    acc = mp.mpc(0)
    coeffs_desc = m[::-1]   # Grad 6 … 0
    quot = []
    for c in coeffs_desc[:-1]:
        acc = acc * L + c
        quot.append(acc)
    quot = quot[::-1]       # aufsteigend, Grad 0 … 5
    Cm = [(quot[k] if k < len(quot) else 0) + mL * invxL[k] for k in range(N + 1)]
    A = [1 / L ** (k + 1) for k in range(N + 1)]            # 1/(L − x)
    g = [A[k] - cs * (B[k] + Cm[k]) for k in range(N + 1)]
    print("n = %d: g-Koeffizienten (|g_k|ρ^k bei k = 0, N/2, N: %s, %s, %s) (%.0f s)" % (n, mp.nstr(abs(g[0]), 3), mp.nstr(abs(g[N // 2]) * rho ** (N // 2), 3), mp.nstr(abs(g[N]) * rho ** N, 3), time.time() - t0))
    # ---- s_reg ≈ J̃_{2n}⁻¹ (ρ-skaliert g)
    J2 = acb_mat(J_tilde("kubisch", N, bits=bits))
    ctx.prec = bits
    gv = acb_mat(N + 1, 1)
    for k in range(N + 1):
        v = g[k] * rho ** k
        gv[k, 0] = acb(arb(mp.nstr(mp.re(v), nd)), arb(mp.nstr(mp.im(v), nd)))
    sreg = J2.solve(gv, algorithm="approx")
    sr = [acb2mp(sreg[k, 0], nd) for k in range(n + 1)]
    print("   s_reg (Galerkin 2n) fertig (%.0f s); |s_reg,k|ρ^k bei k = 0, n/2, n: %s %s %s" % (time.time() - t0, mp.nstr(abs(sr[0]), 3), mp.nstr(abs(sr[n // 2]), 3), mp.nstr(abs(sr[n]), 3)))
    # ---- J̃_n, b, exakte Lösung, LU ohne Pivot (mpmath)
    Jn = acb_mat(J_tilde("kubisch", n, bits=bits))
    ctx.prec = bits
    bvec = [rho ** k / L ** (k + 1) for k in range(n + 1)]
    Jm = [[acb2mp(Jn[i, j], nd) for j in range(n + 1)] for i in range(n + 1)]
    bv = acb_mat(n + 1, 1)
    for k in range(n + 1):
        bv[k, 0] = acb(arb(mp.nstr(mp.re(bvec[k]), nd)), arb(mp.nstr(mp.im(bvec[k]), nd)))
    ue = Jn.solve(bv, algorithm="approx")
    uex = [acb2mp(ue[k, 0], nd) for k in range(n + 1)]
    # Doolittle: nur L wird gebraucht (unteres Dreieck, Diagonale 1)
    Aw = [row[:] for row in Jm]
    Lf = [[mp.mpc(0)] * (n + 1) for _ in range(n + 1)]
    for k in range(n + 1):
        Lf[k][k] = mp.mpc(1)
        piv = Aw[k][k]
        rowk = Aw[k]
        for i in range(k + 1, n + 1):
            lik = Aw[i][k] / piv
            Lf[i][k] = lik
            if lik != 0:
                rowi = Aw[i]
                for j in range(k + 1, n + 1):
                    rowi[j] -= lik * rowk[j]
    print("   LU fertig (%.0f s)" % (time.time() - t0))

    def fwd(y):
        z = [mp.mpc(0)] * (n + 1)
        for i in range(n + 1):
            s_ = y[i]
            Li = Lf[i]
            for j in range(i):
                s_ -= Li[j] * z[j]
            z[i] = s_
        return z
    yb = fwd(bvec)

    def matvec(v):
        return [mp.fsum(Jm[i][j] * v[j] for j in range(n + 1)) for i in range(n + 1)]

    epsL = [(L / rho) ** k for k in range(n + 1)]
    uL = mp.fsum(epsL[k] * uex[k] for k in range(n + 1))

    def auswerten(name, ut):
        r = [bvec[i] - v for i, v in enumerate(matvec(ut))]
        yr = fwd(r)
        rh = [yr[m_] / yb[m_] for m_ in range(n + 1)]
        ab = [abs(x) for x in rh]
        mx = max(range(n + 1), key=lambda i: ab[i])
        prof = []
        for q in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0):
            i = min(n, int(round(q * n)))
            prof.append("%.2f:%s" % (q, mp.nstr(ab[i], 3)))
        eL = mp.fsum(epsL[k] * (uex[k] - ut[k]) for k in range(n + 1))
        # Differenz ρ_{m,n} − ρ_{m,n} geht erst mit n−1; hier Profil
        print("   %s: max_m|ρ_{m,n}| = %s bei m/n = %.3f | Profil |ρ| (m/n): %s | Fehler in L: |e_n|/|u_n(L)| = %s"
              % (name, mp.nstr(ab[mx], 4), mx / n, " ".join(prof), mp.nstr(abs(eL) / abs(uL), 4)), flush=True)
        return rh

    auswerten("P0 (ũ = 0)", [mp.mpc(0)] * (n + 1))
    auswerten("Kontrolle (ũ = u_n exakt)", uex)
    # P1: T_n s = s_reg + c_s/(x − L)
    pol = [cs * invxL[k] * rho ** k for k in range(n + 1)]
    P1 = [sr[k] + pol[k] for k in range(n + 1)]
    auswerten("P1 (T_n s, nur aussen)", P1)
    # P2: T_n s_reg + T_n[(C/d_n) V((x − L)/d_n)] — Koeffizienten per FFT auf |x| = Rc
    cP = []
    K = 0
    fac = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
    Rc = rho * (1 - fac / n)
    M = 6 * n
    t1 = time.time()
    # Reihe von V: Σ c_k (−λ̃ζ)^k/(k+1), c_k = Γ(β − iαk)/(k!Γ(β)); Termzahl adaptiv
    gb = mp.gamma(bP)

    def V(z):
        w = -lamPt * z
        s_ = mp.mpc(0); p = mp.mpc(1); k = 0
        big = mp.mpf(0)
        while True:
            c = mp.gamma(bP - 1j * al * k) / (mp.factorial(k) * gb) * p / (k + 1)
            s_ += c
            ac = abs(c)
            big = max(big, ac)
            if k > 30 and ac < mp.mpf(10) ** (-nd + 5) * (1 + big):
                break
            p *= w
            k += 1
            if k > 20000:
                raise RuntimeError("V-Reihe konvergiert nicht")
        return s_
    vals = []
    for j in range(M):
        x = Rc * mp.expjpi(mp.mpf(2 * j) / M)
        vals.append(Ck / dn * V((x - L) / dn))
    print("   V auf %d Punkten (Rc = ρ(1 − %.0f/n)) ausgewertet (%.0f s); max|V-Teil| = %s" % (M, fac, time.time() - t1, mp.nstr(max(abs(v) for v in vals), 4)))
    # Koeffizienten a_k = (1/M) Σ v_j x_j^{−k}; ρ-skaliert: a_k ρ^k
    hk = []
    for k in range(n + 1):
        s_ = mp.fsum(vals[j] * mp.expjpi(-mp.mpf(2 * j * k) / M) for j in range(M))
        hk.append(s_ / M * (rho / Rc) ** k)
    P2 = [sr[k] + hk[k] for k in range(n + 1)]
    auswerten("P2 (T_n s_reg + T_n[(C/d_n)V])", P2)
    print("fertig (%.0f s)" % (time.time() - t0))


if __name__ == "__main__":
    main()
