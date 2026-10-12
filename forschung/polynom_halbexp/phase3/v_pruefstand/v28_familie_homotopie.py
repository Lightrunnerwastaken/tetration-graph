"""Strang V — (S_γ) für alles, was H1/H2 braucht (Auftrag Leitung 23:20, Punkt 2), Kneser, Linse wie v12:
  (a) Sekante 𝒥_P (v ↦ v∘P + M_P v, M_P = (f∘P − f∘f)/(P − f)) bei P = p_N;
  (b) λ-Familie aus H1: P(λ) = (1−λ)f_K + Q(λ), Q ∈ P_N löst T_N[P∘P − eˣ] = 0 (Newton ab Q = λp_N); dort Tangente J_N(P) und Sekante 𝒥_P;
  (c) σ-Homotopie aus H2: P ∈ P_{N+1} löst T_{N+1}[P∘P − eˣ] = (1−s)c_{N+1}x^{N+1}, P(0) = pad(p_N), P(1) = p_{N+1}; dort Tangente J_{N+1}(P).
Gemessen jeweils ‖Ĵ⁻¹‖ der konformen Fassung T_N^w[(v∘P + m v)∘k] (m = P'∘P bzw. M_P), Vergleich f_K.
Alle Reihen modulo x^{n+1} in arb (arb_poly mit expliziter Abschneidung).
Fassung b (23:35): Kompositionen mit f_K NICHT über Σ_k a_k P^k (das Umentwickeln der Taylor-Reihe um P(0) ≈ 0.5 verstärkt die Fehler der
a_k für k > 500 (relativ 1e-21) um bis zu 1e144 — Fassung a ergab ab N = 400 Unsinn), sondern über die Taylor-Entwicklung um f(x):
  f∘P = eˣ + Σ_{n=1}^{5} g_n η^n/n!,  f'∘P = Σ_{n=0}^{4} g_{n+1} η^n/n!,  M_P = (f∘P − f∘f)/(P − f) = Σ_{n=1}^{5} g_n η^{n−1}/n!,
mit η = P − f (mod x^{n+1} = P − T_n f), g_1 = m = eˣ/f', g_{k+1} = g_k'/f' (alles exakte Reihen aus T_{n+8} f); Restglied ~ η^6, η ~ 1e-10.
Aufruf: python v28_familie_homotopie.py N S [teile: sek,lam,sig] [lambdas]
"""
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

ARGS = sys.argv[1:]
N = int(ARGS[0])
S = ARGS[1]
TEILE = (ARGS[2] if len(ARGS) > 2 else "sek,lam,sig").split(",")
LAMS = [float(x) for x in (ARGS[3] if len(ARGS) > 3 else "0.25,0.5,0.75").split(",")]
sys.argv = ["v20", str(N), S]
import v20_pN_S as V20  # noqa: E402
V12 = V20.V12

sys.stdout.reconfigure(encoding="utf-8")
KMAX = 4000


def tr(p, n):
    return p.truncate(n + 1)


def koeff(p, n):
    return [p[i] if i <= p.degree() else arb(0) for i in range(n + 1)]


def l1(v, rho, n):
    return sum(abs(float(v[i].mid())) * float(rho.mid()) ** i for i in range(n + 1))


def ableitungen_f_um_f(a, n, nmax=6):
    """g_k = f^{(k)}∘f mod x^{n+1}, k = 1..nmax: g_1 = eˣ/f', g_{k+1} = g_k'/f' (Reihen bis Grad n + nmax, dann abgeschnitten)."""
    ne = n + nmax + 1
    f = arb_poly(a[:ne + 2])
    fp = f.derivative()
    ex = arb_poly([arb(1) / arb.fac_ui(k) for k in range(ne + 1)])
    g = [None, reihen_div(ex, fp, ne)]
    for k in range(1, nmax):
        g.append(reihen_div(g[k].derivative(), fp, ne - k))
    return [None] + [tr(gk, n) for gk in g[1:]]


def reihen_div(zaehler, nenner, n):
    """zaehler/nenner mod x^{n+1} (nenner(0) ≠ 0)."""
    zc = koeff(zaehler, n)
    nc = koeff(nenner, n)
    q = [arb(0)] * (n + 1)
    for i in range(n + 1):
        s_ = zc[i]
        for j in range(i):
            s_ -= q[j] * nc[i - j]
        q[i] = s_ / nc[0]
    return arb_poly(q)


def f_um_f(G, eta, n, ex_poly):
    """(f∘P, f'∘P, M_P) mod x^{n+1} aus g_k = G[k] und η = P − f."""
    fP = ex_poly
    fpP = G[1]
    MP = G[1]
    ep = arb_poly([1])
    fak = arb(1)
    for k in range(1, 6):
        ep = tr(ep * eta, n)                          # η^k
        fak = fak * k                                 # k!
        fP = fP + tr(G[k] * ep, n) * (1 / fak)
        if k + 1 < len(G) and k <= 4:
            fpP = fpP + tr(G[k + 1] * ep, n) * (1 / fak)
        MP = MP + tr(G[k + 1] * ep, n) * (1 / (fak * (k + 1))) if k + 1 < len(G) else MP
    return tr(fP, n), tr(fpP, n), tr(MP, n)


def poly_komp(qc, P, n):
    """Q∘P mod x^{n+1} für ein Polynom Q mit Koeffizientenliste qc (Horner)."""
    r = arb_poly([0])
    for c in reversed(qc):
        r = tr(r * P, n) + c
    return r


def norm_w(acoef, mcoef, n, F, bits):
    V12.modell_koeff = lambda n_: (acoef, mcoef)
    J = V12.J_hat(n, F, bits)
    ctx.prec = bits
    X = V20.zu_np(V20.inverse(J))
    return float(np.linalg.norm(X, 2))


def newton(residuum, jacobi, x0, n, rho, it_max=8, tol=1e-45):
    x = x0
    for it in range(it_max):
        R = residuum(x)
        J = jacobi(x)
        B = arb_mat(n + 1, 1)
        for i in range(n + 1):
            B[i, 0] = -R[i]
        d = J.solve(B, algorithm="approx")
        dx = arb_poly([d[i, 0] for i in range(n + 1)])
        x = x + dx
        nd = l1(koeff(dx, n), rho, n)
        if nd < tol:
            return x, it + 1, nd
    return x, it_max, nd


def jacobi_tangente(P, mP, n):
    """J_N(P)δ = T_n[δ∘P + (P'∘P)δ] als Matrix in Monom-Koordinaten (Spalten x^j)."""
    J = arb_mat(n + 1, n + 1)
    pw = arb_poly([1])
    mc = koeff(mP, n)
    for j in range(n + 1):
        if j > 0:
            pw = tr(pw * P, n)
        c = koeff(pw, n)
        for i in range(n + 1):
            J[i, j] = c[i] + (mc[i - j] if i >= j else 0)
    return J


def main():
    F, it, d = V12.theodorsen()
    import os
    bits = int(os.environ.get("V28_BITS", int(1.2 * (N + 1)) + 400))
    ctx.prec = bits
    rho = abs(V12.L_arb())
    a = [arb(x) for x in V20.AK[:KMAX + 1]]
    pN = [arb(x) for x in V20.PN[str(N)]]
    n = N
    f = arb_poly(a[:n + 1])
    ex = [arb(1)]
    for k in range(1, n + 3):
        ex.append(ex[-1] / k)
    # Referenz f_K
    _, mf = V20.koeff_fK(n)
    nf = norm_w(a[:n + 1], mf, n, F, bits)
    print(f"Kneser N = {N}, Linse S = {S}: Referenz ‖(Ĵ^fK)⁻¹‖ = {nf:.6f}  [{bits} Bit]", flush=True)
    t0 = time.time()
    G = ableitungen_f_um_f(a, n)
    ex_poly = arb_poly(ex[:n + 1])
    fT = arb_poly(a[:n + 1])
    # Kontrolle: m aus g_1 gegen v20 (m·f' = eˣ)
    dm = l1([G[1][i] - mf[i] for i in range(n + 1)], rho, n)
    print(f"Kontrolle g_1 = m: Abweichung {dm:.1e}", flush=True)
    if "sek" in TEILE:
        t0 = time.time()
        ctx.prec = bits
        P = arb_poly(pN)
        eta = tr(P - fT, n)
        fP, fpP, MP = f_um_f(G, eta, n, ex_poly)
        ns = norm_w(pN, koeff(MP, n), n, F, bits)
        # Gegenprobe: Tangente p_N'∘p_N exakt vs (f'∘P + η'∘P)
        print(f"(a) Sekante bei P = p_N: ‖(Ĵ_sek)⁻¹‖ = {ns:.6f}; ‖η‖ = {l1(koeff(eta, n), rho, n):.2e}  [{time.time() - t0:.0f} s]", flush=True)
    if "lam" in TEILE:
        for lam in LAMS:
            t0 = time.time()
            ctx.prec = bits
            L_ = arb(lam)

            def teile(Q):
                P = tr(fT * (1 - L_) + Q, n)
                eta = tr(Q - fT * L_, n)
                fP, fpP, MP = f_um_f(G, eta, n, ex_poly)
                qc = koeff(Q, n)
                QP = poly_komp(qc, P, n)
                dQ = Q.derivative()
                QpP = poly_komp(koeff(dQ, n - 1) if n >= 1 else [arb(0)], P, n)
                PP = tr(fP * (1 - L_) + QP, n)
                mP = tr(fpP * (1 - L_) + QpP, n)
                return P, PP, mP, MP

            def resid(Q):
                _, PP, _, _ = teile(Q)
                return [(PP[i] if i <= PP.degree() else arb(0)) - ex[i] for i in range(n + 1)]

            def jac(Q):
                P, _, mP, _ = teile(Q)
                return jacobi_tangente(P, mP, n)

            Q0 = arb_poly([L_ * x for x in pN])
            Q, its, nd = newton(resid, jac, Q0, n, rho)
            P, PP, mP, MP = teile(Q)
            nt = norm_w(koeff(P, n), koeff(mP, n), n, F, bits)
            ns = norm_w(koeff(P, n), koeff(MP, n), n, F, bits)
            abst = l1([Q[i] - L_ * pN[i] if i <= Q.degree() else -L_ * pN[i] for i in range(n + 1)], rho, n)
            print(f"(b) λ = {lam:.2f}: Newton {its} Schritte (letzter Schritt {nd:.1e}), ‖Q − λp_N‖ = {abst:.2e}; ‖(Ĵ_tan(P))⁻¹‖ = {nt:.6f}, "
                  f"‖(Ĵ_sek(P))⁻¹‖ = {ns:.6f}  [{time.time() - t0:.0f} s]", flush=True)
    if "sig" in TEILE:
        M_ = n + 1
        ctx.prec = bits
        P0 = arb_poly([arb(x.mid()) for x in pN])
        if os.environ.get("V28_VERFEINERN", "") == "1":
            # p_N auf ~0.4N Stellen nachschärfen (Newton für T_N[P∘P − eˣ] = 0): sonst ist c_(N+1) von der verstärkten Rundung beherrscht
            def resid_n(P):
                r = arb_poly([0])
                for cc in reversed(koeff(P, n)):
                    r = tr(r * P, n) + cc
                return [(r[i] if i <= r.degree() else arb(0)) - ex[i] for i in range(n + 1)]

            def jac_n(P):
                dP = P.derivative()
                mm = arb_poly([0])
                for cc in reversed(koeff(dP, n - 1)):
                    mm = tr(mm * P, n) + cc
                return jacobi_tangente(P, mm, n)

            P0, its0, nd0 = newton(resid_n, jac_n, P0, n, rho, it_max=6, tol=10.0 ** (-(0.45 * n + 40)))
            P0 = arb_poly([arb(P0[i].mid()) for i in range(n + 1)])
            print(f"    p_N nachgeschärft: {its0} Schritte, letzter {nd0:.1e}", flush=True)
        PP0 = tr(P0.compose(P0), M_) if False else None
        # P0∘P0 mod x^{M+1} (Polynom-Komposition, Horner)
        r = arb_poly([0])
        for cc in reversed(koeff(P0, n)):
            r = tr(r * P0, M_) + cc
        cM = r[M_] - ex[M_]
        print(f"(c) σ-Homotopie N → N+1: c_(N+1) = {float(cM.mid()):.3e}", flush=True)
        for s in (0.0, 0.5, 1.0):
            t0 = time.time()
            ctx.prec = bits
            sA = arb(s)

            def resid_s(P):
                r = arb_poly([0])
                for cc in reversed(koeff(P, M_)):
                    r = tr(r * P, M_) + cc
                out = [(r[i] if i <= r.degree() else arb(0)) - ex[i] for i in range(M_ + 1)]
                out[M_] = out[M_] - (1 - sA) * cM
                return out

            def jac_s(P):
                dP = P.derivative()
                m = arb_poly([0])
                for cc in reversed(koeff(dP, M_ - 1) if dP.degree() >= 0 else [arb(0)]):
                    m = tr(m * P, M_) + cc
                return jacobi_tangente(P, m, M_)

            P, its, nd = newton(resid_s, jac_s, P0, M_, rho)
            dP = P.derivative()
            m = arb_poly([0])
            for cc in reversed(koeff(dP, M_ - 1)):
                m = tr(m * P, M_) + cc
            nt = norm_w(koeff(P, M_), koeff(m, M_), M_, F, int(1.2 * M_) + 400)
            print(f"    s = {s:.1f}: Newton {its} Schritte (letzter {nd:.1e}); ‖(Ĵ_tan,N+1(P(s)))⁻¹‖ = {nt:.6f}  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
