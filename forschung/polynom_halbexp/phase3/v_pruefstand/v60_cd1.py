"""Strang V — v60: unabhängige Nachrechnung von Ks Christoffel–Darboux-Identität (CD1) (k_kern/BERICHT §40.1; Auftrag Leitung 13:20).
  (z − t)·K_n(z,t) = q_{n+1}(z)·λ_n(t) − t·ν_n(z) + Ξ_n(z,t),
  K_n(z,t) = e(z)ᵀJ_n⁻¹ê(t), e(z) = (z^j), ê(t) = (t^{−k}); λ_n(t) = e_nᵀJ_n⁻¹ê(t); ν_n = J_n⁻¹·1; q_{n+1}(z) = z^{n+1} − (J_n⁻¹T_n f^{n+1})(z);
  Ξ_n(z,t) = J_n⁻¹T_n[(f − x)·K_n(f(x), t)](z).
Eigener Code (ungeskalierter Rahmen, J_n direkt aus den Spalten T_n[f^j + m x^j], m = f'∘f bzw. für Kneser m = eˣ/f'), Inverse per arb_mat.solve.
Ξ zweimal: (a) Matrixform (Komposition als Reihe), (b) Konturform (1/2πi)∮ K_n(z,s)(f(s) − s)K_n(f(s),t) ds/s, Trapezregel auf |s| = r_s
(r_s = 0.6ρ, M = 6n + 64 Punkte; für Kneser innerhalb des Konvergenzkreises). Zusätzlich Ks Tabellenwerte am Steg t = −ρ: m·ν_n(t), B(t) = m·q_{n+1}(t)λ_n(t)/t.
Präzisionskontrolle: Wiederholung mit +256 Bit.
Aufruf: python v60_cd1.py kubisch|kneser n1,n2,… [P]
"""
import cmath
import json
import math
import sys
import time
from pathlib import Path

from flint import arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parent.parent / "gemeinsam"
MODELL = sys.argv[1]
NS = [int(x) for x in sys.argv[2].split(",")]
L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"


def koeff(N):
    """Taylor-Koeffizienten von f und m = f'∘f bis Grad N (arb)."""
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
        return [f[k] if k <= f.degree() else arb(0) for k in range(N + 1)], [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
    z = json.loads((GEM / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))["a"]
    g = json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]
    a = [arb((z[k] if k <= 640 else g[k])[0] if isinstance(z[0], list) else (z[k] if k <= 640 else g[k])) for k in range(N + 2)]
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
    return a[:N + 1], mk


def rechne(n, P, punkte):
    ctx.prec = P
    fk, mk = koeff(n + 1)
    f = arb_poly(fk[:n + 2])
    J = arb_mat(n + 1, n + 1)
    col = arb_poly([1])
    for j in range(n + 1):
        for k in range(n + 1):
            J[k, j] = (col[k] if k <= col.degree() else arb(0)) + (mk[k - j] if k >= j else arb(0))
        col = (col * f).truncate(n + 2)
    fn1 = [col[k] if k <= col.degree() else arb(0) for k in range(n + 1)]   # T_n f^{n+1}
    Jc = acb_mat(J)
    Id = acb_mat(n + 1, n + 1)
    for k in range(n + 1):
        Id[k, k] = acb(1)
    X = Jc.solve(Id, algorithm="approx")                                     # J_n⁻¹ (Mittelpunkt), einmal
    XT = X.transpose()
    one = acb_mat(n + 1, 1)
    one[0, 0] = acb(1)
    nu = X * one
    F1 = acb_mat(n + 1, 1)
    for k in range(n + 1):
        F1[k, 0] = acb(fn1[k])
    qv = X * F1
    Jt = Jc.transpose()
    en = acb_mat(n + 1, 1)
    en[n, 0] = acb(1)
    lam_row = XT * en                              # e_nᵀJ⁻¹ als Spalte
    fa = acb_poly([acb(x) for x in fk[:n + 2]])
    fmx = fa - acb_poly([0, 1])

    def ev(vec, z):                                                          # Σ vec_j z^j
        s_ = acb(0)
        for j in range(n, -1, -1):
            s_ = s_ * z + vec[j, 0]
        return s_

    def ehat(t):
        v = acb_mat(n + 1, 1)
        tk = acb(1)
        for k in range(n + 1):
            v[k, 0] = tk
            tk = tk / t
        return v

    erg = []
    rho = abs(acb(arb(L_RE), arb(L_IM)))
    for z, t in punkte:
        et = ehat(t)
        Kv = X * et                                # Koeffizienten von K_n(·, t)
        K = ev(Kv, z)
        lam = sum((lam_row[k, 0] * et[k, 0] for k in range(n + 1)), acb(0))
        q = z ** (n + 1) - ev(qv, z)
        nuz = ev(nu, z)
        # (a) Ξ in Matrixform: T_n[(f − x)·K_n(f(x), t)]
        comp = acb_poly([0])
        for j in range(n, -1, -1):
            comp = (comp * fa + acb_poly([Kv[j, 0]])).truncate(n + 1)
        g = (fmx * comp).truncate(n + 1)
        G = acb_mat(n + 1, 1)
        for k in range(n + 1):
            G[k, 0] = g[k] if k <= g.degree() else acb(0)
        Xi = ev(X * G, z)
        lhs = (z - t) * K
        rhs = q * lam - t * nuz + Xi
        # (b) Ξ in Konturform, Trapezregel auf |s| = r_s: (1/M)Σ K_n(z,s)(f(s) − s)K_n(f(s),t), K_n(z,s) = Σ_k r_z[k]s^{−k}, r_z = J^{−T}e(z)
        ez = acb_mat(n + 1, 1)
        zk = acb(1)
        for k in range(n + 1):
            ez[k, 0] = zk
            zk = zk * z
        rz = XT * ez
        M = 6 * n + 64
        rs = arb(0.6) * rho
        Xi_b = acb(0)
        for jj in range(M):
            s = acb(rs) * acb(0, 2 * arb.pi() * jj / M).exp()
            fs = fa(s)
            Xi_b += ev(rz, 1 / s) * (fs - s) * ev(Kv, fs)
        Xi_b = Xi_b / M
        erg.append(dict(z=z, t=t, K=K, lam=lam, q=q, nu=nuz, Xi=Xi, Xi_b=Xi_b, lhs=lhs, rhs=rhs))
    return erg, Jc, nu, qv, lam_row, fa, fk, mk


def main():
    P0 = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    rho = math.hypot(float(L_RE[:20]), float(L_IM[:20]))
    pkt = [(0.5 * rho * cmath.exp(0.7j), 0.8 * rho * cmath.exp(2.6j)),          # beide in D
           (1.03 * rho * cmath.exp(math.radians(175) * 1j), -1.03 * rho),       # Steg, ausserhalb D
           (1.1 * rho * cmath.exp(2.0j), 1.05 * rho * cmath.exp(-2.5j)),        # beide ausserhalb D
           (-rho, -rho),                                                         # Diagonale auf ∂D (Steg)
           (0.3 + 0.2j, 1.2 * rho * cmath.exp(0.3j))]                            # Γout-Seite ausserhalb D
    for n in NS:
        P = P0 or int(1.5 * n) + 300
        t0 = time.time()
        res = {}
        for PP in (P, P + 256):
            ctx.prec = PP
            punkte = [(acb(z.real, z.imag), acb(t.real, t.imag)) for z, t in pkt]
            res[PP] = rechne(n, PP, punkte)
        e1, Jc, nu, qv, lam_row, fa, fk, mk = res[P]
        e2 = res[P + 256][0]
        print(f"[{MODELL}] n = {n}, P = {P} Bit ({time.time() - t0:.0f}s)", flush=True)
        a3 = e1[3]
        t3 = a3["t"]
        if MODELL == "kubisch":
            f_ = arb_poly(fk[:4])
            mt = acb_poly([acb(x) for x in f_.derivative().coeffs()])(acb_poly([acb(x) for x in fk[:4]])(t3))
        else:
            mt = acb(1)                                                      # Kneser: Werte ohne den Faktor m(t)
            print("  (Kneser: Steg-Werte ohne Faktor m(t))", flush=True)
        if True:
            print(f"  Steg t = −ρ: m·ν_n(t) = {complex(float((mt * a3['nu']).real.mid()), float((mt * a3['nu']).imag.mid())):.6f}, "
                  f"B(t) = m·q_(n+1)(t)λ_n(t)/t = {complex(float((mt * a3['q'] * a3['lam'] / t3).real.mid()), float((mt * a3['q'] * a3['lam'] / t3).imag.mid())):.5f}, "
                  f"m·Ξ_n(t,t)/t = {complex(float((mt * a3['Xi'] / t3).real.mid()), float((mt * a3['Xi'] / t3).imag.mid())):.5f}", flush=True)
        for a_, b_ in zip(e1, e2):
            d = abs(a_["lhs"] - a_["rhs"])
            gr = max(float(abs(a_["lhs"]).mid()), float(abs(a_["q"] * a_["lam"]).mid()), float(abs(a_["t"] * a_["nu"]).mid()), float(abs(a_["Xi"]).mid()))
            st = float(abs(a_["K"] - b_["K"]).mid()) / max(float(abs(a_["K"]).mid()), 1e-300)
            z, t = complex(float(a_["z"].real.mid()), float(a_["z"].imag.mid())), complex(float(a_["t"].real.mid()), float(a_["t"].imag.mid()))
            print(f"  z = {abs(z) / rho:.3f}ρ∠{math.degrees(cmath.phase(z)):.0f}°, t = {abs(t) / rho:.3f}ρ∠{math.degrees(cmath.phase(t)):.0f}°: "
                  f"|LHS − RHS| = {float(d.mid()):.1e}, Grösse der Terme {gr:.2e}, rel. {float(d.mid()) / gr:.1e}; Präzisionsstabilität K {st:.0e}; "
                  f"|K| = {float(abs(a_['K']).mid()):.3e}, |Ξ| = {float(abs(a_['Xi']).mid()):.3e}, Ξ Matrix/Kontur − 1 = {float(abs(a_['Xi_b'] / a_['Xi'] - 1).mid()):.0e}", flush=True)


if __name__ == "__main__":
    main()
