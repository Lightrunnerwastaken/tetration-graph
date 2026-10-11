"""Strang O — Kernmodul: Modelle f, Matrizen J_N, M_N, K_N, Y_N = (I + K_N)^{-1} in ℓ²_ρ (ρ = |L|), alles in arb.

Bezeichnungen (D §11.1): J_N = P_N C_f P_N + M_N, M_N = untere Dreiecks-Toeplitz von m = f'∘f,
K_N = M_N^{-1} P_N C_f P_N = P_N K P_N, Kδ = (δ∘f)/m,  Y_N := (I + K_N)^{-1} = J_N^{-1} M_N.
Umskaliert: Ã = D A D^{-1}, D = diag(ρ^k)  (Orthonormalbasis x^k/ρ^k von H²(|x| < ρ)).

Modelle (alle reell, Fixpunkte L, L̄ von exp, ρ = |L|):
  'kubisch'  : f(x) = x + (x² − 2Re(L)x + |L|²)(α + βx), f'(L) = μ (Standard μ = √L)  [D §9.3]
  'moeb_ell' : elliptische Möbius-Drehung um L, L̄ mit g(0) = 1 (Strang E „Möbius-exp“), |f'(L)| = 1
  'moeb_lox' : loxodromische Möbius-Abbildung f = M^{-1}(μ M), M = (x−L)/(x−L̄)  (f'(L) = μ, f'(L̄) = 1/μ) [D d24]
  'affin'    : f(x) = L + μ(x − L)  (nur Fixpunkt L; P_N invariant — kein Abschneide-Effekt)
"""
import sys

import numpy as np
from flint import arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ctx

L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"
Lc = complex(0.31813150520476413, 1.3372357014306895)
RHO = abs(Lc)


def Lacb():
    return acb(arb(L_RE), arb(L_IM))


# ---------------------------------------------------------------- Modelle (Taylor-Koeffizienten bis Grad n)

def kubik_koeff(mu=None):
    """Koeffizienten (komplexe floats, für numpy) und exakte arb-Koeffizienten des kubischen Modells (reell)."""
    L = Lacb()
    if mu is None:
        mu = L.sqrt()
    rL = (mu - 1) / (L - L.conjugate())
    beta = rL.imag / L.imag
    alpha = rL.real - beta * L.real
    R2 = L.real ** 2 + L.imag ** 2
    c = [R2 * alpha, 1 + R2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]
    return c          # arb-Liste, reell


def radial_koeff(alpha=None):
    """Strang Ms radiales Quadrat f(x) = x + α(x−L)(x−L̄) (auf |x| = ρ radial: f = s(φ)x, s = 1 + 2αρ(cos φ − cos θ)).
    Standard α = 0.31382/ρ (|f'(L)| = |√L| = 1.1717)."""
    A, B = arb(L_RE), arb(L_IM)
    r2 = A * A + B * B
    if alpha is None:
        alpha = arb("0.31382") / r2.sqrt()
    alpha = arb(alpha)
    return [alpha * r2, 1 - 2 * alpha * A, alpha]


def moeb_ell_par():
    """Elliptische Drehung um L, L̄ mit Winkel c, so dass g = f∘f: g(0) = 1. f = (αx+β)/(γx+δ) (arb)."""
    A, B = arb(L_RE), arb(L_IM)
    c = ((1 - A) / B).atan() - (-A / B).atan()
    h = (c / 2).tan()
    al, be = B - A * h, h * (A * A + B * B)
    ga, de = -h, B + h * A
    return al, be, ga, de, c


def moeb_lox_par(mu=None):
    """f = M^{-1}(μ M(x)), M = (x − L)/(x − L̄): f(x) = (a x + b)/(c x + d) (acb, NICHT reell für arg μ ≠ 0)."""
    L = Lacb()
    Lb = L.conjugate()
    if mu is None:
        mu = L.sqrt()
    # M = (x−L)/(x−L̄) ; y = M^{-1}(w) = (L − L̄ w)/(1 − w); w = μ (x−L)/(x−L̄)
    # y = (L(x−L̄) − L̄ μ (x−L)) / ((x−L̄) − μ(x−L)) = ((L − L̄μ) x + (−L L̄ + L̄ μ L)) / ((1−μ) x + (−L̄ + μ L))
    a = L - Lb * mu
    b = -L * Lb + Lb * mu * L
    c = 1 - mu
    d = -Lb + mu * L
    return a, b, c, d


def rational_taylor(a, b, c, d, n, komplex=False):
    """Taylor-Koeffizienten von (a x + b)/(c x + d) bis Grad n."""
    geo = [(-c / d) ** k / d for k in range(n + 1)]
    return [b * geo[0]] + [b * geo[k] + a * geo[k - 1] for k in range(1, n + 1)]


def modell(name, N, mu=None):
    """Gibt (fk, mk, komplex) zurück: Taylor-Koeff. von f und m = f'∘f bis Grad N (arb oder acb)."""
    if name in ("kubisch", "radial"):
        c = kubik_koeff(mu) if name == "kubisch" else radial_koeff(mu)
        f = arb_poly(c)
        fp = f.derivative()
        # m = f'(f(x)) exakt (Polynom Grad 6)
        m = arb_poly([0])
        pw = arb_poly([1])
        for k in range(fp.degree() + 1):
            m += fp[k] * pw
            pw = pw * f
        fk = [f[k] if k <= f.degree() else arb(0) for k in range(N + 1)]
        mk = [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
        return fk, mk, False
    if name in ("kneser", "kneser_z"):
        # Knesers halbe e-Funktion aus gemeinsam/kneser_halbexp_koeff_1000.json (abs. Fehler ~4e-47), m = e^x / f'(x) (exakt: f'∘f · f' = e^x)
        # 'kneser_z': Ps zertifizierte Koeffizienten (Satz P-V, Fehler·ρⁿ ≤ 2.9e-127, Mittelpunkte mit 140 Stellen), nur N ≤ 639
        import json, os
        fn = "kneser_halbexp_koeff_zert_640.json" if name == "kneser_z" else "kneser_halbexp_koeff_1000.json"
        d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gemeinsam", fn), encoding="utf-8"))
        roh = d["a"][:N + 2]
        a = [arb(x[0] if isinstance(x, list) else x) for x in roh]
        fk = a[:N + 1]
        fp = [a[k + 1] * (k + 1) for k in range(N + 1)]
        ex = [arb(1)]
        for k in range(1, N + 1):
            ex.append(ex[-1] / k)
        mk = [arb(0)] * (N + 1)          # m·f' = e^x  ⇒  m_k = (e_k − Σ_{j<k} m_j fp_{k−j}) / fp_0
        for k in range(N + 1):
            s_ = ex[k]
            for j in range(k):
                s_ -= mk[j] * fp[k - j]
            mk[k] = s_ / fp[0]
        return fk, mk, False
    if name == "affin":
        L = Lacb()
        mu_ = L.sqrt() if mu is None else mu
        fk = [L * (1 - mu_), mu_] + [acb(0)] * (N - 1)
        mk = [mu_] + [acb(0)] * N
        return fk, mk, True
    if name in ("moeb_ell", "moeb_lox", "moeb_sym"):
        if name == "moeb_ell":
            a, b, c, d, _ = moeb_ell_par()
            komplex = False
        elif name == "moeb_sym":
            tau = (arb(mu) / 2).tan()
            a, b, c, d = arb(1), tau, -tau, arb(1)
            komplex = False
        else:
            a, b, c, d = moeb_lox_par(mu)
            komplex = True
        fk = rational_taylor(a, b, c, d, N)
        # m = f'(f(x)) = (ad − bc)/(c f(x) + d)²
        det = a * d - b * c
        P = acb_poly if komplex else arb_poly
        fpoly = P(fk)
        a0 = fk[0]
        q0 = c * a0 + d
        hser = fpoly - P([a0])
        qinv, hp = P([0]), P([1])
        for nn in range(N + 1):
            qinv += hp * ((-c) ** nn / q0 ** (nn + 1))
            hp = (hp * hser).truncate(N + 1)
        qinv = qinv.truncate(N + 1)
        m = (qinv * qinv).truncate(N + 1) * det
        mk = [m[k] for k in range(N + 1)]
        return fk, mk, komplex
    raise ValueError(name)


# ---------------------------------------------------------------- Matrizen

def matrizen(name, N, bits=None, mu=None):
    """Umskalierte J̃_N, M̃_N (arb_mat oder acb_mat) für Modell name."""
    if bits is None:
        bits = int(1.0 * N) + 250
    ctx.prec = bits
    fk, mk, komplex = modell(name, N, mu)
    P = acb_poly if komplex else arb_poly
    Mat = acb_mat if komplex else arb_mat
    num = acb if komplex else arb
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt() if name != "moeb_sym" else arb(1)
    f = P(fk)
    pw = [P([1])]
    for _ in range(N):
        pw.append((pw[-1] * f).truncate(N + 1))
    rp = [rho ** k for k in range(N + 1)]
    J = Mat(N + 1, N + 1)
    M = Mat(N + 1, N + 1)
    for j in range(N + 1):
        col = pw[j]
        for k in range(N + 1):
            s = rp[k] / rp[j]
            ck = col[k] if k <= col.degree() else num(0)
            mm = mk[k - j] if k >= j else num(0)
            J[k, j] = (ck + mm) * s
            M[k, j] = mm * s
    return J, M, komplex


def Y_matrix(name, N, bits=None, mu=None):
    """Y_N = (I + K_N)^{-1} = J_N^{-1} M_N, umskaliert, als numpy-Matrix (float/complex)."""
    J, M, komplex = matrizen(name, N, bits, mu)
    Y = J.solve(M, algorithm="approx")
    n = N + 1
    if komplex:
        A = np.array([[complex(Y[i, j].real.mid()) + 1j * float(Y[i, j].imag.mid()) for j in range(n)] for i in range(n)])
    else:
        A = np.array([[float(Y[i, j].mid()) for j in range(n)] for i in range(n)])
    return A


def to_np(Mx, komplex):
    n, m = Mx.nrows(), Mx.ncols()
    if komplex:
        return np.array([[float(Mx[i, j].real.mid()) + 1j * float(Mx[i, j].imag.mid()) for j in range(m)] for i in range(n)])
    return np.array([[float(Mx[i, j].mid()) for j in range(m)] for i in range(n)])


# ---------------------------------------------------------------- numpy-Auswertung des Modells (für Symbole)

def f_numpy(name, mu=None):
    ctx.prec = 200
    if name in ("kubisch", "radial"):
        c = [float(x) for x in (kubik_koeff(mu) if name == "kubisch" else radial_koeff(mu))]
        p = np.polynomial.Polynomial(c)
        return p, p.deriv()
    if name == "moeb_ell":
        a, b, c, d, _ = moeb_ell_par()
        a, b, c, d = (float(t) for t in (a, b, c, d))
    elif name == "moeb_sym":
        tau = np.tan(float(mu) / 2)
        a, b, c, d = 1.0, tau, -tau, 1.0
    elif name == "moeb_lox":
        a, b, c, d = (complex(float(t.real), float(t.imag)) for t in moeb_lox_par(mu))
    elif name == "affin":
        mu_ = np.sqrt(Lc) if mu is None else complex(mu)
        return (lambda x: Lc + mu_ * (x - Lc)), (lambda x: mu_ + 0 * x)
    det = a * d - b * c
    return (lambda x: (a * x + b) / (c * x + d)), (lambda x: det / (c * x + d) ** 2)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    for nm in ("kubisch", "moeb_ell", "moeb_lox", "affin"):
        f, fp = f_numpy(nm)
        print(nm, "f(L) − L =", abs(f(Lc) - Lc), " f'(L) =", fp(Lc), " f(0) =", f(0))


def Ystrich_matrix(name, N, bits=None, mu=None):
    """Dualbild bei a₀ = f(0): Y'_N = S (I − Y_N) S^{-1} = S J_N^{-1} U_N in H²(D(a₀, R₀)), R₀ = |L − a₀|.
    S = Taylor-Verschiebung 0 → a₀, U_N[k][j] = [x^k] (f − a₀)^j  (D3: J = U S + M). Rückgabe numpy (umskaliert mit R₀^k)."""
    if bits is None:
        bits = int(1.2 * N) + 250
    ctx.prec = bits
    fk, mk, komplex = modell(name, N, mu)
    P = acb_poly if komplex else arb_poly
    Mat = acb_mat if komplex else arb_mat
    num = acb if komplex else arb
    a0 = fk[0]
    f = P(fk)
    u = f - P([a0])
    n = N + 1
    pw = [P([1])]
    upw = [P([1])]
    for _ in range(N):
        pw.append((pw[-1] * f).truncate(n))
        upw.append((upw[-1] * u).truncate(n))
    J = Mat(n, n)
    U = Mat(n, n)
    S = Mat(n, n)
    from math import comb
    for j in range(n):
        for k in range(n):
            ck = pw[j][k] if k <= pw[j].degree() else num(0)
            mm = mk[k - j] if k >= j else num(0)
            J[k, j] = ck + mm
            U[k, j] = upw[j][k] if k <= upw[j].degree() else num(0)
            S[j, k] = num(comb(k, j)) * a0 ** (k - j) if k >= j else num(0)
    X = J.solve(U, algorithm="approx")
    Yp = S * X
    L = Lacb()
    R0 = abs(L - a0) if komplex else ((arb(L_RE) - a0) ** 2 + arb(L_IM) ** 2).sqrt()
    rp = [R0 ** k for k in range(n)]
    for i in range(n):
        for j in range(n):
            Yp[i, j] = Yp[i, j] * rp[i] / rp[j]
    return to_np(Yp, komplex)
