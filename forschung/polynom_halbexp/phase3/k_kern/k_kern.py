"""Strang K — Kernmodul (Phase 3).

Modelle und Abschnittsmatrizen, übernommen (Kopie, nicht importiert) aus o_operator/o_kern.py (Strang O),
damit in fremden Ordnern kein __pycache__ entsteht. Bezeichnungen:
  J_N[k][j] = [x^k] f^j + [x^{k-j}] m,  m = f'∘f;   umskaliert  J̃ = D J D^{-1},  D = diag(ρ^k),  ρ = |L|.
  ẽ_t = ((ρ/t)^k)_{k≤N}  (umskalierter abgeschnittener Cauchy-Kern bei t),  G_N(t) = ‖J̃^{-1} ẽ_t‖₂  (Majorante O §33.1).
"""
import json
import os

from flint import arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, ctx

L_RE = "0.31813150520476413531265425158766451720351761387139986692237"
L_IM = "1.33723570143068940890116214319371061253950213846051241687"
Lc = complex(0.31813150520476413, 1.3372357014306895)
RHO = abs(Lc)
HIER = os.path.dirname(os.path.abspath(__file__))
GEMEINSAM = os.path.join(HIER, "..", "..", "gemeinsam")


def Lacb():
    return acb(arb(L_RE), arb(L_IM))


def kubik_koeff(mu=None):
    L = Lacb()
    if mu is None:
        mu = L.sqrt()
    rL = (mu - 1) / (L - L.conjugate())
    beta = rL.imag / L.imag
    alpha = rL.real - beta * L.real
    R2 = L.real ** 2 + L.imag ** 2
    return [R2 * alpha, 1 + R2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]


def radial_koeff(alpha=None):
    A, B = arb(L_RE), arb(L_IM)
    r2 = A * A + B * B
    if alpha is None:
        alpha = arb("0.31382") / r2.sqrt()
    alpha = arb(alpha)
    return [alpha * r2, 1 - 2 * alpha * A, alpha]


def jet_koeff(eps):
    """Gs Familie jet_ε (in x): f_ε = f_kubisch + (ε/ρ³)(x² − 2Re(L)x + ρ²)²  — L, L̄ und μ bleiben exakt fest."""
    c = kubik_koeff()
    A, B = arb(L_RE), arb(L_IM)
    r2 = A * A + B * B
    rho = r2.sqrt()
    q = arb_poly([r2, -2 * A, 1])
    p = arb_poly(c) + (q * q) * (arb(eps) / rho ** 3)
    return [p[k] for k in range(p.degree() + 1)]


def jet4_koeff(eps):
    """Störung, die in L, L̄ von 4. Ordnung verschwindet: f_ε = f_kubisch + (ε/ρ⁷)(x² − 2Re(L)x + ρ²)⁴ (L, μ, f'', f''' in L fest)."""
    c = kubik_koeff()
    A, B = arb(L_RE), arb(L_IM)
    r2 = A * A + B * B
    rho = r2.sqrt()
    q = arb_poly([r2, -2 * A, 1])
    p = arb_poly(c) + (q * q * q * q) * (arb(eps) / rho ** 7)
    return [p[k] for k in range(p.degree() + 1)]


def modell(name, N, mu=None):
    """(fk, mk): Taylor-Koeffizienten von f und m = f'∘f bis Grad N (arb)."""
    if name in ("kubisch", "radial") or name.startswith("jet"):
        if name.startswith("jet4_"):
            c = jet4_koeff(name[5:])
        elif name.startswith("jet"):
            c = jet_koeff(name[3:])
        else:
            c = kubik_koeff(mu) if name == "kubisch" else radial_koeff(mu)
        f = arb_poly(c)
        fp = f.derivative()
        m = arb_poly([0])
        pw = arb_poly([1])
        for k in range(fp.degree() + 1):
            m += fp[k] * pw
            pw = pw * f
        fk = [f[k] if k <= f.degree() else arb(0) for k in range(N + 1)]
        mk = [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
        return fk, mk
    if name in ("kneser", "kneser_z"):
        fn = "kneser_halbexp_koeff_zert_640.json" if name == "kneser_z" else "kneser_halbexp_koeff_1000.json"
        d = json.load(open(os.path.join(GEMEINSAM, fn), encoding="utf-8"))
        roh = d["a"][:N + 2]
        a = [arb(x[0] if isinstance(x, list) else x) for x in roh]
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
    raise ValueError(name)


def J_tilde(name, N, bits=None, mu=None):
    """Umskalierte J̃_N als arb_mat (Präzision bits)."""
    if bits is None:
        bits = int(1.0 * N) + 250
    ctx.prec = bits
    fk, mk = modell(name, N, mu)
    f = arb_poly(fk)
    pw = [arb_poly([1])]
    for _ in range(N):
        pw.append((pw[-1] * f).truncate(N + 1))
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    rp = [rho ** k for k in range(N + 1)]
    J = arb_mat(N + 1, N + 1)
    for j in range(N + 1):
        col = pw[j]
        for k in range(N + 1):
            ck = col[k] if k <= col.degree() else arb(0)
            mm = mk[k - j] if k >= j else arb(0)
            J[k, j] = (ck + mm) * rp[k] / rp[j]
    return J


def to_acb_mat(A):
    n, m = A.nrows(), A.ncols()
    B = acb_mat(n, m)
    for i in range(n):
        for j in range(m):
            B[i, j] = acb(A[i, j])
    return B


def cauchy_rhs(N, ts):
    """Matrix mit Spalten ẽ_t = ((ρ/t)^k)_{k≤N} für komplexe t (Python complex oder acb)."""
    rho = (arb(L_RE) ** 2 + arb(L_IM) ** 2).sqrt()
    B = acb_mat(N + 1, len(ts))
    for c, t in enumerate(ts):
        tt = t if isinstance(t, acb) else acb(t.real, t.imag)
        q = acb(rho) / tt
        p = acb(1)
        for k in range(N + 1):
            B[k, c] = p
            p = p * q
    return B


def spaltennormen(X):
    out = []
    for c in range(X.ncols()):
        s = arb(0)
        for k in range(X.nrows()):
            s += abs(X[k, c]) ** 2
        out.append(float(s.sqrt().mid()))
    return out
