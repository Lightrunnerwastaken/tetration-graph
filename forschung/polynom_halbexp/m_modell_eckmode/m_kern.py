"""Strang M — Kernfunktionen für Modelle f (Polynom, reell oder komplex): J_N(f), Eckmoden, Pol-Antworten.

J_N δ = T_N[δ∘f + m·δ],  m = f'∘f,  δ ∈ P_N.   Matrix: J[k][j] = [x^k] f^j + [x^{k−j}] m.
Gewichte: ℓ²_ρ (ρ = Radius der Kreisscheibe D, Standard |L|), d. h. δ̃_k = δ_k ρ^k.
"""
import math
import sys

import numpy as np
from flint import acb, acb_mat, acb_poly, arb, arb_mat, arb_poly, ctx

L_KOMPLEX = complex(0.31813150520476413, 1.3372357014306895)
RHO = abs(L_KOMPLEX)


def ist_reell(coef):
    return all(abs(complex(c).imag) == 0 for c in coef)


def _zahl(c, reell):
    if reell:
        return arb(repr(float(complex(c).real))) if not isinstance(c, (arb, str)) else arb(c)
    c = complex(c) if not isinstance(c, acb) else c
    return acb(c.real, c.imag) if not isinstance(c, acb) else c


def potenzen(coef, M, reell):
    """[f^0, ..., f^M] modulo x^{M+1}."""
    P_ = arb_poly if reell else acb_poly
    f = P_([_zahl(c, reell) for c in coef])
    out = [P_([1])]
    for _ in range(M):
        out.append((out[-1] * f).truncate(M + 1))
    return out, f


def jmat(coef, N, reell=None):
    """J_N(f) als arb_mat/acb_mat ((N+1)x(N+1)) für ein Polynom f mit Koeffizienten coef (exakt, f polynomial)."""
    if reell is None:
        reell = ist_reell(coef)
    P, f = potenzen(coef, N, reell)
    fp = f.derivative()
    # m = f'∘f mod x^{N+1}
    m = (arb_poly if reell else acb_poly)([0])
    fpc = [fp[i] for i in range(fp.degree() + 1)]
    for i, c in enumerate(fpc):
        m += c * P[i]
    M_ = arb_mat if reell else acb_mat
    J = M_(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            v = P[j][k]
            if j <= k:
                v = v + m[k - j]
            J[k, j] = v
    return J


def skaliere(Xnp, rho, zeilen, spalten):
    return Xnp * (rho ** np.asarray(zeilen, dtype=float))[:, None] / (rho ** np.asarray(spalten, dtype=float))[None, :]


def mat_np(X, reell):
    n, m = X.nrows(), X.ncols()
    if reell:
        return np.array([[float(X[i, j].mid()) for j in range(m)] for i in range(n)])
    return np.array([[complex(X[i, j].mid()) for j in range(m)] for i in range(n)])


def inverse_spalten(J, N, B=None, reell=True):
    """Letzte B Spalten (oder alle) von J⁻¹ als numpy, plus Spaltenindizes."""
    n = N + 1
    if B is None or B >= n:
        B = n
    M_ = arb_mat if reell else acb_mat
    E = M_(n, B)
    for c in range(B):
        E[n - B + c, c] = 1
    X = J.solve(E, algorithm="approx")
    return mat_np(X, reell), np.arange(n - B, n)


def ecksingulaer(coef, N, rho=RHO, B=None, prec=None, k=4):
    """Grösste Singulärwerte der ρ-skalierten Inversen (ℓ²_ρ) + Singulärvektoren."""
    reell = ist_reell(coef)
    ctx.prec = prec or int(0.9 * N) + 300
    J = jmat(coef, N, reell)
    X, cols = inverse_spalten(J, N, B, reell)
    G = skaliere(X, rho, np.arange(N + 1), cols)
    U, S, Vt = np.linalg.svd(G, full_matrices=False)
    return S[:k], U[:, :k], Vt[:k, :], cols


def auf_kreis(vec, M=8192):
    """Werte des Polynoms Σ vec_k e^{ikφ} auf dem Kreis (φ_j = 2πj/M)."""
    n = len(vec)
    return np.fft.ifft(np.concatenate([vec, np.zeros(M - n)])) * M


def quadratisch(alpha, L=L_KOMPLEX):
    """Radiales reelles Quadrat-Modell f(x) = x + α(x−L)(x−L̄): Koeffizienten [c0, c1, c2]."""
    a, r2 = L.real, abs(L) ** 2
    return [alpha * r2, 1 - 2 * alpha * a, alpha]


def kubisch(mu, L=L_KOMPLEX):
    """Ds kubisches Modell f(x) = x + (x² − 2Re(L)x + |L|²)(α + βx) mit f'(L) = μ."""
    rL = (mu - 1) / (L - L.conjugate())
    beta = rL.imag / L.imag
    alpha = rL.real - beta * L.real
    Lr = abs(L)
    return [Lr ** 2 * alpha, 1 + Lr ** 2 * beta - 2 * L.real * alpha, alpha - 2 * L.real * beta, beta]


def polywert(coef, x):
    return sum(c * x ** i for i, c in enumerate(coef))


def polyabl(coef, x):
    return sum(i * c * x ** (i - 1) for i, c in enumerate(coef) if i)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    c = kubisch(np.sqrt(L_KOMPLEX))
    print("kubisch", c, "f(L)-L", abs(polywert(c, L_KOMPLEX) - L_KOMPLEX), "f'(L)", polyabl(c, L_KOMPLEX))
    q = quadratisch(0.2283)
    print("quadratisch", q, "f(L)-L", abs(polywert(q, L_KOMPLEX) - L_KOMPLEX), "f'(L)", polyabl(q, L_KOMPLEX),
          abs(polyabl(q, L_KOMPLEX)))
