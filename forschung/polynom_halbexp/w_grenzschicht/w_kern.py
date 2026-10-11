"""Strang W — Kernfunktionen: Modelle (kubisch, Kneser über Π_N f_K), J_N in arb, Inverse/Singulärvektoren (|L|-umskaliert).

Konventionen (wie Strang D): J_N[k][j] = [x^k] f^j + [x^{k-j}] (f'∘f)  (Linearisierung δ ↦ T_N[δ∘f + (f'∘f)δ]),
umskaliert J̃ = D J D⁻¹, D = diag(ρ^k), ρ = |L|.  SVD J̃⁻¹ = Σ σ_i u_i v_iᵀ: u_i Ausgabe (Koeffizienten δ_k ρ^k),
v_i Eingabe (Gleichungen). Adjungierte Eckmode als Laurent-Polynom: V_i(x) = Σ_k v_{i,k} (ρ/x)^k.
"""
import json
import os
import sys

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

HIER = os.path.dirname(os.path.abspath(__file__))
GEMEINSAM = os.path.join(HIER, "..", "gemeinsam")
L = complex(0.31813150520476413, 1.3372357014306895)
RHO = abs(L)


def kubik_koeff(mu=None, phi_deg=None, Lfix=None):
    """Reelles kubisches Modell f(x) = x + (x² − 2Re(Λ)x + |Λ|²)(α + βx), Fixpunkte Λ, Λ̄, f'(Λ) = μ.
    Standard: Λ = L, μ = √L.  Rückgabe: Koeffizientenliste [c0, c1, c2, c3] (float) und Λ."""
    Lf = L if Lfix is None else Lfix
    if phi_deg is not None:
        Lf = abs(Lf) * np.exp(1j * np.deg2rad(phi_deg))
    if mu is None:
        mu = np.sqrt(L)
    rL = (mu - 1) / (Lf - Lf.conjugate())
    beta = rL.imag / Lf.imag
    alpha = rL.real - beta * Lf.real
    R = abs(Lf)
    coef = [R ** 2 * alpha, 1 + R ** 2 * beta - 2 * Lf.real * alpha, alpha - 2 * Lf.real * beta, beta]
    return coef, Lf


def potenzen(a, M):
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(M):
        P.append((P[-1] * f).truncate(M + 1))
    return P


def J_mat(a, N, P=None):
    """J_N (arb_mat) für das Polynom/Keim mit Koeffizienten a (Liste arb, Länge ≥ N+1 wird auf N+1 gekürzt)."""
    a = list(a[:N + 1]) + [arb(0)] * max(0, N + 1 - len(a))
    if P is None:
        P = potenzen(a, N)
    fpf = arb_poly([0])
    for m in range(1, N + 1):
        fpf += m * a[m] * P[m - 1]
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        Pk = None
        for j in range(N + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return J


def J_mat_m(a, mser, N, P=None):
    """J_N mit vorgegebenem Multiplikator m (Liste arb, Taylor-Koeffizienten) statt f'∘f (für A_f bei Kneser)."""
    a = list(a[:N + 1]) + [arb(0)] * max(0, N + 1 - len(a))
    if P is None:
        P = potenzen(a, N)
    J = arb_mat(N + 1, N + 1)
    for k in range(N + 1):
        for j in range(N + 1):
            J[k, j] = P[j][k] + (mser[k - j] if j <= k else 0)
    return J


def kubik_a(N, coef):
    return [arb(float(c)) for c in coef] + [arb(0)] * (N + 1 - len(coef))


def inverse_umskaliert(J, N, rho=RHO):
    """Volle Inverse, |L|-umskaliert, als float64-Matrix (Einträge sind moderat)."""
    n = N + 1
    Ji = J.solve(arb_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)]), algorithm="approx")
    k = np.arange(n)
    sc = rho ** k
    G = np.array([[float(Ji[i, j].mid()) for j in range(n)] for i in range(n)])
    return G * sc[:, None] / sc[None, :]


def letzte_spalten(J, N, B, extra=None, rho=RHO):
    """J⁻¹ auf die letzten B Einheitsvektoren (+ optionale Zusatzspalten, Liste von float-Vektoren in ROHEN Koeffizienten);
    Rückgabe: G (umskaliert, n×B), Zusatzantworten (roh)."""
    n = N + 1
    extra = extra or []
    E = arb_mat(n, B + len(extra))
    for c in range(B):
        E[N - B + 1 + c, c] = 1
    for i, v in enumerate(extra):
        for kk in range(n):
            E[kk, B + i] = v[kk]
    X = J.solve(E, algorithm="approx")
    k = np.arange(n)
    cols = np.arange(N - B + 1, n)
    G = np.array([[float(X[i, c].mid()) for c in range(B)] for i in range(n)]) * (rho ** k)[:, None] / (rho ** cols)[None, :]
    Z = [np.array([float(X[i, B + e].mid()) for i in range(n)]) for e in range(len(extra))]
    return G, Z


def V_auswerten(v, x, rho=RHO):
    """V(x) = Σ_k v_k (ρ/x)^k (Horner), x komplex (Skalar oder Array)."""
    w = rho / np.asarray(x, dtype=complex)
    s = np.zeros_like(w)
    for c in v[::-1]:
        s = s * w + c
    return s


def U_auswerten(u, x, rho=RHO):
    """U(x) = Σ_k u_k (x/ρ)^k."""
    w = np.asarray(x, dtype=complex) / rho
    s = np.zeros_like(w)
    for c in u[::-1]:
        s = s * w + c
    return s
