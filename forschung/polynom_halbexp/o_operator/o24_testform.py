"""O24 — Test der verallgemeinerten Akkretivität (Runde 6, Weg zu einer expliziten Schranke):
  Re⟨Tδ, Bδ⟩_p ≥ λ‖Bδ‖²_p  für alle Polynome δ ?
T δ = δ + (δ∘f)/m auf ∂D (= J_∞δ/m), B δ := δ auf Γin, (δ∘f)/m auf Γout (Randwerte von δ auf ∂Ω, zurückgezogen);
Gewicht p auf ∂D ('eins' oder 'm2' = |m|², dann ‖Tδ‖_p = ‖J_∞δ‖). Cross-Term X_p(δ) = Re∫ p (δ∘f) δ̄ / m über ganz ∂D.
Ausgabe je N: λ = min Re⟨Tδ,Bδ⟩_p/‖Bδ‖²_p,  c_T = min ‖Tδ‖_p/‖δ‖_{H²},  c_B = min ‖Bδ‖_p/‖δ‖_{H²},  Vorhersage 1/√(1/λ²+(1+1/λ)²)-Schranke.
Gram-Matrizen in arb (flint) aus Doppelgenauigkeits-Werten von f, m an den Quadraturpunkten (konsistent als „exakte“ Daten),
Eigenwerte in mpmath. Reelles Modell ⇒ Integration über φ ∈ [0, π] und 2·Re.
Aufruf: python o24_testform.py modell N1,N2,... gewicht [nq]
"""
import json
import os
import sys

import mpmath as mp
import numpy as np
from flint import acb, acb_mat, arb, ctx

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
Ns = [int(s) for s in sys.argv[2].split(",")]
gew = sys.argv[3] if len(sys.argv) > 3 else "eins"
nq = int(sys.argv[4]) if len(sys.argv) > 4 else 600
R = K.RHO
thL = np.angle(K.Lc)


def f_und_m(name):
    if name == "kneser":
        d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gemeinsam",
                                        "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
        a = np.array([float(x) for x in d["a"]])
        k = np.arange(len(a))

        def f(x):
            return np.polynomial.polynomial.polyval(x, a)

        def fp(x):
            return np.polynomial.polynomial.polyval(x, a[1:] * k[1:])

        return f, (lambda x: np.exp(x) / fp(x))
    f, fp = K.f_numpy(name)
    return f, (lambda x: fp(f(x)))


def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    s = (x + 1) / 2
    t = (1 - np.cos(np.pi * s)) / 2
    dt = np.pi / 2 * np.sin(np.pi * s)
    return a + (b - a) * t, (b - a) * dt * w / 2


f, mfun = f_und_m(name)
Nmax = max(Ns)
ctx.prec = int(8 * Nmax) + 200
mp.mp.dps = int(2.6 * Nmax) + 60
teile = []
for teil, (a_, b_) in (("out", (0.0, thL)), ("in", (thL, np.pi))):
    phi, w = gl(a_, b_, nq)
    x = R * np.exp(1j * phi)
    fx = f(x)
    m = mfun(x)
    p = np.ones_like(phi) if gew == "eins" else np.abs(m) ** 2
    if teil == "out":
        assert np.all(np.abs(fx) > R * (1 - 1e-12)), "Γout: |f| > ρ verletzt"
    else:
        assert np.all(np.abs(fx) < R * (1 + 1e-12)), "Γin: |f| < ρ verletzt"
    # Matrizen (Zeilen = Punkte, Spalten = Grad k): D = (x/ρ)^k, F = (f(x)/ρ)^k / m, gewichtet mit sqrt(p w / π)
    nP = len(phi)
    D = acb_mat(nP, Nmax + 1)
    D0 = acb_mat(nP, Nmax + 1)
    F = acb_mat(nP, Nmax + 1)
    for j in range(nP):
        sw = arb(float(np.sqrt(p[j] * w[j] / np.pi)))
        s0 = arb(float(np.sqrt(w[j] / np.pi)))
        xr = acb(complex(x[j] / R))
        fr = acb(complex(fx[j] / R))
        mi = acb(complex(1 / m[j]))
        px, pf, p0 = sw, sw * mi, s0
        for kk in range(Nmax + 1):
            D[j, kk] = px
            D0[j, kk] = p0
            F[j, kk] = pf
            px = px * xr
            p0 = p0 * xr
            pf = pf * fr
    teile.append((teil, D, F, D0))


def gram(A, B):
    """Re(A^H B) (reell, wegen Symmetrie φ ↔ −φ: Integral über [0, π] mal 2 ergibt 2·Re, Normierung dφ/2π ⇒ Faktor 1/π schon in sw)."""
    C = A.conjugate().transpose() * B
    n = C.nrows()
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = mp.mpf(C[i, j].real.mid().str(mp.mp.dps, radius=False))
    return M


H_all, GB_all, GD_all, GT_all = None, None, None, None
for teil, D, F, D0 in teile:
    T_ = D + F
    B_ = D if teil == "in" else F
    blk = {"H": gram(T_, B_), "GB": gram(B_, B_), "GD": gram(D0, D0), "GT": gram(T_, T_)}
    if H_all is None:
        H_all, GB_all, GD_all, GT_all = blk["H"], blk["GB"], blk["GD"], blk["GT"]
    else:
        H_all += blk["H"]; GB_all += blk["GB"]; GD_all += blk["GD"]; GT_all += blk["GT"]


def gen_min(A, Bm, n, beide=False):
    A = A[:n, :n]
    A = (A + A.T) / 2
    Bm = Bm[:n, :n]
    Lc = mp.cholesky(Bm)
    Li = mp.inverse(Lc)
    M = Li * A * Li.T
    M = (M + M.T) / 2
    ev = mp.eigsy(M, eigvals_only=True)
    return (min(ev), max(ev)) if beide else min(ev)


for N in Ns:
    n = N + 1
    lam = gen_min(H_all, GB_all, n)
    cT2 = gen_min(GT_all, GD_all, n)
    cB2 = gen_min(GB_all, GD_all, n)
    xmin, xmax = gen_min(H_all - GB_all, GB_all, n, beide=True)
    lamf = float(lam)
    schr = 1 / np.sqrt(1 / lamf ** 2 + (1 + 1 / lamf) ** 2) if lamf > 0 else float("nan")
    print(f"{name} N={N} p={gew} nq={nq}: λ = {mp.nstr(lam, 8)};  c_T = {mp.nstr(mp.sqrt(cT2), 8)};  "
          f"c_B = {mp.nstr(mp.sqrt(cB2), 8)};  Schranke aus λ: {schr:.4f};  X/‖B‖² ∈ [{mp.nstr(xmin, 6)}, {mp.nstr(xmax, 6)}]", flush=True)
