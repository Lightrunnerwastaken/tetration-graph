"""O29 — Konturverschiebung auf die „Mittellinien“: X als Integral von W(ξ)·u(z+τ)·conj(u(z)) mit kleiner reeller Verschiebung τ.

Unterer Term (Γin, Im ζ = 0):  X_b = Re∫ A(ξ) u(g(ξ)) conj(u(ξ)) dξ,  A = p·R  (p = |m|², R = (g'/f')^{1/2}/m).
  Analytische Fortsetzung: conj(u(ξ)) = ũ(ξ), ũ(ζ) := conj(u(ζ̄));  p(ζ) := m(x(ζ))·conj(m(x(ζ̄))).
  Kontur ζ = c(ξ) = ξ − iε(ξ) mit Im g(c) = ε (symmetrischer Punkt): dann g(c) = z + τ, z = ξ + iε, τ = Re g(c) − ξ reell,
  X_b = Re∫ W_b(ξ) u(z+τ) conj(u(z)) dξ,  W_b = A(c(ξ))·c'(ξ).
Oberer Term (Γout, Im ζ = π):  X_t = Re∫ p R̄ u(ξ+iπ) conj(u(g(ξ+iπ))) dξ; Kontur s = ξ + iε (nach oben in die Linse),
  Im g(ξ − iε + iπ) = π + ε  ⇒  g(ξ−iε+iπ) = z + τ, z = ξ + i(π+ε);  X_t = Re∫ W_t(ξ) u(z) conj(u(z+τ)) dξ.
Schranke: X ≥ ∫Re W |u|² − ∫|W||Δu||u| ≥ −¼∫ (|W|²τ²/Re W)·|u'|²  (Δu = u(z+τ) − u(z)), falls Re W > 0.
Ausgabe: ε, τ, W, Re W, ρ := |W|²τ²/(4 Re W) längs beider Kurven.
f ausserhalb D: nahe Γin über f(x) = log f(eˣ) (Zweig stetig), in der Linse über f(x) = e^{f⁻¹(x)} (Newton).
"""
import json
import os
import sys

import numpy as np

import o26_streifen as G

sys.stdout.reconfigure(encoding="utf-8")
L, Lb, th0 = G.L, G.Lb, G.th0
d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gemeinsam",
                                "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
a = np.array([float(x) for x in d["a"]])
kk = np.arange(len(a))
RHO = abs(L)


def f_in(x):
    return np.polynomial.polynomial.polyval(x, a)


def fp_in(x):
    return np.polynomial.polynomial.polyval(x, a[1:] * kk[1:])


def finv(y, x0):
    x = x0
    for _ in range(60):
        dx = (f_in(x) - y) / fp_in(x)
        x = x - dx
        if abs(dx) < 1e-15:
            break
    return x


def f_und_fp(x, wo, ref=None):
    """wo: 'D' (|x| ≤ ρ), 'aussen_in' (jenseits Γin), 'linse' (in f(D)∖D). ref: Referenzwert von f für den Zweig."""
    if wo == "D":
        return f_in(x), fp_in(x)
    if wo == "aussen_in":
        ex = np.exp(x)
        fe = f_in(ex)
        fx = np.log(fe)
        if ref is not None:
            fx = fx + 2j * np.pi * np.round((ref - fx).imag / (2 * np.pi))
        return fx, fp_in(ex) * ex / fe
    if wo == "linse":
        y = finv(x, ref)            # ref = Startwert für f⁻¹(x) (in D)
        fx = np.exp(y)
        return fx, fx / fp_in(y), y
    raise ValueError(wo)


M = lambda x: (x - L) / (x - Lb)
Mp = lambda x: (L - Lb) / (x - Lb) ** 2


def xz(z):
    w = np.exp(z + 1j * th0)
    return (L - Lb * w) / (1 - w)


def g_und_R(zz, wo, gref, Rref, fref=None):
    """g(ζ), g'(ζ), R(ζ), m(x(ζ)) mit stetigen Zweigen (gref, Rref Referenzwerte)."""
    x = xz(zz)
    if wo == "linse":
        fx, fpx, _ = f_und_fp(x, "linse", fref)
    else:
        fx, fpx = f_und_fp(x, wo, fref)
    gz = np.log(M(fx)) - 1j * th0
    gz = gz + 2j * np.pi * np.round((gref - gz).imag / (2 * np.pi))
    gp = (Mp(fx) * fpx / M(fx)) * (M(x) / Mp(x))
    m = np.exp(x) / fpx
    Rw = np.sqrt(gp / fpx) / m
    if abs(Rw - Rref) > abs(-Rw - Rref):
        Rw = -Rw
    return gz, gp, Rw, m, fx


def haupt():
    xis = np.linspace(-12, 12, 481)
    mu = G.mu
    Rend = mu ** -1.5
    print("UNTEN (Γin): ξ, ε, τ, W_b, Re W_b, |W|²/ReW, ρ = |W|²τ²/(4ReW)")
    res_b = []
    eps_prev = G.ell.imag / 2
    for xi in xis:
        # Referenzen auf der Linie selbst
        x0 = xz(xi + 0j)
        f0 = f_in(x0)
        g0 = np.log(M(f0)) - 1j * th0
        lb = g0 - xi
        g0 = xi + lb.real + 1j * ((lb.imag + np.pi) % (2 * np.pi) - np.pi)
        R0 = G.daten(xi, False)[2]
        if R0.real < 0:
            R0 = -R0
        # Newton in ε für Im g(ξ − iε) = ε
        eps = eps_prev
        for _ in range(50):
            gz, gp, Rw, m, fx = g_und_R(xi - 1j * eps, "aussen_in", g0, R0, f0)
            F = gz.imag - eps
            dF = -gp.real - 1           # d/dε Im g(ξ − iε) = Im(−i g') = −Re g'
            de = -F / dF
            eps += de
            if abs(de) < 1e-14:
                break
        eps_prev = eps
        gz, gp, Rw, m, fx = g_und_R(xi - 1j * eps, "aussen_in", g0, R0, f0)
        tau = gz.real - xi
        xc = xz(xi + 1j * eps)                         # konj. Punkt in D
        p_c = m * np.conj(np.exp(xc) / fp_in(xc))      # p(ζ) = m(x(ζ))·conj(m(x(ζ̄)))
        res_b.append((xi, eps, tau, p_c * Rw))
    # c'(ξ) = 1 − iε'(ξ)
    E = np.array([r[1] for r in res_b])
    dE = np.gradient(E, xis)
    Wb = np.array([r[3] for r in res_b]) * (1 - 1j * dE)
    TAUb = np.array([r[2] for r in res_b])
    rho_b = np.abs(Wb) ** 2 * TAUb ** 2 / (4 * Wb.real)
    for i in range(0, len(xis), 20):
        print(f"{xis[i]:+6.2f}  ε={E[i]:.4f}  τ={TAUb[i]:+.4f}  W={Wb[i].real:+.4f}{Wb[i].imag:+.4f}i  ReW={Wb[i].real:.4f}  "
              f"|W|²/ReW={abs(Wb[i]) ** 2 / Wb[i].real:.4f}  ρ={rho_b[i]:.5f}")
    print(f"unten: min Re W = {Wb.real.min():.4f}, max ρ = {rho_b.max():.5f}, ε ∈ [{E.min():.4f}, {E.max():.4f}]")
    np.savez("o29_unten.npz", xi=xis, eps=E, tau=TAUb, W=Wb)

    # ------------------------------------------------------------------ OBEN (Γout, Im ζ = π)
    print("OBEN (Γout): ξ, ε, τ, W_t, Re W_t, |W|²/ReW, ρ")
    res_t = []
    eps_prev = G.ell.imag / 2
    for xi in xis:
        x0 = xz(xi + 1j * np.pi)
        f0 = f_in(x0)
        g0 = np.log(M(f0)) - 1j * th0
        lt = g0 - (xi + 1j * np.pi)
        g0 = xi + 1j * np.pi + lt.real + 1j * ((lt.imag + np.pi) % (2 * np.pi) - np.pi)
        R0 = G.daten(xi, True)[2]
        if R0.real < 0:
            R0 = -R0
        eps = eps_prev
        for _ in range(50):
            zl = xi + 1j * (np.pi - eps)
            gz, gp, Rw, m, fx = g_und_R(zl, "D", g0, R0)
            F = gz.imag - (np.pi + eps)
            dF = -gp.real - 1            # d/dε Im g(ξ + i(π − ε)) = −Re g'
            de = -F / dF
            eps += de
            if abs(de) < 1e-14:
                break
        eps_prev = eps
        zl = xi + 1j * (np.pi - eps)
        gz, gp, Rw, m_low, fx = g_und_R(zl, "D", g0, R0)
        tau = gz.real - xi
        zu = xi + 1j * (np.pi + eps)            # Punkt in der Linse
        xu = xz(zu)
        y0 = xz(zl - tau)                         # Startwert: f(x(ζ_low − τ)) ≈ x(z)
        y = finv(xu, y0)
        assert abs(f_in(y) - xu) < 1e-12 and abs(y) < RHO * (1 + 1e-9), (xi, y, xu)
        fpu = np.exp(y) / fp_in(y)                # f'(x) = e^y / f'(y), y = f⁻¹(x)
        m_up = np.exp(xu) / fpu
        P = m_up * np.conj(m_low)                 # p(s) = m(x(s+iπ))·conj(m(x(s̄+iπ)))
        Rt_conj = np.conj(Rw)                     # R̃(s) = conj(R(s̄)) mit s̄ + iπ = zl
        res_t.append((xi, eps, tau, P * Rt_conj, abs(fx_check := 0)))
    Et = np.array([r[1] for r in res_t])
    dEt = np.gradient(Et, xis)
    Wt = np.array([r[3] for r in res_t]) * (1 + 1j * dEt)
    TAUt = np.array([r[2] for r in res_t])
    rho_t = np.abs(Wt) ** 2 * TAUt ** 2 / (4 * Wt.real)
    for i in range(0, len(xis), 20):
        print(f"{xis[i]:+6.2f}  ε={Et[i]:.4f}  τ={TAUt[i]:+.4f}  W={Wt[i].real:+.4f}{Wt[i].imag:+.4f}i  ReW={Wt[i].real:.4f}  "
              f"|W|²/ReW={abs(Wt[i]) ** 2 / Wt[i].real:.4f}  ρ={rho_t[i]:.5f}")
    print(f"oben: min Re W = {Wt.real.min():.4f}, max ρ = {rho_t.max():.5f}, ε ∈ [{Et.min():.4f}, {Et.max():.4f}]")
    np.savez("o29_oben.npz", xi=xis, eps=Et, tau=TAUt, W=Wt)


if __name__ == "__main__":
    haupt()
