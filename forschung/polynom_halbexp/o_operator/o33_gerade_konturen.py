"""O33 — Satz O11 mit WAAGRECHTEN Konturen (einfacher für die Strenge): unten c(ξ) = ξ − iε_b, oben s = ξ + iε_t.
τ wird komplex: unten τ = g(ξ − iε_b) − (ξ + iε_b), oben τ = g(ξ + i(π − ε_t)) − (ξ + i(π + ε_t)).
Schranke: X ≥ −Σ∫ ρ |u'(z + rτ)|² (r ∈ [0,1]) mit ρ = |W|²|τ|²/(4 Re W); Schur-Test wie o32.
Ausgabe je (ε_b, ε_t): min Re W, max ρ, Schur-Λ (r = 0, ½, 1) und c(1 − Λ).
"""
import sys

import numpy as np

import o26_streifen as G
import o29_mittellinie as K29   # Funktionen f_in, fp_in, finv, g_und_R, xz, M ... (rechnet beim Import die Kurven; Ausgabe unterdrückt)

sys.stdout.reconfigure(encoding="utf-8")
L, th0 = G.L, G.th0
xi = np.linspace(-14, 14, 561)
hx = xi[1] - xi[0]


def unten(eps):
    W = np.zeros(len(xi), complex)
    tau = np.zeros(len(xi), complex)
    for i, x_ in enumerate(xi):
        x0 = K29.xz(x_ + 0j)
        f0 = K29.f_in(x0)
        g0 = np.log(K29.M(f0)) - 1j * th0
        lb = g0 - x_
        g0 = x_ + lb.real + 1j * ((lb.imag + np.pi) % (2 * np.pi) - np.pi)
        R0 = G.daten(x_, False)[2]
        if R0.real < 0:
            R0 = -R0
        gz, gp, Rw, m, fx = K29.g_und_R(x_ - 1j * eps, "aussen_in", g0, R0, f0)
        xc = K29.xz(x_ + 1j * eps)
        P = m * np.conj(np.exp(xc) / K29.fp_in(xc))
        W[i] = P * Rw
        tau[i] = gz - (x_ + 1j * eps)
    return W, tau


def oben(eps):
    W = np.zeros(len(xi), complex)
    tau = np.zeros(len(xi), complex)
    for i, x_ in enumerate(xi):
        x0 = K29.xz(x_ + 1j * np.pi)
        f0 = K29.f_in(x0)
        g0 = np.log(K29.M(f0)) - 1j * th0
        lt = g0 - (x_ + 1j * np.pi)
        g0 = x_ + 1j * np.pi + lt.real + 1j * ((lt.imag + np.pi) % (2 * np.pi) - np.pi)
        R0 = G.daten(x_, True)[2]
        if R0.real < 0:
            R0 = -R0
        zl = x_ + 1j * (np.pi - eps)
        gz, gp, Rw, m_low, fx = K29.g_und_R(zl, "D", g0, R0)
        zu = x_ + 1j * (np.pi + eps)
        xu = K29.xz(zu)
        tau_re = gz.real - x_
        y = K29.finv(xu, K29.xz(zl - tau_re))
        assert abs(K29.f_in(y) - xu) < 1e-11, (x_, eps)
        m_up = np.exp(xu) / (np.exp(y) / K29.fp_in(y))
        W[i] = m_up * np.conj(m_low) * np.conj(Rw)
        tau[i] = gz - zu
    return W, tau


# Rand (wie o32)
dat_b = G.daten(xi + 0j, False)
dat_t = G.daten(xi + 0j, True)
pb = dat_b[3]
lt_, Rt_, pt_ = dat_t[1], dat_t[2], dat_t[3]
gtop = xi + 1j * np.pi + lt_
dg = np.gradient(gtop, xi)
w_top = pt_ * np.abs(Rt_) ** 2 / np.abs(dg)
E = 200.0
ext = np.arange(14 + hx, E, hx)
bz = np.concatenate([-ext[::-1], xi, ext]) + 0j
bw = np.concatenate([np.full(len(ext), pb[0]), pb, np.full(len(ext), pb[-1])])
tz = np.concatenate([gtop[0] + (-ext[::-1] - xi[0]), gtop, gtop[-1] + (ext - xi[-1])])
tw = np.concatenate([np.full(len(ext), w_top[0]), w_top, np.full(len(ext), w_top[-1])])
tds = np.concatenate([np.full(len(ext), hx), np.abs(dg) * hx, np.full(len(ext), hx)])
Z = np.concatenate([bz, tz])
Wr = np.concatenate([bw, tw])
dS = np.concatenate([np.full(len(bz), hx), tds])
xs = np.concatenate([-ext[::-1], xi, ext])


def schur(curves):
    out = []
    for r in (0.0, 0.5, 1.0):
        Asup = 0.0
        pts = []
        for (z, tau, rho) in curves:
            zc = np.concatenate([z[0] + (-ext[::-1] - xi[0]), z, z[-1] + (ext - xi[-1])])
            tc = np.concatenate([np.full(len(ext), tau[0]), tau, np.full(len(ext), tau[-1])])
            rc = np.concatenate([np.full(len(ext), rho[0]), rho, np.full(len(ext), rho[-1])])
            zr = zc + r * tc
            pts.append((zr, rc))
            for j in np.where(np.abs(xs) <= 30)[0]:
                Asup = max(Asup, np.sqrt(rc[j]) / (2 * np.pi) * np.sum(dS / (np.abs(Z - zr[j]) ** 2 * np.sqrt(Wr))))
        Bsup = 0.0
        for jz in np.where(np.abs(Z.real) <= 30)[0]:
            s_ = sum(np.sum(np.sqrt(rc) * hx / np.abs(Z[jz] - zr) ** 2) for zr, rc in pts)
            Bsup = max(Bsup, s_ / (2 * np.pi * np.sqrt(Wr[jz])))
        out.append(Asup * Bsup)
    return out


mi, mo = 0.7319047278299879, 1.1724160252312663


def cval(lam):
    s = 1 / lam
    b = np.linspace(0, s, 100001)
    return 1 / np.sqrt(((s ** 2 - b ** 2) / mi ** 2 + (1 + b) ** 2 / mo ** 2).max())


eb_list = [float(v) for v in sys.argv[1].split(",")] if len(sys.argv) > 1 else [0.3343, 0.40, 0.45]
et_list = [float(v) for v in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.3343, 0.30]
for eb in eb_list:
    Wb, tb = unten(eb)
    rb = np.abs(Wb) ** 2 * np.abs(tb) ** 2 / (4 * Wb.real)
    for et in et_list:
        Wt, tt = oben(et)
        rt = np.abs(Wt) ** 2 * np.abs(tt) ** 2 / (4 * Wt.real)
        if Wb.real.min() <= 0 or Wt.real.min() <= 0:
            print(f"ε_b={eb}, ε_t={et}: Re W ≤ 0 (unten {Wb.real.min():.3f}, oben {Wt.real.min():.3f})")
            continue
        lam_r = schur([(xi + 1j * eb, tb, rb), (xi + 1j * (np.pi + et), tt, rt)])
        Lam = max(lam_r)
        print(f"ε_b={eb:.4f}, ε_t={et:.4f}: min ReW unten/oben {Wb.real.min():.3f}/{Wt.real.min():.3f}, max |τ| {np.abs(tb).max():.3f}/{np.abs(tt).max():.3f}, "
              f"max ρ {rb.max():.4f}/{rt.max():.4f};  Schur Λ_r = {lam_r[0]:.4f}/{lam_r[1]:.4f}/{lam_r[2]:.4f}  ⇒ λ ≥ {1 - Lam:.4f}, c ≥ {cval(1 - Lam):.4f}", flush=True)
