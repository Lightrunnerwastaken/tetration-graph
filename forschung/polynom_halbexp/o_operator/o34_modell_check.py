"""O34 — Satz O11 (Konturverschiebung, waagrechte Konturen) für ein Polynom-Modell (radial / kubisch) als Prüfstein gegen Ms Satz M6.
Unten: c = ξ − iε_b, z = ξ + iε_b, τ = g(c) − z, W = p(c)R(c), p(ζ) = m(x(ζ))·conj(m(x(ζ̄))); oben: z = ξ + i(π+ε_t), ζ_l = ξ + i(π−ε_t),
τ = g(ζ_l) − z, W = m(x(z))·conj(m(x(ζ_l)))·conj(R(ζ_l)). f polynomial ⇒ alles ganz; m = f'∘f.
Ausgabe: min Re W, max ρ, Schur-Λ, λ = 1 − Λ und c(λ) (Satz O9) — Vergleich: M6 liefert σ_min ≥ 0.4249 (radial).
Aufruf: python o34_modell_check.py modell eps_b eps_t
"""
import sys

import numpy as np

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
name = sys.argv[1]
eb = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
et = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
f, fp = K.f_numpy(name)
L = K.Lc
Lb = L.conjugate()
RHO = abs(L)
m = lambda x: fp(f(x))
M = lambda x: (x - L) / (x - Lb)
Mp = lambda x: (L - Lb) / (x - Lb) ** 2
th0 = np.angle(M(-RHO))
mu = fp(L)
ell = np.log(mu)


def xz(z):
    w = np.exp(z + 1j * th0)
    return (L - Lb * w) / (1 - w)


def g_all(z, gref, Rref):
    x = xz(z)
    fx = f(x)
    gz = np.log(M(fx)) - 1j * th0
    gz = gz + 2j * np.pi * np.round((gref - gz).imag / (2 * np.pi))
    gp = (Mp(fx) * fp(x) / M(fx)) * (M(x) / Mp(x))
    Rw = np.sqrt(gp / fp(x)) / m(x)
    Rw = np.where(np.abs(Rw - Rref) > np.abs(-Rw - Rref), -Rw, Rw)
    return gz, gp, Rw, x


xi = np.linspace(-16, 16, 641)
hx = xi[1] - xi[0]
# Referenzen auf den Kanten (Zweige stetig, Re R > 0)
def kante(eta):
    z = xi + 1j * eta
    x = xz(z)
    gz = np.log(M(f(x))) - 1j * th0
    d = gz - z
    gz = z + d.real + 1j * ((d.imag + np.pi) % (2 * np.pi) - np.pi)
    gp = (Mp(f(x)) * fp(x) / M(f(x))) * (M(x) / Mp(x))
    Rw = np.sqrt(gp / fp(x)) / m(x)
    Rw = np.where(Rw.real < 0, -Rw, Rw)
    return gz, gp, Rw, x


gb0, gpb0, Rb0, xb0 = kante(0.0)
gt0, gpt0, Rt0, xt0 = kante(np.pi)
# unten
gzc, gpc, Rc, xc = g_all(xi - 1j * eb, gb0, Rb0)
xzz = xz(xi + 1j * eb)
Wb = m(xc) * np.conj(m(xzz)) * Rc
tb = gzc - (xi + 1j * eb)
# oben
gzl, gpl, Rl, xl = g_all(xi + 1j * (np.pi - et), gt0, Rt0)
xu = xz(xi + 1j * (np.pi + et))
Wt = m(xu) * np.conj(m(xl)) * np.conj(Rl)
tt = gzl - (xi + 1j * (np.pi + et))
rb = np.abs(Wb) ** 2 * np.abs(tb) ** 2 / (4 * Wb.real)
rt = np.abs(Wt) ** 2 * np.abs(tt) ** 2 / (4 * Wt.real)
print(f"{name}: μ = {mu:.5f}, ε_b = {eb}, ε_t = {et}: min Re W unten/oben = {Wb.real.min():.4f}/{Wt.real.min():.4f}, "
      f"max|τ| {np.abs(tb).max():.3f}/{np.abs(tt).max():.3f}, max ρ {rb.max():.4f}/{rt.max():.4f}")
# Rand: Boden (p = |m|²), Deckel g(ξ+iπ) mit Gewicht |m|²|R|²/|g'|
pb = np.abs(m(xb0)) ** 2
w_top = np.abs(m(xt0)) ** 2 * np.abs(Rt0) ** 2 / np.abs(gpt0)
E = 200.0
ext = np.arange(16 + hx, E, hx)
bz = np.concatenate([-ext[::-1], xi, ext]) + 0j
bw = np.concatenate([np.full(len(ext), pb[0]), pb, np.full(len(ext), pb[-1])])
tz = np.concatenate([gt0[0] + (-ext[::-1] - xi[0]), gt0, gt0[-1] + (ext - xi[-1])])
tw = np.concatenate([np.full(len(ext), w_top[0]), w_top, np.full(len(ext), w_top[-1])])
tds = np.concatenate([np.full(len(ext), hx), np.abs(gpt0) * hx, np.full(len(ext), hx)])
Z = np.concatenate([bz, tz])
Wr = np.concatenate([bw, tw])
dS = np.concatenate([np.full(len(bz), hx), tds])
xs = np.concatenate([-ext[::-1], xi, ext])
curves = [(xi + 1j * eb, tb, rb), (xi + 1j * (np.pi + et), tt, rt)]
Lams = []
for r in (0.0, 0.5, 1.0):
    pts = []
    Asup = 0.0
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
    Lams.append(Asup * Bsup)
Lam = max(Lams)
lam = 1 - Lam
# Satz O9
mi = np.abs(m(xb0)).min()
mo = np.abs(m(xt0)).min()
s = 1 / lam
bgrid = np.linspace(0, s, 100001)
cval = 1 / np.sqrt(((s ** 2 - bgrid ** 2) / mi ** 2 + (1 + bgrid) ** 2 / mo ** 2).max())
print(f"Schur Λ_r = {Lams[0]:.4f}/{Lams[1]:.4f}/{Lams[2]:.4f} ⇒ λ ≥ {lam:.4f};  m_in = {mi:.4f}, m_out = {mo:.4f} ⇒ c_∞ ≥ {cval:.4f}")
