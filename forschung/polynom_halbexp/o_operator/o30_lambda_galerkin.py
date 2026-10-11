"""O30 — Galerkin-Schätzung von Λ := sup_u D(u)/‖Bu‖²,  D(u) = Σ_{unten,oben} ∫ |W|²/(4 Re W) · |u(z+τ) − u(z)|² dξ
(Mittellinien aus o29). Satz-Kette: X(u) ≥ −D(u) ≥ −Λ‖Bu‖²  ⇒  (FA) mit λ = 1 − Λ.
Galerkin liefert eine UNTERE Schätzung von Λ (sup über Teilraum); die strenge OBERE Schranke kommt aus o31.
Basis wie o27 (Wellenpakete exp(−(ζ−c)²/(2s²) + iσζ)), ‖Bu‖² exakt auf Γin und Linsendeckel (über g auf Γout).
"""
import sys

import numpy as np

import o26_streifen as G

sys.stdout.reconfigure(encoding="utf-8")
A = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
dc = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
s = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
Xi = float(sys.argv[4]) if len(sys.argv) > 4 else 8.0
dsig = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5

ub0 = np.load("o29_unten.npz")
ob0 = np.load("o29_oben.npz")
xi0 = ub0["xi"]
h = xi0[1] - xi0[0]
XE = float(sys.argv[6]) if len(sys.argv) > 6 else 40.0
xi = np.arange(-XE, XE + h / 2, h)


def ext(arr):
    # konstant (Modellwert) fortsetzen ausserhalb [−12, 12]
    return np.interp(xi, xi0, arr.real) + (1j * np.interp(xi, xi0, arr.imag) if np.iscomplexobj(arr) else 0)


ub = {k: ext(ub0[k]) for k in ("eps", "tau", "W")}
ob = {k: ext(ob0[k]) for k in ("eps", "tau", "W")}
# B-Norm-Daten: bis |ξ| ≤ 18 aus o26, sonst Modellwerte
xc = np.clip(xi, -18, 18)
dat_b = G.daten(xc + 0j, False)
dat_t = G.daten(xc + 0j, True)
lb, Rb, pb = dat_b[1], dat_b[2], dat_b[3]
lt, Rt, pt = dat_t[1], dat_t[2], dat_t[3]
lt = lt + 0 * xi
zb = xi + 0j
zt = xi + 1j * np.pi
gt = zt + lt
zb_m = xi + 1j * ub["eps"]
zt_m = xi + 1j * (np.pi + ob["eps"])
rb = np.abs(ub["W"]) ** 2 / (4 * ub["W"].real)
rt = np.abs(ob["W"]) ** 2 / (4 * ob["W"].real)
print(f"Kurven: max |W|²/(4ReW): unten {rb.max():.3f}, oben {rt.max():.3f}; Integration ξ ∈ [−{XE}, {XE}]")
cs = np.arange(-A, A + 1e-9, dc)
sg = np.arange(-Xi, Xi + 1e-9, dsig)
basis = [(c, q) for c in cs for q in sg]
nb = len(basis)


def phi(z, c, q):
    return np.exp(-(z - c) ** 2 / (2 * s * s) + 1j * q * z)


def mat(z):
    return np.array([phi(z, c, q) for c, q in basis]).T


Ub = mat(zb)
Ugt = mat(gt)
Bt = Rt[:, None] * Ugt
GB = (Ub.conj().T @ (pb[:, None] * Ub) + Bt.conj().T @ (pt[:, None] * Bt)) * h
Db_ = mat(zb_m + ub["tau"]) - mat(zb_m)
Dt_ = mat(zt_m + ob["tau"]) - mat(zt_m)
DD = (Db_.conj().T @ (rb[:, None] * Db_) + Dt_.conj().T @ (rt[:, None] * Dt_)) * h
GB = (GB + GB.conj().T) / 2
DD = (DD + DD.conj().T) / 2
d = np.sqrt(np.real(np.diag(GB)))
GBn = GB / np.outer(d, d)
DDn = DD / np.outer(d, d)
ev, V = np.linalg.eigh(GBn)
keep = ev > 1e-11 * ev.max()
P = V[:, keep] / np.sqrt(ev[keep])
Mx = P.conj().T @ DDn @ P
lam = np.linalg.eigvalsh((Mx + Mx.conj().T) / 2)
print(f"Basis {nb} (behalten {keep.sum()}): Λ_Galerkin = {lam[-1]:.5f} (zweiter {lam[-2]:.5f})  ⇒  λ = 1 − Λ ≈ {1 - lam[-1]:.4f}")

# Probe: Rayleigh-Quotient des Maximalvektors direkt nachrechnen (feineres Gitter, ohne Gram-Matrizen)
lam2, Wv = np.linalg.eigh((Mx + Mx.conj().T) / 2)
v = P @ Wv[:, -1]
coef = v / d
xi_f = np.linspace(xi[0], xi[-1], 2 * (len(xi) - 1) + 1)
hf = xi_f[1] - xi_f[0]


def ufun(z):
    out = np.zeros(z.shape, complex)
    for (c, q), cf in zip(basis, coef):
        out += cf * phi(z, c, q)
    return out


db = G.daten(np.clip(xi_f, -18, 18) + 0j, False)
dt_ = G.daten(np.clip(xi_f, -18, 18) + 0j, True)
Bn = np.sum(db[3] * np.abs(ufun(xi_f + 0j)) ** 2) * hf + np.sum(dt_[3] * np.abs(dt_[2]) ** 2 * np.abs(ufun(xi_f + 1j * np.pi + dt_[1] + 0 * xi_f)) ** 2) * hf
epsb = np.interp(xi_f, xi, ub["eps"].real); taub = np.interp(xi_f, xi, ub["tau"].real); rbf = np.interp(xi_f, xi, rb)
epst = np.interp(xi_f, xi, ob["eps"].real); taut = np.interp(xi_f, xi, ob["tau"].real); rtf = np.interp(xi_f, xi, rt)
zbm = xi_f + 1j * epsb
ztm = xi_f + 1j * (np.pi + epst)
Dn = np.sum(rbf * np.abs(ufun(zbm + taub) - ufun(zbm)) ** 2) * hf + np.sum(rtf * np.abs(ufun(ztm + taut) - ufun(ztm)) ** 2) * hf
print(f"Probe Maximalvektor: D/‖B‖² direkt = {Dn / Bn:.5f}  (Matrix-Eigenwert {lam2[-1]:.5f})")
