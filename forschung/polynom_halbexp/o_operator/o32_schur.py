"""O32 — Schur-Test für Λ: D(u) ≤ ∫₀¹ D_r(u) dr,  D_r(u) = Σ_Kurven ∫ ρ(ξ) |u'(z(ξ) + rτ(ξ))|² dξ,  ρ = |W|²τ²/(4 Re W),
u'(z) = (1/2πi) ∫_{∂Ω̃} u(ζ) dζ/(ζ − z)²  (Cauchy; ∂Ω̃ = Boden ℝ ∪ Linsendeckel g(s + iπ)).
Kern k(ξ, ζ) = √ρ(ξ) / (2π |ζ − z_r(ξ)|² √w(ζ)),  w = B-Gewicht (Boden p_b, Deckel p_t|R_t|²/|g'|).
Schur: ‖K‖² ≤ sup_ξ ∫ k |dζ| · sup_ζ ∫ k dξ   (mit Testfunktionen h ≡ 1; optional Gewichtung).
Numerisch (Vorstufe der strengen Rechnung): Rand und Kurven auf ξ ∈ [−12, 12], Enden konstant fortgesetzt (Modellwerte).
"""
import sys

import numpy as np

import o26_streifen as G

sys.stdout.reconfigure(encoding="utf-8")
ub = np.load("o29_unten.npz")
ob = np.load("o29_oben.npz")
xi = ub["xi"]
hx = xi[1] - xi[0]
# Randdaten
dat_b = G.daten(xi + 0j, False)
dat_t = G.daten(xi + 0j, True)
pb = dat_b[3]
lt, Rt, pt = dat_t[1], dat_t[2], dat_t[3]
gtop = xi + 1j * np.pi + lt                     # Deckelpunkte g(s + iπ), s = xi
dg = np.gradient(gtop, xi)                      # g'(s+iπ) längs des Deckels (|dζ| = |g'| ds)
w_top = pt * np.abs(Rt) ** 2 / np.abs(dg)
# Rand als Liste: Punkte, Bogenmass |dζ|, Gewicht w — mit Verlängerung der Enden (Modellwerte) bis ±E
E = 200.0
ext = np.arange(12 + hx, E, hx)
bx = np.concatenate([-ext[::-1], xi, ext])
bz = bx + 0j
bw = np.concatenate([np.full(len(ext), pb[0]), pb, np.full(len(ext), pb[-1])])
bds = np.full(len(bx), hx)
tz = np.concatenate([gtop[0] + (-ext[::-1] - xi[0]), gtop, gtop[-1] + (ext - xi[-1])])
tw = np.concatenate([np.full(len(ext), w_top[0]), w_top, np.full(len(ext), w_top[-1])])
tds = np.concatenate([np.full(len(ext), hx), np.abs(dg) * hx, np.full(len(ext), hx)])
Z = np.concatenate([bz, tz])
Wr = np.concatenate([bw, tw])
dS = np.concatenate([bds, tds])
print(f"Randgewichte: Boden w ∈ [{bw.min():.3f}, {bw.max():.3f}], Deckel w ∈ [{tw.min():.3f}, {tw.max():.3f}]")


def kurve(dat, oben):
    W = dat["W"]
    tau = dat["tau"]
    eps = dat["eps"]
    rho = np.abs(W) ** 2 * tau ** 2 / (4 * W.real)
    z = xi + 1j * ((np.pi if oben else 0.0) + eps)
    # Enden konstant verlängern
    zc = np.concatenate([z[0] + (-ext[::-1] - xi[0]), z, z[-1] + (ext - xi[-1])])
    tc = np.concatenate([np.full(len(ext), tau[0]), tau, np.full(len(ext), tau[-1])])
    rc = np.concatenate([np.full(len(ext), rho[0]), rho, np.full(len(ext), rho[-1])])
    return zc, tc, rc


kb = kurve(ub, False)
kt = kurve(ob, True)
xs = np.concatenate([-ext[::-1], xi, ext])
for r in (0.0, 0.5, 1.0):
    Asup, Bsum = 0.0, np.zeros(len(Z))
    parts = []
    for (zc, tc, rc) in (kb, kt):
        zr = zc + r * tc
        # A(ξ) = √ρ/(2π) ∫ |dζ|/(|ζ − z|² √w)   (nur Punkte mit |ξ| ≤ 30 auswerten, Rest Modell)
        sel = np.abs(xs) <= 30
        Avals = []
        for j in np.where(sel)[0]:
            dd = np.abs(Z - zr[j]) ** 2
            Avals.append(np.sqrt(rc[j]) / (2 * np.pi) * np.sum(dS / (dd * np.sqrt(Wr))))
        Asup = max(Asup, max(Avals))
        parts.append(max(Avals))
        # B(ζ) = (1/√w) (1/2π) ∫ √ρ dξ/|ζ − z(ξ)|²
        for jz in range(0, len(Z)):
            pass
    # B(ζ) für Randpunkte mit |Re ζ| ≤ 30 (vektorisiert über Kurvenpunkte)
    Bsup = 0.0
    zr_b = kb[0] + r * kb[1]
    zr_t = kt[0] + r * kt[1]
    selz = np.abs(Z.real) <= 30
    for jz in np.where(selz)[0]:
        s_ = np.sum(np.sqrt(kb[2]) * hx / np.abs(Z[jz] - zr_b) ** 2) + np.sum(np.sqrt(kt[2]) * hx / np.abs(Z[jz] - zr_t) ** 2)
        Bsup = max(Bsup, s_ / (2 * np.pi * np.sqrt(Wr[jz])))
    print(f"r = {r}: sup_ξ A = {Asup:.4f} (unten {parts[0]:.4f}, oben {parts[1]:.4f}),  sup_ζ B = {Bsup:.4f},  "
          f"Schur Λ_r ≤ {Asup * Bsup:.4f}")
