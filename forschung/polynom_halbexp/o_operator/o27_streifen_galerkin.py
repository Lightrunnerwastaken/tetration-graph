"""O27 — Galerkin für (FA) in der Streifenkoordinate ζ = log M(x) − iθ₀ (o26): D = {0 < Im ζ < π}.
δ = u·J^{−1/2} (J = dx/dζ, H²-isometrisch), Tδ ↔ u + R·u∘g, R = (g'/f')^{1/2}/m, B-Daten: u (unten), R·u∘g (oben = Γout).
  G_B = ∫_unten p|u|² dξ + ∫_oben p|R u∘g|² dξ,   H = Re∫_unten p (u + R u∘g) ū + Re∫_oben p (u + R u∘g) conj(R u∘g),
  λ_Gal = min H/G_B über span{Wellenpakete φ_{c,σ}(ζ) = exp(−(ζ − c)²/(2s²) + iσζ)}  (ganze Funktionen ⇒ in E²(Ω̃)).
Quadratur: Trapezregel in ξ (Integranden analytisch, gaußsch abfallend). Normierung der Basis gegen Auslöschung.
Aufruf: python o27_streifen_galerkin.py modell A dc s Xi dsig
"""
import sys

import numpy as np
import mpmath as mp

sys.argv += [] if len(sys.argv) > 1 else ["kneser"]
import o26_streifen as G   # rechnet beim Import die Tabelle (Ausgabe wird unten unterdrückt)

name = sys.argv[1]
A = float(sys.argv[2]) if len(sys.argv) > 2 else 12.0
dc = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
s = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
Xi = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
dsig = float(sys.argv[6]) if len(sys.argv) > 6 else 0.5

xi = np.arange(-A - 8 * s, A + 8 * s + 1e-9, 0.02)
dat_b = [G.daten(t, False) for t in xi]
dat_t = [G.daten(t, True) for t in xi]
lb = np.array([d[1] for d in dat_b]); Rb = np.array([d[2] for d in dat_b]); pb = np.array([d[3] for d in dat_b])
lt = np.array([d[1] for d in dat_t]); Rt = np.array([d[2] for d in dat_t]); pt = np.array([d[3] for d in dat_t])
# Wurzelzweig stetig ab dem L-Ende (μ^{-3/2})
for Rarr in (Rb, Rt):
    prev = G.mu ** -1.5
    for i in range(len(Rarr)):
        if abs(Rarr[i] - prev) > abs(-Rarr[i] - prev):
            Rarr[i] = -Rarr[i]
        prev = Rarr[i]
zb = xi + 0j
zt = xi + 1j * np.pi
gb = zb + lb
gt = zt + lt
cs = np.arange(-A, A + 1e-9, dc)
sg = np.arange(-Xi, Xi + 1e-9, dsig)
basis = [(c, q) for c in cs for q in sg]
nb = len(basis)
print(f"{name}: Basis {nb} Pakete (c ∈ [−{A},{A}] Δc={dc}, s={s}, σ ∈ [−{Xi},{Xi}] Δσ={dsig}), Gitter {len(xi)}", flush=True)


def phi(z, c, q):
    return np.exp(-(z - c) ** 2 / (2 * s * s) + 1j * q * z)


Ub = np.array([phi(zb, c, q) for c, q in basis]).T      # (Punkte, Basis)
Ugb = np.array([phi(gb, c, q) for c, q in basis]).T
Ut = np.array([phi(zt, c, q) for c, q in basis]).T
Ugt = np.array([phi(gt, c, q) for c, q in basis]).T
w = 0.02
Bb = Ub
Bt = Rt[:, None] * Ugt
Tb = Ub + Rb[:, None] * Ugb
Tt = Ut + Rt[:, None] * Ugt
GB = (Bb.conj().T @ (pb[:, None] * Bb) + Bt.conj().T @ (pt[:, None] * Bt)) * w
HH = (Bb.conj().T @ (pb[:, None] * Tb) + Bt.conj().T @ (pt[:, None] * Tt)) * w
HH = (HH + HH.conj().T) / 2
GB = (GB + GB.conj().T) / 2
d = np.sqrt(np.real(np.diag(GB)))
GBn = GB / np.outer(d, d)
HHn = HH / np.outer(d, d)
ev, V = np.linalg.eigh(GBn)
keep = ev > 1e-11 * ev.max()
print(f"Gram: Eigenwerte {ev.min():.2e} … {ev.max():.2e}, behalten {keep.sum()}/{nb}", flush=True)
P = V[:, keep] / np.sqrt(ev[keep])
M = P.conj().T @ HHn @ P
lam, W = np.linalg.eigh((M + M.conj().T) / 2)
print(f"λ_Gal = {lam[0]:.6f}  (zweiter {lam[1]:.6f}),  max {lam[-1]:.4f}")
v = P @ W[:, 0]
coef = v / d
# Lage des Minimalvektors: Gewicht je Zentrum c und je Frequenz σ
amp = np.abs(coef) ** 2
byc = {}
bys = {}
for (c, q), a_ in zip(basis, amp):
    byc[c] = byc.get(c, 0) + a_
    bys[q] = bys.get(q, 0) + a_
tot = sum(byc.values())
print("Gewicht nach Zentrum c:", " ".join(f"{c:+.0f}:{byc[c] / tot:.2f}" for c in cs if byc[c] / tot > 0.02))
print("Gewicht nach Frequenz σ:", " ".join(f"{q:+.1f}:{bys[q] / tot:.2f}" for q in sg if bys[q] / tot > 0.02))
