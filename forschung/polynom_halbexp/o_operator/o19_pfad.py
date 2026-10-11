"""O19 — Pfadbedingungen (S) für die Index-Homotopie von Satz O7: f_t = (1−t) f_M + t f_Z, t ∈ [0,1],
f_M = Möbius-exp (Strang E, invertierbar nach Satz O2), f_Z = Zielmodell (kubisch bzw. radial, μ = f_Z'(L)).
Alle f_t fixieren L, L̄ (konvexe Kombination). Geprüft auf feinen Gittern (numerisch, kein Einschluss):
 (S1) f_t holomorph um D̄ ∪ f_t(D̄): Pol 3.9527 von f_M nicht in f_t(D̄) ∪ D̄ (min Abstand);
 (S2) f_t univalent auf D̄: f_t' ≠ 0 auf D̄ (min |f_t'|) und f_t(∂D) einfache Kurve (min Abstand nicht-benachbarter Punkte / Bogenlänge);
 (S3) Γin/Γout: |f_t(ρe^{iφ})| − ρ wechselt das Vorzeichen genau in ±arg L;
 (S4) m_t = f_t'∘f_t ≠ 0 auf D̄ (min |m_t| auf ∂D und im Innern), ∂Ω_t = Γin ∪ f_t(Γout) einfach;
 (S5) lokales Symbol: min_σ |1 + μ_t^{−3/2+iσ}| > 0, 0 < arg μ_t < π.
Aufruf: python o19_pfad.py kubisch|radial
"""
import sys

import numpy as np

import o_kern as K

sys.stdout.reconfigure(encoding="utf-8")
ziel = sys.argv[1] if len(sys.argv) > 1 else "kubisch"
L = K.Lc
R = K.RHO
fM, fMp = K.f_numpy("moeb_ell")
fZ, fZp = K.f_numpy(ziel)
pol = 3.952682509
phi = np.linspace(0, 2 * np.pi, 20001)[:-1]
x = R * np.exp(1j * phi)
# Innengitter
rr, pp = np.meshgrid(np.linspace(0, 1, 120), np.linspace(0, 2 * np.pi, 480, endpoint=False))
xin = (R * rr * np.exp(1j * pp)).ravel()
thL = np.angle(L)
sig = np.linspace(-80, 80, 32001)
print(f"Pfad Möbius-exp → {ziel}:")
for t in np.linspace(0, 1, 11):
    f = lambda z: (1 - t) * fM(z) + t * fZ(z)
    fp = lambda z: (1 - t) * fMp(z) + t * fZp(z)
    fx, fin = f(x), f(xin)
    # (S1)
    d_pol = min(np.abs(np.concatenate([fx, fin]) - pol).min(), np.abs(np.concatenate([x, xin]) - pol).min()) if t < 1 else np.inf
    # (S2)
    min_fp = min(np.abs(fp(x)).min(), np.abs(fp(xin)).min())
    # Einfachheit der Randkurve: Selbstabstand (Punkte mit Indexabstand > 50)
    sub = fx[::20]
    D = np.abs(sub[:, None] - sub[None, :])
    n = len(sub)
    idx = np.arange(n)
    gap = np.minimum(np.abs(idx[:, None] - idx[None, :]), n - np.abs(idx[:, None] - idx[None, :]))
    selbst = D[gap > 25].min()
    # (S3)
    s = np.sign(np.abs(fx) - R)
    wechsel = phi[np.nonzero(s != np.roll(s, 1))[0]]
    w_deg = np.degrees(wechsel)
    # (S4)
    m_b = np.abs(fp(fx)).min()
    m_i = np.abs(fp(fin)).min()
    # (S5)
    mu = fp(L)
    zz = np.exp((-1.5 + 1j * sig) * np.log(mu))
    sym = np.abs(1 + zz).min()
    print(f" t={t:.1f}: μ_t={mu:.4f} (|μ|={abs(mu):.4f}, arg={np.angle(mu):.4f}); Pol-Abstand {d_pol:.3f}; min|f'| {min_fp:.3f}; "
          f"Selbstabstand f(∂D) {selbst:.3f}; Vorzeichenwechsel |f|−ρ bei {np.array2string(w_deg, precision=2)}°; "
          f"min|m| Rand/innen {m_b:.3f}/{m_i:.3f}; min|1+μ^(s−1)| {sym:.4f}", flush=True)
print(f"(arg L = {np.degrees(thL):.3f}°)")
