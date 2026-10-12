"""mk_kneser — Kneser-Halbexponentialfunktion f_K in der normierten Variablen z = x/ρ (ρ = |L|), Gleitkomma-Fassung.
F(z) := f_K(ρz)/ρ = (1/ρ) Σ b_k z^k,  b_k := a_k ρ^k  (a_k aus gemeinsam/kneser_halbexp_koeff_4000.json, Strang A; NICHT zertifiziert).
F(F(z)) = e^{ρz}/ρ;  m(z) = F'(F(z)) = e^{ρz}/F'(z);  F(e^{iθ}) = e^{iθ}, F'(e^{iθ}) = √L, θ = arg L.
Fortsetzung ausserhalb D̄ (nur über die Funktionalgleichung, keine weiteren Daten):
  nahe Γin (eˣ ∈ D):   F(z) = log(ρ F(w))/ρ,  F'(z) = F'(w)·w/F(w),  w = e^{ρz}/ρ
  nahe Γout (Linse):   F(z) = e^{ρ h(z)}/ρ,   F'(z) = e^{ρ h(z)}/F'(h(z)),  h = F⁻¹-Zweig mit Werten in D
Benutzt nur die Strang-A-Daten (kein Motor).
"""
import json
import os

import mpmath as mp
import numpy as np

HIER = os.path.dirname(os.path.abspath(__file__))
GEM = os.path.join(HIER, "..", "..", "..", "gemeinsam")  # V-Kopie: Pfad angepasst


def lade_b(K=None, dps=60):
    """b_k = a_k ρ^k als float-Array (k = 0..K) und L (mpmath)."""
    mp.mp.dps = dps
    d = json.load(open(os.path.join(GEM, "kneser_halbexp_koeff_4000.json"), encoding="utf-8"))
    s = json.load(open(os.path.join(GEM, "kneser_singularitaet_L.json"), encoding="utf-8"))
    L = mp.mpc(mp.mpf(s["L"][0]), mp.mpf(s["L"][1]))
    rho = abs(L)
    a = d["a"] if K is None else d["a"][:K + 1]
    b = np.array([float(mp.mpf(x) * rho ** k) for k, x in enumerate(a)])
    return b, L


B, L_MP = lade_b()
L = complex(L_MP)
RHO = abs(L)
TH = float(mp.arg(L_MP))
MU = complex(mp.sqrt(L_MP))
EITH = np.exp(1j * TH)
KK = np.arange(len(B))
B1 = (KK * B)[1:]          # Koeffizienten von ρF'(z) = Σ k b_k z^{k−1}


def _potenzsumme(c, z, block=256):
    """Σ c_k z^k vektorisiert (np.polyval läuft als Python-Schleife über alle Koeffizienten — viel zu langsam)."""
    z = np.asarray(z, complex)
    flach = z.ravel()
    out = np.empty(flach.shape, complex)
    k = np.arange(len(c))
    for i in range(0, len(flach), block):
        zz = flach[i:i + block]
        lz = np.log(np.where(zz == 0, 1, zz))
        P = np.exp(np.outer(lz, k))
        P[zz == 0, 1:] = 0
        P[zz == 0, 0] = 1
        out[i:i + block] = P @ c
    return out.reshape(z.shape)


def F_innen(z):
    return _potenzsumme(B, z) / RHO


def Fp_innen(z):
    return _potenzsumme(B1, z) / RHO


def h_newton(w, v0=None, it=40):
    """h(w): F(h) = w, h ∈ D. Startwert: lineare Näherung an der nächsten Ecke."""
    w = np.asarray(w, complex)
    if v0 is None:
        oben = w.imag >= 0
        v0 = np.where(oben, EITH + (w - EITH) / MU, np.conj(EITH) + (w - np.conj(EITH)) / np.conj(MU))
    v = np.array(v0, complex)
    for _ in range(it):
        dv = (_potenzsumme(B[:801], v) / RHO - w) / (_potenzsumme(B1[:800], v) / RHO)   # Startwerte: 800 Glieder (Rest < 1e-16)
        v = v - dv
        if np.all(np.abs(dv) < 1e-15):
            break
    return v


def F(z, art="auto"):
    """F an beliebigen z nahe D̄: innen Taylor, aussen nahe Γin per log, aussen nahe Γout per e^{ρh}."""
    z = np.asarray(z, complex)
    out = np.empty_like(z)
    innen = np.abs(z) <= 1
    out[innen] = F_innen(z[innen])
    aus = ~innen
    if np.any(aus):
        za = z[aus]
        gin = za.real < np.cos(TH)            # Γin-Seite (cos φ < cos θ)
        res = np.empty_like(za)
        if np.any(gin):
            w = np.exp(RHO * za[gin]) / RHO
            res[gin] = np.log(RHO * F_innen(w)) / RHO
        if np.any(~gin):
            v = h_newton(za[~gin])
            res[~gin] = np.exp(RHO * v) / RHO
        out[aus] = res
    return out


def Fp(z):
    z = np.asarray(z, complex)
    out = np.empty_like(z)
    innen = np.abs(z) <= 1
    out[innen] = Fp_innen(z[innen])
    aus = ~innen
    if np.any(aus):
        za = z[aus]
        gin = za.real < np.cos(TH)
        res = np.empty_like(za)
        if np.any(gin):
            w = np.exp(RHO * za[gin]) / RHO
            res[gin] = Fp_innen(w) * w * RHO / (RHO * F_innen(w))
        if np.any(~gin):
            v = h_newton(za[~gin])
            res[~gin] = np.exp(RHO * v) / Fp_innen(v)
        out[aus] = res
    return out


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    print("ρ =", RHO, "θ =", TH, "μ = √L =", MU, "K =", len(B) - 1)
    print("F(e^{iθ}) − e^{iθ} =", F_innen(EITH) - EITH, " F'(e^{iθ}) − √L =", Fp_innen(EITH) - MU)
    # Funktionalgleichung im Innern prüfen
    zz = 0.7 * np.exp(1j * np.linspace(0, 6, 7))
    print("max|F(F(z)) − e^{ρz}/ρ| (|z| = 0.7):", np.abs(F(F(zz)) - np.exp(RHO * zz) / RHO).max())
    # Stetigkeit der Fortsetzung über den Kreis
    for ph in (2.0, 2.8, np.pi, -2.2):
        z1, z2 = 0.999999 * np.exp(1j * ph), 1.000001 * np.exp(1j * ph)
        print(f"φ = {ph:.2f}: Sprung F {abs(F(np.array([z1]))[0]-F(np.array([z2]))[0]):.1e}, F' {abs(Fp(np.array([z1]))[0]-Fp(np.array([z2]))[0]):.1e}")
    for ph in (0.0, 0.6, 1.2):
        z1, z2 = 0.999999 * np.exp(1j * ph), 1.000001 * np.exp(1j * ph)
        print(f"φ = {ph:.2f} (Γout): Sprung F {abs(F(np.array([z1]))[0]-F(np.array([z2]))[0]):.1e}, F' {abs(Fp(np.array([z1]))[0]-Fp(np.array([z2]))[0]):.1e}")
    phs = np.linspace(TH + 1e-3, 2 * np.pi - TH - 1e-3, 2001)
    print("max |F| auf Γin:", np.abs(F_innen(np.exp(1j * phs))).max(), " min |m| auf Γin:",
          np.abs(Fp_innen(F_innen(np.exp(1j * phs)))).min())
    pho = np.linspace(-TH + 1e-3, TH - 1e-3, 2001)
    v = h_newton(np.exp(1j * pho))
    print("max |h| auf Γout:", np.abs(v).max(), " Rest |F(h) − w|:", np.abs(F_innen(v) - np.exp(1j * pho)).max())
