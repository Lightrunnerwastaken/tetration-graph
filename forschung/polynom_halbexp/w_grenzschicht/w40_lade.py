"""Lädt das Leitordnungsmodell aus w40 (npz) und liefert H(t) (unskaliert auf ρ = 1, t = x/ρ) sowie das Γin-Profil Λ_Modell(τ)."""
import numpy as np

import w_kern as W


def lade(K=20, NC=3, cut=100, ref=True):
    d = np.load(W.HIER + ("/w40_modell_ref.npz" if ref else f"/w40_modell_K{K}_n{NC}_c{cut}.npz"))
    sol, src, gam, P, K1, nc = d["sol"], d["src"], float(d["gam"]), complex(d["P"]), complex(d["K1"]), int(d["NC"])
    Lt = W.L / W.RHO

    def ecke(t):
        z1 = np.exp(1j * gam) * (t - Lt) / (t - P)
        z2 = np.conj(np.exp(1j * gam) * (np.conj(t) - Lt) / (np.conj(t) - P))
        cols = []
        for n in range(1, nc + 1):
            w1 = np.exp(n * K1 * np.log(z1))
            w2 = np.exp(n * np.conj(K1) * np.log(z2))
            cols += [w1.real + w2.real, -w1.imag + w2.imag]
        return np.array(cols).T

    def H(t):
        t = np.atleast_1d(np.asarray(t, complex))
        M = np.hstack([np.log(np.abs(t[:, None] - src[None, :])), np.ones((len(t), 1)), ecke(t)])
        return -np.log(np.abs(t)) + M @ sol

    def Lam(tau_deg):
        """Profil auf Γin: H(e^{i(arg L + τ)}), τ in Grad (≥ 0)."""
        return H(np.exp(1j * (np.angle(Lt) + np.deg2rad(np.atleast_1d(tau_deg)))))
    return H, Lam, d["rB"], d["phi_grid"]
