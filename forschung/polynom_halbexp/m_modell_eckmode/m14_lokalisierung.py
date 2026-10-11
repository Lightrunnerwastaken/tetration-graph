"""M14 — Lokalisierung der Akkretivität von 𝒯_P = I + 𝒦_P auf L²(∂D) (Poisson-Fortsetzung), radiales Modell, Exp-Lage.
Zerlegung der Eins χ_L² + χ_L̄² + χ_m² = 1 (logarithmischer Übergang zwischen η₁ und η₂ = η₁e^Λ um ±θ).
IMS: Re<𝒯u,u> = Σ_j Re<𝒯χ_ju, χ_ju> + Re<(𝒦_P − Σ_j χ_j𝒦_Pχ_j)u, u>.
Gemessen (Fourier −K..K, Gitter): λ_min Re(I + χ_j𝒦χ_j) für Ecke und Mitte, Norm des Kommutatorterms R = Re(𝒦 − Σχ𝒦χ).
Aufruf: python m14_lokalisierung.py eta1 eta2 K [G]
"""
import sys
import time

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
eta1, eta2, K = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
G = int(sys.argv[4]) if len(sys.argv) > 4 else 60000
eps, thg = 0.31382, 76.618
th = np.deg2rad(thg)
c0 = 1 - 2 * eps * np.cos(th)
f = lambda z: eps + c0 * z + eps * z * z
fp = lambda z: c0 + 2 * eps * z
h = lambda w: (-c0 + np.sqrt(c0 * c0 - 4 * eps * (eps - w) + 0j)) / (2 * eps)

# Gitter mit Verdichtung an ±θ: Mittelpunktsregel in einer Variablen t mit φ = θ ± (Abstand) ...
# einfache Variante: gleichmässig, aber fein
phi = (np.arange(G) + 0.5) / G * 2 * np.pi - np.pi
wq = 2 * np.pi / G
z = np.exp(1j * phi)
gin = np.cos(phi) < np.cos(th)
pi_ = np.where(gin, f(z), h(z))
w = np.where(gin, 1 / fp(f(z)), fp(z))


def chi_ecke(phi, c):
    d = np.abs(np.angle(np.exp(1j * (phi - c))))
    t = (np.log(np.maximum(d, 1e-300)) - np.log(eta1)) / np.log(eta2 / eta1)
    t = np.clip(t, 0, 1)
    return np.cos(np.pi * t / 2)


cL = chi_ecke(phi, th)
cLq = chi_ecke(phi, -th)
cm = np.sqrt(np.clip(1 - cL ** 2 - cLq ** 2, 0, 1))
k = np.arange(-K, K + 1)
E = np.exp(1j * np.outer(phi, k))                        # e^{ikφ}
# Poisson-Fortsetzung von e^{ikψ} am Punkt p: p^k (k ≥ 0), conj(p)^{|k|} (k < 0)
Pp = np.where(k[None, :] >= 0, pi_[:, None] ** np.maximum(k, 0)[None, :],
              np.conj(pi_)[:, None] ** np.maximum(-k, 0)[None, :])


def form(chi_left, chi_right):
    """Matrix B_{lk} = <χ_left 𝒦_P (χ_right e^{ik·}), e^{il·}> — braucht P(χ_right e^{ik·}) am Punkt π(φ):
    χ_right ist keine trigonometrische Funktion; Näherung: P(χ e^{ik·})(π(φ)) ≈ Σ_j ĉ_j P(e^{i(k+j)·})(π(φ)) mit den
    Fourier-Koeffizienten ĉ von χ (abgeschnitten bei |j| ≤ J)."""
    J = 4 * K
    jj = np.arange(-J, J + 1)
    chat = (np.exp(-1j * np.outer(jj, phi)) @ chi_right) * wq / (2 * np.pi)
    kk = np.arange(-K - J, K + J + 1)
    Pbig = np.where(kk[None, :] >= 0, pi_[:, None] ** np.maximum(kk, 0)[None, :],
                    np.conj(pi_)[:, None] ** np.maximum(-kk, 0)[None, :])
    # Spalte k: Σ_j chat_j Pbig[:, k+j]
    cols = np.zeros((G, len(k)), complex)
    for i, kv in enumerate(k):
        cols[:, i] = Pbig[:, (kv + jj) - kk[0]] @ chat
    vals = (chi_left * w)[:, None] * cols
    return E.conj().T @ (vals * wq / (2 * np.pi))


t0 = time.time()
# volle Form (χ ≡ 1)
Bfull = E.conj().T @ ((w[:, None] * Pp) * wq / (2 * np.pi))
I = np.eye(len(k))
lam_full = np.linalg.eigvalsh(I + 0.5 * (Bfull + Bfull.conj().T))[0]
print(f"η₁={eta1}, η₂={eta2}, K={K}: λ_min Re(I+𝒦) gesamt = {lam_full:.6f} ({time.time()-t0:.0f}s)", flush=True)
Bs = {}
for name, ch in (("L", cL), ("L̄", cLq), ("m", cm)):
    B = form(ch, ch)
    Bs[name] = B
    lam = np.linalg.eigvalsh(I + 0.5 * (B + B.conj().T))[0]
    print(f"  Stück {name}: λ_min Re(I + χ𝒦χ) = {lam:.6f}", flush=True)
R = Bfull - sum(Bs.values())
Rh = 0.5 * (R + R.conj().T)
print(f"  Kommutatorterm: ‖Re(𝒦 − Σχ𝒦χ)‖ = {np.linalg.norm(Rh, 2):.4f}, λ_min = {np.linalg.eigvalsh(Rh)[0]:.4f} ({time.time()-t0:.0f}s)",
      flush=True)
