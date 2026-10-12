"""Strang V — Ks Plan A3 (i)–(iii): Eck-Korrektur D_N = Ĵ_N⁻¹ − P_N Ĵ_{N'}⁻¹ P_N (Widom-Differenz, Definition wie k05), eigener Code (v12/v20).

(i)  Universalität: Singulärwerte (Paare L/L̄?), Lage der führenden Singulärvektoren im Index ξ = j/N und im Winkel, Modelle kubisch / Kneser.
(ii) Linsenform nahe L: dieselben Grössen für andere τ₁, p, S (Aufrufparameter wie v12).
(iii) „Kante fern von L trivial“: ‖T_x(χ) D_N T_w(χ_k)‖ mit glatter Winkel-Abschneidung χ = 0 in |θ ∓ θ_L| < δ, = 1 ab 2δ (cos²-Rampe);
     T_x(χ): Toeplitz-Matrix des Symbols χ(θ) auf den x-Koeffizienten (Kreis |x| = ρ), T_w(χ_k): Symbol χ(arg k(e^{iψ})) auf den w-Koeffizienten.
Zeilen von D_N = Unbekannte (x-Index j), Spalten = w-Gleichungen. N' = 3N (Ks Bedingung N' ≥ 3N).
Aufruf: python v26_widom.py MODELL N1,N2,... [S] [tau1_grad] [p]     MODELL kubisch | kneser (f_K aus kneser_halbexp_koeff_4000.json)
"""
import math
import sys
import time

import numpy as np
from flint import arb, arb_mat, ctx

ARGS = sys.argv[1:]
MODELL = ARGS[0]
Ns = [int(x) for x in ARGS[1].split(",")]
S = ARGS[2] if len(ARGS) > 2 else "0.0166"
TAU1 = ARGS[3] if len(ARGS) > 3 else "12"
PEXP = ARGS[4] if len(ARGS) > 4 else "1.3"
sys.argv = ["v20", "100", S]
import v20_pN_S as V20  # noqa: E402
V12 = V20.V12
V12.S = float(S)
V12.TAU1 = math.radians(float(TAU1))
V12.PEXP = float(PEXP)
V12.MODELL = "kubisch" if MODELL == "kubisch" else "kneser_z"

sys.stdout.reconfigure(encoding="utf-8")
thL = V12.thL


def chi(theta, delta):
    """glatte Abschneidung: 0 nahe ±θ_L (|Abstand| < δ), 1 ab 2δ."""
    out = np.ones_like(theta)
    for c in (thL, -thL):
        dd = np.abs(np.angle(np.exp(1j * (theta - c))))
        r = np.clip((dd - delta) / delta, 0, 1)
        out = out * np.where(dd < delta, 0.0, np.where(dd > 2 * delta, 1.0, np.sin(np.pi * r / 2) ** 2))
    return out


def toeplitz(symb, n):
    M = len(symb)
    c = np.fft.fft(symb) / M                      # c[m] = Fourier-Koeffizient m (m mod M)
    T = np.empty((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            T[i, j] = c[(i - j) % M]
    return T


def profil(v, n):
    """Masse |v_j|² im Indexband ξ ∈ (0.6, 1], Index der grössten Masse."""
    m = np.abs(v) ** 2
    m = m / m.sum()
    j = np.arange(n)
    return float(m[j > 0.6 * (n - 1)].sum()), int(np.argmax(m)) / (n - 1)


def winkel(vec, n, kpts=None):
    """Winkel (Grad) des Maximums von |Σ v_j z^j| auf dem Einheitskreis; für w-Vektoren über arg k(e^{iψ}) in x-Winkel umgerechnet."""
    M = 4096
    vals = np.fft.ifft(np.concatenate([vec, np.zeros(M - n)])) * M      # Σ v_j e^{ijψ}
    i = int(np.argmax(np.abs(vals)))
    ang = 2 * np.pi * i / M
    if kpts is not None:
        ang = kpts[i]
    return math.degrees(math.atan2(math.sin(ang), math.cos(ang)))


def main():
    F, it, d = V12.theodorsen()
    M = 4096
    psi = 2 * np.pi * np.arange(M) / M
    big = 4 * len(F)
    coef = np.zeros(big, dtype=complex)
    coef[:len(F)] = F
    vals = np.fft.ifft(coef) * big
    idx = (np.arange(M) * (big // M)).astype(int)
    argk = psi + np.imag(vals[idx])                 # arg k(e^{iψ}) = ψ + Im F(e^{iψ})
    print(f"[{MODELL}, Linse S = {S}, τ₁ = {TAU1}°, p = {PEXP}]", flush=True)
    for N in Ns:
        t0 = time.time()
        Nr = 3 * N
        bits = int(1.2 * Nr) + 400
        ctx.prec = bits
        if MODELL == "kneser":
            kk = V20.koeff_fK(Nr)
            V12.modell_koeff = lambda n_, kk=kk: (kk[0][:n_ + 1], kk[1][:n_ + 1])
        Jr = V12.J_hat(Nr, F, bits)
        ctx.prec = bits
        Xr = V20.zu_np(V20.inverse(Jr))
        n = N + 1
        Jn = arb_mat(n, n)
        for i in range(n):
            for j in range(n):
                Jn[i, j] = Jr[i, j]                   # geschachtelt: Ĵ_N = linker oberer Block von Ĵ_{N'}
        Xn = V20.zu_np(V20.inverse(Jn))
        D = Xn - Xr[:n, :n]
        U, sv, Vh = np.linalg.svd(D)
        pr_l, jm_l = profil(U[:, 0], n)
        pr_r, jm_r = profil(Vh[0, :], n)
        wl = winkel(U[:, 0], n)
        wr = winkel(Vh[0, :].conj(), n, argk)
        txt = []
        for dg in (10, 20, 40):
            delta = math.radians(dg)
            Tx = toeplitz(chi(psi, delta), n)
            Tw = toeplitz(chi(argk, delta), n)
            txt.append(f"δ = {dg}°: {np.linalg.norm(Tx @ D @ Tw, 2):.3e}")
        print(f"N = {N:4d} (N' = {Nr}): ‖Ĵ_N⁻¹‖ = {np.linalg.norm(Xn, 2):.5f}, ‖D_N‖ = {sv[0]:.4f}; Singulärwerte {sv[1]:.4f}, {sv[2]:.4f}, {sv[3]:.4f}, "
              f"{sv[4]:.4f}; führend links (x): Masse ξ>0.6 = {pr_l:.2f}, Spitze ξ = {jm_l:.2f}, Winkel {wl:+.1f}°; rechts (w): Masse {pr_r:.2f}, "
              f"Spitze {jm_r:.2f}, Winkel {wr:+.1f}° (arg L = {math.degrees(thL):.1f}°); fern von L komprimiert: {'; '.join(txt)}  [{time.time() - t0:.0f} s]",
              flush=True)


if __name__ == "__main__":
    main()
