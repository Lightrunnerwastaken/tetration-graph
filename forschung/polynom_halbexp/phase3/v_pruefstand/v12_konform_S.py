"""Strang V — unabhängige Nachrechnung von (S_γ): sup_N ‖Ĵ_N⁻¹‖ < ∞ (Ks Lemma K1′/K1″), eigener Code.

Aussage (wie K, k_kern/BERICHT §5.1): k: 𝔻 → Int γ, k(w) = ρ·w·E(w), E = e^F, F analytisch, k(0) = 0; γ = Γout ∪ Linse über Γin
(Polarradius ρe^{s(φ)}). Unbekannte ξ_j = x_jρ^j (ℓ²), Gleichungen in w: Ĵ_Nξ = T_N^w[(J_∞x)∘k] (ℓ²). Gemessen ‖Ĵ_N⁻¹‖ = 1/σ_min(Ĵ_N).

Unabhängig von K gerechnet:
  * Modell kubisch selbst hergeleitet (reelles Kubik mit f(L) = L, f'(L) = √L; 4×4-System in arb); Kneser = Ps zertifizierte Mitten.
  * konforme Abbildung: eigene Theodorsen-Iteration (FFT M = 2^16), danach PRÜFUNG der Randbedingung |k| = R(arg k) an 8192 Punkten
    (float, ausserhalb von D bzw. auf ∂D — keine Auswertung im Innern) und Startwert-unabhängige Konvergenz.
  * Ĵ_N DIREKT als Spalten T_N^w[(f∘k/ρ)^j + (m∘k)·(wE)^j] (nicht als Produkt C_k·J̃_N wie K); alle Reihen mit arb_poly und expliziter
    Abschneidung (kein arb_series ⇒ keine ctx.cap-Falle); F_0…F_N als exakte Zahlen (float → arb exakt), nur F_n mit n ≤ N gehen ein.
  * Präzision 1.2N + 400 Bit (Einträge bis ~1.65^N); Inverse per arb_mat.solve(approx); Kontrolle: Radien, ‖Ĵ X − I‖ an Stichproben,
    Wiederholung mit +250 Bit (Option).
Aufruf: python v12_konform_S.py MODELL N1,N2,... [S] [tau1_grad] [p] [form] [bits_plus]
  MODELL kubisch | kneser_z;  form "K" (Ks Profil S(1 − e^{−(τ/τ₁)^p})) oder "sin" (S·sin(πτ'/τ_end)^p, anderes Profil)
  oder "kerbe" (Ks Profil mit Gauss-Kerbe auf s_min bei τ₀, Breite w; Umgebungsvariable V12_KERBE = "τ₀,w,s_min")
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from flint import acb, arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parents[1] / "gemeinsam"
MODELL = sys.argv[1]
Ns = [int(x) for x in sys.argv[2].split(",")]
S = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0166
TAU1 = math.radians(float(sys.argv[4])) if len(sys.argv) > 4 else math.radians(12.0)
PEXP = float(sys.argv[5]) if len(sys.argv) > 5 else 1.3
FORM = sys.argv[6] if len(sys.argv) > 6 else "K"
BITS_PLUS = int(sys.argv[7]) if len(sys.argv) > 7 else 0
import os
MFFT = int(os.environ.get("V12_MFFT", 2 ** 16))           # Abbildungsgenauigkeit variieren (Fallen-Test)
STOER = float(os.environ.get("V12_STOER", 0.0))          # zufällige Störung der F_n (Robustheit gegen Abbildungsfehler)


# ------------------------------------------------------------------ Geometrie und konforme Abbildung (float, nur auf dem Rand)
Lf = complex(0.31813150520476413, 1.3372357014306895)
rho_f = abs(Lf)
thL = math.atan2(Lf.imag, Lf.real)
tau_end = 2 * math.pi - 2 * thL          # Länge von Γin (von L über −ρ nach L̄)


def s_profil(phi):
    phi = np.mod(phi, 2 * np.pi)
    t = phi - thL                          # Winkel ab arg L in Γin-Richtung
    inn = (t > 0) & (t < tau_end)
    tp = np.where(inn, np.minimum(t, tau_end - t), 0.0)
    if FORM == "K":
        s = S * (1 - np.exp(-(tp / TAU1) ** PEXP))
    elif FORM == "doppel":
        # Ks Linse über Γin und EINBUCHTUNG über Γout: s = −A(1 − e^{−c·τ'/A}), τ' = Winkelabstand zu ±θ auf der Γout-Seite (V-Idee VE2)
        A_, c_ = (float(x) for x in os.environ.get("V12_EIN", "0.6,0.85").split(","))
        base = S * (1 - np.exp(-(tp / TAU1) ** PEXP))
        ph_ = np.mod(phi, 2 * np.pi)
        tout = np.where(ph_ < np.pi, thL - ph_, ph_ - (2 * np.pi - thL))
        aus = (ph_ < thL) | (ph_ > 2 * np.pi - thL)
        ein = -A_ * (1 - np.exp(-c_ * np.maximum(tout, 0) / A_))
        return np.where(inn, base, np.where(aus, ein, 0.0))
    elif FORM == "kerbe":
        # Ks Profil mit Kerbe: in |τ − τ₀| ≲ w auf s_min abgesenkt (beidseitig symmetrisch, τ = Abstand zu L bzw. L̄); V12_KERBE = "τ₀,w,s_min" (Grad, Grad, Höhe)
        t0_, w_, smin_ = (float(x) for x in os.environ.get("V12_KERBE", "15,4,0.003").split(","))
        base = S * (1 - np.exp(-(tp / TAU1) ** PEXP))
        g_ = np.exp(-((np.degrees(tp) - t0_) / w_) ** 2)
        s = base - (base - np.minimum(base, smin_)) * g_
    else:
        s = S * np.sin(np.pi * np.clip(t, 0, tau_end) / tau_end) ** PEXP
    return np.where(inn, s, 0.0)


def theodorsen():
    psi = 2 * np.pi * np.arange(MFFT) / MFFT
    n = np.fft.fftfreq(MFFT, d=1.0 / MFFT)
    theta = psi.copy()
    daempf = float(os.environ.get("V12_DAEMPF", "1.0"))          # 1.0 = ungedämpft (Standard); < 1 für stark nicht-kreisförmige Gebiete
    for it in range(int(os.environ.get("V12_ITMAX", "200"))):
        g = s_profil(theta)
        gh = np.fft.fft(g)
        conj = np.real(np.fft.ifft(-1j * np.sign(n) * gh))
        neu = psi + conj
        d = np.max(np.abs(neu - theta))
        theta = daempf * neu + (1 - daempf) * theta
        if d < 1e-15:
            break
    g = s_profil(theta)
    gh = np.fft.fft(g) / MFFT
    F = np.zeros(MFFT // 2, dtype=complex)
    F[0] = gh[0].real
    F[1:] = 2 * gh[1:MFFT // 2]
    return F, it, d


def pruefe_abbildung(F, npkt=8192):
    """Randbedingung log|k| − log ρ = s(arg k) an npkt Punkten von ∂𝔻 (Auswertung der Reihe per FFT, nur auf dem Rand)."""
    m = len(F)
    big = 4 * m
    coef = np.zeros(big, dtype=complex)
    coef[:m] = F
    vals = np.fft.ifft(coef) * big           # F(e^{iψ}) auf big Punkten
    psi = 2 * np.pi * np.arange(big) / big
    logk = np.log(rho_f) + vals + 1j * psi
    k = np.exp(logk)
    res = np.abs(np.real(logk) - np.log(rho_f) - s_profil(np.angle(k)))
    # Monotonie des Randwinkels (Injektivität auf dem Rand)
    arg = np.unwrap(np.angle(k))
    mono = np.min(np.diff(arg))
    return float(np.max(res)), float(np.min(np.abs(k)) / rho_f), float(np.max(np.abs(k)) / rho_f), float(mono), float(np.min(np.abs(k - Lf)))


# ------------------------------------------------------------------ Modelle (arb)
def L_arb():
    return (-acb(-1).lambertw()).conjugate()


def modell_koeff(N):
    """(a_0..a_{N}, m_0..m_{N}) als arb: f und m = f'∘f."""
    L = L_arb()
    if MODELL == "kubisch":
        mu = L.sqrt()
        # f(x) = c0 + c1 x + c2 x² + c3 x³ reell, f(L) = L, f'(L) = μ  ⇒  4 reelle Gleichungen
        A = arb_mat(4, 4)
        bvec = arb_mat(4, 1)
        pw = [acb(1), L, L * L, L * L * L]
        dpw = [acb(0), acb(1), 2 * L, 3 * L * L]
        for j in range(4):
            A[0, j] = pw[j].real
            A[1, j] = pw[j].imag
            A[2, j] = dpw[j].real
            A[3, j] = dpw[j].imag
        bvec[0, 0], bvec[1, 0], bvec[2, 0], bvec[3, 0] = L.real, L.imag, mu.real, mu.imag
        c = A.solve(bvec)
        a = [c[j, 0] for j in range(4)] + [arb(0)] * max(0, N - 3)
        a = a[:N + 1]
        f = arb_poly(a[:4])
        fp = f.derivative()
        m = arb_poly([0])
        pwr = arb_poly([1])
        for k in range(fp.degree() + 1):
            m += fp[k] * pwr
            pwr = (pwr * f).truncate(N + 1)
        mk = [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
        return a, mk
    if MODELL == "kneser_z":
        d = json.loads((GEM / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))
        assert N + 1 < len(d["a"]), "Kneser-Daten nur bis 640"
        a = [arb(d["a"][k][0]) for k in range(N + 2)]
        fp = [a[k + 1] * (k + 1) for k in range(N + 1)]
        ex = [arb(1)]
        for k in range(1, N + 1):
            ex.append(ex[-1] / k)
        mk = [arb(0)] * (N + 1)                  # m·f' = e^x
        for k in range(N + 1):
            s_ = ex[k]
            for j in range(k):
                s_ -= mk[j] * fp[k - j]
            mk[k] = s_ / fp[0]
        return a[:N + 1], mk
    raise ValueError(MODELL)


# ------------------------------------------------------------------ Ĵ_N direkt (reell: Linse symmetrisch ⇒ F reell)
def J_hat(N, F, bits):
    ctx.prec = bits
    L = L_arb()
    rho = abs(L)
    assert np.max(np.abs(F[:N + 1].imag)) < 1e-13, "F nicht reell"
    Fr = [arb(float(F[k].real)) for k in range(N + 1)]           # exakte Zahlen ⇒ k exakt analytisch, fest (unabhängig von N)
    E = [arb(0)] * (N + 1)                                     # E = exp(F): E_n = (1/n)Σ k F_k E_{n−k}
    E[0] = Fr[0].exp()
    for n in range(1, N + 1):
        s_ = arb(0)
        for k in range(1, n + 1):
            s_ += k * Fr[k] * E[n - k]
        E[n] = s_ / n
    wE = arb_poly([arb(0)] + E[:N])
    Q = [arb_poly([1])]
    for j in range(1, N + 1):
        Q.append((Q[-1] * wE).truncate(N + 1))
    a, mk = modell_koeff(N)
    fk, mkk = arb_poly([0]), arb_poly([0])                     # f∘k/ρ = Σ a_nρ^{n−1}Q_n,  m∘k = Σ m_nρ^nQ_n
    rn = arb(1)
    for n in range(N + 1):
        fk += (a[n] * rn / rho) * Q[n]
        mkk += (mk[n] * rn) * Q[n]
        rn *= rho
    fk = fk.truncate(N + 1)
    mkk = mkk.truncate(N + 1)
    J = arb_mat(N + 1, N + 1)
    R = arb_poly([1])
    for j in range(N + 1):
        if j > 0:
            R = (R * fk).truncate(N + 1)
        MQ = (mkk * Q[j]).truncate(N + 1)
        Rc, Mc = R.coeffs(), MQ.coeffs()
        for i in range(N + 1):
            J[i, j] = (Rc[i] if i < len(Rc) else 0) + (Mc[i] if i < len(Mc) else 0)
    return J


def norm_inverse(J):
    n = J.nrows()
    I = arb_mat(n, n)
    for i in range(n):
        I[i, i] = 1
    X = J.solve(I, algorithm="approx")
    Xf = np.array([[float(X[i, j].mid()) for j in range(n)] for i in range(n)])
    sv = np.linalg.svd(Xf, compute_uv=False)
    res = 0.0                                                   # Residuum ‖JX − I‖ an Stichproben (arb)
    for c in (0, n // 3, n - 1):
        for i in range(0, n, max(1, n // 40)):
            s_ = arb(0)
            for k in range(n):
                s_ += J[i, k] * X[k, c]
            res = max(res, abs(float((s_ - (1 if i == c else 0)).mid())))
    return sv, res


def main():
    t0 = time.time()
    F, it, d = theodorsen()
    if STOER > 0:
        rng = np.random.default_rng(12345)
        F = F.copy()
        F[1:] += STOER * rng.standard_normal(len(F) - 1) / np.arange(1, len(F)) ** 1.5     # glatte zufällige Störung, reell
        print(f"  Störung der Abbildung: F_n += {STOER}·N(0,1)/n^1.5 (n ≥ 1)", flush=True)
    res, kmin, kmax, mono, distL = pruefe_abbildung(F)
    print(f"[{MODELL}, M = {MFFT}, Linse {FORM} S = {S}, τ₁ = {math.degrees(TAU1):.1f}°, p = {PEXP}] Theodorsen: {it} Iterationen, letzte Änderung {d:.1e}; "
          f"F_0 = {F[0].real:.5f}, max|F_n| (n ≥ 2000) = {np.max(np.abs(F[2000:])):.1e}", flush=True)
    print(f"  Prüfung der Abbildung: max|log|k| − log ρ − s(arg k)| = {res:.2e}, min|k|/ρ = {kmin:.7f}, max|k|/ρ = {kmax:.6f}, "
          f"min Δarg = {mono:.2e} (> 0 ⇒ Randparametrisierung injektiv), Abstand zu L = {distL:.2e}", flush=True)
    for N in Ns:
        t1 = time.time()
        bits = int(1.2 * N) + 400 + BITS_PLUS
        J = J_hat(N, F, bits)
        t2 = time.time()
        sv, resid = norm_inverse(J)
        print(f"N = {N:5d}  ‖Ĵ⁻¹‖ = {sv[0]:.6f}  (weitere: {sv[1]:.4f}, {sv[2]:.4f})  Residuum ‖ĴX − I‖ (Stichprobe) = {resid:.1e}  "
              f"[{bits} Bit, Aufbau {t2 - t1:.0f} s, Inverse {time.time() - t2:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
