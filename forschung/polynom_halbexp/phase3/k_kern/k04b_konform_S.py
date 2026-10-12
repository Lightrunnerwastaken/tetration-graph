"""k04 — Test von Lemma K1′: sind die Abschnitte in der KONFORMEN Gleichungskoordinate von D̄ ∪ Linse polynomial stabil?

Ĵ_N := C_k^{(N)} J̃_N,  C_k^{(N)}[i][j] = [w^{i−j}] E(w)^j,  k(w) = ρ w E(w),  E = e^{F},  k: 𝔻 → Int γ konform (k(0) = 0),
γ = Γout ∪ Linse über Γin (Polarradius ρ·e^{s(α)}, s = 0 auf Γout).  Gemessen: ‖Ĵ_N^{-1}‖ = 1/σ_min(Ĵ_N)  (ℓ²_ρ(x) ← ℓ²(w)),
zum Vergleich ‖J̃_N^{-1}‖ (ℓ²_ρ) und die Lage des schlechtesten Vektors.

Konstruktion (deterministisch, ohne Datenfit):
  1. Linsenprofil s(τ) = S·(1 − exp(−(τ'/τ₁)^{p})), τ' = min(τ, τ_end − τ) Winkel ab arg L bzw. arg L̄ auf Γin, S = 0.0166, τ₁ = 12°, p = 1.3
     (liegt über Ws Profil Λ/0.7 ausser knapp bei 4°; in L tangential mit Exponent 1.3 > 1).
  2. Theodorsen: α(θ) = θ + H[s(α(θ))] per FFT (M Punkte), F(w) = Σ F_n w^n mit F_0 = Mittel, F_n = 2·ĝ_n, g = s∘α.
  3. F als EXAKTES Polynom vom Grad N (Gleitkomma-Koeffizienten als exakte Zahlen) ⇒ k exakt analytisch; E = e^F und E^j in arb.
     (Nicht strukturierte Rundungsfehler in C würden mit ‖J̃‖ ~ e^{0.47N} verstärkt — darum alles in arb.)
Aufruf: python k04_konform.py kubisch 200,400,600 [bits_faktor]
"""
import math
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from flint import arb, arb_mat, arb_poly, ctx  # noqa: E402
from k_kern import J_tilde, Lc, RHO  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
argL = math.atan2(Lc.imag, Lc.real)
TAU_END = 2 * math.pi - 2 * argL          # Winkellänge von Γin
S_GLOBAL = 0.0166
TAU1_GLOBAL = 12.0
P_GLOBAL = 1.3


def s_profil(alpha, S=None, tau1=None, p=None):
    S = S_GLOBAL if S is None else S
    tau1 = math.radians(TAU1_GLOBAL) if tau1 is None else tau1
    p = P_GLOBAL if p is None else p
    a = np.mod(alpha, 2 * math.pi)
    a = np.where(a > math.pi, 2 * math.pi - a, a)          # symmetrisch: |Winkel| in [0, π]
    tau = a - argL                                          # ab arg L Richtung π
    tp = np.minimum(tau, 2 * math.pi - 2 * argL - tau)      # Abstand zur nächsten Ecke (L oder L̄) — bei a ≤ π genügt tau
    tp = np.where(tau > 0, tau, 0.0)
    return np.where(tau > 0, S * (1 - np.exp(-(tp / tau1) ** p)), 0.0)


def theodorsen(M=2 ** 15, it=60, **kw):
    th = 2 * math.pi * np.arange(M) / M
    n = np.fft.fftfreq(M, d=1.0 / M)
    alpha = th.copy()
    for _ in range(it):
        g = s_profil(alpha, **kw)
        gh = np.fft.fft(g)
        Hg = np.real(np.fft.ifft(-1j * np.sign(n) * gh))
        alpha_neu = th + Hg
        d = np.max(np.abs(alpha_neu - alpha))
        alpha = alpha_neu
        if d < 1e-15:
            break
    g = s_profil(alpha, **kw)
    gh = np.fft.fft(g) / M
    F = np.zeros(M // 2)
    F[0] = gh[0].real
    F[1:] = 2 * gh[1:M // 2].real          # reell wegen Symmetrie
    # Kontrolle: Randkurve
    w = np.exp(1j * th)
    return F, d


def C_matrix(F, N):
    """C[i][j] = [w^{i−j}] E^j, E = exp(Σ_{n≤N} F_n w^n) — exakt in arb."""
    Fk = [arb(float(F[n])) if n < len(F) else arb(0) for n in range(N + 1)]
    E = [arb(0)] * (N + 1)
    E[0] = Fk[0].exp()
    for m in range(1, N + 1):
        acc = arb(0)
        for k in range(1, m + 1):
            acc += k * Fk[k] * E[m - k]
        E[m] = acc / m
    Ep = arb_poly(E)
    C = arb_mat(N + 1, N + 1)
    P = arb_poly([1])
    for j in range(N + 1):
        for i in range(j, N + 1):
            d = i - j
            C[i, j] = P[d] if d <= P.degree() else arb(0)
        P = (P * Ep)
        P = arb_poly([P[i] for i in range(min(P.degree(), N) + 1)])
    return C


def float_mat(X):
    n = X.nrows()
    return np.array([[float(X[i, j].mid()) for j in range(n)] for i in range(n)])


def main():
    name = sys.argv[1]
    Ns = [int(x) for x in sys.argv[2].split(",")]
    bf = float(sys.argv[3]) if len(sys.argv) > 3 else 1.2
    global S_GLOBAL, TAU1_GLOBAL, P_GLOBAL
    S_GLOBAL = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0166
    TAU1_GLOBAL = float(sys.argv[5]) if len(sys.argv) > 5 else 12.0
    P_GLOBAL = float(sys.argv[6]) if len(sys.argv) > 6 else 1.3
    print('τ₁ =', TAU1_GLOBAL, '°, p =', P_GLOBAL)
    print("Linsenhöhe S =", S_GLOBAL, " (Eckrate c ≈ 0.012, Schwelle c/0.95 ≈ 0.0126)")
    F, d = theodorsen()
    # Kontrolle der Abbildung
    th = np.linspace(0, 2 * math.pi, 4001)
    kw = RHO * np.exp(1j * th) * np.exp(np.polynomial.polynomial.polyval(np.exp(1j * th), F[:4096]))
    print("Theodorsen: letzte Änderung %.1e, F_0 = %.5f, max|F_n| (n ≥ 2000) = %.1e" % (d, F[0], np.max(np.abs(F[2000:4096]))))
    print("Randkurve k(∂𝔻): min |k|/ρ = %.6f, max |k|/ρ = %.6f, Abstand zu L = %.2e" %
          (np.min(np.abs(kw)) / RHO, np.max(np.abs(kw)) / RHO, np.min(np.abs(kw - Lc))))
    for N in Ns:
        t0 = time.time()
        bits = int(bf * N) + 300
        J = J_tilde(name, N, bits=bits)
        ctx.prec = bits
        C = C_matrix(F, N)
        Jh = C * J
        t1 = time.time()
        I = arb_mat(N + 1, N + 1)
        for i in range(N + 1):
            I[i, i] = arb(1)
        X = Jh.solve(I, algorithm="approx")
        Y = J.solve(I, algorithm="approx")
        t2 = time.time()
        Xf, Yf = float_mat(X), float_mat(Y)
        U, S, Vt = np.linalg.svd(Xf)
        UJ, sJ, VJt = np.linalg.svd(Yf)
        # Eck-Anteil: ‖Ĵ^{-T} u_i‖ = ‖C^{-T} J̃^{-T} u_i‖ für die grössten Moden von J̃^{-1}
        eck = [np.linalg.norm(Xf.T @ UJ[:, i]) for i in range(3)]
        u1 = U[:, 0]                  # Ausgangsvektor (Unbekannte, ℓ²_ρ) der grössten Verstärkung von Ĵ^{-1}
        kk = np.arange(N + 1)
        schwer = np.sum(kk * u1 ** 2) / np.sum(u1 ** 2)
        # Winkel, an dem der Ausgangsvektor als Polynom auf |x| = ρ am grössten ist
        phi = np.linspace(0, math.pi, 2001)
        vals = np.abs(np.polynomial.polynomial.polyval(np.exp(1j * phi), u1))
        print("N = %4d  ‖Ĵ⁻¹‖ = %.4g  (2. %.4g, 3. %.4g)   ‖J̃⁻¹‖ = %.4g   ‖Ĵ⁻¹‖/√N = %.4f   Ausgangsvektor: Index-Schwerpunkt %.2f N, "
              "Maximum bei %.1f°  [%.0f s + %.0f s, %d Bit]\n"
              "          Eckmoden von J̃⁻¹: σ = %.3g, %.3g, %.3g  →  ‖Ĵ^{-T}u_i‖ = %.3g, %.3g, %.3g"
              % (N, S[0], S[1], S[2], sJ[0], S[0] / math.sqrt(N), schwer / N, math.degrees(phi[np.argmax(vals)]),
                 t1 - t0, t2 - t1, bits, sJ[0], sJ[1], sJ[2], eck[0], eck[1], eck[2]), flush=True)


if __name__ == "__main__":
    main()
