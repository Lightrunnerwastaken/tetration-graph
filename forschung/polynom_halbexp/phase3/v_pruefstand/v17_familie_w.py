"""Strang V — Plan A3: Familie in w-Form. (S_γ) für die Sekantenoperatoren 𝒥_P v = v∘P + M_P v, P = f + η, kubisches Modell.

Für kubisches f ist M_P = (f∘P − f∘f)/(P − f) = c₁ + c₂(P + f) + c₃(P² + Pf + f²) EXAKT ein Polynom in P, f (keine Reihendivision, P-Warnung 9.2(6)).
Ĵ^P_N := T_N^w[(𝒥_Px)∘k] mit derselben konformen Abbildung k wie v12 (Linse S, Ks Profilform). Gemessen ‖(Ĵ^P_N)⁻¹‖ gegen ε für
  η = ε·g(x/ρ),  g ∈ {„glatt“: w², „zufaellig“: Σ_{k≤N} r_k w^k/√(N+1) mit r_k ~ N(0,1) (fester Seed), „hoch“: dieselbe Summe nur k ≥ 0.9N,
       „Lhoch“: Σ_{k≥0.9N} cos(kθ_L)w^k/√(0.1N), „Lalle“: Σ_{k≥1} cos(kθ_L)k^{−1.5}w^k/Σk^{−1.5} (auf L gerichtet wie p_N − f_K)}.
Frage (Ks K1″(a) für 𝒥_P): bleibt ‖(Ĵ^P_N)⁻¹‖ für ε ~ R_N (≈ 1e-10…1e-13) gleichmässig beschränkt, und wie hängt die Störung von ε und N ab?
Aufruf: python v17_familie_w.py N1,N2 S eps1,eps2,... art
"""
import math
import os
import sys
import time

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

ARGS = sys.argv[1:]
Ns = [int(x) for x in ARGS[0].split(",")]
S = ARGS[1]
EPS = [float(x) for x in ARGS[2].split(",")]
ART = ARGS[3]
sys.argv = ["v12", "kubisch", "100", S]
import v12_konform_S as V12  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")


def J_hat_P(N, F, bits, eps, art):
    ctx.prec = bits
    L = V12.L_arb()
    rho = abs(L)
    Fr = [arb(float(F[k].real)) for k in range(N + 1)]
    E = [arb(0)] * (N + 1)
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
    a, _ = V12.modell_koeff(N)
    c = a[:4]
    # f∘k und η∘k (beide als Reihe in w, nicht durch ρ geteilt)
    fk = arb_poly([0])
    rn = arb(1)
    for n in range(4):
        fk += (c[n] * rn) * Q[n]
        rn *= rho
    fk = fk.truncate(N + 1)
    rng = np.random.default_rng(2026)
    if art == "glatt":
        g = {2: 1.0}
    elif art.startswith("band"):
        # band<A>_<B>: Zufallskoeffizienten (fester Seed) nur für A·N ≤ k < B·N, auf sup_{|w|=1}|g| = 1 normiert
        A_, B_ = (float(x) for x in art[4:].split("_"))
        r = rng.standard_normal(N + 1)
        g = {k: float(r[k]) for k in range(N + 1) if A_ * N <= k < B_ * N}
        cc = np.zeros(N + 1)
        for k, v in g.items():
            cc[k] = v
        sup = np.max(np.abs(np.fft.fft(cc, 16 * (N + 1))))
        g = {k: v / sup for k, v in g.items()}
    elif art in ("Lhoch", "Lalle"):
        # Störung in Richtung L (wie p_N − f_K: Koeffizienten ~ cos(kθ_L), auf D̄ bei x = L, L̄ am grössten, dort ist die Linse 0)
        thL = math.atan2(1.3372357014306895, 0.31813150520476413)
        if art == "Lhoch":
            ks = range(int(0.9 * N), N + 1)
            g = {k: math.cos(k * thL) / math.sqrt(len(ks)) for k in ks}
        else:
            nrm = sum(k ** -1.5 for k in range(1, N + 1))
            g = {k: math.cos(k * thL) * k ** -1.5 / nrm for k in range(1, N + 1)}
        cc = np.zeros(N + 1)                       # auf sup_{|w|=1}|g| = 1 normieren (Gitter 16N)
        for k, v in g.items():
            cc[k] = v
        sup = np.max(np.abs(np.fft.fft(cc, 16 * (N + 1))))
        g = {k: v / sup for k, v in g.items()}
    else:
        r = rng.standard_normal(N + 1) / math.sqrt(N + 1)
        g = {k: float(r[k]) for k in range(N + 1) if (art == "zufaellig" or k >= int(0.9 * N))}
    ek = arb_poly([0])
    for k, gk in g.items():                       # η(x) = ε·g(x/ρ) ⇒ η∘k = ε Σ g_k (k/ρ)^k = ε Σ g_k Q_k
        ek += arb(gk) * Q[k]
    ek = (ek * arb(eps)).truncate(N + 1)
    Pk = (fk + ek)
    teil = os.environ.get("V17_TEIL", "")             # "K": nur Komposition stören, "M": nur Multiplikator stören (Diagnose)
    Pk_M = fk if teil == "K" else Pk
    if teil == "M":
        Pk, Pk_M = fk, fk + ek
    MP = (arb_poly([c[1]]) + c[2] * (Pk_M + fk) + c[3] * ((Pk_M * Pk_M).truncate(N + 1) + (Pk_M * fk).truncate(N + 1) + (fk * fk).truncate(N + 1))).truncate(N + 1)
    Pr = (Pk * (1 / rho)).truncate(N + 1)
    J = arb_mat(N + 1, N + 1)
    R = arb_poly([1])
    for j in range(N + 1):
        if j > 0:
            R = (R * Pr).truncate(N + 1)
        MQ = (MP * Q[j]).truncate(N + 1)
        Rc, Mc = R.coeffs(), MQ.coeffs()
        for i in range(N + 1):
            J[i, j] = (Rc[i] if i < len(Rc) else 0) + (Mc[i] if i < len(Mc) else 0)
    return J


def main():
    F, it, d = V12.theodorsen()
    print(f"Linse S = {S}, Störung {ART}", flush=True)
    for N in Ns:
        for eps in [0.0] + EPS:
            t0 = time.time()
            bits = int(1.2 * N) + 400
            J = J_hat_P(N, F, bits, eps, ART)
            sv, res = V12.norm_inverse(J)
            print(f"N = {N:5d}  ε = {eps:8.1e}  ‖(Ĵ^P)⁻¹‖ = {sv[0]:.6f}  (weitere {sv[1]:.4f})  Residuum {res:.1e}  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
