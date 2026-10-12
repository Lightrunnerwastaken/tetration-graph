"""exp.polynom_phase3 — Kernzahlen aus Phase 3 des Forschungsteams (forschung/polynom_halbexp/phase3, 11.–12.10.2026) unabhängig nachgerechnet.

 Geprüft wird, was sich mit vertretbarem Aufwand selbst nachrechnen lässt:
   1. (S_γ): Schreibt man die abgeschnittenen Gleichungen in der konformen Koordinate w von D̄ ∪ Linse über Γin (k: 𝔻 → Int γ, k(0) = 0,
      Linsenprofil S(1 − e^{−(τ/12°)^{1.3}}), S = 0.0166), so bleibt ‖Ĵ_N⁻¹‖ konstant ≈ 1.2512 (kubisches Modell mit f(L) = L, f'(L) = √L);
   2. Kontrolle ohne Linse (S = 0): ‖J̃_N⁻¹‖ wächst (Eckmode), bei N = 200 ≈ 1.852;
   3. dieselbe Konstante mit Knesers zertifizierten Koeffizienten (≈ 1.2503);
   4. Tail-Schranke (K2) |a_k|ρ^k ≤ 461.36·k^{−6} ist mit den zertifizierten a_1 … a_640 verträglich;
   5. Integraldarstellung der Eckprofile: F_β(u) = Σ Γ(β − iαk)(−u)^k/(k!Γ(β)) = Γ(β)^{−1}∫₀^∞ t^{β−1}e^{−t − u t^{−iα}}dt, α = L/(4π).
 Konstruktion von Ĵ_N wie in phase3/v_pruefstand/v12_konform_S.py (Strang V), hier kompakt neu geschrieben.
"""
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from flint import acb, arb, arb_mat, arb_poly, ctx

from gemeinsam import Lauf, plt

lauf = Lauf("exp.polynom_phase3", __file__)
G = Path(__file__).resolve().parent.parent / "forschung" / "polynom_halbexp" / "gemeinsam"
ZERT = json.loads((G / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))["a"]

Lf = complex(0.31813150520476413, 1.3372357014306895)
RHO = abs(Lf)
TH = math.atan2(Lf.imag, Lf.real)
TAU_END = 2 * math.pi - 2 * TH
MFFT = 2 ** 16


def profil(phi, S):
    t = np.mod(phi, 2 * np.pi) - TH
    inn = (t > 0) & (t < TAU_END)
    tp = np.where(inn, np.minimum(t, TAU_END - t), 0.0)
    return np.where(inn, S * (1 - np.exp(-(tp / math.radians(12.0)) ** 1.3)), 0.0)


def theodorsen(S):
    """log k(w) = log ρ + log w + F(w); Randbedingung Re F = s(arg k) auf |w| = 1."""
    psi = 2 * np.pi * np.arange(MFFT) / MFFT
    n = np.fft.fftfreq(MFFT, d=1.0 / MFFT)
    theta = psi.copy()
    for _ in range(300):
        neu = psi + np.real(np.fft.ifft(-1j * np.sign(n) * np.fft.fft(profil(theta, S))))
        d = np.max(np.abs(neu - theta))
        theta = neu
        if d < 1e-15:
            break
    gh = np.fft.fft(profil(theta, S)) / MFFT
    F = np.zeros(MFFT // 2, dtype=complex)
    F[0] = gh[0].real
    F[1:] = 2 * gh[1:MFFT // 2]
    return F


def L_arb():
    return (-acb(-1).lambertw()).conjugate()


def koeff(modell, N):
    L = L_arb()
    if modell == "kubisch":
        mu = L.sqrt()
        A, b = arb_mat(4, 4), arb_mat(4, 1)
        pw, dpw = [acb(1), L, L * L, L * L * L], [acb(0), acb(1), 2 * L, 3 * L * L]
        for j in range(4):
            A[0, j], A[1, j], A[2, j], A[3, j] = pw[j].real, pw[j].imag, dpw[j].real, dpw[j].imag
        b[0, 0], b[1, 0], b[2, 0], b[3, 0] = L.real, L.imag, mu.real, mu.imag
        c = A.solve(b)
        a = ([c[j, 0] for j in range(4)] + [arb(0)] * N)[:N + 1]
        f = arb_poly(a[:4])
        fp, m, pwr = f.derivative(), arb_poly([0]), arb_poly([1])
        for k in range(fp.degree() + 1):
            m += fp[k] * pwr
            pwr = (pwr * f).truncate(N + 1)
        return a, [m[k] if k <= m.degree() else arb(0) for k in range(N + 1)]
    a = [arb(ZERT[k][0]) for k in range(N + 2)]                  # Kneser: m = f'∘f aus m·f' = eˣ
    fp = [a[k + 1] * (k + 1) for k in range(N + 1)]
    ex = [arb(1)]
    for k in range(1, N + 1):
        ex.append(ex[-1] / k)
    mk = [arb(0)] * (N + 1)
    for k in range(N + 1):
        s_ = ex[k]
        for j in range(k):
            s_ -= mk[j] * fp[k - j]
        mk[k] = s_ / fp[0]
    return a[:N + 1], mk


def norm_inv(modell, N, F):
    ctx.prec = int(1.2 * N) + 400
    rho = abs(L_arb())
    Fr = [arb(float(F[k].real)) for k in range(N + 1)]
    E = [Fr[0].exp()] + [arb(0)] * N
    for n in range(1, N + 1):
        E[n] = sum((k * Fr[k] * E[n - k] for k in range(1, n + 1)), arb(0)) / n
    wE = arb_poly([arb(0)] + E[:N])
    Q = [arb_poly([1])]
    for _ in range(N):
        Q.append((Q[-1] * wE).truncate(N + 1))
    a, mk = koeff(modell, N)
    fk, mkk, rn = arb_poly([0]), arb_poly([0]), arb(1)
    for n in range(N + 1):
        fk += (a[n] * rn / rho) * Q[n]
        mkk += (mk[n] * rn) * Q[n]
        rn *= rho
    fk, mkk = fk.truncate(N + 1), mkk.truncate(N + 1)
    J, R = arb_mat(N + 1, N + 1), arb_poly([1])
    for j in range(N + 1):
        if j:
            R = (R * fk).truncate(N + 1)
        Rc, Mc = R.coeffs(), (mkk * Q[j]).truncate(N + 1).coeffs()
        for i in range(N + 1):
            J[i, j] = (Rc[i] if i < len(Rc) else 0) + (Mc[i] if i < len(Mc) else 0)
    I = arb_mat(N + 1, N + 1)
    for i in range(N + 1):
        I[i, i] = 1
    X = J.solve(I, algorithm="approx")
    Xf = np.array([[float(X[i, j].mid()) for j in range(N + 1)] for i in range(N + 1)])
    return float(np.linalg.svd(Xf, compute_uv=False)[0])


F_linse, F_null = theodorsen(0.0166), theodorsen(0.0)
Ns = [100, 150, 200]
mit = {N: norm_inv("kubisch", N, F_linse) for N in Ns}
lauf.pruefe("(S_γ) kubisch: ‖Ĵ_N⁻¹‖ in der konformen Koordinate (Linse S = 0.0166) konstant 1.2512 ± 0.001 für N = 100, 150, 200",
            all(abs(v - 1.2512) < 1e-3 for v in mit.values()), {N: round(v, 6) for N, v in mit.items()})

ohne = {N: norm_inv("kubisch", N, F_null) for N in (100, 200)}
lauf.pruefe("Kontrolle ohne Linse: ‖J̃_N⁻¹‖ wächst (Eckmode), bei N = 200 ≈ 1.852",
            ohne[200] > ohne[100] and abs(ohne[200] - 1.852) < 0.01, {N: round(v, 4) for N, v in ohne.items()})

kn = norm_inv("kneser", 200, F_linse)
lauf.pruefe("(S_γ) mit Knesers zertifizierten Koeffizienten: ‖Ĵ_200⁻¹‖ ≈ 1.2503 (kubisch und Kneser gleich auf 1e-3)",
            abs(kn - 1.2503) < 1e-3 and abs(kn - mit[200]) < 2e-3, {"Kneser": round(kn, 6), "kubisch": round(mit[200], 6)})

mp.mp.dps = 60
rho = mp.mpf(RHO)
quot = [abs(mp.mpf(ZERT[k][0])) * rho ** k * k ** 6 for k in range(1, len(ZERT))]
lauf.pruefe("(K2): die zertifizierten a_1 … a_640 erfüllen |a_k|ρ^k·k⁶ ≤ 461.36", max(quot) <= mp.mpf("461.36"),
            {"max": mp.nstr(max(quot), 6), "bei_k": int(np.argmax([float(q) for q in quot])) + 1})

mp.mp.dps = 40
Lm = -mp.lambertw(-1)
Lm = mp.conj(Lm) if mp.im(Lm) < 0 else Lm
al = Lm / (4 * mp.pi)
abw = {}
for name, beta in (("β_A", mp.mpf("0.5") - 2j * al), ("β_P", mp.mpf("0.5") + 1j * al)):
    for u in (mp.mpf("0.5"), mp.mpf(2)):
        reihe = mp.nsum(lambda k: mp.gamma(beta - 1j * al * k) * (-u) ** k / mp.factorial(k), [0, mp.inf]) / mp.gamma(beta)
        integral = mp.quad(lambda t: t ** (beta - 1) * mp.exp(-t - u * t ** (-1j * al)),       # Stützstellen nahe 0: t^{β−1} ist dort singulär
                           [0, mp.mpf("1e-12"), mp.mpf("1e-6"), mp.mpf("1e-2"), 1, 10, mp.inf], maxdegree=10) / mp.gamma(beta)
        abw[f"{name}, u={mp.nstr(u, 2)}"] = abs(reihe - integral)
lauf.pruefe("Eckprofile: Gammareihe = Integraldarstellung (Abweichung < 1e-15)", max(abw.values()) < mp.mpf("1e-15"),
            {k: mp.nstr(v, 3) for k, v in abw.items()})

lauf.daten = {"S_gamma_kubisch": mit, "ohne_linse": ohne, "S_gamma_kneser_200": kn, "K2_max": mp.nstr(max(quot), 8)}
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(list(mit), list(mit.values()), "o-", label="mit Linse (konforme Koordinate)")
ax.plot(list(ohne), list(ohne.values()), "s--", label="ohne Linse")
ax.set_xlabel("Grad N"); ax.set_ylabel("Norm der Inversen"); ax.legend(fontsize=8)
ax.set_title("Phase 3: Abschnitte in der konformen Koordinate bleiben beschränkt")
lauf.bild(fig)
lauf.zusammenfassung = (f"‖Ĵ_N⁻¹‖ = {mit[200]:.6f} (kubisch, N = 100–200 konstant), Kneser {kn:.6f}; ohne Linse wächst sie "
                        f"({ohne[100]:.3f} → {ohne[200]:.3f}); (K2) mit den zertifizierten Koeffizienten verträglich; Eckprofil-Integral bestätigt.")
lauf.schreibe()
