"""Strang V — Familie in w-Form, Neumann-Grösse: q_N(P) := ‖Ĵ_N⁻¹ (Ĵ^P_N − Ĵ_N)‖₂ für P = f + η, kubisches Modell (Aufbau wie v17).

Eine Familien-Erweiterung von (S_γ) über eine Neumann-Reihe braucht q_N(P) < 1 gleichmässig für alle P der Newton-Kugel. v17 hat nur σ_max gemessen
(dort ist der Effekt bei zufälliger Störung quadratisch in ε und wächst mit N — Hinweis auf Vergrösserung von η auf der Linse, |x| > ρ).
Hier direkt: D = Ĵ^P − Ĵ in arb (exakte Differenz), links q_L = ‖Ĵ⁻¹D‖₂ und rechts q_R = ‖DĴ⁻¹‖₂ (über Ĵᵀ Z = Dᵀ), je per solve(approx)
und float-SVD. (Ĵ^P)⁻¹ = Ĵ⁻¹(I + DĴ⁻¹)⁻¹ ⇒ q_R < 1 genügt für ‖(Ĵ^P)⁻¹‖ ≤ ‖Ĵ⁻¹‖/(1 − q_R). Ĵ hat riesige Einträge (|f∘k| > ρ auf Γout,
die Spalten (f∘k/ρ)^j wachsen wie (max|f∘k|/ρ)^j), daher kann q_L ≫ 1 sein, ohne dass (Ĵ^P)⁻¹ sich ändert. Dazu ‖D‖₂ und
die wahre Änderung der Inversen Δ_inv = ‖(Ĵ^P)⁻¹ − Ĵ⁻¹‖₂ sowie die Glieder der Reihe (Ĵ^P)⁻¹ = Σ_k (−X)^k Ĵ⁻¹, X = Ĵ⁻¹D:
t_k = ‖X^k Ĵ⁻¹‖₂ für k = 1, 2, 3 (Produkte in arb, erst dann float — float-Produkte würden auslöschen), und
das Wachstum von η∘k auf der Linse: max_{|w|=1}|η∘k(w)| / ε (gemessen an 4096 Punkten aus den Reihenkoeffizienten).
Aufruf: python v18_familie_neumann.py N1,N2,... S eps1,eps2,... art
"""
import math
import sys
import time

import numpy as np
from flint import arb_mat, ctx

ARGS = sys.argv[1:]
Ns = [int(x) for x in ARGS[0].split(",")]
S = ARGS[1]
EPS = [float(x) for x in ARGS[2].split(",")]
ART = ARGS[3]
sys.argv = ["v17", "400", S, "1e-10", ART]
import v17_familie_w as V17  # noqa: E402
import v12_konform_S as V12  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
import os  # noqa: E402
KURZ = os.environ.get("V18_KURZ", "") == "1"          # nur Δ_inv (schnell)


def zu_np(M):
    n, m = M.nrows(), M.ncols()
    A = np.empty((n, m))
    for i in range(n):
        for j in range(m):
            A[i, j] = float(M[i, j].mid())
    return A


def main():
    F, it, d = V12.theodorsen()
    print(f"Linse S = {S}, Störung {ART}", flush=True)
    for N in Ns:
        bits = int(1.2 * N) + 400
        ctx.prec = bits
        J0 = V17.J_hat_P(N, F, bits, 0.0, ART)
        Id = arb_mat(N + 1, N + 1)
        for i in range(N + 1):
            Id[i, i] = 1
        J0i = J0.solve(Id, algorithm="approx")
        J0in = zu_np(J0i)
        c1 = zu_np(J0)[:, 1]                       # Spalte 1 = T_N[f∘k/ρ + (m∘k)·(k/ρ)] — nur zur Grössenangabe
        ww = np.exp(2j * np.pi * np.arange(4096) / 4096)
        print(f"N = {N}: max_|w|=1 |Ĵ e₁| = {np.max(np.abs(np.polynomial.polynomial.polyval(ww, c1))):.4f}, ‖Ĵ‖₂ ≈ {np.linalg.norm(zu_np(J0), 2):.3e}", flush=True)
        for eps in EPS:
            t0 = time.time()
            J1 = V17.J_hat_P(N, F, bits, eps, ART)
            ctx.prec = bits
            if KURZ:
                J1i = J1.solve(Id, algorithm="approx")
                dinv = np.linalg.norm(zu_np(J1i - J0i), 2)
                print(f"N = {N:5d}  ε = {eps:8.1e}  Δ_inv = {dinv:.4e} (/ε = {dinv / eps:.4e})  [{time.time() - t0:.0f} s]", flush=True)
                continue
            D = J1 - J0
            X = J0.solve(D, algorithm="approx")
            Z = J0.transpose().solve(D.transpose(), algorithm="approx")
            J1i = J1.solve(Id, algorithm="approx")
            dinv = np.linalg.norm(zu_np(J1i - J0i), 2)
            Xn, Dn, Zn = zu_np(X), zu_np(D), zu_np(Z)
            XJ = X * J0i
            X2J = X * XJ
            X3J = X * X2J
            t1, t2, t3 = (np.linalg.norm(zu_np(M_), 2) for M_ in (XJ, X2J, X3J))
            q = np.linalg.norm(Xn, 2)
            qr = np.linalg.norm(Zn, 2)
            dn = np.linalg.norm(Dn, 2)
            # Spalte 0 von D ist T_N[(M_P − m)∘k] (Q_0 = 1, (P/ρ)^0 = 1): enthält η∘k linear; Grösse von η∘k auf |w| = 1 über die Spalte j = 0
            col = Dn[:, 0]
            ww = np.exp(2j * np.pi * np.arange(4096) / 4096)
            mx = np.max(np.abs(np.polynomial.polynomial.polyval(ww, col)))
            print(f"N = {N:5d}  ε = {eps:8.1e}  Δ_inv = {dinv:.4e} (/ε = {dinv / eps:.4e})  t₁ = {t1:.4e}  t₂ = {t2:.4e}  t₃ = {t3:.4e}  q_R = ‖ΔĴĴ⁻¹‖ = {qr:.4e} (/ε = {qr / eps:.4e})  q_L = ‖Ĵ⁻¹ΔĴ‖ = {q:.4e}  ‖ΔĴ‖ = {dn:.4e}  "
                  f"max|ΔĴ e₀| auf |w|=1 = {mx:.3e} (/ε = {mx / eps:.3e})  [{time.time() - t0:.0f} s]", flush=True)


if __name__ == "__main__":
    main()
