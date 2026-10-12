"""Strang V — Reserve-Fenster, Teil (a): Untergrenze der Linse = n-te-Wurzel-Profil Λ(τ) der Eckmoden auf Γin, eigene Messung.

G_N(t) = ‖J̃_N⁻¹ẽ_t‖ (Os Majorante) auf dem Kreis |t| = ρ bei t = ρe^{i(θ_L+τ)}; Ansatz ln G_N = a(τ) + κ(τ)·ln N + Λ(τ)·N, angepasst über
mehrere N (alle Eckmoden zugleich, keine SVD). Λ(τ) ist die Wachstumsrate, die eine Linse der Höhe s(τ) mindestens kompensieren muss
(G(ρe^{s}e^{iτ}) ≈ e^{N(Λ(τ) − s)} wegen der Steg-Indizes k ≈ N; gemessen κ_S = 1.00 in v12).
J̃_N wie in v12 (k(w) = ρw, also E = 1): Spalten T_N[(f(ρw)/ρ)^j + m(ρw)w^j], arb_poly mit expliziter Abschneidung, 1.2N + 400 Bit;
Lösung J̃X = B (Spalten ẽ_t) per arb_mat.solve(approx); Kontrolle des Residuums an Stichproben.
Aufruf: python v14_lambda.py MODELL N1,N2,...
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from flint import acb, acb_mat, arb, arb_mat, ctx

ARGS = sys.argv[1:]
MODELL = ARGS[0]
Ns = [int(x) for x in ARGS[1].split(",")]
sys.argv = ["v12", MODELL, "100", "0.0"]
import v12_konform_S as V12  # noqa: E402  (nur Modell und Matrixaufbau)

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
L = complex(0.31813150520476413, 1.3372357014306895)
rho = abs(L)
thL = math.atan2(L.imag, L.real)
TAUS = [0.5, 1, 2, 3, 4, 6, 8, 11, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 103.38]


def main():
    erg = {"modell": MODELL, "taus": TAUS, "N": {}}
    for N in Ns:
        t0 = time.time()
        bits = int(1.2 * N) + 400
        F0 = np.zeros(N + 1, dtype=complex)
        J = V12.J_hat(N, F0, bits)                       # = J̃_N (ρ-skaliert)
        nt = len(TAUS)
        B = arb_mat(N + 1, 2 * nt)                       # J̃ reell ⇒ Real- und Imaginärteil der rechten Seiten getrennt lösen
        for c_, tau in enumerate(TAUS):
            ph = thL + math.radians(tau)
            q = acb(arb(math.cos(-ph)), arb(math.sin(-ph)))      # ρ/t = e^{−iφ} auf dem Kreis
            p = acb(1)
            for k in range(N + 1):
                B[k, c_] = p.real
                B[k, nt + c_] = p.imag
                p = p * q
        X = J.solve(B, algorithm="approx")
        G = []
        for c_ in range(nt):
            s_ = arb(0)
            for k in range(N + 1):
                s_ += X[k, c_] ** 2 + X[k, nt + c_] ** 2
            G.append(float(s_.sqrt().mid()))
        res = 0.0
        for i in range(0, N + 1, max(1, N // 30)):
            s_ = arb(0)
            for k in range(N + 1):
                s_ += J[i, k] * X[k, nt - 1]
            res = max(res, abs(float((s_ - B[i, nt - 1]).mid())))
        erg["N"][N] = G
        print(f"N = {N}: ln G = " + " ".join(f"{math.log(g):.4f}" for g in G) + f"   [Residuum {res:.1e}, {bits} Bit, {time.time() - t0:.0f} s]", flush=True)
        (HIER / f"v14_lambda_{MODELL}.json").write_text(json.dumps(erg, indent=0), encoding="utf-8")
    if len(Ns) >= 3:
        Nv = np.array(Ns, dtype=float)
        A = np.vstack([np.ones_like(Nv), np.log(Nv), Nv]).T
        print("\nAnpassung ln G = a + κ ln N + Λ N:")
        for c_, tau in enumerate(TAUS):
            y = np.array([math.log(erg["N"][N][c_]) for N in Ns])
            sol, *_ = np.linalg.lstsq(A, y, rcond=None)
            r2 = (y[-1] - y[-2]) / (Nv[-1] - Nv[-2])
            print(f"  τ = {tau:7.2f}°   Λ = {sol[2]:.5f}   κ = {sol[1]:+.3f}   (Rohrate letztes Paar {r2:.5f})", flush=True)


if __name__ == "__main__":
    main()
