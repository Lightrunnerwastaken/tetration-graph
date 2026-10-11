"""Strang P, Prüfung 22 (Runde 3, Auftrag 2): Ist ein KOHÄRENTES Zweig-Zertifikat möglich? (numerisch, nicht rigoros)

Kohärente Norm: Unbekannte δ ∈ P_M mit ‖δ‖_G := sup_{G}|δ|, G = D̄ ∪ P(D̄) (P = p_M; = K' für P ≈ f_K),
Residuen mit sup auf D̄. Der kritische Faktor ist die Inversen-Norm
    γ_M := ‖J_M(P)⁻¹‖_{D̄ → G}  ≤  sqrt(M+1) · max_k ‖J_M⁻¹ x^k‖_G / ρ^k          (Parseval + Dreiecksungleichung)
und sup_{P(D̄)}|q| = sup_{∂D}|q∘P| (Maximumprinzip). Alle übrigen Grössen (Taylor-Abschneidung ~ log M, Markov/Cauchy auf G,
Ableitung längs des Pfades) sind polynomial in M. Wächst γ_M polynomial, trägt der kohärente Weg; wächst er wie 1.6^M, nicht.
Vergleich: dieselbe Grösse mit G = D̄ allein und die ℓ¹_ρ-Norm (C: |J⁻¹|_ρ ≈ 2.5).
Aufruf: python p22_kohaerent_test.py M1 M2 ...
"""
import json
import math
import sys
import time
from pathlib import Path

from flint import acb, acb_mat, arb, arb_mat, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
RHO = 1.37455701
MS = [int(a) for a in sys.argv[1:]] or [20, 40, 60, 80]


def lade(M):
    for name in ["polynome_bis_194.json", "polynome_195_bis_300.json"]:
        d = json.load(open(GEM / name, encoding="utf-8"))
        if str(M) in d:
            return d[str(M)]
    raise KeyError(M)


def jacobi(a, M):
    """J[k][j] = [x^k] P^j + [x^{k-j}](P'∘P), wie p18 (Punktwerte)."""
    from flint import arb_poly
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(M):
        P.append((P[-1] * f).truncate(M + 1))
    fpf = arb_poly([0])
    for m in range(1, M + 1):
        fpf += m * a[m] * P[m - 1]
    J = arb_mat(M + 1, M + 1)
    for k in range(M + 1):
        for j in range(M + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return J


def horner(a, z):
    s = acb(0)
    for c in reversed(a):
        s = s * z + c
    return s


def main():
    erg = {}
    for M in MS:
        t0 = time.time()
        ctx.prec = int(3.4 * (M * math.log10(1.7) + 45))
        a = [arb(x) for x in lade(M)]
        J = jacobi(a, M)
        eins = arb_mat([[1 if i == j else 0 for j in range(M + 1)] for i in range(M + 1)])
        Y = J.solve(eins)                       # J⁻¹ (Spalte k = J⁻¹ x^k)
        K = 6 * M + 60
        zs = [acb(arb(2 * i + 1) / K).exp_pi_i() * RHO for i in range(K)]
        ws = [horner(a, z) for z in zs]
        Z1 = acb_mat([[z ** j for j in range(M + 1)] for z in zs])
        Z2 = acb_mat([[w ** j for j in range(M + 1)] for w in ws])
        Ya = acb_mat(Y)
        V1, V2 = Z1 * Ya, Z2 * Ya
        g_D, g_P, g_l1, kmax = 0.0, 0.0, 0.0, None
        zeilen = []
        for k in range(M + 1):
            s1 = max(abs(complex(V1[i, k])) for i in range(K))
            s2 = max(abs(complex(V2[i, k])) for i in range(K))
            l1 = sum(abs(float(Y[j, k])) * RHO ** j for j in range(M + 1))
            sk = RHO ** k
            if max(s1, s2) / sk > max(g_D, g_P):
                kmax = k
            g_D, g_P, g_l1 = max(g_D, s1 / sk), max(g_P, s2 / sk), max(g_l1, l1 / sk)
            if k in (0, 1, 2, M // 4, M // 2, 3 * M // 4, M - 1, M):
                zeilen.append((k, s1 / sk, s2 / sk, l1 / sk))
        maxP = max(abs(complex(w)) for w in ws)
        print(f"M={M:4d} prec={ctx.prec} K={K} max|P| auf ∂D = {maxP:.4f}  ({time.time() - t0:.1f}s)")
        print(f"   max_k sup_∂D|J⁻¹x^k|/ρ^k = {g_D:.4g}   max_k sup_P(∂D)|J⁻¹x^k|/ρ^k = {g_P:.4g} (k* = {kmax})"
              f"   max_k ℓ¹_ρ/ρ^k = {g_l1:.4g}   γ_M ≤ {math.sqrt(M + 1) * max(g_D, g_P):.4g}")
        for k, s1, s2, l1 in zeilen:
            print(f"      k={k:4d}: D̄ {s1:.4g}   P(D̄) {s2:.4g}   ℓ¹ {l1:.4g}")
        erg[M] = dict(g_D=g_D, g_P=g_P, g_l1=g_l1, kmax=kmax, maxP=maxP, zeilen=zeilen)
        sys.stdout.flush()
    (HIER / "p22_kohaerent.json").write_text(json.dumps(erg, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
