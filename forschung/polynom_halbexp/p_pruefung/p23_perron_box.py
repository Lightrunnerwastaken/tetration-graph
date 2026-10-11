"""Strang P, Prüfung 23: Box-FORM für das Zweig-Zertifikat — Perron-Profil statt ρ^{-k}-Boxen (numerisch).

Im parametrischen Krawczyk-Test längs des Prädiktors z̃(s) = P₀ + sΔ₁ (Δ₁ = −c J(P₀)⁻¹x^M) lautet die Enthaltensbedingung
für den Term (I − Y J(Z))(Z − z̃) mit Z − z̃ = Box(r):  |S|·|A|·r ≤ θ·r,  A := J(P₀)⁻¹·∂_t J(P₀ + tΔ₁)|_{t=0}.
Mit ρ^{-k}-Boxen ist θ = ‖A‖_w (gewichtete Zeilensumme) — wächst ~1.6^M (Grund der vielen s-Teilschritte in p18).
Mit der Perron-Box (r = Perron-Vektor von |A|) ist θ = Spektralradius ϱ(|A|) — die bestmögliche Box-Form.
Ausgabe je M: ‖A‖_w, ϱ(|A|), |c|, ‖Δ₁‖_ρ, p18-Teilschritte zum Vergleich.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from flint import arb, arb_mat, arb_poly, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
RHO = 1.37455701
MS = [int(a) for a in sys.argv[1:]] or [20, 30, 40, 44, 50, 60]


def lade(N):
    for name in ["polynome_bis_194.json", "polynome_195_bis_300.json"]:
        d = json.load(open(GEM / name, encoding="utf-8"))
        if str(N) in d:
            return d[str(N)]
    raise KeyError(N)


def F_J(a, M):
    f = arb_poly(a)
    P = [arb_poly([1])]
    for _ in range(M):
        P.append((P[-1] * f).truncate(M + 1))
    g, fpf = arb_poly([0]), arb_poly([0])
    for m in range(M + 1):
        g += a[m] * P[m]
        if m:
            fpf += m * a[m] * P[m - 1]
    F = [g[k] - 1 / arb(math.factorial(k)) for k in range(M + 1)]
    J = arb_mat(M + 1, M + 1)
    for k in range(M + 1):
        for j in range(M + 1):
            J[k, j] = P[j][k] + (fpf[k - j] if j <= k else 0)
    return F, J


def main():
    raus = {}
    for M in MS:
        ctx.prec = 4 * M + 400
        a0 = [arb(x) for x in lade(M - 1)] + [arb(0)]
        F0, J0 = F_J(a0, M)
        c = F0[M]                                       # Kopf-Rest c_M (Rest ist ~ 1e-60 aus der Datei)
        eM = arb_mat([[1 if k == M else 0] for k in range(M + 1)])
        t = J0.solve(eM)
        D1 = [-c * t[k, 0] for k in range(M + 1)]
        h = arb(10) ** (-40)
        ap = [a0[k] + h * D1[k] for k in range(M + 1)]
        am = [a0[k] - h * D1[k] for k in range(M + 1)]
        _, Jp = F_J(ap, M)
        _, Jm = F_J(am, M)
        dJ = (Jp - Jm) * (1 / (2 * h))
        A = J0.solve(dJ)
        An = np.array([[abs(float(A[i, j])) for j in range(M + 1)] for i in range(M + 1)])
        w = np.array([RHO ** (-k) for k in range(M + 1)])
        norm_w = float(np.max((An @ w) / w))
        # Perron: Potenziteration auf |A|
        v = np.ones(M + 1)
        lam = 0.0
        for _ in range(3000):
            u = An @ v
            lam_neu = float(np.max(u / np.maximum(v, 1e-300)))
            v = u / np.max(u)
            if abs(lam_neu - lam) < 1e-12 * lam_neu:
                break
            lam = lam_neu
        # Collatz–Wielandt: obere Schranke max_k (|A|v)_k / v_k
        cw = float(np.max((An @ v) / np.maximum(v, 1e-300)))
        d1 = sum(abs(float(D1[k])) * RHO ** k for k in range(M + 1))
        prof = [float(v[k] / w[k]) for k in range(0, M + 1, max(1, M // 8))]
        print(f"M={M:3d}  |c|={abs(float(c)):.3e}  ‖Δ₁‖_ρ={d1:.3e}  ‖A‖_w={norm_w:.4g}  ϱ(|A|)≈{lam:.4g} (CW {cw:.4g})"
              f"  Verhältnis {norm_w / max(cw, 1e-300):.3g}")
        print(f"      Perron-Profil / ρ^-k (k = 0, M/8, …): " + " ".join(f"{x:.2g}" for x in prof))
        sys.stdout.flush()
        raus[M] = dict(c=abs(float(c)), d1=d1, norm_w=norm_w, perron=cw, profil=prof)
    (HIER / "p23_perron.json").write_text(json.dumps(raus, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
