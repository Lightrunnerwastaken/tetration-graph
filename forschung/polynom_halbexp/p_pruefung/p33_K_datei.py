"""Strang P, Prüfung 33 (Meilenstein 6): Datei für M/O — (K1) Kugeln b_k = a_kρ^k, k ≤ 640; (K2) |b_k| ≤ C·k^{−6} ∀k ≥ 1; (K3).

(K2): b_k = (1/2π)∫ f_K(ρe^{iθ})e^{−ikθ}dθ; f_K ∈ C⁵(∂D), ∂_θ⁵f absolut stetig (stetig in θ_L, glatt sonst, ∂_θ⁶f ∈ L¹),
sechsmal partiell integrieren: |b_k| ≤ k^{−6}·‖∂_θ⁶f‖_{L¹(0,2π)}/(2π) =: C·k^{−6}. ‖·‖_{L¹} = 2·(rechts + links + Reste) aus p32.
Restsummen: Σ_{k>K} k^j|b_k| ≤ C·K^{j−5}/(5−j) (Integralvergleich, j = 0, 1, 2).
Aufruf: python p33_K_datei.py
"""
import json
from pathlib import Path

from flint import acb, arb, ctx

HIER = Path(__file__).resolve().parent
ZIEL = HIER.parent / "gemeinsam" / "kneser_halbexp_b_zert.json"


def main():
    ctx.prec = 600
    rho = abs(acb(-1).lambertw())
    d = json.loads((HIER.parent / "gemeinsam" / "kneser_halbexp_koeff_zert_640.json").read_text(encoding="utf-8"))
    a = [arb(m) + arb(0, 1) * arb(r) for m, r, _ in d["a"]]
    K = len(a) - 1
    b = [a[k] * rho ** k for k in range(K + 1)]
    teile = {}
    for s in ("rechts", "links"):
        z = json.loads((HIER / f"p32_rand_{s}.json").read_text(encoding="utf-8"))
        teile[s] = (arb(z["integral_obere"]) * arb("1.000001"), arb(z["rest"]) * arb("1.000001"), z)
    L1 = 2 * sum((ii + rr for ii, rr, _ in teile.values()), arb(0))
    C = (L1 / (2 * arb.pi())).upper()
    C = arb(float(C) * (1 + 1e-12))                      # obere Schranke als Dezimalzahl
    C_str = f"{float(C.upper()):.6g}"
    C = arb(C_str)
    assert C >= L1 / (2 * arb.pi())
    p = 6
    Ka = arb(K)
    rest0 = C * Ka ** (1 - p) / (p - 1)                  # Σ_{k>K}|b_k|
    rest1 = C * Ka ** (2 - p) / (p - 2)                  # Σ_{k>K} k|b_k|
    rest2 = C * Ka ** (3 - p) / (p - 3)                  # Σ_{k>K} k²|b_k|
    s400_0 = sum((abs(b[k]) for k in range(401, K + 1)), arb(0)) + rest0
    s400_1 = sum((k * abs(b[k]) for k in range(401, K + 1)), arb(0)) + rest1
    rad400 = sum((b[k].rad() for k in range(0, 401)), arb(0))
    S1 = sum((k * abs(b[k]) for k in range(1, K + 1)), arb(0)) + rest1      # Σ k|b_k| (≥ ρ sup|f'|)
    S2 = sum((k * (k - 1) * abs(b[k]) for k in range(2, K + 1)), arb(0)) + rest2
    print(f"‖∂_θ⁶f‖_L¹ ≤ {L1.str(6)} ⇒ C = {C_str} (p = 6)")
    print(f"Σ_(k≤400) Radius = {rad400.str(3)};  Σ_(k>400)|b_k| ≤ {s400_0.str(3)};  Σ_(k>400) k|b_k| ≤ {s400_1.str(3)}")
    print(f"Σ_(k>640): |b_k| ≤ {rest0.str(3)}, k|b_k| ≤ {rest1.str(3)}, k²|b_k| ≤ {rest2.str(3)};  Σk|b_k| ≤ {S1.str(6)}, Σk(k−1)|b_k| ≤ {S2.str(6)}")
    assert rad400 < arb("1e-12") and s400_0 < arb("1e-12") and s400_1 < arb("1e-9")
    aus = {
        "beschreibung": "ZERTIFIZIERT (Strang P, Satz P-V + p32/p33): b_k = a_k ρ^k (ρ = |L|) von Knesers f_K. "
                        "(K1) Kugeln b_mitte ± b_radius für 0 ≤ k ≤ K. (K2) |b_k| ≤ C·k^(−p) für ALLE k ≥ 1 (insb. k > K): sechsmal partiell "
                        "integrieren auf |x| = ρ, C = ‖∂_θ⁶ f_K(ρe^{iθ})‖_L¹/(2π), zertifiziert in p32 (Zellen + exakte obere Struktur bei L, p32a). "
                        "(K3) f_K(L) = L, f_K'(L) = √L als Randwerte: aus sexp_K(z) = P(e^{L(z+θ_K(z))}) auf Im z ≥ 0.36 (p32a, ARCHITEKTUR §9.8) folgt "
                        "f_K(x) = P(√L·w·R(q)), w = P⁻¹(x), R(0) = 1, also f_K → L und f_K' → √L·R(0) = √L für x → L in D̄.",
        "K": K, "p": str(p), "C": C_str,
        "b_mitte": [x.mid().str(60, radius=False) for x in b],
        "b_radius": [x.rad().str(5, radius=False) for x in b],
        "restsummen": {"sum_b_k_ueber_400": s400_0.str(5), "sum_k_b_k_ueber_400": s400_1.str(5),
                       "sum_b_k_ueber_K": rest0.str(5), "sum_k_b_k_ueber_K": rest1.str(5), "sum_k2_b_k_ueber_K": rest2.str(5),
                       "radius_summe_bis_400": rad400.str(5), "S1_sum_k_b_k": S1.str(8), "S2_sum_k_km1_b_k": S2.str(8)},
        "K3": {"f_K(L)": "L", "f_K'(L)": "sqrt(L) = e^(L/2)", "fundstelle": "p_pruefung/p32a_theta_K.py (Struktursatz), "
               "p32_rand_ableitungen.py (Strukturformel, Werte nahe L), ARCHITEKTUR §9.8"},
        "L1_theta6": L1.str(8),
        "teilintegrale": {s: {"zellen": ii.str(8), "rest": rr.str(3), "anzahl_zellen": len(z["zellen"])} for s, (ii, rr, z) in teile.items()},
    }
    ZIEL.write_text(json.dumps(aus, indent=0, ensure_ascii=False), encoding="utf-8")
    print(f"geschrieben: {ZIEL}")


if __name__ == "__main__":
    main()
