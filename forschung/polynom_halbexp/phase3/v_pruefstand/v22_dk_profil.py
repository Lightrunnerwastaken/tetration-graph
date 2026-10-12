"""Strang V — Profil d_k = |p_N,k − a_k|ρ^k des echten p_N − f_K (p_N: gemeinsam/pN_N100_bis_1000.json, a_k: kneser_halbexp_koeff_4000.json).
Ausgabe je N: Summe und Maximum in Gradbändern ξ = k/N, oberste Koeffizienten. Zweck: Bedingung (b) für gewichtete Normen (Auftrag Leitung 22:52)."""
import json
import sys
from pathlib import Path

from mpmath import mp, mpf, lambertw, fabs

sys.stdout.reconfigure(encoding="utf-8")
mp.dps = 80
GEM = Path(__file__).resolve().parents[2] / "gemeinsam"
PN = json.loads((GEM / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"]
AK = json.loads((GEM / "kneser_halbexp_koeff_4000.json").read_text(encoding="utf-8"))["a"]
rho = fabs(lambertw(-1))
erg = {}
for Ns in sorted(PN, key=int):
    N = int(Ns)
    d = [float(fabs(mpf(PN[Ns][k]) - mpf(AK[k])) * rho ** k) for k in range(N + 1)]
    erg[N] = d
    bands = [(0, 0.5), (0.5, 0.8), (0.8, 0.9), (0.9, 0.95), (0.95, 0.99), (0.99, 1.01)]
    txt = "  ".join(f"[{a:.2f},{b:.2f}): Σ {sum(d[int(a * N):min(N + 1, int(b * N + 0.999))]):.2e}" for a, b in bands)
    print(f"N = {N:4d}: Σ = {sum(d):.3e}; {txt}; d_N = {d[N]:.2e}, d_(N−1) = {d[N - 1]:.2e}, d_(N−5) = {d[N - 5]:.2e}", flush=True)
(Path(__file__).resolve().parent / "v22_dk_profil.json").write_text(json.dumps({str(k): v for k, v in erg.items()}), encoding="utf-8")
