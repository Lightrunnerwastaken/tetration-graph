"""O20 — Lemma O7a: Normen des lokalen Eckoperators T_L: E²(S) → H²(H), F ↦ F|_H + F(μ·)/μ.
ω(σ) = |1 + μ^{−3/2+iσ}|²·(1 + e^{−2πσ})/(1 + e^{−2(π+β)σ}), β = arg μ;  ‖T_L⁻¹‖ = (inf ω)^{−1/2}, ‖T_L‖ = (sup ω)^{1/2}."""
import sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237", "1.33723570143068940890116214319371061253950213846051241687")
for name, mu in (("√L (Kneser/kubisch)", mp.sqrt(L)), ("radial", mp.mpc(1, "0.6105986212909795")), ("Möbius-exp", mp.expj("0.7051039595"))):
    beta = mp.arg(mu)
    def om(s):
        z = mp.exp((-mp.mpf(3) / 2 + 1j * s) * mp.log(mu))
        w = (1 + mp.exp(-2 * mp.pi * s)) / (1 + mp.exp(-2 * (mp.pi + beta) * s))
        return abs(1 + z) ** 2 * w
    v = [om(x) for x in mp.linspace(-60, 60, 24001)]
    print(f"{name}: β = {mp.nstr(beta, 6)}, inf ω = {mp.nstr(min(v), 6)}, sup ω = {mp.nstr(max(v), 6)}, ‖T_L⁻¹‖ = {mp.nstr(1/mp.sqrt(min(v)), 6)}, "
          f"‖T_L‖ = {mp.nstr(mp.sqrt(max(v)), 6)}, |μ|^−3 = {mp.nstr(abs(mu)**-3, 6)}")
