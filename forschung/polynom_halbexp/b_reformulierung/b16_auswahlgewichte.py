"""b16: Auswahlgewichte. Für k > N sei Λ_N(k) := ℓ(Π_N z^k) = ℓ(J_N⁻¹ T_N[f_K^k]), ℓ = Modus-(+1)-Koeffizient von θ_P im Korridor.
Dann (linear, exakt): ℓ(Π_N u − u) = Σ_{k>N} u_k Λ_N(k) − ℓ(Schwanz u)  ⇒  c_N ≈ Σ_{k>N} a_k Λ_N(k).
Frage: Wie sieht Λ_N(k)·L^k (bzw. ·conj) als Funktion von k aus? Lokal bei k ≈ N oder auf Skala N verteilt?"""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
Kc = [mp.mpf(c) for c in K]
def ell(coef, n=1):
    th = [horner(coef, z) * gw for (w, z), gw in zip(pts, gewicht)]
    e = [mp.exp(-2j * mp.pi * n * w) for w, _ in pts]
    return mp.fsum(t * x for t, x in zip(th, e)) / M
res = {}
for N in (40, 80, 120):
    KM = 4 * N + 200
    pw, g, Ji = baue(N, J=KM)
    Lam = {}
    for k in range(N + 1, KM + 1):
        x = Ji * arb_mat([[pw[k][i]] for i in range(N + 1)])
        Lam[k] = ell([mp.mpf(x[i, 0].mid().str(50, radius=False)) for i in range(N + 1)])
    cN = mp.fsum(Kc[k] * Lam[k] for k in Lam)
    print(f"\nN={N}: Σ_k a_k Λ_N(k) = {mp.nstr(cN, 6)}  (Vergleich b10: c_N)")
    for k in [N + 1, N + 2, N + 5, N + 10, N + 20, int(1.5 * N), 2 * N, 3 * N, 4 * N]:
        v = Lam[k] * mp.conj(L) ** k      # Gewicht relativ zu |L|^-k-Koeffizienten bei L̄ (Modus +1 ↔ L̄)
        print(f"   k={k}: Λ_N(k)·L̄^k = {mp.nstr(v, 5)}   |·| = {mp.nstr(abs(v), 4)}   |Λ·L^k| = {mp.nstr(abs(Lam[k] * L ** k), 4)}")
    res[N] = {k: (mp.nstr(mp.re(v), 20), mp.nstr(mp.im(v), 20)) for k, v in Lam.items()}
json.dump(res, open("b16_gewichte.json", "w"))
