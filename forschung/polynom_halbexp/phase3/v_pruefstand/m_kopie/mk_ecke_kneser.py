"""mk_ecke_kneser — Eckschranken für Kneser in Streifenkoordinaten (streng relativ zu (K1)+(K2) über S2 = sup_{D̄}|F''|).
Für t ≤ log q (q klein), beide Kanten, normierte Variable z = x/ρ, Ecke L := e^{iθ}, μ = √L = F'(L):
  s(q) := 2 sinθ·q/(1 − q) ≥ |z − L| (Sehnenlänge; Unterkante M = e^{iθ}q, Oberkante M = −e^{iθ}q).
  Unterkante (Γin): F(z) = L + μδ + R, |R| ≤ S2|δ|²/2 (Sehne ⊂ D̄)
     |G − t − ℓ| ≤ gG := −log(1 − S2 s/(2|μ|)) − log(1 − (|μ−1|s + S2 s²/2)/(2 sinθ − s)),
     W_in = 1/m = F'(z)e^{−ρz} (aus F∘F = e^{ρz}/ρ):  |W_in − 1/μ| ≤ gWin := s e^{−ρ(cosθ − s)}(S2 + ρ|μ|).
  Oberkante (Γout): h(y) als einziger Fixpunkt von Φ(v) = L + (y − L − R(v))/μ in der Kugel B = {|v − c| ≤ κs}, c = L + (y − L)/μ
     (B ⊂ D̄, Φ(B) ⊂ B, Kontraktion — alles geprüft), |v − L| ≤ s(1/|μ| + κ) =: dv;
     |H − (t + iπ − ℓ)| ≤ gH := −log(1 − S2 dv/(2|μ|)) − log(1 − (s + dv)/(2 sinθ − dv)),  |W_out − μ| ≤ gWout := S2 s.
  Zusätzlich für M25: Re H' ≥ 1 − eH' und r = |ξ'(t+iπ)|/|ξ'(H)| ≤ rmax.
"""
from flint import arb


def eckschranken(K, q, kappa=arb("0.3")):
    th, rho, mu, S2 = K.th, K.rho, K.mu, K.S2
    amu = abs(mu)
    b = mu.arg()
    s = 2 * th.sin() * q / (1 - q)
    # Unterkante
    u1 = S2 * s / (2 * amu)
    u2 = (abs(mu - 1) * s + S2 * s * s / 2) / (2 * th.sin() - s)
    assert bool(u1 < 1) and bool(u2 < 1)
    gG = -(1 - u1).log() - (1 - u2).log()
    gWin = s * (-rho * (th.cos() - s)).exp() * (S2 + rho * amu)
    # Oberkante: Fixpunktargument
    ah = 2 * (s / 2).asin()                      # Bogenwinkel α zur Sehne s (obere Schranke für Γout-Punkte mit |y−L| ≤ s)
    # 1 − |c| ≥ (s/|μ|)[sin(α/2 + b) − s/(2|μ|)] für alle Sehnen ≤ s: Minimum über α' ∈ [0, α] — sin(α'/2 + b) ≥ sin(b) für α' klein
    tiefe = (s / amu) * (b.sin() - s / (2 * amu))  # konservativ (sin(α'/2 + b) ≥ sin b für α' ≤ π − 2b)
    eps_r = kappa * s
    a1 = 1 / amu + kappa
    ok_B = bool(eps_r <= tiefe)
    ok_inv = bool(S2 * s * s * a1 * a1 / (2 * amu) <= eps_r)
    ok_kon = bool(S2 * s * a1 / amu < 1)
    assert ok_B and ok_inv and ok_kon, (ok_B, ok_inv, ok_kon)
    dv = s * a1
    w1 = S2 * dv / (2 * amu)
    w2 = (s + dv) / (2 * th.sin() - dv)
    assert bool(w1 < 1) and bool(w2 < 1)
    gH = -(1 - w1).log() - (1 - w2).log()
    gWout = S2 * s
    # H' = [y'/(y−L)]·[(y−L)/((v−L)F'(v))]·[1 − (v−L)/(v−L̄)]  mit y'/(y−L) = 1/(1−M)
    f1 = q / (1 - q)                               # |y'/(y−L) − 1|
    num = amu - S2 * dv / 2                        # |μ + R/δ_v| ≥ …
    den = amu + S2 * dv                            # |F'(v)| ≤ …
    f2lo = num / den
    f3 = dv / (2 * th.sin() - dv)
    reHp_lo = (1 - f1) * f2lo * (1 - f3) - 0       # untere Schranke |·|; Re ≥ |·|·cos(arg) — konservativ über Betrag/Arg:
    # Argumentschranke: arg jedes Faktors ≤ arcsin(Abweichung) ⇒ Re H' ≥ |H'|·cos(Σ arcsin)
    a_sum = (f1 / (1 - f1)).asin() if bool(f1 < arb("0.5")) else arb(1)
    a_sum += ((S2 * dv / 2 + S2 * dv) / (amu - S2 * dv)).asin()
    a_sum += f3.asin()
    reHp = reHp_lo * a_sum.cos()
    # r: e^{t − Re H}·|M_H − 1|²/|M_{t+iπ} − 1|², t − Re H ∈ Re ℓ ± gH, |M_H| ≤ q e^{Re ℓ + gH}
    qH = q * (mu.log().real + gH).exp()
    rmax = (mu.log().real + gH).exp() * (1 + qH) ** 2 / (1 - q) ** 2
    return dict(s=s, gG=gG, gWin=gWin, gH=gH, gWout=gWout, dv=dv, reHp=reHp, rmax=rmax, b=b, amu=amu)


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    from mk_kneser_arb import KneserArb
    K = KneserArb()
    for T in (3, 4, 5, 6, 8, 13):
        q = arb(-T).exp()
        e = eckschranken(K, q)
        print(f"T = {T}: s ≤ {float(e['s'].upper()):.3e}, |G−t−ℓ| ≤ {float(e['gG'].upper()):.3e}, |W_in−1/μ| ≤ {float(e['gWin'].upper()):.3e}, "
              f"|H−(t+iπ−ℓ)| ≤ {float(e['gH'].upper()):.3e}, |W_out−μ| ≤ {float(e['gWout'].upper()):.3e}, Re H' ≥ {float(e['reHp'].lower()):.5f}, "
              f"r ≤ {float(e['rmax'].upper()):.5f}")
