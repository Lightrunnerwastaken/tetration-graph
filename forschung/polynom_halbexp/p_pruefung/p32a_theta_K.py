"""Strang P, Prüfung 32a (Meilenstein 6, Punkt 2): exakte obere Struktur von sexp_K — θ_K mit Fehlerschranke (rigoros, arb).

Behauptung (bewiesen durch diese Rechnung + Satz P-V):
  Es gibt eine 1-periodische, auf Im z ≥ y₁ holomorphe und für Im z → ∞ beschränkte Funktion θ_K mit
      sexp_K(z) = P(e^{L(z + θ_K(z))})   für alle Im z ≥ y₁ (alle Re z),
  θ_K(z) = Σ_{n≥0} κ_n qⁿ (q = e^{2πiz}) und |θ_K − θ̃| ≤ δθ auf Im z ≥ y₁, also |κ_n − c̃_n| ≤ δθ·e^{2πn y₁}.
Beweis:
 (a) Streifen y₁ ≤ Im z ≤ 3, |Re z| ≤ ½: Überdeckung durch Quadrate Z_b. Φ(ζ) := P(e^{Lζ}) (ganz), ζ̃(z) = z + θ̃(z), Ũ = Φ∘ζ̃.
     Für festes z ∈ Z_b: g(τ) = Φ(ζ̃(z) + τ) − sexp_K(z), |g(0)| = |Ũ(z) − sexp_K(z)| ≤ ε_U, |g'(0)| ≥ m₁, |g''| ≤ M₂ auf |τ| ≤ r_τ
     (m₁, M₂ aus einem Taylor-Modell von Φ um ζ_c = ζ̃(z_c), ζ̃(Z_b) ⊂ D(ζ_c, s₁)). Krawczyk mit Y = 1/g'(0) auf D(0, r_τ), r_τ = 1e-120:
     κ = M₂r_τ/m₁, ε_U/m₁ < (1−κ)r_τ ⇒ genau eine Nullstelle τ*(z), |τ*| ≤ ε_U/((1−κ)m₁); θ_K(z) := θ̃(z) + τ*(z), holomorph (IFT).
     Überlappungen: Eindeutigkeit in D(θ̃(z), r_τ) ⇒ gleich. Periodizität: θ = θ_K(z) löst auch die Gleichung bei z+1
     (P(Lw) = e^{P(w)}, sexp_K(z+1) = e^{sexp_K(z)}), Eindeutigkeit ⇒ θ_K(z+1) = θ_K(z).
 (b) Im z ≥ 3: |w̃| ≤ w_oben (p28), P univalent auf D(0, 0.1) (|P' − 1| ≤ κ' < 1), Rouché ⇒ W(z) := P⁻¹(sexp_K(z)) ∈ D(0, 0.1)
     eindeutig, W(z+1) = L·W(z); E := W e^{−Lz} ist 1-periodisch, |E| ≤ C|q|^{−0.213} ⇒ hebbar in q = 0;
     |E − Ẽ| ≤ ε_U/(1−κ')·e^{0.159 + 1.337·3} auf |q| = e^{−6π} ⇒ (Maximumprinzip) auf der ganzen Scheibe; inf|Ẽ| > 0 ⇒ θ_K := log E / L.
 (c) Identitätssatz: sexp_K(z) = P(e^{L(z+θ_K)}) auf |Re z| ≤ ½ ⇒ auf ganz Im z > y₁ (beide Seiten holomorph).
Aufruf: python p32a_theta_K.py [y1] [prec]
Ausgabe: p32a_theta_K.json (delta_theta, y1, kappa_max, Anzahl Quadrate).
"""
import json
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

import kc_rig
from kc_kneser import Kneser, amax

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
Y1 = sys.argv[1] if len(sys.argv) > 1 else "0.36"
PREC = int(sys.argv[2]) if len(sys.argv) > 2 else 700


def ob(x):
    return x.abs_upper()


def main():
    t0 = time.time()
    K = Kneser(PREC)
    R = K.R
    L = R.L
    c = K.c
    epsU = K.eps                          # ε ≥ ε_U (gilt für Ũ auf Im ≥ 0.35)
    y1 = arb(Y1)
    assert y1 > arb("0.35")
    # ---------------- (a) Streifen y1 ≤ Im ≤ 3
    s = arb("0.02")                       # halbe Kantenlänge der Quadrate
    nx = 26                               # Re ∈ [−0.52, 0.52] (Überlapp über ±½ für die periodische Fortsetzung)
    ny = int(round((3.0 - float(y1.mid())) / 0.04)) + 1
    delta_a = arb(0)
    kappa_max = arb(0)
    anzahl = 0
    r_tau = arb("1e-120")
    I = acb(0, 1)
    tpiI = 2 * acb.pi() * I
    dc = [n * cn for n, cn in enumerate(c)]
    for i in range(nx):
        xc = arb(-0.50 + 0.04 * i)
        for j in range(ny):
            yc = y1 + s + arb(0.04) * j
            zc = acb(xc, yc)
            Zb = acb(arb(xc, s), arb(yc, s))
            # θ̃' über Z_b (Reihe in q, harmlos) ⇒ ζ̃(Z_b) ⊂ D(ζ_c, s1), ζ̃(z) = z + θ̃(z)
            qb = (tpiI * Zb).exp()
            dth = R._horner(dc, qb) * tpiI
            thc = R.theta_ball(c, zc)
            zeta_c = zc + acb(thc.real.mid(), thc.imag.mid())
            s1 = s * arb(2).sqrt() * (1 + ob(dth)) + ob(R.theta_ball(c, zc) - acb(thc.real.mid(), thc.imag.mid()))
            r = 2 * s1
            a, tau = R.U_taylor([acb(0)], zeta_c, 16, r, arb("0.2"))      # Taylor-Modell von Φ(ζ) = P(e^{Lζ})
            m1 = abs(a[1]).lower() - sum((j * ob(a[j]) * s1 ** (j - 1) for j in range(2, 17)), arb(0)) - tau / (r - s1)
            M2 = sum((j * (j - 1) * ob(a[j]) * (s1 + r_tau) ** (j - 2) for j in range(2, 17)), arb(0))                 + 2 * tau / (r - s1 - r_tau) ** 2
            if not (m1 > 0):
                raise RuntimeError(f"|Φ'| nicht von 0 weg bei z = {zc.str(5)}")
            kap = M2 * r_tau / m1
            if not (kap < 1 and epsU / m1 < (1 - kap) * r_tau):
                raise RuntimeError(f"Krawczyk τ versagt bei z = {zc.str(5)}")
            kappa_max = amax(kappa_max, kap)
            bound = epsU / ((1 - kap) * m1)
            delta_a = amax(delta_a, bound)
            anzahl += 1
        if i % 5 == 0:
            print(f"  Spalte {i}/{nx}: δθ(a) ≤ {delta_a.str(3)}, κ ≤ {kappa_max.str(3)}  ({time.time() - t0:.0f}s)", flush=True)
    # ---------------- (b) Im ≥ 3
    r = arb("0.1")
    kap1 = sum((n * ob(R.p[n]) * r ** (n - 1) for n in range(2, R.n0 + 1)), arb(0)) \
        + R.C * (2 * r) ** (R.n0 + 1) / r                # |P' − 1| auf D(0, 0.1) (Rest per Cauchy auf Radius r)
    kap2 = sum((ob(R.p[n]) * r ** (n - 1) for n in range(2, R.n0 + 1)), arb(0)) + R.C * r ** R.n0
    assert kap1 < 1
    w_oben = K.w_oben
    # Rouché: |P(w) − L| ≥ 0.1(1 − κ₂) auf |w| = 0.1; |sexp_K − L| ≤ |w̃|(1 + κ₂) + ε auf Im ≥ 3; |L·W| < 0.1
    rouche = w_oben * (1 + kap2) + epsU < r * (1 - kap2)
    LW_ok = R.absL * (w_oben + epsU / (1 - kap1)) < r
    assert rouche and LW_ok
    fakt = ob((abs(L.real) * arb("0.5") + L.imag * 3).exp())   # sup e^{−Re(Lz)}, |Re z| ≤ ½, Im z = 3
    dE = epsU / (1 - kap1) * fakt * arb("1.0001")
    # inf |Ẽ| = inf |e^{Lθ̃}| auf |q| ≤ e^{−6π}
    tpi = 2 * arb.pi()
    qmax = (-tpi * 3).exp()
    rest = sum((ob(cn) * qmax ** n for n, cn in enumerate(c) if n >= 1), arb(0))
    Et_inf = ((L * c[0]).real - R.absL * rest).exp()
    assert dE < Et_inf / 2
    delta_b = 2 * dE / (Et_inf * R.absL)
    delta = amax(delta_a, delta_b)
    print(f"(a) Streifen [{Y1}, 3]: {anzahl} Quadrate, κ ≤ {kappa_max.str(3)}, δθ ≤ {delta_a.str(3)}")
    print(f"(b) Im ≥ 3: |P'−1| ≤ {kap1.str(3)} auf D(0,0.1), Rouché ok, |E − Ẽ| ≤ {dE.str(3)}, inf|Ẽ| ≥ {Et_inf.str(4)}, δθ ≤ {delta_b.str(3)}")
    print(f"⇒ |θ_K − θ̃| ≤ δθ = {delta.str(3)} auf Im z ≥ {Y1}  ({time.time() - t0:.0f}s)")
    (HIER / "p32a_theta_K.json").write_text(json.dumps({
        "beschreibung": "θ_K (Knesers exakte obere Struktur) gegen θ̃ aus p26_iter6: |θ_K − θ̃| ≤ delta_theta auf Im z ≥ y1; "
                        "|κ_n − c̃_n| ≤ delta_theta·e^{2πn·y1}. Beweis im Kopf von p32a_theta_K.py.",
        "y1": Y1, "delta_theta": delta.str(5, radius=False), "delta_a": delta_a.str(5), "delta_b": delta_b.str(5),
        "kappa_max": kappa_max.str(5), "quadrate": anzahl, "prec": PREC}, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
