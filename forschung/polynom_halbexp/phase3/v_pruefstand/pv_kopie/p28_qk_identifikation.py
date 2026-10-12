"""Strang P, Satz P-V Schritt 4: quasikonforme Glättung (Beltrami) + Identifikation (Paulsen–Cowgill) — rigoros (arb).

Konstruktion (Lipschitz-Abschneider, Ebene, Fundamentalstreifen −½ ≤ Re z < ½):
  χ(x): 0 für x ≤ ½−ρ, linear auf [½−ρ, ½−ρ/3], 1 danach;   ψ(y): 1 für y ≤ H, linear fallend auf [H, H+ρ/2], 0 danach.
  W₀(z) = χ(Re z)·d(z − 1)   (Re z ∈ [½−ρ, ½]),  W₀ = 0 sonst  ⇒  W₀(z+1) − W₀(z) = d(z + W₀(z)) auf der Naht (exakt).
  Y(z) = W₀(z) + e(z + W₀(z)) (periodisch: Identität e(ζ) − e(ζ+1+d(ζ)) = d(ζ)),  V₀ = ψ(Im z)·Y.
  F₀ = T̃∘(id + W₀) für |Im z| ≤ H,  F₀ = Ũ∘(id + V₀) für Im z ≥ H − ρ/2 (dort ψ = 1 ⇒ beide gleich), gespiegelt unten.
  F₀ ist Lipschitz, erfüllt F₀(z+1) = e^{F₀(z)} EXAKT, Beltrami-Koeffizient μ₀ = ∂̄F₀/∂F₀ mit |μ₀| ≤ k, Träger kompakt.
  Messbarer Riemannscher Abbildungssatz auf dem Zylinder (q = e^{2πiz}): φ qk, φ(z+1) = φ(z)+1, μ_φ = μ₀;  F := F₀∘φ⁻¹ analytisch,
  exakte Lösung, reell, Grenzwerte L/L̄.  Schranke |φ − id| ≤ (1/2π)·ε/(1−ε), ε ≤ (1/π)‖h‖₄·‖1/(|ζ||ζ−q|)‖_{4/3},
  ‖h‖₄ ≤ k|A|^{1/4}/(1 − k‖S‖₄),  ‖S‖₄ ≤ 12 (Bañuelos–Wang; jede endliche Schranke genügt, da k winzig).
Identifikation: Paulsen–Cowgill Prop. 2 (Fassung Paper VI Thm 11.11) — Holomorphie auf Re z > −2 per Funktionalgleichung,
  dafür Prüfung A: F ∉ (−∞, 0] auf dem Streifen, Prüfung B: F ∉ (0, 1] für 0 < Re ζ < ½ (Taylor-Modelle auf Scheiben).
Aufruf: python p28_qk_identifikation.py DATEN.json [prec]
"""
import json
import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

import kc_rig

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
DATEN = HIER / sys.argv[1]
PREC = int(sys.argv[2]) if len(sys.argv) > 2 else 600


def amax(a, b):
    return (a + b + abs(a - b)) / 2


def main():
    t0 = time.time()
    ctx.prec = PREC
    Z = json.loads((HIER / "p27_zert.json").read_text(encoding="utf-8"))
    H, RHO = arb(Z["H"]), arb(Z["rho"])
    e_max, d_max = arb(Z["e_max"]).abs_upper(), arb(Z["d_max"]).abs_upper()
    marg = arb("0.05")
    e1 = e_max / marg          # |e'| auf dem inneren Band (Cauchy, Rand 0.05)
    d1 = d_max / marg
    chi1 = 1 / (RHO * 2 / 3)   # |χ'|
    psi1 = 2 / RHO             # |ψ'|
    # Größen
    dbarW = chi1 / 2 * d_max
    dW = chi1 / 2 * d_max + d1
    Ymax = d_max + e_max
    dbarY = (1 + e1) * dbarW
    dY = dW * (1 + e1) + e1
    dbarV = psi1 / 2 * Ymax + dbarY
    dV = psi1 / 2 * Ymax + dY
    k1 = dbarW / (1 - dW)
    k2 = dbarV / (1 - dV)
    k = amax(k1, k2)
    print(f"e_max = {e_max.str(4)}, d_max = {d_max.str(4)};  |μ₀| ≤ k = {k.str(4)}")
    # Beltrami auf dem Zylinder: q = e^{2πiz}, Träger in |q| ∈ [e^{−2π(H+ρ/2)}, e^{2π(H+ρ/2)}]
    tpi = 2 * arb.pi()
    rmin = (-tpi * (H + RHO / 2)).exp()
    rmax = (tpi * (H + RHO / 2)).exp()
    A = arb.pi() * rmax ** 2
    S4 = arb(12)
    h4 = k * A ** arb("0.25") / (1 - k * S4)
    rA = (A / arb.pi()).sqrt()
    L43 = (3 * arb.pi() * rA ** (arb(2) / 3)) ** arb("0.75")      # ‖1/|ζ−q|‖_{L^{4/3}(A)} ≤ (∫_{Scheibe gleicher Fläche})^{3/4}
    eps = h4 * L43 / rmin / arb.pi()
    assert eps < arb("0.5")
    phi_fehler = eps / (1 - eps) / tpi
    print(f"‖h‖₄ ≤ {h4.str(3)},  ε ≤ {eps.str(3)},  |φ − id| ≤ {phi_fehler.str(4)}  ({time.time() - t0:.0f}s)")
    raus = dict(k=k.str(6), eps=eps.str(6), phi_fehler=phi_fehler.str(6), e1=e1.str(4), d1=d1.str(4))
    # ---------------- Identifikation: Prüfungen A und B mit Taylor-Modellen niedriger Ordnung
    R = kc_rig.Rig(PREC, n0=200, m_extra=200)
    d = json.loads(DATEN.read_text(encoding="utf-8"))
    t = [acb(arb(arb(a).mid())) for a, b in d["t"]]
    c = [acb(arb(arb(a).mid()), arb(arb(b).mid())) for a, b in d["c"]]
    rr = arb("0.025")
    schritt = 0.035                  # Gitter; Scheibenradius 0.025 ≥ halbe Diagonale 0.0247
    # Fehler von F gegenüber T̃ bzw. Ũ (grob, mit Lipschitz-Schranke 10 für |∇F₀|): vernachlässigbar gegen 1e-3
    fehl = 10 * (d_max + phi_fehler + Ymax)
    ok_A, ok_B = True, True
    n = 0
    ys = [0.0]
    y = 0.0
    while y < 3.0:
        y += schritt
        ys.append(y)
    xs = [-0.5 + schritt * i for i in range(int(1.0 / schritt) + 1)]
    xs = sorted(set(xs + [0.0125 * i for i in range(1, 6)]))          # feiner nahe Re = 0
    dt = [nn * t[nn] for nn in range(1, len(t))]
    ddt = [nn * dt[nn] for nn in range(1, len(dt))]
    # Prüfung B nahe 0 analytisch: G(0) = 1, G'(0) ≥ g1, |G''| ≤ M2 auf D(0, 0.1) ⇒ kein ζ mit 0 < Re ζ, |ζ| < √2·g1/M2, G(ζ) ∈ (0,1]
    g1 = R.T_ball(dt, acb(0)).real - 10 * (d_max + phi_fehler)
    M2 = R.T_ball(ddt, kc_rig.disc(acb(0), arb("0.1"))).abs_upper() + 1
    rB = -amax(-arb("0.1"), -(arb(2).sqrt() * g1 / M2))
    print(f"B nahe 0: G'(0) ≥ {g1.str(5)}, |G''| ≤ {M2.str(4)} ⇒ analytisch erledigt für |ζ| < {rB.str(4)}")
    min_abst = None
    for yv in ys:
        for xv in xs:
            zc = acb(arb(xv), arb(yv))
            if yv <= float(H.mid()):
                val = R.T_ball(t, zc)                                   # Punktwert (enger Ball)
                var = R.T_ball(dt, kc_rig.disc(zc, rr)).abs_upper() * rr + val.real.rad() + val.imag.rad()
                val = acb(val.real.mid(), val.imag.mid())
            else:
                a, tau = R.U_taylor(c, zc, 10, rr, arb("0.12"))
                val = a[0]
                var = sum((abs(a[j]) * rr ** j for j in range(1, 11)), arb(0)) + tau
            var = var + fehl
            # A: Abstand zu (−∞, 0]:  Im ≠ 0 oder Re > 0
            im_lo = abs(val.imag) - var
            re_lo = val.real - var
            okA = bool(im_lo > 0) or bool(re_lo > 0)
            # B (nur 0 < Re < ½): nicht in (0, 1]: Im ≠ 0 oder Re > 1
            okB = True
            if 0 < xv < 0.5 and not (abs(zc) < rB - rr):          # Scheibe ganz in D(0, r_B): analytisch erledigt
                okB = bool(im_lo > 0) or bool(val.real - var > 1)
            if not okA:
                ok_A = False
                print(f"  Prüfung A offen bei {xv:+.2f}{yv:+.2f}i: F ≈ {val.str(5)}, ±{var.str(3)}")
            if not okB:
                ok_B = False
                print(f"  Prüfung B offen bei {xv:+.2f}{yv:+.2f}i: F ≈ {val.str(5)}, ±{var.str(3)}")
            n += 1
    print(f"Prüfungen A/B auf {n} Scheiben (0 ≤ Im ≤ 3): A {'erfüllt' if ok_A else 'OFFEN'}, B {'erfüllt' if ok_B else 'OFFEN'}  ({time.time() - t0:.0f}s)")
    # Im z ≥ 3: |w| = |L^{z+θ̃}| ≤ e^{0.318·½ − 1.337·3}·|e^{L θ̃}|, |P(w) − L| ≤ 2|w| für |w| ≤ 0.1 ⇒ Im F > 0
    Lc = R.L
    thmax = sum((abs(ck) * (-tpi * 3 * kk).exp() for kk, ck in enumerate(c) if kk >= 1), arb(0))
    wmax = (Lc.real * arb("0.5") - Lc.imag * 3 + (Lc * c[0]).real).exp() * (abs(Lc) * thmax).exp()
    assert wmax < arb("0.1")
    Bp = 1 + sum((abs(R.p[nn]) * arb("0.1") ** (nn - 1) for nn in range(2, R.n0 + 1)), arb(0)) + R.C * arb("0.1") ** R.n0
    abw = Bp * wmax + fehl
    assert abw < Lc.imag
    print(f"Im z ≥ 3: |w| ≤ {wmax.str(3)} ⇒ |F − L| ≤ {abw.str(3)} < Im L = {Lc.imag.str(5)} ⇒ F ∉ ℝ")
    raus.update(pruefung_A=ok_A, pruefung_B=ok_B, w_oben=wmax.str(4))
    (HIER / "p28_zert.json").write_text(json.dumps(raus, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
