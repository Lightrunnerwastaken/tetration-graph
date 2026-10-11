"""M17 — Endstück (Ecke L) streng: ε_L ≥ ‖1_{t≤−T1}(K̃ − A)1_{t≤−T1}‖ (Streifen, unitäres Bild, beide Kanten).
Unterkante (Γin): M = e^{iθ}q, Oberkante (Γout): M = −e^{iθ}q, q = e^t. Analytische Funktionen von M (|M| ≤ R):
  g1 = Λ(ξ(M)) − ℓ,  g2 = 1/f'(f(ξ(M))) − 1/μ,  g3 = Λ(h(ξ(M))) − ℓ,  g4 = f'(ξ(M)) − μ,
  Λ(x) = log[(c0 + ε(x+L))/(c0 + ε(x+L̄))]; G(t) = t + Λ(ξ), H(t) = t + iπ − Λ(h(ξ)).
Schwarz-Lemma: |g_i(M)| ≤ γ_i|M|, γ_i = max_{|M|=R}|g_i|/R (Analytizität auf |M| ≤ R per Kästchenüberdeckung geprüft).
Schranke (Schur, gleichmässige Majoranten in Δ = t' − t, geschlossene Integrale ∫e^{∓α/2}/(cosh α − cos y)dα = π/sin(y/2)).
Aufruf: python m17_ende.py R "T1,T1,..."
"""
import sys

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 100
R = arb(sys.argv[1]) if len(sys.argv) > 1 else arb("0.3")
T1s = [arb(x) for x in (sys.argv[2] if len(sys.argv) > 2 else "3,4,5").split(",")]

pi = arb.pi()
eps = arb("0.31382")
th = arb("76.618") * pi / 180
c0 = 1 - 2 * eps * th.cos()
L = acb(th.cos(), th.sin())
Lb = L.conjugate()
mu = acb(1, 2 * eps * th.sin())
ell = mu.log()
b = mu.arg()
amu = abs(mu)


def f(z):
    return eps + c0 * z + eps * z * z


def fp(z):
    return c0 + 2 * eps * z


def xi(M):
    return (Lb * M - L) / (M - 1)


def rad_h(x):
    return c0 * c0 - 4 * eps * (eps - x)


def h(x):
    return (-c0 + rad_h(x).sqrt()) / (2 * eps)


def ratio(x):
    return (c0 + eps * (x + L)) / (c0 + eps * (x + Lb))


def ok_halbachse(z):
    """True, wenn die Kugel z die abgeschlossene negative reelle Achse sicher meidet."""
    return bool(z.real > 0) or bool(z.imag > 0) or bool(z.imag < 0)


# 1) Analytizität auf |M| ≤ R: Quadratgitter-Kästchen über [−R, R]²; relevant, wenn das Kästchen die Scheibe trifft
nbox = 60
Rf = float(R.upper())
okA = True
schlecht = 0
for i in range(nbox):
    for j in range(nbox):
        xc = -Rf + 2 * Rf * (i + 0.5) / nbox
        yc = -Rf + 2 * Rf * (j + 0.5) / nbox
        rr = Rf / nbox
        dx, dy = max(abs(xc) - rr, 0.0), max(abs(yc) - rr, 0.0)
        if dx * dx + dy * dy > Rf * Rf * (1 - 1e-9):
            continue                      # Kästchen liegt ausserhalb der Scheibe
        M = acb(arb(xc, rr * (1 + 1e-12)), arb(yc, rr * (1 + 1e-12)))
        x = xi(M)
        v = h(x)
        good = ok_halbachse(ratio(x)) and ok_halbachse(rad_h(x)) and ok_halbachse(ratio(v)) and bool(abs(fp(f(x))) > 0)
        if not good:
            okA = False
            schlecht += 1
print(f"Analytizität auf |M| ≤ {Rf}: {'geprüft' if okA else 'NICHT geprüft'} ({schlecht} kritische Kästchen)")

# 2) Maxima auf |M| = R
nk = 720
mx = [arb(0)] * 4
for k in range(nk):
    a0 = 2 * pi * k / nk
    a1 = 2 * pi * (k + 1) / nk
    am = (a0 + a1) / 2
    phi = arb(am.mid(), ((a1 - a0) / 2).upper())
    M = R * acb(phi.cos(), phi.sin())
    x = xi(M)
    v = h(x)
    gs = [ratio(x).log() - ell, 1 / fp(f(x)) - 1 / mu, ratio(v).log() - ell, fp(x) - mu]
    mx = [arb(max(m.upper(), abs(g).upper())) for m, g in zip(mx, gs)]
gam = [m / R for m in mx]
print("γ1..γ4 =", [str(g.upper())[:10] for g in gam])


def J(y_lo, y_hi, shift_sign, du, kern):
    """Obere Schranke für ∫ e^{−Δ/2} m(Δ) dΔ, m Majorante von P0/P1/|∇P0|/|∇P1| bei Höhe y ∈ [y_lo, y_hi]
    (in (0,π)), horizontale Lage Re ℓ·shift_sign ± du. kern ∈ {P0, P1, G0, G1}.
    cosh(α − τ) ≥ e^{−du} cosh α. Für P0/G0: Nenner cosh(α−τ) − cos y ≥ e^{−du}(cosh α − e^{du} cos y_lo') usw."""
    ed = du.exp()
    # Faktor aus e^{−Δ/2} = e^{−α/2} e^{−(Re ℓ·s + τ)/2} ≤ e^{−α/2} e^{−s Re ℓ/2} e^{du/2}
    vor = (-shift_sign * ell.real / 2).exp() * (du / 2).exp()
    if kern in ("P0", "G0"):
        # cos y maximal bei y_lo (falls cos y > 0) — sichere Wahl: cmax = max(cos y) über [y_lo, y_hi]
        cmax = y_lo.cos() if bool(y_lo < pi / 2) else arb(-1) * (-(y_hi.cos()))
        cmax = arb(max(y_lo.cos().upper(), y_hi.cos().upper()))
        ceff = ed * cmax if bool(cmax > 0) else cmax / ed       # Nenner ≥ e^{−du}(cosh α − ceff)
    else:   # P1, G1: Nenner cosh(α−τ) + cos y; cmin = min cos y
        cmin = arb(min(y_lo.cos().lower(), y_hi.cos().lower()))
        ceff = -(cmin / ed if bool(cmin > 0) else cmin * ed)    # Nenner ≥ e^{−du}(cosh α − ceff)
    assert bool(ceff < 1)
    yeff = ceff.acos()                    # cosh α − cos(yeff)
    integral = pi / (yeff / 2).sin()      # ∫ e^{−α/2}/(cosh α − cos yeff) dα
    if kern in ("P0", "P1"):
        smax = arb(1) if bool(y_lo < pi / 2) and bool(y_hi > pi / 2) else arb(max(y_lo.sin().upper(), y_hi.sin().upper()))
        return vor * ed * smax / (2 * pi) * integral
    return vor * ed / (2 * pi) * integral


for T1 in T1s:
    q1 = (-T1).exp()
    g1, g2, g3, g4 = gam
    # Unterkante: Punkt G(t), Höhe b ± g1 q1, Lage t + Re ℓ ± g1 q1
    du0 = g1 * q1
    yl0, yh0 = b - du0, b + du0
    # Oberkante: Punkt H(t) = t + iπ − Λ(h(ξ)), Höhe π − b ± g3 q1, Lage t − Re ℓ ± g3 q1
    du1 = g3 * q1
    yl1, yh1 = pi - b - du1, pi - b + du1
    # Gewichte ρ: |ρ0/ρ_e' − 1| ≤ 2 q1/(1 − q1)
    dr = 2 * q1 / (1 - q1)
    d_in = dr * (1 / amu + g2 * q1) + g2 * q1           # |(ρ0/ρ')W_in − μ⁻¹|
    d_out = dr * (amu + g4 * q1) + g4 * q1              # |(ρ1/ρ')W_out − μ|
    Jm = {}
    for e, (yl, yh, sgn, du, dvor, wmod, gq) in enumerate(
            [(yl0, yh0, 1, du0, d_in, 1 / amu, g1 * q1), (yl1, yh1, -1, du1, d_out, amu, g3 * q1)]):
        for ep, (kP, kG) in enumerate([("P0", "G0"), ("P1", "G1")]):
            Jm[(e, ep)] = dvor * J(yl, yh, sgn, du, kP) + wmod * gq * J(yl, yh, sgn, du, kG)
    zeilen = [Jm[(0, 0)] + Jm[(0, 1)], Jm[(1, 0)] + Jm[(1, 1)]]
    spalten = [Jm[(0, 0)] + Jm[(1, 0)], Jm[(0, 1)] + Jm[(1, 1)]]
    Rz = arb(max(z.upper() for z in zeilen))
    Cs = arb(max(z.upper() for z in spalten))
    epsL = (Rz * Cs).sqrt()
    print(f"T1 = {T1}: q1 = {float(q1.upper()):.3e}, Blöcke J = {[float(Jm[k].upper()) for k in sorted(Jm)]}, "
          f"ε_L ≤ {float(epsL.upper()):.5f}")
