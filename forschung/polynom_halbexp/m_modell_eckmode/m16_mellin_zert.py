"""M16 — Zertifikat (arb) der Eckkonstante des Modells: c_mod = inf_σ λ_min(I + Re A(σ)) (L²-Fassung, 2×2-Symbol).
A(σ) = [[μ⁻¹|μ|^s F0(b), μ⁻¹|μ|^s F1(b)], [μ|μ|^{−s} F0(π−b), μ|μ|^{−s} F1(π−b)]],  s = −1/2 + iσ, b = arg μ,
F0(y) = sin(s(π−y))/sin(sπ), F1(y) = sin(s y)/sin(sπ)  (Streifen-Poissonkern, Fourier mit Gewicht e^{(t−t')/2}).
Verfahren: Intervall-Bisektion auf [−S, S] mit σ als arb-Kugel; Schwänze |σ| > S über ‖A‖_F ≤ 2√2 e^{−|σ|b}(|μ|^{-3}+|μ|^{3})^{1/2}.
Aufruf: python m16_mellin_zert.py [Ziel] [S] [mu_wahl: radial|sqrtL]
"""
import sys
import time

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
ctx.prec = 80
ziel = arb(sys.argv[1]) if len(sys.argv) > 1 else arb("0.99")
S = int(sys.argv[2]) if len(sys.argv) > 2 else 30
wahl = sys.argv[3] if len(sys.argv) > 3 else "radial"

pi = arb.pi()
if wahl == "radial":
    eps = arb("0.31382")
    th = arb("76.618") * pi / 180
    mu = acb(1, 2 * eps * th.sin())
else:   # μ = √L, L = −W(−1) obere Fixstelle von exp
    L = -acb(-1).lambertw()
    L = acb(L.real, abs(L.imag))
    mu = L.sqrt()
amu = abs(mu)
b = mu.arg()
lam = amu.log()
print(f"μ = {mu}, |μ| = {amu}, b = {b}")


def lam_min(sig):
    """Kugel-Einschluss von λ_min(I + Re A(σ)) für σ als arb-Kugel."""
    s = acb(arb(-0.5), sig)
    sp = (s * pi).sin()
    F0b = (s * (pi - b)).sin() / sp
    F1b = (s * b).sin() / sp
    F0c = (s * b).sin() / sp          # F0(π−b)
    F1c = (s * (pi - b)).sin() / sp   # F1(π−b)
    p1 = (s * lam).exp() / mu         # μ⁻¹|μ|^s
    p2 = (-s * lam).exp() * mu        # μ|μ|^{−s}
    a11, a12, a21, a22 = p1 * F0b, p1 * F1b, p2 * F0c, p2 * F1c
    # Re(I + A) = [[1 + Re a11, (a12 + conj a21)/2], [conj(...), 1 + Re a22]]
    d1 = 1 + a11.real
    d2 = 1 + a22.real
    c = (a12 + a21.conjugate()) / 2
    q = abs((d1 - d2) / 2) ** 2 + abs(c) ** 2          # ≥ 0; obere Schranke genügt
    wq = arb(q.upper()).sqrt().upper()
    return arb(((d1 + d2) / 2).lower()) - arb(wq)        # exakte Unterschranke (Kugel mit Radius 0 + Rundung)


t0 = time.time()
# Schwanzschranke
sch = 2 * arb(2).sqrt() * (-S * b).exp() * (amu ** -3 + amu ** 3).sqrt()
print(f"Schwanz |σ| ≥ {S}: λ_min ≥ 1 − {sch.upper()}  (also ≥ Ziel: {bool(1 - sch > ziel)})")
# Bisektion
stack = [(arb(-S), arb(S))]
n_ok = 0
lo_min = None
fail = []
while stack:
    a, c = stack.pop()
    mid = (a + c) / 2
    rad = (c - a) / 2
    sig = arb(mid.mid(), rad.upper() + mid.rad())
    v = lam_min(sig)
    if v > ziel:
        n_ok += 1
        lb = v.lower()
        if lo_min is None or lb < lo_min:
            lo_min = lb
        continue
    if rad < arb("1e-9"):
        fail.append((a, c, v))
        continue
    stack.append((a, mid))
    stack.append((mid, c))
print(f"Intervalle zertifiziert: {n_ok}, gescheitert: {len(fail)}; kleinste Unterschranke {lo_min} ({time.time()-t0:.0f}s)")
for a, c, v in fail[:5]:
    print("  FEHLT:", a, c, v)
# Punktwerte zum Vergleich
best = min((lam_min(arb(x / 100)), x / 100) for x in range(-1000, 1000))
print(f"Rasterminimum (Schritt 0.01): {best[0]} bei σ = {best[1]}")
