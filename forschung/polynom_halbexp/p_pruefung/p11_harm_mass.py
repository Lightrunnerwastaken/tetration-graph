"""Strang P, Prüfung 11: Konstanten aus C §12.9/§13.4 nachrechnen.
(1) min g/ω(D-Bogen) auf ∂(W + d) (d = 0.1, ohne |ζ − L| < 0.25), C: ≥ 0.0649.
(2) Schlitz-Grenzwert ω_ζ(D-Bogen ab |s| ≥ s1)/g(ζ) für ζ = L(1+u), u → 0, gegen 1/(κ0 s1).
Halbebenen-Bild τ = W(ζ): D-Bogen ↦ s < 0, A-Bogen ↦ s > 0; ω_τ((−∞,0)) = arg τ / π."""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
L = mp.findroot(lambda z: mp.exp(z) - z, mp.mpc(0.318, 1.337))
rho, a0 = abs(L), mp.mpf("0.49856328794111443467961909249313329400247")
R0 = abs(L - a0)
Lb = mp.conj(L)
beta = mp.arg((a0 + R0 - L) / (a0 + R0 - Lb))
omega = mp.arg((-rho - L) / (-rho - Lb)) - beta


def W(z):
    w = (z - L) / (z - Lb) * mp.expj(-beta)
    ang = mp.arg(w)
    if ang < -mp.mpf(10) ** -20:
        ang += 2 * mp.pi
    return abs(w) ** (mp.pi / omega) * mp.expj(ang * mp.pi / omega)


t0 = W(mp.mpf(10) ** 25)
kappa0 = 2 * mp.pi * mp.im(t0) / abs(t0) ** 2


def g(z):
    t = W(z)
    return mp.log(abs((t - mp.conj(t0)) / (t - t0)))


print("τ0 =", mp.nstr(t0, 8), " κ0 =", mp.nstr(kappa0, 8))
# (1) Kontur ∂(W + d): Punkte ausserhalb W̄ mit Abstand d
d = mp.mpf("0.1")
best = None
for k in range(2000):
    th = 2 * mp.pi * k / 2000
    for (c, r) in [(mp.mpf(0), rho), (a0, R0)]:
        z = c + (r + d) * mp.expj(th)
        if abs(z) < rho + d - mp.mpf(10) ** -12 or abs(z - a0) < R0 + d - mp.mpf(10) ** -12:
            continue
        if abs(z - L) < 0.25 or abs(z - Lb) < 0.25:
            continue
        t = W(z)
        om = mp.arg(t) / mp.pi
        q = g(z) / om
        if best is None or q < best[0]:
            best = (q, z)
print(f"(1) min g/ω(D-Bogen) auf ∂(W+0.1) ohne |ζ−L|<0.25: {mp.nstr(best[0], 5)} bei ζ = {mp.nstr(best[1], 5)}")
# (2) Schlitz
for s1 in [mp.mpf("0.0019"), mp.mpf("0.0109"), mp.mpf("0.024")]:
    for u in [mp.mpf("1e-4"), mp.mpf("1e-3"), mp.mpf("1e-2")]:
        z = L * (1 + u)
        t = W(z)
        om = (mp.pi - mp.arg(t + s1)) / mp.pi       # ω von (−∞, −s1]
        om = (mp.arg(t + s1)) / mp.pi                # (−∞,−s1) liegt links: Winkel von τ aus gesehen
        print(f"(2) s1 = {mp.nstr(s1, 3)}, u = {mp.nstr(u, 2)}: ω/g = {mp.nstr(om / g(z), 6)}   1/(κ0 s1) = {mp.nstr(1 / (kappa0 * s1), 6)}")
