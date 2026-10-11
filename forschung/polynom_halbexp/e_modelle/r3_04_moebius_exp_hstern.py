"""Runde 3, Schritt 4: „Möbius-exp“ als Prüfstein für Strang C — Drehung um L, L̄ (Fixpunkte von exp), c so dass g(0) = 1.
Gleiche W-Geometrie wie bei exp (a₀ = 0.478 statt 0.4986), aber KEINE Eckmode. (H*) im Band und Green-Gesetz am Pol von f_M."""
import sys, time
import mpmath as mp
from flint import acb
from moebius_lin import setup_allg, pol_antwort, auswerten

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
Ns = [int(x) for x in sys.argv[1].split(",")]
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L
A, B = mp.re(L), mp.im(L)
c = mp.atan((1 - A) / B) + mp.atan(A / B)
S0 = setup_allg(A, B, c, 4, 300)
a0 = mp.mpf(S0["a0"].mid().str(30, radius=False))
zp = mp.mpf(S0["pol"].mid().str(30, radius=False))
rL, R0 = abs(L), abs(L - a0)
Mf = lambda z: (z - L) / (z - mp.conj(L))
d1, d2 = mp.arg(Mf(-rL)), mp.arg(Mf(a0 + R0))
lo, hi = sorted([d1, d2])
beta, om = (lo, hi - lo) if lo < 0 < hi else (hi, lo + 2 * mp.pi - hi)


def Wk(z):
    m = Mf(z) * mp.expj(-beta)
    a = mp.arg(m)
    if a < 0:
        a += 2 * mp.pi
    return abs(m) ** (mp.pi / om) * mp.expj(a * mp.pi / om)


Winf = Wk(mp.mpf(10) ** 12)
green = lambda z: mp.log(abs((Wk(z) - mp.conj(Winf)) / (Wk(z) - Winf)))
alpha = om / mp.pi
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
print(f"Möbius-exp: c = {mp.nstr(c, 8)}, a₀ = {mp.nstr(a0, 8)}, R₀ = {mp.nstr(R0, 6)}, Pol ζ_p = {mp.nstr(zp, 6)}, Aussenwinkel "
      f"ω = {mp.nstr(om * 180 / mp.pi, 5)}°, α = {mp.nstr(alpha, 4)}; e^(−g(ζ_p)) = {mp.nstr(mp.exp(-green(zp)), 6)}, "
      f"Max-Formel {mp.nstr(max(rL / zp, R0 / (zp - a0)), 6)}")
K = 256
zs = [acb(mp.nstr(mp.cos(2 * mp.pi * j / K) * rL, 25), mp.nstr(mp.sin(2 * mp.pi * j / K) * rL, 25)) for j in range(K)]   # |z| = |L|


def maxD(coef, zeta=None):
    best = mp.mpf(0)
    for z in zs:
        v = auswerten(coef, z) - (1 / (zeta - z) if zeta is not None else 0)
        best = max(best, abs(mp.mpc(v.real.mid().str(20, radius=False), v.imag.mid().str(20, radius=False))))
    return best


orte = {"Ecke L, Winkelhalbierende": (L, mp.expj(psi), True), "links −|L|": (-rL, mp.mpc(-1, 0), False),
        "rechts a₀+R₀": (a0 + R0, mp.mpc(1, 0), False), "∂D oben 120°": (rL * mp.expjpi(mp.mpf(2) / 3), mp.expjpi(mp.mpf(2) / 3), False)}
for N in Ns:
    t0 = time.time()
    S = setup_allg(A, B, c, N, 1000)
    Ji = S["J"].inv()
    zeile = []
    for name, (b, nrm, ecke) in orte.items():
        w = []
        for t in (mp.mpf("0.5"), mp.mpf(2), mp.mpf(8)):
            d = t / (N ** alpha if ecke else N)
            zeta = b + d * nrm
            coef = pol_antwort(S, acb(mp.nstr(mp.re(zeta), 25), mp.nstr(mp.im(zeta), 25)), Ji)
            w.append(mp.nstr(d * maxD(coef), 3))
        zeile.append(f"{name}: {'/'.join(w)}")
    zpa = acb(mp.nstr(zp, 25))
    coef = pol_antwort(S, zpa, Ji)
    E = maxD(coef, zpa)
    print(f"N={N} ({time.time() - t0:.0f}s): " + " | ".join(zeile), flush=True)
    print(f"      Pol: max_(|z|=|L|)|E_N(ζ_p,·)| = {mp.nstr(E, 4)}, ·e^(Ng) = {mp.nstr(E * mp.exp(N * green(zp)), 4)}", flush=True)
