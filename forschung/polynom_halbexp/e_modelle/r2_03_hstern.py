"""Runde 2, Schritt 3: (H*) und Gesetz R im Möbius-Drehungsmodell (symmetrisch, Fixpunkte ±i).
Ḡ' = D̄(0,1) ∪ D̄(a₀,R₀), a₀ = τ, R₀ = √(1+τ²); Ecken bei ±i mit Aussenwinkel ω = π − c/2 (Cayley: Sektor), α = ω/π.
(1) dist(ζ,Ḡ')·max_{|z|≤1}|Π_N k_ζ(z)| für ζ im Abstand t/N (Ecke: t/N^α) vor dem Rand;
(2) am Pol ζ_p = 1/τ: max|E_N(ζ_p,·)| gegen ρ^N (Gesetz R) und e^{−N g(ζ_p)} (Maximumprinzip, Green-Funktion explizit)."""
import sys, time
import mpmath as mp
from flint import arb, acb, ctx
from moebius_lin import setup, pol_antwort, auswerten

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 30
c = mp.mpf(sys.argv[1])
Ns = [int(x) for x in sys.argv[2].split(",")]
tau = mp.tan(c / 2)
R0 = mp.sqrt(1 + tau ** 2)
om = mp.pi - c / 2
alpha = om / mp.pi
rho = mp.sin(c / 2) / mp.cos(c)
eg = mp.sin(mp.pi * c / (2 * om)) / mp.sin(mp.pi * (mp.pi - 2 * c) / (2 * om))
print(f"c = {c}: τ = {mp.nstr(tau, 6)}, R₀ = {mp.nstr(R0, 6)}, ω = {mp.nstr(om * 180 / mp.pi, 5)}°, α = {mp.nstr(alpha, 4)}; "
      f"Gesetz R am Pol: ρ = {mp.nstr(rho, 6)}, Maximumprinzip e^(−g(ζ_p)) = {mp.nstr(eg, 6)}")
K = 256
zs = [acb(mp.nstr(mp.cos(2 * mp.pi * j / K), 25), mp.nstr(mp.sin(2 * mp.pi * j / K), 25)) for j in range(K)]


def maxD(coef, zeta):
    best = mp.mpf(0)
    for z in zs:
        v = auswerten(coef, z)
        best = max(best, abs(mp.mpc(v.real.mid().str(20, radius=False), v.imag.mid().str(20, radius=False))))
    return best


orte = {
    "Ecke i, Winkelhalbierende": (mp.mpc(0, 1), mp.expj(mp.pi / 2 + c / 4), True),
    "links −1": (mp.mpc(-1, 0), mp.mpc(-1, 0), False),
    "rechts a₀+R₀": (tau + R0, mp.mpc(1, 0), False),
    "D(a₀) oben 60°": (tau + R0 * mp.expj(mp.pi / 3), mp.expj(mp.pi / 3), False),
    "D oben 120°": (mp.expj(2 * mp.pi / 3), mp.expj(2 * mp.pi / 3), False),
}
for N in Ns:
    t0 = time.time()
    S = setup(c, N, 900)
    Ji = S["J"].inv()
    zeile = []
    for name, (b, nrm, ecke) in orte.items():
        werte = []
        for t in (mp.mpf("0.5"), mp.mpf(2), mp.mpf(8)):
            d = t / (N ** alpha if ecke else N)
            zeta = b + d * nrm
            coef = pol_antwort(S, acb(mp.nstr(mp.re(zeta), 25), mp.nstr(mp.im(zeta), 25)), Ji)
            werte.append(mp.nstr(d * maxD(coef, zeta), 3))
        zeile.append(f"{name}: {'/'.join(werte)}")
    # Pol
    zp = 1 / tau
    coef = pol_antwort(S, acb(mp.nstr(zp, 25)), Ji)
    Emax = mp.mpf(0)
    for z in zs:
        v = auswerten(coef, z) - 1 / (acb(mp.nstr(zp, 25)) - z)
        Emax = max(Emax, abs(mp.mpc(v.real.mid().str(20, radius=False), v.imag.mid().str(20, radius=False))))
    print(f"N={N} ({time.time() - t0:.0f}s): " + " | ".join(zeile), flush=True)
    print(f"      Pol: max_D̄|E_N(ζ_p,·)| = {mp.nstr(Emax, 4)},  /ρ^N = {mp.nstr(Emax / rho ** N, 4)},  /e^(−Ng) = {mp.nstr(Emax / eg ** N, 4)}", flush=True)
