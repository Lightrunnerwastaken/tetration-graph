"""b25 (Runde 2): Test von X* (Pol-Antwort im Band um G'): Q_N(ζ) := dist(ζ, Ḡ')·max_{z∈D̄}|Π_N k_ζ(z)| für ζ knapp ausserhalb
von G' (Abstand t/N bzw. t·N^{-0.883}) an verschiedenen Randpunkten: Ecke L (Winkelhalbierende und schräg), linker Rand von D,
rechter Rand von D(a0,|L−a0|), oberer Rand von D(a0) nahe 1+1.1i, und Pol auf dem Schnitt. X* behauptet Q_N ≤ C N^κ (κ < 4.69)."""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 100
aL = abs(L); a0 = mp.mpf(K[0]); rA = abs(L - a0)
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 480) for j in range(480)]
kreis += [L * mp.expjpi(mp.mpf(t) / 1800) for t in range(-60, 61)] + [mp.conj(L) * mp.expjpi(mp.mpf(t) / 1800) for t in range(-60, 61)]
def dist_G(z):
    d1 = abs(z) - aL; d2 = abs(z - a0) - rA
    return min(d1, d2) if (d1 > 0 and d2 > 0) else mp.mpf(0)
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
orte = {
  "Ecke/Halbierende": lambda d: L + d * mp.expj(psi),
  "Ecke/radial":      lambda d: L * (1 + d / aL),
  "Ecke/flach(15°)":  lambda d: L + d * mp.expj(mp.arg(L - a0) - mp.pi / 2 + mp.radians(15)),
  "links (−|L|)":      lambda d: -(aL + d),
  "D(a0) rechts":      lambda d: a0 + rA + d,
  "D(a0) oben 60°":    lambda d: a0 + (rA + d) * mp.expjpi(mp.mpf(1) / 3),
  "D oben 120°":       lambda d: (aL + d) * mp.expjpi(mp.mpf(2) / 3),
}
for N in (60, 100, 140):
    pw, g, Ji = baue(N, J=N)
    for name, ort in orte.items():
        zeile = []
        for t in (0.5, 2, 8):
            d = mp.mpf(t) / N
            zeta = ort(d)
            dd = dist_G(zeta)
            P = loese(Ji, rhs_pol(N, pw, g, zeta))
            m = max(abs(horner(P, z)) for z in kreis)
            zeile.append(f"t={t}: dist·max|Π k| = {mp.nstr(dd * m, 3)}")
        print(f"N={N} {name:18s}: " + " | ".join(zeile), flush=True)
