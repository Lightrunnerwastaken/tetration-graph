"""Strang V — wie weit weicht das echte p_N auf der Linse von f_K ab? (Geometrische Voraussetzung für K1″(a) mit 𝒥_P: P(Int γ) nahe f(Int γ).)
f_K auf der Linse über die Fortsetzung F(z) = log(ρF(w))/ρ, w = e^{ρz}/ρ (mk_kneser_arb, zertifizierte Koeffizienten + (K2)-Schwanz, Kugeln),
p_N exakt (60 Stellen) an denselben Punkten. Linse s(τ) = S(1 − e^{−(τ/τ₁)^p}) über Γin (τ = Winkelabstand zu arg L), τ₁ = 12°, p = 1.3.
Aufruf: python v29_pN_auf_linse.py S N1,N2,..."""
import json
import math
import sys
from pathlib import Path

from flint import acb, arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER / "m_kopie"))
from mk_kneser_arb import KneserArb  # noqa: E402

S = float(sys.argv[1])
Ns = [int(x) for x in sys.argv[2].split(",")]
KN = KneserArb(str(HIER / "kneser_halbexp_b_zert_V.json"))
ctx.prec = 200
PN = json.loads((HIER.parents[1] / "gemeinsam" / "pN_N100_bis_1000.json").read_text(encoding="utf-8"))["p"]
thL = float(KN.th.mid())
tau_end = 2 * math.pi - 2 * thL
for N in Ns:
    ctx.prec = 600
    p = [arb(arb(x).mid()) for x in PN[str(N)]]     # Mitten als exakte Zahlen (p_N wie gegeben)
    ctx.prec = 200
    worst = (0.0, None)
    zeilen = []
    for tg in (2, 5, 10, 20, 40, 60, 90, 103.4):
        t = math.radians(tg)
        s = S * (1 - math.exp(-(t / math.radians(12)) ** 1.3))
        phi = thL + t
        z = acb(arb(math.exp(s) * math.cos(phi)), arb(math.exp(s) * math.sin(phi)))     # z = x/ρ
        xb = z * KN.rho
        x = acb(xb.real.mid(), xb.imag.mid())          # exakter Punkt (sonst Umwicklung in Horner)
        Fz = KN.F_gin(x / KN.rho)[0]                     # f_K(x)/ρ
        ctx.prec = 1400
        v = acb(0)
        for c in reversed(p):
            v = v * x + c
        ctx.prec = 200
        d = abs(v / KN.rho - Fz) * KN.rho
        assert float(d.rad()) < 0.1 * float(d.mid()) + 1e-11, (tg, d)
        dv = float(d.upper())
        zeilen.append(f"{tg:g}°: {dv:.1e}")
        if dv > worst[0]:
            worst = (dv, tg)
    print(f"S = {S}, N = {N}: |p_N − f_K| auf der Linse: " + ", ".join(zeilen) + f"  (max {worst[0]:.1e} bei {worst[1]}°)", flush=True)
