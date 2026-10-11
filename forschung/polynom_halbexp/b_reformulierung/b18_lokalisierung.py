"""b18: Lokalisierung des Konsistenzfehlers r_N = T_N[Schwanz_N(f_K)∘f_K] auf dem Kreis |x| = |L|:
Γin = {|f_K(x)| < |L|} (linker Bogen L → −|L| → L̄), Γout = Rest. Ist r_N auf Γin klein und nur auf Γout gross?"""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b11_projektion_moden.py", encoding="utf-8").read().split("tests = {}")[0])
Kc = [mp.mpf(c) for c in K]
aL = abs(L); argL = mp.arg(L)
for N in (40, 80, 120):
    pw, g, Ji = baue(N, J=700)
    r = [arb(0)] * (N + 1)
    for j in range(N + 1, 701):
        cj = arb(K[j])
        for k in range(N + 1):
            r[k] += cj * pw[j][k]
    rv = [mp.mpf(x.mid().str(40, radius=False)) for x in r]
    zeile = []
    for phi_deg in (0, 30, 60, 70, 76.6, 80, 90, 120, 150, 180):
        x = aL * mp.expjpi(mp.mpf(phi_deg) / 180)
        zeile.append(f"{phi_deg}°:{mp.nstr(abs(horner(rv, x)), 2)}")
    print(f"N={N}: |r_N(|L|e^(iφ))|  (Γout: |φ| < 76.6°, Γin: |φ| > 76.6°):  " + "  ".join(zeile), flush=True)
