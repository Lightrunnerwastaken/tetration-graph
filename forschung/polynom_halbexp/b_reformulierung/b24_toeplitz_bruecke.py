"""b24 (Runde 2): Brücke zu Strang D. Ist J_N⁻¹ T_N F ≈ P_N T(σ) F für F = F_ζ := (Schwanz_N k_ζ)∘f ?
Exakt: Π_N k_ζ − T_N k_ζ = J_N⁻¹ T_N F_ζ,  F_ζ(x) = (f(x)/ζ)^{N+1}/(ζ − f(x))  (C's Formel).
Ds Grenzoperator: G = T(σ), σ = f'(x)e^{−x}·1_{Γin}(x) auf |x| = |L|, Γin = {|f(x)| < |L|}.
Vergleich der beiden Polynome auf |z| = |L| und an einzelnen Punkten (z = 0.5, 1, L)."""
import json, sys
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 50
aL = abs(L)
a0 = mp.mpf(K[0])
Kc = [mp.mpf(c) for c in K]
dK = [k * Kc[k] for k in range(1, len(Kc))]
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
MF = 2048
xs = [aL * mp.expjpi(2 * mp.mpf(j) / MF) for j in range(MF)]
fx = [horner(Kc, x) for x in xs]
sig = [horner(dK, x) * mp.exp(-x) if abs(fv) < aL else mp.mpf(0) for x, fv in zip(xs, fx)]
print("Anteil Γin:", sum(1 for s in sig if s != 0) / MF)


def toeplitz_poly(N, zeta):
    Fv = [(fv / zeta) ** (N + 1) / (zeta - fv) for fv in fx]
    prod = [s * F for s, F in zip(sig, Fv)]
    c = []
    for k in range(N + 1):
        ck = mp.fsum(p * mp.expjpi(-2 * mp.mpf(k * j) / MF) for j, p in enumerate(prod)) / MF
        c.append(ck / aL ** k)
    return c


kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 360) for j in range(360)] + [L]
for N in (60, 100):
    pw, g, Ji = baue(N, J=N)
    for u in (mp.mpf("0.01"), mp.mpf("0.05"), mp.mpf("0.15")):
        zeta = L + u * mp.expj(psi)
        mp.mp.dps = 100
        P = loese(Ji, rhs_pol(N, pw, g, zeta))
        TNk = [zeta ** (-(k + 1)) for k in range(N + 1)]
        exakt = [P[k] - TNk[k] for k in range(N + 1)]
        mp.mp.dps = 50
        Tp = toeplitz_poly(N, zeta)
        d = max(abs(horner(exakt, z) - horner(Tp, z)) for z in kreis)
        m1 = max(abs(horner(exakt, z)) for z in kreis); m2 = max(abs(horner(Tp, z)) for z in kreis)
        print(f"N={N}, u={mp.nstr(u,2)}: max|J⁻¹T_N F_ζ| = {mp.nstr(m1,3)}, max|P_N T(σ)F_ζ| = {mp.nstr(m2,3)}, max|Differenz| = {mp.nstr(d,3)};"
              f"  bei z=0.5: {mp.nstr(horner(exakt, 0.5),4)} vs {mp.nstr(horner(Tp, 0.5),4)}", flush=True)
