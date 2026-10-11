"""b23 (Runde 2): Test der Aussage X (Pol-Antwort).  E_N(ζ, z) := Π_N k_ζ(z) − k_ζ(z),  M_N(ζ) := max_{|z|=|L|} |E_N(ζ,z)|
(= sup über D̄, Maximumprinzip, da k_ζ für ζ ∉ D̄ auf D̄ holomorph ist).
(a) ζ auf dem Schnitt ζ = L + u e^{iψ}, ψ = Winkelhalbierende der Aussenecke von G' (87.15°);  (b) ζ auf einer äusseren Kontur.
(c) Schranke I_N := (1/π) ∫ |Sprung(ζ_u)| M_N(ζ_u) du (L und L̄) gegen den wahren Wert sup_D |Π_N f_K − f_K| (b14)."""
import json, sys, math
import mpmath as mp
from flint import arb, arb_poly, arb_mat, ctx
sys.stdout.reconfigure(encoding="utf-8")
exec(open("b22_polmoden.py", encoding="utf-8").read().split("mp.mp.dps = 100\nres = {}")[0])
mp.mp.dps = 100
aL = abs(L)
a0 = mp.mpf(K[0])
psi = (mp.arg(L) + mp.pi / 2 + mp.arg(L - a0) - mp.pi / 2) / 2 + 0     # Mitte von (arg L − 90°, arg(L−a0) + 90°)? s. u.
# Aussenkeil: Richtungen (arg(L−a0) − 90°, arg L + 90°) = (7.68°, 166.62°); Halbierende:
psi = ((mp.arg(L - a0) - mp.pi / 2) + (mp.arg(L) + mp.pi / 2)) / 2
print("ψ =", mp.nstr(mp.degrees(psi), 6), "°")
NZ = 720
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / NZ) for j in range(NZ)]
# Punkte nahe L dichter
kreis += [L * mp.expjpi(mp.mpf(t) / 1800) for t in range(-60, 61)]
kreis += [mp.conj(L) * mp.expjpi(mp.mpf(t) / 1800) for t in range(-60, 61)]

D1 = mp.mpc("-0.032003287962208", "-0.011244050502060")    # Strang A: f_K ∋ D1 (1 − z/L)^{1+γ} (Hauptzweig, radialer Schnitt)
def sprung(u):
    # Schnitt entlang ψ: stetige Fortsetzung von D aus; arg(1 − z/L) = π + arg t − arg L mit arg t = ψ bzw. ψ − 2π.
    t_arg1 = psi; t_arg2 = psi - 2 * mp.pi
    w1 = mp.log(u / aL) + 1j * (mp.pi + t_arg1 - mp.arg(L))
    w2 = mp.log(u / aL) + 1j * (mp.pi + t_arg2 - mp.arg(L))
    return abs(D1 * (mp.exp((1 + gam) * w1) - mp.exp((1 + gam) * w2)))
us = [mp.mpf(10) ** (-e / 6) for e in range(21, 2, -1)]     # 3e-4 … 0.56
erg = {}
for N in (60, 100, 140, 180):
    pw, g, Ji = baue(N, J=N)
    Ms = {}
    for u in us:
        zeta = L + u * mp.expj(psi)
        P = loese(Ji, rhs_pol(N, pw, g, zeta))
        Ms[u] = max(abs(horner(P, z) - 1 / (zeta - z)) for z in kreis)
    # äussere Kontur: Rand von G' um 0.2 vergrössert (ohne die Schnittstücke) — Kreisbögen
    aussen = []
    for j in range(48):
        phi = 2 * mp.pi * j / 48
        z1 = (aL + mp.mpf("0.2")) * mp.expj(phi); z2 = a0 + (abs(L - a0) + mp.mpf("0.2")) * mp.expj(phi)
        aussen.append(z1 if abs(z1 - a0) > abs(L - a0) + mp.mpf("0.2") else z2)
    Mout = []
    for zeta in aussen:
        P = loese(Ji, rhs_pol(N, pw, g, zeta))
        Mout.append(max(abs(horner(P, z) - 1 / (zeta - z)) for z in kreis[::4]))
    # Schranke über den Schnitt (Trapez in log u), beide Schnitte (L, L̄) ⇒ Faktor 2, Vorfaktor 1/(2π)
    lu = [mp.log(u) for u in us]
    f_int = [sprung(u) * Ms[u] * u for u in us]
    I = sum((f_int[i] + f_int[i + 1]) / 2 * (lu[i + 1] - lu[i]) for i in range(len(us) - 1)) * 2 / (2 * mp.pi)
    erg[N] = (Ms, max(Mout), I)
    print(f"N={N}: u·M_N(L+u e^iψ) für u = " + ", ".join(f"{mp.nstr(u,2)}:{mp.nstr(u*Ms[u],3)}" for u in us), flush=True)
    print(f"      max M_N auf äusserer Kontur (Abstand 0.2 von G') = {mp.nstr(max(Mout), 3)};  Schnitt-Schranke I_N = {mp.nstr(I, 3)}", flush=True)

json.dump({N: {"M": {str(u): str(m) for u, m in e[0].items()}, "Mout": str(e[1]), "I": str(e[2])} for N, e in erg.items()}, open("b23_polfehler.json", "w"))
