"""b42 (Runde 4): Konvergiert x₁ = Π_N f_K auch auf dem Halbmond Ḡ'∖D̄ (äusserer Bogen von D(a₀,R₀))?
f_K dort exakt über f_K(z) = exp(f_K(Log z)) (Lemma 8.4(ii)), Log z ∈ D, Knesers Reihe (4000 Koeffizienten).
Vergleich: ε_N^D = sup_{|z|=|L|}|x₁ − f_K| und ε_N^G = sup_{äusserer Bogen ohne D(L,0.02)}|x₁ − f_K|, dazu die Stelle des Maximums.
Daten: gemeinsam/lin_PiN_fK_bis_1000.json (Strang A)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
d4 = json.load(open(G + "kneser_halbexp_koeff_4000.json"))
lin = json.load(open(G + "lin_PiN_fK_bis_1000.json"))
Ns = [int(n) for n in sys.argv[1:]] or [100, 150, 200, 250, 300, 400, 500]
mp.mp.dps = max(60, max(Ns) // 3 + 30)
Kc = [mp.mpf(c) for c in (d4["a"] if isinstance(d4, dict) else d4)]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
aL = abs(L); a0 = Kc[0]; R0 = abs(L - a0); phiL = mp.arg(L - a0)


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


bogen = []
for j in range(1, 200):
    ph = -phiL + 2 * phiL * j / 200
    z = a0 + R0 * mp.expj(ph)
    if abs(z) > aL and abs(z - L) > 0.02 and abs(z - mp.conj(L)) > 0.02:
        bogen.append((ph, z))
fG = [mp.exp(horner(Kc, mp.log(z))) for _, z in bogen]
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 240) for j in range(240)]
fD = [horner(Kc, z) for z in kreis]
for N in Ns:
    if str(N) not in lin:
        continue
    b = [mp.mpf(x) for x in lin[str(N)]["b"]]
    eD = max(abs(horner(b, z) - v) for z, v in zip(kreis, fD))
    werte = [(abs(horner(b, z) - v), ph) for (ph, z), v in zip(bogen, fG)]
    eG, phm = max(werte)
    rechts = abs(horner(b, a0 + R0) - mp.exp(horner(Kc, mp.log(a0 + R0))))
    print(f"N={N}: ε^D = {mp.nstr(eD, 3)}, ε^G (äusserer Bogen) = {mp.nstr(eG, 3)} bei φ = {mp.nstr(phm * 180 / mp.pi, 4)}°,"
          f" am rechten Ende a₀+R₀: {mp.nstr(rechts, 3)}", flush=True)
