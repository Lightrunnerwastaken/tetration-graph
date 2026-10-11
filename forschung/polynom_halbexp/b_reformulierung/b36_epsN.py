"""b36 (Runde 3): ε_N := sup_{|z|=|L|} |Π_N f_K − f_K| aus gemeinsam/lin_PiN_fK_bis_1000.json und Knesers Reihe (4000 Koeff.),
dazu der Quotient sup|x₁ − p_N| / ε_N² mit den b33-Werten (x₁ = Π_N f_K)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 40
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
d4 = json.load(open(G + "kneser_halbexp_koeff_4000.json"))
Kc = [mp.mpf(c) for c in (d4["a"] if isinstance(d4, dict) else d4)]
lin = json.load(open(G + "lin_PiN_fK_bis_1000.json"))
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
aL = abs(L)
kreis = [aL * mp.expjpi(2 * mp.mpf(j) / 360) for j in range(360)] + [L]


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


fK = {i: horner(Kc, z) for i, z in enumerate(kreis)}
b33 = {150: 8.64e-22, 200: 6.55e-23, 232: 1.78e-23, 300: 1.83e-24}
for N in (150, 200, 225, 250, 300):
    if str(N) not in lin:
        continue
    b = [mp.mpf(x) for x in lin[str(N)]["b"]]
    eps = max(abs(horner(b, z) - fK[i]) for i, z in enumerate(kreis))
    zus = f", sup|x₁ − p_N|/ε_N² = {b33[N] / float(eps)**2:.3f}" if N in b33 else ""
    print(f"N={N}: ε_N = {mp.nstr(eps, 4)}{zus}", flush=True)
# N = 232 liegt nicht in der Datei: interpoliere log ε linear in log N zwischen 225 und 250 (nur zur Orientierung)
