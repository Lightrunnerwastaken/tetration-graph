"""b13: Inneres Profil bei L: e_N(z) = p_N(z) − f_K(z) für z = L(1 − t/N) (radial nach innen), t = 0 … 50.
Skalierung prüfen: N^{β}·|e_N| mit β = 5.447 (Skala 1/N) bzw. β/θ = 4.876 (Ecke von G')."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 50
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_1000.json"))["a"]]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s
ts = [0, 0.5, 1, 2, 4, 8, 16, 32]
print("|e_N(L(1−t/N))|·N^4.876 :")
for N in (40, 80, 120, 160, 194):
    row = []
    for t in ts:
        z = L * (1 - mp.mpf(t) / N)
        e = horner(P[N], z) - horner(K, z)
        row.append(mp.nstr(abs(e) * mp.mpf(N) ** 4.876, 3))
    print(f"N={N}:", row)
print("\n|e_N| bei festem Abstand r = |z − L| (radial) — Abklingexponent in N:")
for r in (0.002, 0.01, 0.05, 0.2):
    vals = []
    for N in (80, 120, 160, 194):
        z = L * (1 - mp.mpf(r) / abs(L))
        vals.append(abs(horner(P[N], z) - horner(K, z)))
    ex = mp.log(vals[-1] / vals[0]) / mp.log(mp.mpf(194) / 80)
    print(f"r={r}: |e_N| (N=80,120,160,194) = {[mp.nstr(v,3) for v in vals]}, Exponent ≈ {mp.nstr(ex, 3)}")
