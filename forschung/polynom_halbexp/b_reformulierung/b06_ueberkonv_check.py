"""b06: Überkonvergenz prüfen: p_N gegen f_K bei Punkten ausserhalb D; q_N (Taylor) zum Vergleich (Kneser bis k=1000)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
KK = json.load(open(G + "kneser_halbexp_koeff_1000.json"))
K = [mp.mpf(c) for c in KK["a"]]
fe = [mp.mpf(c) for c in KK["fehler_abs"]]
print("Fehler der Kneser-Koeffizienten k=100,300,600,1000:", [mp.nstr(fe[k], 3) for k in (100, 300, 600, 1000)])
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s
def fK(z):
    if abs(z) <= 1.05:
        return horner(K[:400], z)
    return mp.exp(fK(mp.log(z)))
# Taylorteilsummen und p_N bei x=1.6 (reell) und z = 1.5 i
for z in (mp.mpf("1.6"), mp.mpf("1.45"), mp.mpc(0, 1.5), mp.mpc(-1.5, 0)):
    ref = fK(z) if mp.re(z) > 0 else None
    if ref is None:
        ref = horner(P[194], z)
    print("z =", mp.nstr(z, 4), " Referenz", mp.nstr(ref, 12))
    for N in (40, 80, 120, 194):
        print(f"   N={N}: |p_N - ref| = {mp.nstr(abs(horner(P[N], z) - ref), 3)}   |q_N - ref| = {mp.nstr(abs(horner(K[:N+1], z) - ref), 3)}")
    best = min((abs(horner(K[:n + 1], z) - ref), n) for n in range(20, 1000, 20))
    print("   beste Taylorteilsumme bis 1000:", mp.nstr(best[0], 3), "bei n =", best[1])
