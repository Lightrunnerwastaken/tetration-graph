"""b01: Struktur des Fehlers d_N = p_N - q_N (q_N = Knesers Taylor-Polynom) und von E_N = p_N∘p_N - exp.
Nur Daten aus gemeinsam/ (kein Kneser-Motor)."""
import json, sys
import mpmath as mp
sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
G = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/gemeinsam/"
P = {int(k): [mp.mpf(c) for c in v] for k, v in json.load(open(G + "polynome_bis_194.json")).items()}
K = [mp.mpf(c) for c in json.load(open(G + "kneser_halbexp_koeff_80.json"))]
L = -mp.lambertw(-1); L = mp.conj(L) if mp.im(L) < 0 else L
aL = abs(L)

def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s

print("|L| =", mp.nstr(aL, 12))
# 1. Koeffizientenprofil von d_N, skaliert mit |L|^k
for N in (20, 40, 60, 70):
    d = [P[N][k] - K[k] for k in range(N + 1)]
    prof = [mp.nstr(abs(d[k]) * aL ** k, 3) for k in range(0, N + 1, max(1, N // 10))]
    print(f"N={N}: |d_k| |L|^k für k=0,{N//10},..:", prof)
    print("     letzte 4 k:", [mp.nstr(d[k] * aL ** k, 4) for k in range(N - 3, N + 1)])
# Kneser-Koeffizienten skaliert
print("Kneser a_k |L|^k k^6.45 (k=10..80 step 10):", [mp.nstr(K[k] * aL ** k * mp.mpf(k) ** 6.45, 4) for k in range(10, 81, 10)])

# 2. Abstand p_N - f_K auf Kreisen |x| = r (f_K über die 81 Koeffizienten; gut für r <= ~1.1)
for r in (0.5, 0.8, 1.0, 1.1):
    row = []
    for N in (10, 20, 30, 40, 50, 60):
        m = 0
        for j in range(64):
            z = r * mp.expjpi(2 * mp.mpf(j) / 64)
            m = max(m, abs(horner(P[N], z) - horner(K, z)))
        row.append(mp.nstr(m, 3))
    print(f"r={r}: max|p_N - f_K| N=10..60:", row)

# 3. E_N = p∘p - exp auf Kreisen: Wo ist die abgeschnittene Gleichung gut erfüllt?
for N in (20, 40, 80):
    row = []
    for r in (0.3, 0.5, 0.7, 0.9, 1.1, 1.3):
        m = 0
        for j in range(64):
            z = r * mp.expjpi(2 * mp.mpf(j) / 64)
            m = max(m, abs(horner(P[N], horner(P[N], z)) - mp.exp(z)))
        row.append(mp.nstr(m, 3))
    print(f"N={N}: max|E_N| auf r=0.3..1.3:", row)
