"""X01 — Alternativweg „Positivität / Vorzeichen / Monotonie in N“: trägt irgendeine Vorzeichenstruktur?

Prüft (numerisch, Daten aus gemeinsam/):
 (a) Vorzeichen der Kneser-Koeffizienten a_k (erste negative Stelle, Anteil, Gewicht Σ_{a_k<0}|a_k|ρ^k);
 (b) Vorzeichen der Koeffizienten von p_N (gleich wie Kneser bis auf den Kopf?), von p_N − f_K und von p_{N+1} − p_N;
 (c) Monotonie in N: a_0^(N), p_N(1), max_{|x|=ρ}|p_N|, ℓ¹_ρ-Norm von p_N − q_N;
 (d) „Maximum auf der positiven Achse“: max_φ |p_N(r e^{iφ})| = p_N(r)?  (wäre für (i) aus README hilfreich);
 (e) Hankel-Minoren der Kneser-Koeffizienten (Totalpositivität à la Stieltjes?).
"""
import json
import os
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
mp.mp.dps = 60
HIER = os.path.dirname(os.path.abspath(__file__))
G = os.path.join(HIER, "..", "gemeinsam")
L = mp.mpc("0.31813150520476413531265425158766451720351761387139986692237",
           "1.33723570143068940890116214319371061253950213846051241687")
RHO = abs(L)

K = [mp.mpf(s) for s in json.load(open(os.path.join(G, "kneser_halbexp_koeff_4000.json")))["a"]]
P = {int(k): [mp.mpf(s) for s in v] for k, v in json.load(open(os.path.join(G, "polynome_bis_194.json"))).items()}
P.update({int(k): [mp.mpf(s) for s in v] for k, v in json.load(open(os.path.join(G, "polynome_195_bis_300.json"))).items()})
PG = {int(k): [mp.mpf(s) for s in v] for k, v in json.load(open(os.path.join(G, "pN_N100_bis_1000.json")))["p"].items()}


def horner(c, z):
    s = mp.mpc(0)
    for x in reversed(c):
        s = s * z + x
    return s


# (a) Kneser-Vorzeichen
neg = [k for k in range(len(K)) if K[k] < 0]
print("(a) Kneser: erste negative Koeffizienten k =", neg[:12])
for kmax in (50, 200, 1000, 4000):
    nn = sum(1 for k in neg if k <= kmax)
    w = sum(abs(K[k]) * RHO ** k for k in neg if k <= kmax)
    tot = sum(abs(K[k]) * RHO ** k for k in range(min(kmax + 1, len(K))))
    print(f"   k ≤ {kmax}: {nn} negative ({nn / (kmax + 1):.3f}), Gewicht Σ_neg |a_k|ρ^k = {mp.nstr(w, 5)} von {mp.nstr(tot, 6)}")
# Abstand der Vorzeichenwechsel (Periode 2π/arg L = 4.70 erwartet)
wech = [k for k in range(20, 400) if (K[k] > 0) != (K[k + 1] > 0)]
d = [b - a for a, b in zip(wech[:-1], wech[1:])]
print("   mittlerer Abstand der Vorzeichenwechsel (k 20…400):", mp.nstr(mp.mpf(sum(d)) / len(d), 4),
      "  π/arg L =", mp.nstr(mp.pi / mp.arg(L), 4))

# (b) Vorzeichen von p_N, p_N − f_K, p_{N+1} − p_N
print("\n(b) Vorzeichen")
for N in (60, 120, 194, 300):
    c = P[N]
    ab = [k for k in range(N + 1) if (c[k] > 0) != (K[k] > 0)]
    e = [c[k] - K[k] for k in range(N + 1)]
    se = "".join("+" if x > 0 else "-" for x in e[:60])
    dd = [P[N + 1][k] - (c[k] if k <= N else 0) for k in range(N + 2)] if N + 1 in P else None
    print(f"   N={N}: Vorzeichen p_N ≠ Kneser bei k = {ab[:10]}{'…' if len(ab) > 10 else ''} (Anzahl {len(ab)}, erstes {ab[0] if ab else '-'}; k/N = {ab[0] / N if ab else 0:.3f})")
    print(f"        sign(p_N − f_K)_k, k<60: {se}")
    if dd:
        sd = "".join("+" if x > 0 else "-" for x in dd[:60])
        print(f"        sign(p_N+1 − p_N)_k, k<60: {sd}")

# (c) Monotonie in N
print("\n(c) Monotonie in N")
T_halb = None
f1 = sum(K[:4001])  # f_K(1) aus Taylor (Fehler ~ (1/ρ)^4000 vernachlässigbar)
NS = sorted(P)
a0 = [P[N][0] - K[0] for N in NS]
p1 = [horner(P[N], 1).real - f1 for N in NS]
wa0 = [NS[i] for i in range(1, len(NS)) if (a0[i] > 0) != (a0[i - 1] > 0)]
wp1 = [NS[i] for i in range(1, len(NS)) if (p1[i] > 0) != (p1[i - 1] > 0)]
print("   Vorzeichenwechsel von a_0^(N) − a_0^K bei N =", wa0)
print("   Vorzeichenwechsel von p_N(1) − f_K(1)   bei N =", wp1)
mon = [NS[i] for i in range(1, len(NS)) if abs(a0[i]) > abs(a0[i - 1])]
print("   |a_0^(N) − a_0^K| nicht fallend bei N =", mon[:30], "…" if len(mon) > 30 else "")

# (d) Maximum auf dem Kreis |x| = r
print("\n(d) max_φ |p_N(r e^{iφ})| gegen p_N(r)")
M = 720
for N in (50, 120, 194, 300):
    c = P[N]
    out = []
    for fr in (0.5, 0.9, 1.0):
        r = fr * RHO
        vals = [abs(horner(c, r * mp.expjpi(2 * mp.mpf(j) / M))) for j in range(M)]
        jm = max(range(M), key=lambda j: vals[j])
        out.append(f"r={fr}ρ: max {mp.nstr(vals[jm], 8)} bei φ={360 * jm / M:.1f}°, p_N(r)={mp.nstr(horner(c, r).real, 8)}")
    print(f"   N={N}: " + " | ".join(out))
# auch für die grossen N
for N in (500, 1000):
    c = PG[N]
    r = RHO
    vals = [abs(horner(c, r * mp.expjpi(2 * mp.mpf(j) / M))) for j in range(M)]
    jm = max(range(M), key=lambda j: vals[j])
    print(f"   N={N}: r=ρ: max {mp.nstr(vals[jm], 8)} bei φ={360 * jm / M:.1f}°, p_N(ρ) = {mp.nstr(horner(c, r).real, 8)}")

# (e) Hankel-Minoren
print("\n(e) Hankel-Determinanten det[a_{i+j}]_{i,j<n} (Kneser) und verschobene det[a_{i+j+1}]")
for n in range(1, 9):
    H0 = mp.matrix(n, n)
    H1 = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            H0[i, j] = K[i + j]
            H1[i, j] = K[i + j + 1]
    print(f"   n={n}: {mp.nstr(mp.det(H0), 6)}   {mp.nstr(mp.det(H1), 6)}")
