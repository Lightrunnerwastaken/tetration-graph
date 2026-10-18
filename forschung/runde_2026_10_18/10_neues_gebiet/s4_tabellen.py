"""s4: (a) Glattheitsexponent von Knesers Erzeuger am Fixpunkt für verschiedene Basen b > η:
        J_K - J_L ~ (w-L)^{α},  α = 1 + 2πi/log λ,  λ = L ln b;   J_K ist genau C^{⌊Re α⌋} (nicht C^{⌊Re α⌋+1}), falls a_1 ≠ 0.
    (b) Resonanzsuche zwischen Fixpunkten von exp: λ_j = λ_k^m (λ = L, da exp' = exp), |j|,|k| ≤ 40, 2 ≤ |m| ≤ 12 und m = -1.
Aufruf: python s4_tabellen.py
"""
import mpmath as mp
mp.mp.dps = 30

print("(a) Basis b, L, Re α = 1 + 2π·Im(log λ)/|log λ|², Glattheit C^k")
for b in ["1.45", "1.5", "1.6", "2", "e", "3", "5", "10", "100", "1e6"]:
    bb = mp.e if b == "e" else mp.mpf(b)
    lb = mp.log(bb)
    L = -mp.lambertw(-lb) / lb
    if mp.im(L) < 0:
        L = mp.conj(L)
    lam = L * lb
    al = 1 + 2j * mp.pi / mp.log(lam)
    print(f"  b={b:>5}: L={mp.nstr(L, 10):>28}  |λ|={mp.nstr(abs(lam), 6):>8}  α={mp.nstr(al, 8):>24}  C^{int(mp.floor(mp.re(al)))}")

print("\n(b) Fixpunkte L_k = -W_k(-1) von exp; kleinster relativer Abstand |L_k^m - L_j|/|L_j| (j ≠ k)")
Ls = {k: -mp.lambertw(-1, k) for k in range(-40, 41)}
best = []
for k, Lk in Ls.items():
    for m in list(range(-12, -1)) + [-1] + list(range(2, 13)):
        v = Lk ** m
        for j, Lj in Ls.items():
            if j == k:
                continue
            best.append((abs(v - Lj) / abs(Lj), k, m, j))
best.sort()
for d, k, m, j in best[:6]:
    print(f"  L_{k}^{m} ≈ L_{j}: rel. Abstand {mp.nstr(d, 5)}")
