"""Test der Kette und exakte Möbius-Formel.

Für φ = Möbius-Kreisdiffeo (Parameter a, |a| < 1) ist θ∘φ Randwert einer in der Kreisscheibe holomorphen Funktion,
also hat θ_T genau die negativen Moden von p:  θ̂_T(−n) = i aⁿ/(2πn)  und  N(T) = Σ|θ̂_T(−n)|² = Li₂(|a|²)/(4π²),
Dirichlet-Fläche-Defekt E(T) = 2 Σ n|θ̂_T(−n)|² = −ln(1−|a|²)/(2π²).
Aufruf: python s1_moebius.py DPS K
"""
import sys
import time

import mpmath as mp

import e4_lib as E

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
k = int(sys.argv[2]) if len(sys.argv) > 2 else 5
E.init(dps)
MOD = 8


def zeige(name, r):
    print(f"--- {name}")
    for key in ("a", "b", "c", "d"):
        print(f"  {key}: {mp.nstr(r[key], 12)}   |.|={mp.nstr(abs(r[key]), 3)}")
    for j in range(2, 5):
        print(f"  V{j}: |.|={mp.nstr(abs(r[f'V{j}']), 3)}")


t = time.time()
r0 = E.kenngroessen(("trig", []), k, MOD)
print(f"Kneser (p=0), k={k}, dps={dps}, {time.time()-t:.1f}s")
zeige("Kneser", r0)
print("  neg. Moden |θ̂(−m)|:", [mp.nstr(abs(r0[f'neg{m}']), 3) for m in range(1, MOD + 1)])

for a in (mp.mpf("0.05"), mp.mpf("0.3"), mp.mpc("0.2", "0.25")):
    t = time.time()
    r = E.kenngroessen(("moeb", a), k, MOD)
    zeige(f"Möbius a={a}  ({time.time()-t:.1f}s)", r)
    fehl = 0
    for n in range(1, MOD + 1):
        soll = 1j * a ** n / (2 * mp.pi * n)
        fehl = max(fehl, abs(r[f"neg{n}"] - soll))
    N = mp.fsum(abs(r[f"neg{n}"]) ** 2 for n in range(1, MOD + 1))
    print(f"  max_n≤{MOD} |θ̂_T(−n) − i aⁿ/(2πn)| = {mp.nstr(fehl, 3)}")
    print(f"  θ̂_T(−1) = {mp.nstr(r['neg1'], 15)}  soll {mp.nstr(1j*a/(2*mp.pi), 15)}")
    print(f"  N (n≤{MOD}) = {mp.nstr(N, 15)},  Li2(|a|²)/(4π²) = {mp.nstr(mp.polylog(2, abs(a)**2)/(4*mp.pi**2), 15)}")
