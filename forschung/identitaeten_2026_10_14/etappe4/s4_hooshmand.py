"""Lineare Tetration (Hooshmand) und Abklingen der negativen Moden.

(1) Lineare Tetration T(x) = e^x auf [0,1]: θ_T = 1 + h(x) − x, 1+θ_T' = h'(x). Defekte der Identitäten, Momente, negative Moden.
(2) Lemma R (Bericht): für reell-analytisches T ist der Anteil negativer Frequenz von θ_T reell-analytisch, also klingen
    θ̂_T(−m) exponentiell ab; für die nur C¹-glatte lineare Tetration dagegen algebraisch.
    Vergleich: T = sexp(x + t sin 4πx) mit t = 0.03.
Aufruf: python s4_hooshmand.py DPS K MMAX
"""
import sys
import time

import mpmath as mp

import e4_lib as E

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
k = int(sys.argv[2]) if len(sys.argv) > 2 else 6
MM = int(sys.argv[3]) if len(sys.argv) > 3 else 24
E.init(dps)

for name, stoer in (("linear (Hooshmand)", ("lin", None)),
                    ("sexp(x + 0.03 sin 4πx)", ("trig", [(2, mp.mpf("0.03"), 0)])),
                    ("Kneser", ("trig", []))):
    t = time.time()
    r = E.kenngroessen(stoer, k, moden=MM, momente=4)
    print(f"=== {name}  ({time.time()-t:.1f}s, dps {dps}, k {k})")
    print(f"   M[θ_T] = {mp.nstr(r['Mtheta'], 15)}")
    for c, txt in (("a", "M[θ²]−M[θ]²"), ("b", "σ-Produkt − 1"), ("c", "∫slog'²/h' − 1 = M[1/(1+θ')]−1"),
                   ("d", "M[log(T'h'(T))]"), ("V3", "M[(θ−Mθ)³]"), ("V4", "M[(θ−Mθ)⁴]")):
        print(f"   {txt:34s} = {mp.nstr(r[c], 12)}")
    print("   |θ̂_T(−m)|, m = 1..:", " ".join(mp.nstr(abs(r[f'neg{m}']), 3) for m in range(1, MM + 1)))
    N = mp.fsum(abs(r[f'neg{m}']) ** 2 for m in range(1, MM + 1))
    Ed = 2 * mp.fsum(m * abs(r[f'neg{m}']) ** 2 for m in range(1, MM + 1))
    print(f"   N_{MM} = {mp.nstr(N, 8)},  E_{MM} = 2Σ m|θ̂(−m)|² = {mp.nstr(Ed, 8)}")
    if name.startswith("linear"):
        ms = [4, 8, 16, MM]
        for m in ms:
            print(f"     m={m:3d}: m²|θ̂(−m)| = {mp.nstr(m**2*abs(r[f'neg{m}']), 5)}  m³|θ̂(−m)| = {mp.nstr(m**3*abs(r[f'neg{m}']), 5)}")
