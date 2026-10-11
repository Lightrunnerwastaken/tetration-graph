"""Export: Polynom-Wurzeln p_N des Möbius-Drehungsmodells als Prüfstein für Beweisgerüste (→ gemeinsam/modell_moebius_*.json).
g(x) = A + B·tan(arctan((x−A)/B) + c), kanonische Wurzel f_M(x) = A + B·tan(arctan((x−A)/B) + c/2).
Fälle: (A,B) = (0,1), c = 0.8 (konvergent, Rate 0.559) und c = 1.0 (Rate 0.887); (A,B) = L, c = 0.70510 (g(0) = 1, Rate 0.388)."""
import json, sys
from pathlib import Path
import mpmath as mp
from flint import arb
from kern import zweig

mp.mp.dps = 80
HIER = Path(__file__).parent
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 80
LR, LI = mp.mpf("0.31813150520476413"), mp.mpf("1.3372357014306895")
faelle = {"drehung_c0.8": (mp.mpf(0), mp.mpf(1), mp.mpf("0.8")),
          "drehung_c1.0": (mp.mpf(0), mp.mpf(1), mp.mpf("1.0")),
          "moebius_exp": (LR, LI, mp.atan((1 - LR) / LI) + mp.atan(LR / LI))}
aus = {"beschreibung": __doc__, "faelle": {}}
for name, (A, B, c) in faelle.items():
    g = lambda x: A + B * mp.tan(mp.atan((x - A) / B) + c)
    fM = lambda x: A + B * mp.tan(mp.atan((x - A) / B) + c / 2)
    gk_mp = mp.taylor(g, 0, NMAX)
    fk = mp.taylor(fM, 0, NMAX)
    loes = zweig(lambda N: [arb(mp.nstr(x, 78)) for x in gk_mp[:N + 1]], NMAX, bits=1100)
    aus["faelle"][name] = {"A": mp.nstr(A, 20), "B": mp.nstr(B, 20), "c": mp.nstr(c, 20),
                           "g_koeff": [mp.nstr(x, 60) for x in gk_mp], "fM_koeff": [mp.nstr(x, 60) for x in fk],
                           "p_N": {str(N): [mp.nstr(x, 60) for x in loes[N]] for N in loes}}
    print(name, "fertig", flush=True)
json.dump(aus, open(HIER.parent / "gemeinsam" / f"modell_moebius_N{NMAX}.json", "w", encoding="utf-8"), indent=0, ensure_ascii=False)
