"""Prüft kern.zweig gegen gemeinsam/polynome_bis_194.json (b = e) und gegen bekannte Fälle."""
import json, sys, time
from pathlib import Path
import mpmath as mp
from kern import zweig, g_exp_basis, horner

mp.mp.dps = 60
HIER = Path(__file__).parent
ref = json.loads((HIER.parent / "gemeinsam" / "polynome_bis_194.json").read_text())
t0 = time.time()
loes = zweig(g_exp_basis(mp.e), 40, bits=400)
print("Zeit", round(time.time() - t0, 1))
fehler = max(abs(loes[N][k] - mp.mpf(ref[str(N)][k])) for N in loes for k in range(N + 1))
print("max Abweichung zu polynome_bis_194 (N<=40):", mp.nstr(fehler, 3))
# dritte Iterierte: N=1 Start 1/3 + x
t0 = time.time()
l3 = zweig(g_exp_basis(mp.e), 20, bits=300, ordnung=3)
print("Zeit3", round(time.time() - t0, 1))
for N in (1, 2, 3, 5, 10, 20):
    print(N, [mp.nstr(c, 12) for c in l3[N][:4]], "p(1)=", mp.nstr(horner(l3[N], mp.mpf(1)), 15))
