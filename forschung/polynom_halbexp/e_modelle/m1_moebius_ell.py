"""Modell 1: elliptische Möbius-Abbildung g(x) = tan(arctan x + c) (Drehung der projektiven Geraden, Fixpunkte ±i, |g'(±i)| = 1).

Möbius-Wurzel f_M(x) = tan(arctan x + c/2) = (x + τ)/(1 − τx), τ = tan(c/2); Pol bei 1/τ, f_M(0) = τ.
Andere Wurzeln: α⁻¹(α + 1/2 + θ(α)), α = arctan(x)/c, singulär bei ±i (Stärke π/c).
Frage: konvergiert p_N gegen f_M? In welchem Gebiet, mit welcher Rate?
"""
import json, sys, time
from pathlib import Path
import mpmath as mp
from kern import zweig, g_moebius_ell, horner

mp.mp.dps = 60
HIER = Path(__file__).parent
(HIER / "daten").mkdir(exist_ok=True)
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 80
CS = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.4, 0.8, 1.2, 1.5]
erg = {}
for c in CS:
    t0 = time.time()
    tau = mp.tan(mp.mpf(c) / 2)
    fM = lambda z: (z + tau) / (1 - tau * z)
    fk = [tau] + [(1 + tau ** 2) * tau ** (k - 1) for k in range(1, NMAX + 1)]
    try:
        loes = zweig(g_moebius_ell(mp.mpf(c)), NMAX, bits=600)
    except RuntimeError as e:
        print("c =", c, "FEHLER", e)
        continue
    zeilen = {}
    for N in sorted(loes):
        a = loes[N]
        d1 = horner(a, mp.mpf("0.5")) - fM(mp.mpf("0.5"))
        dm = horner(a, mp.mpf("-0.5")) - fM(mp.mpf("-0.5"))
        di = abs(horner(a, mp.mpc(0, "0.9")) - fM(mp.mpc(0, "0.9")))
        dk = max(abs(a[k] - fk[k]) for k in range(min(6, N + 1)))
        dsum = mp.fsum(abs(a[k] - fk[k]) for k in range(N + 1))
        zeilen[N] = {"x=0.5": d1, "x=-0.5": dm, "x=0.9i": di, "max_k<=5": dk, "l1": dsum, "a0": a[0] - tau}
    erg[c] = zeilen
    print(f"\n== c = {c}: tau = {mp.nstr(tau, 6)}, Pol 1/tau = {mp.nstr(1 / tau, 6)}, Zeit {time.time() - t0:.0f}s")
    for N in [1, 2, 3, 5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150]:
        if N in zeilen:
            z = zeilen[N]
            print(f"N={N:4d}  a0-tau={mp.nstr(z['a0'], 3):>10}  x=.5:{mp.nstr(z['x=0.5'], 3):>10}  x=-.5:{mp.nstr(z['x=-0.5'], 3):>10}"
                  f"  x=.9i:{mp.nstr(z['x=0.9i'], 3):>10}  max|da_k|,k<=5:{mp.nstr(z['max_k<=5'], 3):>10}")
json.dump({str(c): {str(N): {k: mp.nstr(v, 8) for k, v in z.items()} for N, z in zz.items()} for c, zz in erg.items()},
          open(HIER / "daten" / f"m1_moebius_ell_N{NMAX}.json", "w"), indent=0)
