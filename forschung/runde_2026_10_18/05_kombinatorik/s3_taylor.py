"""Taylor-Koeffizienten von Knesers sexp_e: Singularitäten-Kaskade bei -2, -3, ...

Behauptung (Beweis im Bericht): D(z) = sexp(z) - log(z+2) ist holomorph auf C \\ (-inf, -3].
Also  c_n = (-1)^{n+1}/(n 2^n) + d_n  mit  limsup |d_n|^{1/n} = 1/3, und lokal bei -3
   D ~ log(ell - K'),  ell = log(1/(1-zeta)),  zeta = -z/3,  K' = ln(3 A),  A = sexp'(0) (= sexp'(-1)).
Transfer: d_n ~ (-1/3)^n (1/n) sum_{k>=1} (1/Gamma)^{(k)}(0)/k! * G^{(k)}(ln n),  G(s) = log(s - K').
Daten (nur lesend): etappe3/daten/taylor_130_256.json und taylor_100_128.json (unabhängige Einstellungen).
Aufruf: python s3_taylor.py   Ausgabe daten/s3_taylor.json
"""
import json
from pathlib import Path
import mpmath as mp

mp.mp.dps = 140
HIER = Path(__file__).resolve().parent
E3 = HIER.parents[1] / "identitaeten_2026_10_14" / "etappe3" / "daten"
a = json.load(open(E3 / "taylor_130_256.json"))["zentren"]
b = json.load(open(E3 / "taylor_100_128.json"))["zentren"]
c = [mp.mpf(s) for s in a["0.0"]]
c2 = [mp.mpf(s) for s in b["0.0"]]
A = c[1]
Kp = mp.log(3 * A)
print("A = sexp'(0) =", mp.nstr(A, 25), "  K' = ln(3A) =", mp.nstr(Kp, 20))
# (1/Gamma)^{(k)}(0)/k!
rg = mp.taylor(lambda s: mp.rgamma(s), 0, 8)
out = {"A": str(A), "Kp": str(Kp), "zeilen": []}
print(" n   c_n                       sign  d_n*(-3)^n*n     Transfer(5 Gl.)   Abw.Einst.(d_n rel)")
for n in range(1, 101):
    lead = (-1) ** (n + 1) / (n * mp.mpf(2) ** n)
    d = c[n] - lead
    s = mp.log(n) - Kp
    # G(s)=log(s): G^{(k)}(s)/... = (-1)^{k-1}(k-1)!/s^k
    tr = sum(rg[k] * (-1) ** (k - 1) * mp.factorial(k - 1) / s ** k for k in range(1, 6))
    norm = d * (-3) ** n * n
    rel = None
    if n < len(c2) - 8:
        rel = abs((c2[n] - lead) - d) / abs(d)
    out["zeilen"].append({"n": n, "c": mp.nstr(c[n], 20), "d_norm": mp.nstr(norm, 15), "transfer": mp.nstr(tr, 15),
                          "rel_einst": None if rel is None else mp.nstr(rel, 3)})
    if n <= 12 or n % 10 == 0:
        print(f"{n:3d} {mp.nstr(c[n], 18):>26} {'+' if c[n] > 0 else '-'}   {mp.nstr(norm, 10):>14} {mp.nstr(tr, 10):>14}  {'' if rel is None else mp.nstr(rel, 2)}")
# Vorzeichenmuster
sg = ''.join('+' if x > 0 else '-' for x in c[1:101])
print("Vorzeichen c_1..c_100:", sg)
alt = [n for n in range(1, 101) if (c[n] > 0) != ((-1) ** (n + 1) > 0)]
print("n mit sign(c_n) != (-1)^{n+1}:", alt)
out["vorzeichen"] = sg
out["abweichend_von_alternierend"] = alt
# Wurzeltest |d_n|^{1/n}
for n in (20, 40, 60, 80, 100):
    d = c[n] - (-1) ** (n + 1) / (n * mp.mpf(2) ** n)
    print(f"|d_{n}|^(1/n) = {mp.nstr(abs(d) ** (mp.mpf(1) / n), 10)}   (1/3 = 0.3333333333)")
# Kontrolle Zentrum -0.4: Radius 1.6, Rest Radius 2.6
cm = [mp.mpf(s) for s in a["-0.4"]]
for n in (20, 40, 60, 80):
    d = cm[n] - (-1) ** (n + 1) / (n * mp.mpf("1.6") ** n)
    print(f"Zentrum -0.4: |d_{n}|^(1/n) = {mp.nstr(abs(d) ** (mp.mpf(1) / n), 10)}  (1/2.6 = {mp.nstr(1/mp.mpf('2.6'), 10)}),"
          f"  d*(-2.6)^n*n = {mp.nstr(d * (-mp.mpf('2.6')) ** n * n, 8)},  Transfer: "
          f"{mp.nstr(sum(rg[k] * (-1) ** (k - 1) * mp.factorial(k - 1) / (mp.log(n) - mp.log(mp.mpf('2.6') * A)) ** k for k in range(1, 6)), 8)}")
json.dump(out, open(HIER / "daten" / "s3_taylor.json", "w"), indent=1)
