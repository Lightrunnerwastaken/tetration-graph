"""Strang A, Aufgabe 1: Taylor-Koeffizienten von f_K(z) = T_e(slog_e(z) + 1/2) und von A(z) = slog_e(z) um 0.

Verfahren: Werte auf |z| = r (M Punkte, obere Hälfte gerechnet, untere durch Spiegelung, f_K und A sind reell),
Motor = fatou.gp-Fork (gp-tetration) bei dps 57 (bei dps ≥ 59 versagt der Fork für komplexe sexp-Argumente; gemessen
10.10.2026). DFT mit arb/acb (python-flint). Drei Radien zum Fehlervergleich. Punkte in Paketen zu 512.
Ausgabe: daten/kreis_r*.json (Rohwerte A, f), gemeinsam/kneser_halbexp_koeff_300.json, …_1000.json, kneser_slog_koeff_300.json.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, arb, ctx  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
DATEN = HIER / "daten"
GEM = HIER.parent / "gemeinsam"
DATEN.mkdir(exist_ok=True)
mp.mp.dps = 75
ctx.prec = 260
DPS = 57
M = 4096
RADIEN = ["1.34", "1.30", "1.25", "1.20"]
PAKET = 512
L = -mp.lambertw(-1)
L = mp.conj(L) if mp.im(L) < 0 else L


def werte(r_str, gp):
    """A(z_j), f(z_j) und Prüfrest T(A+1) − e^z für j = 0..M/2 (obere Hälfte), aus Datei falls vorhanden."""
    pfad = DATEN / f"kreis_r{r_str}_M{M}.json"
    if pfad.is_file():
        d = json.loads(pfad.read_text(encoding="utf-8"))
        return [mp.mpc(*x) for x in d["A"]], [mp.mpc(*x) for x in d["f"]], mp.mpf(d["pruefrest"])
    r = mp.mpf(r_str)
    zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M // 2 + 1)]
    zs[0], zs[-1] = mp.mpc(r, 0), mp.mpc(-r, 0)
    s = gp.session("exp(1)")
    A, f, rest = [], [], mp.mpf(0)
    for i in range(0, len(zs), PAKET):
        teil = zs[i:i + PAKET]
        a = s.slog_batch(teil)
        a = [mp.mpc(mp.re(x), 0) if mp.im(z) == 0 else x for x, z in zip(a, teil)]
        w = s.sexp_batch([x + mp.mpf("0.5") for x in a])
        w = [mp.mpc(mp.re(x), 0) if mp.im(z) == 0 else x for x, z in zip(w, teil)]
        e1 = s.sexp_batch([x + 1 for x in a])
        rest = max(rest, max(abs(u - mp.exp(z)) for u, z in zip(e1, teil)))
        A += a
        f += w
        print(f"  r = {r_str}: {i + len(teil)}/{len(zs)} Punkte, Prüfrest bisher {mp.nstr(rest, 3)}", flush=True)
    pfad.write_text(json.dumps({"r": r_str, "M": M, "dps_motor": DPS, "punkte": "z_j = r·exp(2πij/M), j = 0..M/2",
                                "pruefrest": mp.nstr(rest, 5),
                                "A": [[mp.nstr(mp.re(x), 62), mp.nstr(mp.im(x), 62)] for x in A],
                                "f": [[mp.nstr(mp.re(x), 62), mp.nstr(mp.im(x), 62)] for x in f]}), encoding="utf-8")
    return A, f, rest


def koeff(halb, r_str):
    """Taylor-Koeffizienten aus Werten der oberen Hälfte (reelle Funktion): c_k = DFT/M, a_k = c_k / r^k."""
    voll = halb + [mp.conj(x) for x in halb[1:-1][::-1]]
    assert len(voll) == M
    v = [acb(arb(mp.nstr(mp.re(x), 70)), arb(mp.nstr(mp.im(x), 70))) for x in voll]
    d = acb.dft(v)
    r = mp.mpf(r_str)
    return [mp.mpf(d[k].real.mid().str(70, radius=False)) / M / r ** k for k in range(M // 2)]


t0 = time.time()
gp = FatouGP(dps=DPS, fatou_gp="fork")
ergebnis = {}
for r_str in RADIEN:
    A, f, rest = werte(r_str, gp)
    ergebnis[r_str] = {"a": koeff(f, r_str), "A": koeff(A, r_str), "rest": rest}
    print(f"r = {r_str}: Prüfrest max|T(A+1) − e^z| = {mp.nstr(rest, 3)}  ({time.time() - t0:.0f} s)", flush=True)
gp.close()

# Fehlerschätzung: Abstand der Radien 1.34 und 1.30 (beide genauer als 1.25) und Motorfehler ε/r^k
a1, a2, a3 = (ergebnis[r]["a"] for r in RADIEN[:3])
A1, A2 = ergebnis[RADIEN[0]]["A"], ergebnis[RADIEN[1]]["A"]
eps = max(ergebnis[r]["rest"] for r in RADIEN)
alt = [mp.mpf(x) for x in json.loads((GEM / "kneser_halbexp_koeff_80.json").read_text(encoding="utf-8"))]
KMAX = 1200
fehler, rel = [], []
for k in range(KMAX + 1):
    fe = max(abs(a1[k] - a2[k]), 10 * eps / mp.mpf(RADIEN[0]) ** k)   # Motorfehler: grob ε·r^−k (×10 Sicherheit)
    fehler.append(fe)
    rel.append(fe / abs(a1[k]) if a1[k] != 0 else mp.inf)
print("\n k     a_k (r=1.34)                       |a(1.34)−a(1.30)|  |a(1.30)−a(1.25)|  rel.Fehler   |a_k|·|L|^k·(k+1)^6.45")
for k in list(range(0, 11)) + [20, 40, 60, 80, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1200]:
    print(f"{k:5d}  {mp.nstr(a1[k], 30):>38}  {mp.nstr(abs(a1[k] - a2[k]), 3):>10}  {mp.nstr(abs(a2[k] - a3[k]), 3):>10}  "
          f"{mp.nstr(rel[k], 3):>10}  {mp.nstr(abs(a1[k]) * abs(L) ** k * (k + 1) ** mp.mpf('6.447'), 5)}")
d80 = max(abs(a1[k] - alt[k]) / abs(alt[k]) for k in range(81))
print(f"\nVergleich mit kneser_halbexp_koeff_80.json: max relativer Abstand k ≤ 80: {mp.nstr(d80, 3)}")
print(f"Prüfsummen: a_0 − T(−1/2)-Wert alt: {mp.nstr(a1[0] - alt[0], 3)};  Σ_{{k≤1200}} a_k = {mp.nstr(mp.fsum(a1[:KMAX + 1]), 40)}")

beschr = ("Taylor-Koeffizienten a_k von f_K(z) = T_e(slog_e(z) + 1/2) um z = 0 (Kneser, Basis e). Cauchy-DFT auf |z| = 1.34 "
          "(M = 4096, Motor fatou.gp-Fork bei dps 57, Motorfehler ≤ " + mp.nstr(eps, 2) + "), Kontrolle mit |z| = 1.30 und 1.25. "
          "fehler_abs = max(|a(1.34) − a(1.30)|, 10·ε·1.34^−k), rel = fehler_abs/|a_k|. Strang A, " + time.strftime("%d.%m.%Y %H:%M"))
for kmax, name in ((300, "kneser_halbexp_koeff_300.json"), (1000, "kneser_halbexp_koeff_1000.json")):
    (GEM / name).write_text(json.dumps({
        "beschreibung": beschr, "kmax": kmax,
        "a": [mp.nstr(a1[k], 50) for k in range(kmax + 1)],
        "fehler_abs": [mp.nstr(fehler[k], 3) for k in range(kmax + 1)],
        "rel_fehler": [mp.nstr(rel[k], 3) for k in range(kmax + 1)]}, indent=0), encoding="utf-8")
fA = [max(abs(A1[k] - A2[k]), 10 * eps / mp.mpf(RADIEN[0]) ** k) for k in range(301)]
(GEM / "kneser_slog_koeff_300.json").write_text(json.dumps({
    "beschreibung": "Taylor-Koeffizienten von A(z) = slog_e(z) (Kneser, A(0) = −1) um 0, gleiche Rechnung wie kneser_halbexp_koeff_300.json.",
    "A": [mp.nstr(A1[k], 50) for k in range(301)], "fehler_abs": [mp.nstr(x, 3) for x in fA]}, indent=0), encoding="utf-8")
print(f"geschrieben ({time.time() - t0:.0f} s)")
