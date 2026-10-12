"""zert_kette — prüft die Kette der parallel gerechneten Übergänge (Strang X):
  (0) Ps zertifizierte Box von p_193 (p24_zustand.json) ⊂ B_193;
  (N) für N = N1 … N2−1: Übergang N → N+1 ok (Krawczyk streng, Start-Inklusion, Röhre) und K_{N+1} ⊂ B_{N+1};
Ist alles erfüllt, folgt induktiv: p_N ∈ B_N für alle N bis N2, d. h. der σ-Homotopie-Zweig p_1 → p_N2 ist rigoros zertifiziert
(bis 193 durch P: p18 + p24). Ausgabe: Tabelle und Zusammenfassung (kette.json).
Aufruf: python zert_kette.py N1 N2
"""
import json
import sys
from pathlib import Path

from flint import arb, ctx

import zert_kern as Z

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
N1, N2 = int(sys.argv[1]), int(sys.argv[2])
ctx.prec = 2 * N2 + 600


def lade_box(pfad):
    d = json.loads(pfad.read_text(encoding="utf-8"))
    return d, [Z.dekodiere(x) for x in d["box"]] if "box" in d else None


def enthalten(B, K):
    if len(B) != len(K):
        return False, "Länge"
    for k, (b, t) in enumerate(zip(B, K)):
        lo_ok = (t.mid() - t.rad()) - (b.mid() - b.rad())
        hi_ok = (b.mid() + b.rad()) - (t.mid() + t.rad())
        if not (lo_ok > 0 and hi_ok > 0):
            return False, f"Koeffizient {k}"
    return True, ""


ergebnis = {"N1": N1, "N2": N2, "fehler": [], "ok": True}
z = json.loads((HIER / "daten" / "p24_zustand.json").read_text(encoding="utf-8"))
Pbox = [Z.dekodiere(x) for x in z["box"]]
dB, B = lade_box(HIER / "boxen" / f"B_{N1}.json")
ok0, grund = enthalten(B, Pbox)
print(f"(0) P-Box p_{z['N']} ⊂ B_{N1}: {ok0} {grund}")
if z["N"] != N1 or not ok0:
    ergebnis["ok"] = False
    ergebnis["fehler"].append(f"Start: {grund or 'N passt nicht'}")
zeilen = []
for N in range(N1, N2):
    pK = HIER / "erg" / f"K_{N + 1}.json"
    pB = HIER / "boxen" / f"B_{N + 1}.json"
    if not pK.exists():
        ergebnis["ok"] = False
        ergebnis["fehler"].append(f"{N}→{N + 1}: kein Ergebnis")
        continue
    dK, K = lade_box(pK)
    if not dK["ok"] or K is None:
        ergebnis["ok"] = False
        ergebnis["fehler"].append(f"{N}→{N + 1}: Übergang nicht ok")
        continue
    if N + 1 < N2 or pB.exists():
        if not pB.exists():
            ergebnis["ok"] = False
            ergebnis["fehler"].append(f"B_{N + 1} fehlt")
            continue
        _, Bn = lade_box(pB)
        oki, grund = enthalten(Bn, K)
        if not oki:
            ergebnis["ok"] = False
            ergebnis["fehler"].append(f"K_{N + 1} ⊄ B_{N + 1}: {grund}")
    p = dK["versuche"][-1]
    zeilen.append(dict(M=N + 1, NC=p["NC"], reserve=p["reserve"], t3=p["t3"], theta=p["theta"], r=p["r"], t=dK["t_gesamt"]))
for zl in zeilen[:: max(1, len(zeilen) // 25)] + zeilen[-1:]:
    print(zl)
ergebnis["uebergaenge"] = zeilen
(HIER / f"kette_{N1}_{N2}.json").write_text(json.dumps(ergebnis, indent=1), encoding="utf-8")
print(f"KETTE {N1} → {N2}: ok = {ergebnis['ok']}; Fehler: {ergebnis['fehler'][:10]}")
