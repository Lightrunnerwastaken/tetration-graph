"""Strang V — unabhängige Prüfung der Kette 193 → 400 (eigener Code, nicht zert_kette.py).
Prüft: (1) Ps Box von p_193 (p24_zustand.json, md5 = Ps Datei) ⊂ B_193; (2) für jeden Übergang N → N+1: erg/K_{N+1}.json hat ok = True und im letzten
Versuch streng, start_ok, innen, in_roehre = True; (3) K_{N+1} ⊂ B_{N+1} koordinatenweise mit Radien (strikt, in arb mit 1600 Bit);
(4) Lückenlosigkeit 193…400. Ausgabe: Zusammenfassung und kleinste Reserve, grösster Wrapping-Term t3.
"""
import json
import sys
from pathlib import Path

from flint import arb, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
LAUF = HIER / "lauf"
ctx.prec = 1600


def dek(x):
    m, r, e = x
    z = arb(10) ** int(e)
    return arb(int(m)) * z, arb(int(r)) * z      # (Mitte, Radius) exakt als arb


def enthalten(B, K):
    """[K] ⊂ int [B] koordinatenweise: K_mid − K_rad > B_mid − B_rad und K_mid + K_rad < B_mid + B_rad."""
    if len(B) != len(K):
        return False, "Länge"
    for k, ((bm, br), (km, kr)) in enumerate(zip(B, K)):
        if not ((km - kr) - (bm - br) > 0 and (bm + br) - (km + kr) > 0):
            return False, f"Koeffizient {k}"
    return True, ""


def main():
    fehler = []
    z = json.loads((HIER / "daten" / "p24_zustand.json").read_text(encoding="utf-8"))
    assert z["N"] == 193
    P = [dek(x) for x in z["box"]]
    B193 = [dek(x) for x in json.loads((LAUF / "boxen" / "B_193.json").read_text(encoding="utf-8"))["box"]]
    ok, g = enthalten(B193, P)
    if not ok:
        fehler.append(f"P-Box ⊄ B_193: {g}")
    min_res, max_t3, max_theta = 1.0, 0.0, 0.0
    for N in range(193, 400):
        M = N + 1
        d = json.loads((LAUF / "erg" / f"K_{M}.json").read_text(encoding="utf-8"))
        v = d["versuche"][-1]
        flags = [d.get("ok"), v.get("ok"), v.get("streng"), v.get("start_ok"), v.get("innen"), v.get("in_roehre")]
        if not all(flags):
            fehler.append(f"Übergang {N}→{M}: Flags {flags}")
            continue
        min_res = min(min_res, v["reserve"])
        max_t3 = max(max_t3, v["t3"])
        max_theta = max(max_theta, v["theta"])
        K = [dek(x) for x in d["box"]]
        if M <= 399 or (LAUF / "boxen" / f"B_{M}.json").exists():
            pB = LAUF / "boxen" / f"B_{M}.json"
            if pB.exists():
                Bm = [dek(x) for x in json.loads(pB.read_text(encoding="utf-8"))["box"]]
                ok, g = enthalten(Bm, K)
                if not ok:
                    fehler.append(f"K_{M} ⊄ B_{M}: {g}")
            elif M < 400:
                fehler.append(f"B_{M} fehlt")
    print(f"Kette 193 → 400 (unabhängige Prüfung): Fehler {fehler if fehler else 'keine'}; kleinste Reserve {min_res:.7f}, "
          f"grösstes t3 {max_t3:.2e}, grösstes θ {max_theta:.2e}")
    (HIER / "v16_kette.json").write_text(json.dumps({"fehler": fehler, "min_reserve": min_res, "max_t3": max_t3, "max_theta": max_theta},
                                                    indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
