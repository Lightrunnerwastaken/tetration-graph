"""Strang V — (K2) zusammensetzen: Abschnitte aus v06 prüfen und C₆ bilden.
(1) Jeder Abschnitt: zweig_ok, Integral endlich; Abschnitte einer Seite lückenlos (d_b von i = d_a von i+1, von Rand bis δ_min), Rest vorhanden.
(2) Zweigtreue an jeder Abschnittsgrenze: Punkt-Krawczyk für sexp_K(s) = x* am Grenzpunkt ⇒ winzige Scheibe D* mit genau einer Lösung;
    D* ⊂ letzte z-Scheibe des oberen UND erste z-Scheibe des unteren Abschnitts ⇒ derselbe Zweig. Startzweig bei θ = 0 bzw. π: reelles slog
    (Startscheibe symmetrisch zur reellen Achse wird nicht verlangt; der erste Abschnitt startet bei reellem x, Newton reell).
(3) ‖∂_θ⁶f‖_{L¹(0,2π)} ≤ 2·Σ_Seiten(Σ Integrale + Rest)  (Symmetrie f(x̄) = conj f(x));  C₆ = ‖∂_θ⁶f‖_{L¹}/(2π);  |b_k| ≤ C₆k^{−6} für alle k ≥ 1.
Aufruf: python v07_k2_zusammen.py [ausgeschlossene Abschnitte, z. B. rechts_8 links_9]
"""
import json
import math
import sys
from pathlib import Path

from flint import acb, arb

AUSSCHLUSS = set(sys.argv[1:])        # z. B. rechts_8 links_9 (überlappende Varianten)
sys.argv = [sys.argv[0], "test", "40"]
import p32_rand_ableitungen as P32

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent


def lese_scheibe(paar):
    # paar = (zc als acb-String, zr als arb-String)
    s = paar[0].replace("j", "")
    # acb-String "[a +/- e] + [b +/- f]j" bzw. ohne Imaginärteil
    teile = s.split(" + ")
    if len(teile) == 2:
        re, im = arb(teile[0]), arb(teile[1])
    else:
        t = s.split(" - [")
        if len(t) == 2:
            re, im = arb(t[0]), -arb("[" + t[1])
        else:
            re, im = arb(s), arb(0)
    zc = acb(re, im)
    zr = arb(paar[1])
    # Mitte als Punkt, Radius vergrössert um die Unschärfe der Darstellung (konservativ: Scheibe wird KLEINER gemacht)
    unsch = zc.real.rad() + zc.imag.rad()
    return acb(zc.real.mid(), zc.imag.mid()), zr - unsch


def main():
    RD = P32.Rand()
    K = RD.K
    L = RD.L
    rho = RD.R.absL
    thL = L.imag
    I = acb(0, 1)
    gesamt = arb(0)
    alles_ok = True
    bericht = {}
    for seite, s in (("rechts", -1), ("links", 1)):
        teile = []
        for p in sorted(HIER.glob(f"k2_{seite}_*.json")):
            if p.name.endswith("_zellen.json"):
                continue
            d = json.loads(p.read_text(encoding="utf-8"))
            if d["index"] >= 90 or f"{seite}_{d['index']}" in AUSSCHLUSS:          # Testläufe / ausgeschlossene Varianten
                continue
            teile.append(d)
        teile.sort(key=lambda d: -d["d_a"])
        summe = arb(0)
        fehler = []
        for i, d in enumerate(teile):
            if not d["zweig_ok"]:
                fehler.append(f"Abschnitt {d['index']}: Zweigtreue intern offen")
            if not math.isfinite(d["integral_obere"]):
                fehler.append(f"Abschnitt {d['index']}: Integral nicht endlich")
            summe += arb(d["integral_obere"]) * (1 + arb("1e-12"))      # (F2): float-Obergrenze mit Sicherheitsfaktor
            if i + 1 < len(teile):
                d2 = teile[i + 1]
                if abs(d["d_b"] - d2["d_a"]) > 1e-12 * max(1.0, d2["d_a"]):
                    fehler.append(f"Lücke zwischen {d['index']} (d_b = {d['d_b']}) und {d2['index']} (d_a = {d2['d_a']})")
                # Zweigtreue an der Grenze d* = d2["d_a"]
                zc1, zr1 = lese_scheibe(d["letzte_z"])
                zc2, zr2 = lese_scheibe(d2["erste_z"])
                xs = rho * (I * (thL + s * arb(d2["d_a"]))).exp()
                wk, rk = K.slog_krawczyk(xs, zc2, arb(0))
                ok = wk is not None and bool(abs(wk - zc1) + rk < zr1) and bool(abs(wk - zc2) + rk < zr2)
                if not ok:
                    fehler.append(f"Zweigtreue an der Grenze {d['index']}|{d2['index']} offen")
        if not teile or abs(teile[0]["d_a"] - (float(thL.mid()) if seite == "rechts" else math.pi - float(thL.mid()))) > 1e-9:
            fehler.append("Rand fehlt")
        if not teile or not teile[-1]["bis_L"] or "rest_obere" not in teile[-1]:
            fehler.append("Rest nahe L fehlt: " + str(teile[-1].get("rest_fehler") if teile else ""))
            rest = float("inf")
        else:
            rest = teile[-1]["rest_obere"] * (1 + 1e-12)
        summe_seite = summe + arb(rest)
        bericht[seite] = {"abschnitte": len(teile), "integral": float(summe.upper()), "rest": rest, "fehler": fehler}
        alles_ok = alles_ok and not fehler
        gesamt += summe_seite
        print(f"{seite}: {len(teile)} Abschnitte, ∫ ≤ {float(summe.upper()):.6g}, Rest ≤ {rest:.3g}, Fehler: {fehler}")
    L1 = 2 * gesamt
    C6 = L1 / (2 * arb.pi())
    print(f"‖∂_θ⁶f‖_L¹ ≤ {float(L1.upper()):.6g}  ⇒  C₆ ≤ {float(C6.upper()):.6g}  (Ziel ≤ 5e5); alles geprüft: {alles_ok}")
    bericht.update(L1=float(L1.upper()), C6=float(C6.upper()), ok=alles_ok)
    (HIER / "k2_ergebnis.json").write_text(json.dumps(bericht, indent=1, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
