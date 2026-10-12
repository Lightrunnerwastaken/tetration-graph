"""Strang V — (K1)+(K2)+(K3)-Datei für Ms mk_kneser_arb (Format {"K","b_mitte","b_radius","C","p"}), nach Ps p33, mit drei Änderungen:
  * C aus v07 (k2_ergebnis.json, repariertes p32 = v06), nur wenn dort ok = True;
  * Radien enthalten die Dezimalrundung der 140-stelligen Mitten der Ps-Datei (Ms Befund m27) und werden AUFgerundet ausgegeben;
  * Ziel ist der eigene Ordner (phase3/v_pruefstand), nicht gemeinsam/.
Aufruf: python v09_K_datei.py
"""
import json
from pathlib import Path

import sys
from flint import acb, arb, ctx
sys.stdout.reconfigure(encoding="utf-8")

HIER = Path(__file__).resolve().parent
ZIEL = HIER.parent / "kneser_halbexp_b_zert_V.json"
QUELLE = HIER.parents[2] / "gemeinsam" / "kneser_halbexp_koeff_zert_640.json"


def ziffern(s):
    t = s.replace("-", "").replace(".", "").split("e")[0].split("E")[0].lstrip("0")
    return len(t)


def main():
    ctx.prec = 700
    erg = json.loads((HIER / "k2_ergebnis.json").read_text(encoding="utf-8"))
    assert erg["ok"], "v07 meldet offene Punkte"
    rho = abs(acb(-1).lambertw())
    d = json.loads(QUELLE.read_text(encoding="utf-8"))
    b = []
    for n, (mi, er, _) in enumerate(d["a"]):
        m = arb(mi)
        dez = abs(m) * arb(10) ** (-(ziffern(mi) - 2))         # Rundung der Dezimalmitte
        a = arb(m.mid(), (m.rad() + arb(er) * (1 + arb("1e-4")) + dez).upper())
        b.append(a * rho ** n)
    K = len(b) - 1
    C = arb(erg["C6"]) * (1 + arb("1e-6"))
    C_str = "%.6g" % (float(C.upper()) * 1.00001)
    assert arb(C_str) >= arb(erg["C6"])
    p = 6
    Ka = arb(K)
    rest = [arb(C_str) * Ka ** (j + 1 - p) / (p - j - 1) for j in range(3)]     # Σ_{k>K} k^j|b_k|, j = 0, 1, 2
    aus = {
        "beschreibung": "Strang V (Phase 3): b_k = a_k ρ^k von Knesers f_K. (K1) Kugeln für 0 ≤ k ≤ K aus Ps Satz P-V (Radien inkl. "
                        "Dezimalrundung der Mitten, aufgerundet). (K2) |b_k| ≤ C·k^(−p) für ALLE k ≥ 1: sechsmal partiell integrieren auf |x| = ρ, "
                        "C ≥ ‖∂_θ⁶ f_K(ρe^{iθ})‖_L¹/(2π), zertifiziert mit v06 (repariertes p32: Zellen in d = |θ − θ_L|, nan-sicher, Zweigtreue je Zelle) "
                        "und v07 (Abschnittsgrenzen, Rest nach p32.rest_schranke). (K3) f_K(L) = L, f_K'(L) = √L aus p32a (Ps Struktursatz).",
        "K": K, "p": str(p), "C": C_str,
        "b_mitte": [x.mid().str(60, radius=False) for x in b],
        "b_radius": ["%.6e" % (float(x.rad().upper()) * 1.001 + 1e-300) for x in b],
        "restsummen_ueber_K": [r.str(5) for r in rest],
        "k2_ergebnis": erg,
    }
    ZIEL.write_text(json.dumps(aus, indent=0, ensure_ascii=False), encoding="utf-8")
    print(f"C = {C_str}; Σ_(k>{K}) k^j|b_k| ≤ {[r.str(3) for r in rest]}; geschrieben: {ZIEL}")


if __name__ == "__main__":
    main()
