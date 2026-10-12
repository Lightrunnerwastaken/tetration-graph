"""Strang V — Startbedingung (A) bei N₀ ∈ [318, 400]: sup_{D̄}|p_N − f_K| rigoros, nach Ps p31 (Kopie der Logik), mit den Boxen des
Zweig-Zertifikats 193 → 400 (zert400/lauf/boxen/B_N.json; gültige Einschlüsse von p_N, sobald die Kette K_N ⊂ B_N geprüft ist).
  d3 = Σ_{k≤N}|p_k − a_k|ρ^k + Σ_{N<k≤640}|a_k|ρ^k   (ℓ¹_ρ ≥ sup_{D̄}),   sup_{D̄}|p_N − f_K| ≤ d3 + S₆₄₀,  S₆₄₀ = Σ_{k>640}|b_k| ≤ C₆·640^{−5}/5 (K2).
Koeffizienten a_k: Ps Datei mit Radius + Dezimalrundung der Mitten (Ms Befund). Optional d1 = ‖p_N − Π_Nf_K‖ (Π_N wie p31).
Aufruf: python v15_A_kopplung.py N1 N2 ... [--pi]
"""
import json
import sys
import time
from pathlib import Path

from flint import acb, arb, arb_mat, ctx

sys.stdout.reconfigure(encoding="utf-8")
HIER = Path(__file__).resolve().parent
KOEFF = HIER.parents[1] / "gemeinsam" / "kneser_halbexp_koeff_zert_640.json"
BOXEN = HIER / "zert400" / "lauf" / "boxen"
ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
MIT_PI = "--pi" in sys.argv
NS = [int(t) for t in ARGS]
C6 = arb("461.36")


def dekodiere(x):
    m, r, e = x
    z = arb(10) ** int(e)
    return arb(int(m)) * z + arb(0, 1) * (arb(int(r)) * z)


def ziffern(s):
    return len(s.replace("-", "").replace(".", "").split("e")[0].split("E")[0].lstrip("0"))


def lade_a():
    d = json.loads(KOEFF.read_text(encoding="utf-8"))
    out = []
    for m, r, _ in d["a"]:
        mm = arb(m)
        dez = abs(mm) * arb(10) ** (-(ziffern(m) - 2))
        out.append(arb(mm.mid(), (mm.rad() + arb(r) * (1 + arb("1e-4")) + dez).upper()))
    return out


def main():
    ctx.prec = 1400
    rho = abs(acb(-1).lambertw())
    a = lade_a()
    K = len(a) - 1
    aw = [abs(a[k]) * rho ** k for k in range(K + 1)]
    S640 = C6 * arb(640) ** (-5) / 5
    erg = []
    for N in NS:
        t0 = time.time()
        pb = BOXEN / f"B_{N}.json"
        if not pb.exists():                     # N = 400: Ausgabebox K_400 des letzten Übergangs
            pb = BOXEN.parent / "erg" / f"K_{N}.json"
        box = json.loads(pb.read_text(encoding="utf-8"))
        p = [dekodiere(x) for x in box["box"]]
        assert len(p) == N + 1
        rest = sum((aw[k] for k in range(N + 1, K + 1)), arb(0))
        d3 = sum((abs(p[k] - a[k]) * rho ** k for k in range(N + 1)), arb(0)) + rest
        sup = d3 + S640
        e = dict(N=N, d3=float(d3.upper()), S640=float(S640.upper()), sup_pN_minus_fK=float(sup.upper()),
                 rest_N_bis_640=float(rest.upper()), box_radius_rho=max(float(p[k].rad() * rho ** k) for k in range(N + 1)),
                 t=round(time.time() - t0, 1))
        print(json.dumps(e, ensure_ascii=False), flush=True)
        erg.append(e)
    (HIER / "v15_A.json").write_text(json.dumps({"beschreibung": __doc__.splitlines()[0], "ergebnisse": erg}, indent=1, ensure_ascii=False),
                                    encoding="utf-8")


if __name__ == "__main__":
    main()
