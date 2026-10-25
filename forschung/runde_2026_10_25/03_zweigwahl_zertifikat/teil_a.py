"""Teil A in Intervallen (Basis e und 2): Schlichtheit von sexp auf |z| < r (Noshiro-Warschawski) und die Schranken
   sexp z = 1 + A z (1+eps0),  Log sexp z = A z (1+eps1),  |sexp'/A - 1| <= es   auf |z| <= r = 0.1.
Früher (l1_schranken.py, Teil A) nur Basis e und eps-Werte als Gitter-Maximum (720 Punkte, nicht rigoros).
Hier: Maximumprinzip + Überdeckung des Kreises |z| = r durch M Bälle, Auswertung mit acb.

Annahme (E2) über die Motordaten (etappe3/daten/taylor_*.json, Zentrum 0):
   |c_k(Kneser) - c_k(Datei)| <= DELTA für k <= 256   und   |c_k(Kneser)| <= 1 für k > 256.
Aufruf: C:/Users/janis/miniconda3/python.exe teil_a.py     (< 1 min)
"""
import json
from pathlib import Path

from flint import acb, arb, ctx

ctx.prec = 128
HIER = Path(__file__).resolve().parent
E3 = HIER.parent.parent / "identitaeten_2026_10_14" / "etappe3" / "daten"
DATEI = {"e": E3 / "taylor_130_256.json", "2": E3 / "taylor_b2_80_256.json"}
DELTA = "1e-30"
r = arb("0.1")
MBAELLE = 16384


def lade(b):
    d = json.loads(DATEI[b].read_text())
    z = d["zentren"] if "zentren" in d else d
    return [arb(x) for x in z["0.0"]]


def main():
    out = {}
    for b in ["e", "2"]:
        c0 = lade(b)
        N = len(c0) - 1
        A = c0[1]
        # Bälle um die Datenwerte (Annahme E2)
        c = [acb(arb(x.mid(), DELTA)) for x in c0]
        # Rest k > N: |c_k| <= 1  =>  |Rest(z)| <= Σ_{k>N} r^k <= 2 r^{N+1}, |Rest'(z)| <= Σ k r^{k-1} <= 2 (N+1) r^N
        rest0 = 2 * r ** (N + 1)
        rest1 = 2 * (N + 1) * r ** N
        # Noshiro–Warschawski: Σ_{k>=2} k |c_k| r^{k-1} < A  =>  Re sexp' > 0 auf |z| < r  =>  schlicht
        nw = arb(0)
        for k in range(2, N + 1):
            nw += k * (c0[k].abs_upper() + arb(DELTA)) * r ** (k - 1)
        nw += rest1
        # Toleranz: wie gross dürfte DELTA sein?  Σ_{k=2}^N k r^{k-1} <= 1/(1-r)^2 - 1
        delta_max = (A - (nw - rest1)) / (1 / (1 - r) ** 2 - 1)
        e0 = e1 = es = arb(0)
        pi = arb.pi()
        for j in range(MBAELLE):
            # Ball, der den Bogen arg in [2πj/M, 2π(j+1)/M] auf |z| = r enthält
            phi = (2 * j + 1) * pi / MBAELLE
            mid = r * acb(phi.cos(), phi.sin())
            rad = r * pi / MBAELLE * arb("1.0001")
            z = mid + acb(arb(0, rad.mid() + rad.rad()), arb(0, rad.mid() + rad.rad()))
            s = acb(0)
            for ck in reversed(c):
                s = s * z + ck
            s += acb(arb(0, rest0.mid() + rest0.rad()), arb(0, rest0.mid() + rest0.rad()))
            sd = acb(0)
            for k in range(N, 0, -1):
                sd = sd * z + k * c[k]
            sd += acb(arb(0, rest1.mid() + rest1.rad()), arb(0, rest1.mid() + rest1.rad()))
            Az = A * z
            e0 = e0.max(((s - 1) / Az - 1).abs_upper())
            e1 = e1.max((s.log() / Az - 1).abs_upper())
            es = es.max((sd / A - 1).abs_upper())
        # arg(1+eps1) <= arcsin(eps1)  (wird für Im u2 in (-α-0.04, π+α+0.04) gebraucht)
        argmax = e1.asin()
        out[b] = {"A": A.str(20, radius=False), "N": N, "NW_summe_obere": nw.str(8, radius=False),
                  "NW_ok": bool(nw < A), "delta_max": delta_max.str(5, radius=False),
                  "eps0_obere": e0.str(6, radius=False), "eps1_obere": e1.str(6, radius=False),
                  "es_obere": es.str(6, radius=False), "arcsin_eps1": argmax.str(6, radius=False)}
        print(b, out[b])
    (HIER / "daten").mkdir(exist_ok=True)
    (HIER / "daten" / "teil_a.json").write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
