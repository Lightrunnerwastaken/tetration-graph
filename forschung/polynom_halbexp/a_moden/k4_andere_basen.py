"""Strang A, Anfrage Strang E (BRETT 19:51): Kneser-Werte für andere Basen und Drittel-Iterierte.

sexp_b(1/2), sexp_b(−1/2) für b = 2, 10, 1.5; sexp_e(1/3), sexp_e(−2/3); Taylor-Koeffizienten der halben b-Funktion
f_b(z) = sexp_b(slog_b(z) + 1/2) um 0 (k ≤ 40) für b = 2, 10 (Cauchy-DFT auf |z| = 0.8·|L_b|, M = 256, Kontrolle 0.7·|L_b|).
Motor: fatou.gp-Fork bei dps 50. Prüfung: sexp_b(slog_b(z) + 1) = b^z auf dem Kreis.
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from fatou_backend import FatouGP  # noqa: E402

HIER = Path(__file__).resolve().parent
GEM = HIER.parent / "gemeinsam"
mp.mp.dps = 60
DPS = 50
t0 = time.time()
gp = FatouGP(dps=DPS, fatou_gp="fork")
aus = {"beschreibung": "Kneser-Werte (fatou.gp-Fork, dps 50) für Strang E. Strang A, " + time.strftime("%d.%m.%Y %H:%M"),
       "werte": {}, "halbe_funktion_koeff": {}}
BASEN = {"2": mp.mpf(2), "10": mp.mpf(10), "1.5": mp.mpf("1.5"), "exp(1)": mp.e}
for bs, b in BASEN.items():
    s = gp.session(bs)
    pkt = [mp.mpf("0.5"), mp.mpf("-0.5")] + ([mp.mpf(1) / 3, mp.mpf(-2) / 3] if bs == "exp(1)" else [])
    v = s.sexp_batch(pkt)
    # Kontrolle: sexp(x+1) = b^sexp(x)
    v1 = s.sexp_batch([x + 1 for x in pkt])
    ktr = max(abs(w1 - b ** w) for w, w1 in zip(v, v1))
    for x, w in zip(pkt, v):
        aus["werte"][f"sexp_{bs}({mp.nstr(x, 6)})"] = mp.nstr(mp.re(w), 45)
    print(f"b = {bs}: " + ", ".join(f"sexp({mp.nstr(x, 4)}) = {mp.nstr(mp.re(w), 35)}" for x, w in zip(pkt, v))
          + f";  Kontrolle |sexp(x+1) − b^sexp(x)| = {mp.nstr(ktr, 3)}  ({time.time() - t0:.0f} s)", flush=True)
    if bs in ("2", "10"):
        lnb = mp.log(b)
        Lb = -mp.lambertw(-lnb) / lnb
        Lb = mp.conj(Lb) if mp.im(Lb) < 0 else Lb
        assert abs(b ** Lb - Lb) < mp.mpf(10) ** -50
        M = 256
        koeff = {}
        for fak in ("0.8", "0.7"):
            r = mp.mpf(fak) * abs(Lb)
            zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M)]
            A = s.slog_batch(zs)
            f = s.sexp_batch([a + mp.mpf("0.5") for a in A])
            e1 = s.sexp_batch([a + 1 for a in A])
            pr = max(abs(w - b ** z) for w, z in zip(e1, zs))
            koeff[fak] = [mp.re(mp.fsum(f[j] * mp.expjpi(-2 * k * mp.mpf(j) / M) for j in range(M)) / M) / r ** k for k in range(41)]
            print(f"   |L_b| = {mp.nstr(abs(Lb), 10)}, r = {fak}|L_b|: Prüfrest {mp.nstr(pr, 3)}", flush=True)
        fehler = [abs(x - y) for x, y in zip(koeff["0.8"], koeff["0.7"])]
        aus["halbe_funktion_koeff"][bs] = {"L_b": [mp.nstr(mp.re(Lb), 30), mp.nstr(mp.im(Lb), 30)],
                                           "a": [mp.nstr(x, 40) for x in koeff["0.8"]], "fehler_abs": [mp.nstr(x, 3) for x in fehler]}
        print(f"   a_0..a_4 = {[mp.nstr(x, 15) for x in koeff['0.8'][:5]]}; max Fehler (Radienvergleich) {mp.nstr(max(fehler), 3)}; "
              f"Σa_k − sexp(1/2) = {mp.nstr(mp.fsum(koeff['0.8']) - mp.re(v[0]), 3)}", flush=True)
gp.close()
(GEM / "kneser_andere_basen.json").write_text(json.dumps(aus, indent=1, ensure_ascii=False), encoding="utf-8")
print(f"geschrieben ({time.time() - t0:.0f} s)")
