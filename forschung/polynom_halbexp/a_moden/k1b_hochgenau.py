"""Strang A, Aufgabe 1b: f_K-Koeffizienten mit ~70 Stellen (Anfrage Strang D: k ≤ 400, ≥ 60 Stellen).

Problem: der fatou.gp-Fork liefert komplexe sexp nur bis dps 57 (≈ 48 Stellen); slog geht bei dps 80 (≈ 71 Stellen).
Verfahren: w0 = sexp(A(z) + 1/2) bei dps 57, dann ein Newton-Schritt auf A(w) = A(z) + 1/2 mit slog bei dps 80:
    w1 = w0 − (A(w0) − A(z) − 1/2) / A'(w0),   A'(w0) = A'(z) / f'(z)   (aus A(f(z)) = A(z) + 1/2),
A'(z), f'(z) spektral aus den dps-57-Werten auf demselben Kreis (relativer Fehler ≪ 1e-25 genügt, Korrektur ~1e-48).
Prüfungen je Punkt: |A(w1) − A(z) − 1/2| (Newton-Rest) und |A(e^z) − A(z) − 1| (slog-Zweig/Genauigkeit), beide bei dps 80.
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
mp.mp.dps = 100
ctx.prec = 360
M = 8192
RADIEN = sys.argv[1].split(",") if len(sys.argv) > 1 else ["1.34", "1.30"]
PAKET = 128
H = mp.mpf(1) / 2


def zu_acb(x):
    return acb(arb(mp.nstr(mp.re(x), 95)), arb(mp.nstr(mp.im(x), 95)))


def von_acb(x):
    return mp.mpc(mp.mpf(x.real.mid().str(95, radius=False)), mp.mpf(x.imag.mid().str(95, radius=False)))


def voll(halb):
    return halb + [mp.conj(x) for x in halb[1:-1][::-1]]


def dft(werte_voll, inverse=False):
    d = acb.dft([zu_acb(x) for x in werte_voll], inverse=inverse)
    return [von_acb(x) for x in d]


def ableitung_auf_kreis(halb, r):
    """Werte von g'(z_j) aus Werten g(z_j) auf dem Kreis (spektral): c_k = DFT/M = g_k r^k; g'(z) = Σ k g_k z^{k−1}."""
    c = dft(voll(halb))
    c = [x / M for x in c]
    # g'(z_j) = Σ_k k c_k e^{2πijk/M} / z_j  →  inverse DFT von (k c_k) mal M, geteilt durch z_j
    kc = [k * c[k] for k in range(M)]
    w = dft(kc, inverse=True)
    zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M // 2 + 1)]
    return [w[j] * M / zs[j] for j in range(M // 2 + 1)]


def rechne(r_str):
    pfad = DATEN / f"kreis_hoch_r{r_str}_M{M}.json"
    if pfad.is_file():
        d = json.loads(pfad.read_text(encoding="utf-8"))
        return [mp.mpc(*x) for x in d["A"]], [mp.mpc(*x) for x in d["f"]], d
    r = mp.mpf(r_str)
    zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M // 2 + 1)]
    zs[0], zs[-1] = mp.mpc(r, 0), mp.mpc(-r, 0)
    t0 = time.time()
    gp57 = FatouGP(dps=57, fatou_gp="fork")
    s57 = gp57.session("exp(1)")
    A57, f57 = [], []
    for i in range(0, len(zs), PAKET):
        a = s57.slog_batch(zs[i:i + PAKET])
        A57 += a
        f57 += s57.sexp_batch([x + H for x in a])
    gp57.close()
    print(f"  r = {r_str}: dps-57-Werte fertig ({time.time() - t0:.0f} s)", flush=True)
    dA = ableitung_auf_kreis(A57, r)
    df = ableitung_auf_kreis(f57, r)
    gp80 = FatouGP(dps=80, fatou_gp="fork", eval_timeout=300)
    s80 = gp80.session("exp(1)")
    A80, Aw0, Ae = [], [], []
    for i in range(0, len(zs), PAKET):
        tt = time.time()
        A80 += s80.slog_batch(zs[i:i + PAKET])
        t1 = time.time()
        Aw0 += s80.slog_batch(f57[i:i + PAKET])
        t2 = time.time()
        Ae += s80.slog_batch([mp.exp(z) for z in zs[i:i + PAKET]])
        print(f"    Paket {i}: slog(z) {t1-tt:.1f}s, slog(w0) {t2-t1:.1f}s, slog(e^z) {time.time()-t2:.1f}s", flush=True)
    f1 = [w - (aw - a - H) * fd / ad for w, aw, a, fd, ad in zip(f57, Aw0, A80, df, dA)]
    A1 = []
    for i in range(0, len(zs), PAKET):
        A1 += s80.slog_batch(f1[i:i + PAKET])
    gp80.close()
    # reelle Endpunkte exakt reell
    for idx in (0, -1):
        A80[idx], f1[idx] = mp.mpc(mp.re(A80[idx]), 0), mp.mpc(mp.re(f1[idx]), 0)
    newton_vorher = max(abs(aw - a - H) for aw, a in zip(Aw0, A80))
    newton_rest = max(abs(b - a - H) for b, a in zip(A1, A80))
    fg = max(abs(b - a - 1) for b, a in zip(Ae, A80))
    korrektur = max(abs(x - y) for x, y in zip(f1, f57))
    dA57 = max(abs(x - y) for x, y in zip(A80, A57))
    d = {"r": r_str, "M": M, "punkte": "z_j = r·exp(2πij/M), j = 0..M/2", "newton_rest_vorher": mp.nstr(newton_vorher, 4),
         "newton_rest": mp.nstr(newton_rest, 4), "slog_funktionalgl": mp.nstr(fg, 4), "korrektur_f": mp.nstr(korrektur, 4),
         "abstand_A_dps57_dps80": mp.nstr(dA57, 4),
         "A": [[mp.nstr(mp.re(x), 85), mp.nstr(mp.im(x), 85)] for x in A80],
         "f": [[mp.nstr(mp.re(x), 85), mp.nstr(mp.im(x), 85)] for x in f1]}
    pfad.write_text(json.dumps(d), encoding="utf-8")
    print(f"  r = {r_str}: |A(w0)−A−½| vorher {d['newton_rest_vorher']}, nachher {d['newton_rest']}; "
          f"|A(e^z)−A−1| {d['slog_funktionalgl']}; Korrektur |w1−w0| {d['korrektur_f']}; |A80−A57| {d['abstand_A_dps57_dps80']} "
          f"({time.time() - t0:.0f} s)", flush=True)
    return A80, f1, d


def koeff(halb, r_str):
    c = dft(voll(halb))
    r = mp.mpf(r_str)
    return [mp.re(c[k]) / M / r ** k for k in range(M // 2)]


if __name__ == "__main__":
    t0 = time.time()
    erg = {}
    for r_str in RADIEN:
        A, f, info = rechne(r_str)
        erg[r_str] = (koeff(f, r_str), koeff(A, r_str), info)
    a1, A1, i1 = erg[RADIEN[0]]
    a2, A2, i2 = erg[RADIEN[1]]
    eps = max(mp.mpf(i1["newton_rest"]), mp.mpf(i1["slog_funktionalgl"]), mp.mpf(i2["newton_rest"]), mp.mpf(i2["slog_funktionalgl"]))
    L = -mp.lambertw(-1)
    print(f"Motorfehler ε (max Newton-Rest, Funktionalgl.) = {mp.nstr(eps, 3)}")
    print(" k     a_k                                   |a(r1)−a(r2)|   rel     ε·r1^−k/|a_k|")
    KMAX = 1500
    fehler = [max(abs(a1[k] - a2[k]), 10 * eps / mp.mpf(RADIEN[0]) ** k) for k in range(KMAX + 1)]
    for k in [0, 1, 2, 5, 10, 20, 50, 100, 150, 200, 300, 400, 500, 700, 1000, 1200, 1500]:
        print(f"{k:5d}  {mp.nstr(a1[k], 34):>42}  {mp.nstr(abs(a1[k] - a2[k]), 3):>10}  {mp.nstr(abs(a1[k] - a2[k]) / abs(a1[k]), 3):>9}"
              f"  {mp.nstr(10 * eps / mp.mpf(RADIEN[0]) ** k / abs(a1[k]), 3)}")
    alt = json.loads((GEM / "kneser_halbexp_koeff_300.json").read_text(encoding="utf-8"))
    dalt = max(abs(a1[k] - mp.mpf(alt["a"][k])) / abs(a1[k]) for k in range(301))
    print(f"Vergleich mit der dps-57-Fassung (k ≤ 300): max rel. Abstand {mp.nstr(dalt, 3)}")
    stamp = time.strftime("%d.%m.%Y %H:%M")
    beschr = (f"Taylor-Koeffizienten a_k von f_K(z) = T_e(slog_e(z) + 1/2) um 0. Cauchy-DFT auf |z| = {RADIEN[0]} (M = {M}); "
              f"Werte: sexp bei dps 57 + ein Newton-Schritt mit slog bei dps 80 (Motorfehler ε = {mp.nstr(eps, 2)}); Kontrolle |z| = {RADIEN[1]}. "
              f"fehler_abs = max(|a(r1) − a(r2)|, 10·ε·r1^−k). Strang A, {stamp}.")
    for kmax, name in ((300, "kneser_halbexp_koeff_300.json"), (1000, "kneser_halbexp_koeff_1000.json")):
        (GEM / name).write_text(json.dumps({
            "beschreibung": beschr, "kmax": kmax, "a": [mp.nstr(a1[k], 75) for k in range(kmax + 1)],
            "fehler_abs": [mp.nstr(fehler[k], 3) for k in range(kmax + 1)],
            "rel_fehler": [mp.nstr(fehler[k] / abs(a1[k]), 3) for k in range(kmax + 1)]}, indent=0), encoding="utf-8")
    (GEM / "kneser_halbexp_koeff_400.json").write_text(json.dumps([mp.nstr(a1[k], 75) for k in range(401)]), encoding="utf-8")
    fA = [max(abs(A1[k] - A2[k]), 10 * eps / mp.mpf(RADIEN[0]) ** k) for k in range(1001)]
    (GEM / "kneser_slog_koeff_300.json").write_text(json.dumps({
        "beschreibung": "Taylor-Koeffizienten von A(z) = slog_e(z) (Kneser, A(0) = −1) um 0, gleiche Rechnung wie kneser_halbexp_koeff_300.json "
                        "(slog bei dps 80). " + stamp,
        "A": [mp.nstr(A1[k], 75) for k in range(301)], "fehler_abs": [mp.nstr(x, 3) for x in fA[:301]]}, indent=0), encoding="utf-8")
    print(f"geschrieben ({time.time() - t0:.0f} s)")
