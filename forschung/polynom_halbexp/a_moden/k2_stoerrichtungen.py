"""Strang A, Aufgabe 2: Störrichtungen δ_n(z) = T'(A(z) + 1/2)·e^{2πinA(z)}, A = slog_e (Kneser), n = ±1, ±3, ±5.

Taylor-Koeffizienten um 0 aus Werten auf |z| = r (Kreisdaten von k1b: A und f_K bei ~70 Stellen):
    T'(A + 1/2) = f_K'(z)/A'(z)  (spektrale Ableitungen auf dem Kreis),
    untere Hälfte: δ_n(z̄) = conj(δ_{−n}(z)).
Prüfungen: (a) linearisierte Gleichung δ∘f + (f'∘f)·δ = 0 punktweise (Taylor-Reihen, Punkte mit |z|, |f(z)| < |L|) und
koeffizientenweise; (b) Wachstum |δ_n(tL)| für t → 1 gegen |z − L|^{Re(1 − 2πin/L)} (n < 0) bzw. bei L̄ (n > 0);
(c) Vergleich zweier Radien.
Ausgabe: gemeinsam/stoerrichtungen_delta.json (Re/Im der Koeffizienten von δ_1, δ_3, δ_5; δ_{−n} = konjugiert).
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.stdout.reconfigure(encoding="utf-8")
import mpmath as mp  # noqa: E402
from flint import acb, arb, ctx  # noqa: E402

HIER = Path(__file__).resolve().parent
DATEN = HIER / "daten"
GEM = HIER.parent / "gemeinsam"
mp.mp.dps = 90
ctx.prec = 330
KMAX = 1500
NS = [1, 3, 5]


def lade(r_str):
    for name in (f"kreis_hoch_r{r_str}_M8192.json", f"kreis_r{r_str}_M4096.json"):
        p = DATEN / name
        if p.is_file():
            d = json.loads(p.read_text(encoding="utf-8"))
            return [mp.mpc(*x) for x in d["A"]], [mp.mpc(*x) for x in d["f"]], d["M"], name
    raise FileNotFoundError(r_str)


def zu_acb(x):
    return acb(arb(mp.nstr(mp.re(x), 85)), arb(mp.nstr(mp.im(x), 85)))


def von_acb(x):
    return mp.mpc(mp.mpf(x.real.mid().str(85, radius=False)), mp.mpf(x.imag.mid().str(85, radius=False)))


def dft(v, inverse=False):
    return [von_acb(x) for x in acb.dft([zu_acb(x) for x in v], inverse=inverse)]


def voll(halb):
    return halb + [mp.conj(x) for x in halb[1:-1][::-1]]


def kreisdaten(r_str):
    A, f, M, name = lade(r_str)
    r = mp.mpf(r_str)
    zs = [r * mp.expjpi(2 * mp.mpf(j) / M) for j in range(M // 2 + 1)]

    def abl(halb):
        c = [x / M for x in dft(voll(halb))]
        w = dft([k * c[k] for k in range(M)], inverse=True)
        return [w[j] * M / zs[j] for j in range(M // 2 + 1)]

    dA, df = abl(A), abl(f)
    Tp = [a / b for a, b in zip(df, dA)]          # T'(A + 1/2) auf der oberen Hälfte
    koeffs = {}
    for n in NS:
        oben_p = [t * mp.expjpi(2 * n * a) for t, a in zip(Tp, A)]     # δ_n auf der oberen Hälfte
        oben_m = [t * mp.expjpi(-2 * n * a) for t, a in zip(Tp, A)]    # δ_{−n} auf der oberen Hälfte
        werte = oben_p + [mp.conj(x) for x in oben_m[1:-1][::-1]]      # δ_n(z̄) = conj(δ_{−n}(z))
        c = dft(werte)
        koeffs[n] = [c[k] / M / r ** k for k in range(min(KMAX + 1, M))]
    fk = [mp.re(x) / M / r ** k for k, x in enumerate(dft(voll(f))[:KMAX + 1])]
    Ak = [mp.re(x) / M / r ** k for k, x in enumerate(dft(voll(A))[:KMAX + 1])]
    return koeffs, fk, Ak, name, max(abs(x) for x in Tp)


def reihe(c, z, kmax=None):
    s = mp.mpc(0)
    for x in reversed(c[:kmax] if kmax else c):
        s = s * z + x
    return s


def reihe_abl(c, z):
    return reihe([k * c[k] for k in range(1, len(c))], z)


if __name__ == "__main__":
    t0 = time.time()
    L = -mp.lambertw(-1)
    L = mp.conj(L) if mp.im(L) < 0 else L
    KD = {r: kreisdaten(r) for r in ("1.34", "1.30", "1.25", "1.20")}
    RAD = {1: ("1.34", "1.30"), 3: ("1.30", "1.25"), 5: ("1.25", "1.20")}   # Aliasing ~ (r/|L|)^M·M^(4.447n−2) klein halten
    K1 = {n: KD[RAD[n][0]][0][n] for n in NS}
    K2 = {n: KD[RAD[n][1]][0][n] for n in NS}
    f1, A1 = KD["1.34"][1], KD["1.34"][2]
    print(f"Kreisdaten: {[KD[r][3] for r in KD]}; Radienpaare {RAD}  ({time.time() - t0:.0f} s)")
    # (c) Radienvergleich, relativ zur Hüllkurve
    for n in NS:
        d = [abs(K1[n][k] - K2[n][k]) / max(abs(K1[n][k]), mp.mpf(10) ** -200) for k in range(KMAX + 1)]
        print(f"δ_{n}: rel. Abstand der Radien bei k = 0, 50, 100, 200, 300, 500, 800: "
              + ", ".join(mp.nstr(d[k], 2) for k in (0, 50, 100, 200, 300, 500, 800)))
        print(f"     |δ_{n},k|·|L|^k bei k = 0, 10, 50, 100, 200, 300, 500: "
              + ", ".join(mp.nstr(abs(K1[n][k]) * abs(L) ** k, 3) for k in (0, 10, 50, 100, 200, 300, 500)))
    # (a) linearisierte Gleichung punktweise
    print("\n(a) δ(f(z)) + f'(f(z))·δ(z) (relativ zu |δ(f(z))|):")
    pkt = [mp.mpf(0), mp.mpf("0.3"), mp.mpf("0.5"), mp.mpf("-0.5"), mp.mpc("0.3", "0.3"), mp.mpc("-0.4", "0.6")]
    for z in pkt:
        fz = reihe(f1, z)
        zeile = []
        for n in NS:
            for s in (1, -1):
                c = K1[n] if s == 1 else [mp.conj(x) for x in K1[n]]
                lhs = reihe(c, fz) + reihe_abl(f1, fz) * reihe(c, z)
                zeile.append(f"n={s * n:+d}: {mp.nstr(abs(lhs) / abs(reihe(c, fz)), 2)}")
        print(f"  z = {mp.nstr(z, 3)} (|f(z)| = {mp.nstr(abs(fz), 4)}): " + ", ".join(zeile))
    # koeffizientenweise: [δ∘f + (f'∘f)δ]_k, k ≤ 40, mit Reihen bis KMAX (flint-Polynome)
    from flint import acb_poly
    ctx.prec = 400
    K = 60
    fpol = acb_poly([acb(arb(mp.nstr(x, 80))) for x in f1[:KMAX + 1]])
    # Potenzen von f bis Ordnung K: (f)^j mod x^{K+1}
    pot = [acb_poly([1])]
    ftr = acb_poly([acb(arb(mp.nstr(x, 80))) for x in f1[:K + 1]])
    for j in range(KMAX):
        pot.append((pot[-1] * ftr).truncate(K + 1))
    fpf = acb_poly([0])
    for j in range(1, KMAX + 1):
        fpf += j * acb(arb(mp.nstr(f1[j], 80))) * pot[j - 1]
    for n in NS:
        c = [zu_acb(x) for x in K1[n]]
        dof = acb_poly([0])
        for j in range(KMAX + 1):
            dof += c[j] * pot[j]
        rest = (dof + (fpf * acb_poly(c[:K + 1])).truncate(K + 1)).truncate(K + 1)
        norm = max(abs(von_acb(dof[k])) * mp.mpf("0.5") ** k for k in range(K + 1))
        res = max(abs(von_acb(rest[k])) * mp.mpf("0.5") ** k for k in range(K + 1))
        print(f"  koeffizientenweise (k ≤ {K}, Gewicht 0.5^k): δ_{n}: max|Rest| = {mp.nstr(res, 3)}, max|[δ∘f]_k| = {mp.nstr(norm, 3)}")
    ctx.prec = 330
    # (b) Wachstum zum Rand: δ_n(tL) mit A(tL) aus dem Motor (dps 57), T'(A+1/2) = f'(z)/A'(z) aus den Reihen
    from fatou_backend import FatouGP
    gp = FatouGP(dps=57, fatou_gp="fork")
    s = gp.session("exp(1)")
    ts = [mp.mpf(x) for x in ("0.9", "0.95", "0.98", "0.99", "0.995")]
    pts = [t * L for t in ts] + [t * mp.conj(L) for t in ts]
    Av = s.slog_batch(pts)
    gp.close()
    print("\n(b) Wachstum zum Rand: Exponent α in |δ_n(z)| ~ |z − L|^α (aus aufeinanderfolgenden t)")
    for n in NS:
        for s_ in (1, -1):
            nn = s_ * n
            ws = []
            for z, a in zip(pts, Av):
                Tp = reihe_abl(f1, z) / reihe_abl(A1, z)
                ws.append(abs(Tp * mp.expjpi(2 * nn * a)))
            obenL = ws[:len(ts)]
            untenL = ws[len(ts):]
            ex = lambda w: [mp.log(w[i + 1] / w[i]) / mp.log((1 - ts[i + 1]) / (1 - ts[i])) for i in range(len(ts) - 1)]
            soll = mp.re(1 - 2j * mp.pi * nn / L)
            print(f"  n = {nn:+d}: bei L: α ≈ {', '.join(mp.nstr(x, 4) for x in ex(obenL))};  bei L̄: α ≈ "
                  f"{', '.join(mp.nstr(x, 4) for x in ex(untenL))};  Soll Re(1 − 2πin/L) = {mp.nstr(soll, 5)} (an der wachsenden Seite)")
    # Ausgabe
    (GEM / "stoerrichtungen_delta.json").write_text(json.dumps({
        "beschreibung": ("Taylor-Koeffizienten um 0 von δ_n(z) = T'(A(z)+1/2)·exp(2πi n A(z)), A = Kneser-slog_e (A(0) = −1, A(1) = 0), "
                         "n = 1, 3, 5. δ_{−n} hat die konjugierten Koeffizienten. Reelle Richtungen: Re und Im der Koeffizienten. "
                         "δ_n (n > 0) wächst bei L̄ wie |z − L̄|^(1 − 4.447 n), δ_{−n} bei L. Cauchy-DFT (M = 4096, Motor dps 57), "
                         "Kontrolle mit zweitem Radius; Radien je n: " + str(RAD) + ". Strang A, " + time.strftime("%d.%m.%Y %H:%M")),
        "kmax": KMAX,
        "delta": {str(n): {"re": [mp.nstr(mp.re(x), 40) for x in K1[n]], "im": [mp.nstr(mp.im(x), 40) for x in K1[n]],
                           "rel_abstand_radien": [mp.nstr(abs(K1[n][k] - K2[n][k]) / max(abs(K1[n][k]), mp.mpf(10) ** -200), 2)
                                                  for k in range(KMAX + 1)]} for n in NS}}), encoding="utf-8")
    print(f"geschrieben ({time.time() - t0:.0f} s)")
