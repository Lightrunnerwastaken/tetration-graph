"""exp.adressen — Mischbasen-Phase für komplexe fliehende Basen (Test von bw.verm.adressen).

 Gestoppte Kette (Paper III, komplex): w = T_d(x) (Kneser/Paulsen, komplexe Basis), iterieren bis |c_d w| ≥ R (Stufe k),
 dann oberste Schälung über den Exponenten u = c_d w:  L_b(e^u) = (u − 2πi·round(Im u/2π))/c_b,  danach k-mal L_b(v) = Log v / c_b.
 D(x) = L_b^{k+1}(T_d(x+k+1)),  Φ(x) = slog_b(D(x)) − x  (= lim_n F(x+n) − (x+n), F = slog_b ∘ T_d).
 Fassung 1 addierte k + 1 — falsch (Sprünge um 1, wo die Stoppstufe wechselt).
 Prüfungen: Unabhängigkeit von R (Schwänze konvergieren), Φ(x+1) = Φ(x), Sprünge in x, Adressen s_j = round(Im(c_d z_j)/2π).
"""
import mpmath as mp

from gemeinsam import Lauf, plt
from kneser_sicher import kneser_sicher

mp.mp.dps = 40
lauf = Lauf("exp.adressen", __file__)
ZPI = 2j * mp.pi


def gp_ausdruck(z):
    z = mp.mpc(z)
    return f"({mp.nstr(z.real, 30)})+({mp.nstr(z.imag, 30)})*I" if z.imag != 0 else mp.nstr(z.real, 30)


def phase(b, d, xs, R):
    cb, cd = mp.log(b), mp.log(d)
    Td = kneser_sicher("sexp", gp_ausdruck(d), xs, dps=40)
    punkte, stufen, adressen = [], [], []
    for w in Td:
        if w is None:
            punkte.append(None); stufen.append(None); adressen.append(None)
            continue
        w, k, adr = mp.mpc(w), 0, []
        while abs(cd * w) < R:
            adr.append(int(mp.nint(mp.im(cd * w) / (2 * mp.pi))))
            w = mp.exp(cd * w)
            k += 1
            if k > 12:
                break
        u = cd * w
        if mp.re(u) < 0 or abs(u) < R:                     # Bahn fällt zurück oder flieht nicht
            punkte.append(None); stufen.append(None); adressen.append(None)
            continue
        adr.append(int(mp.nint(mp.im(u) / (2 * mp.pi))))
        v = (u - ZPI * adr[-1]) / cb
        for _ in range(k):
            v = mp.log(v) / cb
        punkte.append(v)
        stufen.append(k)
        adressen.append(adr)
    gueltig = [i for i, v in enumerate(punkte) if v is not None]
    sl_g = kneser_sicher("slog", gp_ausdruck(b), [punkte[i] for i in gueltig], dps=40)
    sl = [None] * len(xs)
    for i, w in zip(gueltig, sl_g):
        sl[i] = w
    phi = [None if sl[i] is None else sl[i] - xs[i] for i in range(len(xs))]
    print('   ', len([v for v in phi if v is not None]), 'von', len(xs), 'Punkten mit Wert', flush=True)
    return phi, stufen, adressen, punkte


paare = {"reell b = 3, d = 2 (Kontrolle)": (mp.mpf(3), mp.mpf(2)),
         "b = 4 + 0.5i, d = 3 + 0.3i": (mp.mpc(4, 0.5), mp.mpc(3, 0.3)),
         "b = 8 − 4i, d = 6 + 2i": (mp.mpc(8, -4), mp.mpc(6, 2))}
NX = 24
xs = [mp.mpf(j) / NX for j in range(NX)]
xs1 = [x + 1 for x in xs]
ergebnis = {}
for name, (b, d) in paare.items():
    try:
        phi_a, st_a, adr_a, _ = phase(b, d, xs, mp.mpf(10) ** 12)
        phi_b, st_b, adr_b, _ = phase(b, d, xs, mp.mpf(10) ** 30)
        phi_1, _, _, _ = phase(b, d, xs1, mp.mpf(10) ** 30)
    except Exception as fehler:
        print("  übersprungen:", name, type(fehler).__name__, str(fehler)[:160])
        continue
    ok = [i for i in range(NX) if None not in (phi_a[i], phi_b[i], phi_1[i])]
    if len(ok) < NX // 2:
        print('  zu wenige Werte:', name, len(ok), flush=True)
        ergebnis[name] = {'fehlend': NX - len(ok)}
        continue
    konv = max(abs(phi_a[i] - phi_b[i]) for i in ok)
    period = max(abs(phi_b[i] - phi_1[i]) for i in ok)
    spruenge = max(abs(phi_b[i + 1] - phi_b[i]) for i in range(NX - 1) if i in ok and i + 1 in ok)
    mu = mp.fsum(phi_b[i] for i in ok) / len(ok)
    amplitude = max(abs(phi_b[i] - mu) for i in ok)
    adr_satz = sorted({tuple(adr_b[i]) for i in ok})
    ergebnis[name] = {"phi": phi_b, "ok": ok, "fehlend": NX - len(ok), "konv": konv, "period": period, "sprung": spruenge, "mu": mu, "amp": amplitude, "adressen": adr_satz,
                      "stufen": sorted({st_b[i] for i in ok})}
    print(f"  {name}: fehlend {NX - len(ok)}, μ = {mp.nstr(mu, 10)}, Amplitude {mp.nstr(amplitude, 4)}, R-Abhängigkeit {mp.nstr(konv, 3)}, "
          f"Φ(x+1)−Φ(x) {mp.nstr(period, 3)}, grösster Schritt {mp.nstr(spruenge, 3)}, Adressen {adr_satz[:4]}")

for name, e in list(ergebnis.items()):
    if "phi" not in e:
        lauf.pruefe(f"{name}: Motor liefert genug Werte", False, f"{e['fehlend']} von {NX} fehlen")
        del ergebnis[name]
        continue
    lauf.pruefe(f"{name}: gestoppte Schwänze konvergieren (R = 10¹² vs 10³⁰)", e["konv"] < mp.mpf(10) ** -10, mp.nstr(e["konv"], 3))
    lauf.pruefe(f"{name}: Phase ist 1-periodisch, Φ(x+1) = Φ(x)", e["period"] < mp.mpf(10) ** -10, mp.nstr(e["period"], 3))
    lauf.pruefe(f"{name}: Phase stetig auf dem Gitter (kein Sprung grösser als 10·mittlerer Schritt)",
                e["sprung"] < 10 * max(e["amp"], mp.mpf(10) ** -12) * 2 * mp.pi / NX + mp.mpf(10) ** -12,
                f"μ = {mp.nstr(e['mu'], 12)}, Amplitude {mp.nstr(e['amp'], 4)}, Adressen der Bahnen {e['adressen'][:3]}")
nichtnull = [n for n, e in ergebnis.items() if any(any(s != 0 for s in a) for a in e["adressen"])]
lauf.pruefe("mindestens ein komplexes Paar mit Adressziffern ≠ 0 untersucht", len(nichtnull) >= 1, nichtnull)

lauf.daten = {n: {"fehlend": e["fehlend"], "mu": mp.nstr(e["mu"], 15), "amplitude": mp.nstr(e["amp"], 6), "R_abh": mp.nstr(e["konv"], 3), "periode": mp.nstr(e["period"], 3),
                  "adressen": [list(a) for a in e["adressen"]], "stufen": e["stufen"]} for n, e in ergebnis.items()}

fig, ax = plt.subplots(1, 2, figsize=(10.5, 4))
for name, e in ergebnis.items():
    ax[0].plot([float(xs[i]) for i in e["ok"]], [float(mp.re(e["phi"][i] - e["mu"])) for i in e["ok"]], "o-", ms=3, label=name)
    ax[1].plot([float(xs[i]) for i in e["ok"]], [float(mp.im(e["phi"][i] - e["mu"])) for i in e["ok"]], "o-", ms=3, label=name)
ax[0].set_title("Re Φ(x) − μ"); ax[1].set_title("Im Φ(x) − μ")
for a_ in ax:
    a_.set_xlabel("x"); a_.legend(fontsize=7)
lauf.bild(fig)
lauf.zusammenfassung = "; ".join(f"{n}: μ = {mp.nstr(e['mu'], 10)}, Φ 1-periodisch bis {mp.nstr(e['period'], 2)}, Adressen {e['adressen'][:2]}"
                                 for n, e in ergebnis.items())
lauf.schreibe()
