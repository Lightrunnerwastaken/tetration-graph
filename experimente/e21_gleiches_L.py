"""exp.gleiches_L — Basen mit gemeinsamem Fixpunkt L und Basiswechsel an diesem Fixpunkt (Rekonstruktion, Aug. 2026).

 Familie: f_u(z) = exp(u z) hat L als Fixpunkt ⇔ u = (Log L + 2πik)/L, Multiplikator λ_k = Log L + 2πik.
 Reelle Basen (Hauptzweig, b > 0: Im u = 0; b < 0: Im u = π) teilen ein L genau dann, wenn
   Re(1/L) = 0 (gleiches Vorzeichen) bzw. Re(1/L) = ∓1/(2m) (eine positive, eine negative Basis).
 Basiswechsel: T_d(x) = Ψ_d(λ_d^x w_d) regulär an L, Abel-Funktion α_b = log σ_b / log λ_b an L. Dann
   α_b(T_d(x)) = C1 + C2 x + a1 λ_d^x w_d / log λ_b + O(λ_d^{2x}),  C2 = log λ_d / log λ_b,
   a1 = u_d/(2(λ_d − 1)) − u_b/(2(λ_b − 1)).   Keine bleibende periodische Welle; exakt linear nur für b = d.
"""
import cmath
import math

import mpmath as mp
import numpy as np
from scipy.special import lambertw

from fixpunkt import AmFixpunkt
from gemeinsam import Lauf, plt

mp.mp.dps = 40
lauf = Lauf("exp.gleiches_L", __file__)
ZWEI_PI_I = 2j * mp.pi


def periode(u, n_vor=6000, qmax=12):
    """Verhalten des Turms (Hauptzweig): Periode, −1 Flucht, 0 unklar."""
    z = 1 + 0j
    u = complex(u)
    try:
        for _ in range(n_vor):
            z = cmath.exp(u * z)
            if abs(z) > 1e300:
                return -1
        bahn = [z]
        for _ in range(qmax):
            bahn.append(cmath.exp(u * bahn[-1]))
    except OverflowError:
        return -1
    for q in range(1, qmax + 1):
        if abs(bahn[q] - bahn[0]) < 1e-8 * (1 + abs(bahn[0])):
            return q
    return 0


def name_periode(q):
    return {-1: "Flucht", 0: "unklar"}.get(q, f"Periode {q}")


# ------------------------------------------------------------------ (A) Familie zu einem L
Le = -mp.lambertw(-1, -1)                          # oberer Fixpunkt von e (Log L = L)
fam = [(mp.log(Le) + ZWEI_PI_I * k) / Le for k in range(-3, 4)]
res = max(abs(mp.exp(u * Le) - Le) for u in fam)
lauf.pruefe("Familie: u_k = (Log L + 2πik)/L haben alle denselben Fixpunkt L (L = L_e, k = −3…3)", res < mp.mpf(10) ** -35, mp.nstr(res, 3))
lam_fam = [Le * u for u in fam]
lauf.pruefe("Multiplikatoren der Familie unterscheiden sich um 2πi", max(abs(lam_fam[i + 1] - lam_fam[i] - ZWEI_PI_I) for i in range(6)) < 1e-35)

# ------------------------------------------------------------------ (B) reelle Basen mit gemeinsamem L
reell = []
# positiv–positiv: L = ±i, b = e^{π/2 + 2πk}
for k in range(-2, 3):
    reell.append({"art": "positiv, L = i", "b": mp.e ** (mp.pi / 2 + 2 * mp.pi * k), "u": mp.pi / 2 + 2 * mp.pi * k, "L": mp.mpc(0, 1)})
# negativ–negativ: L = i·W(π)/π, b = −exp((π/2 + 2πk)·π/W(π))
y = mp.lambertw(mp.pi).real / mp.pi
for k in range(-1, 2):
    lnr = (mp.pi / 2 + 2 * mp.pi * k) / y
    reell.append({"art": "negativ, L = i·W(π)/π", "b": -mp.e ** lnr, "u": mp.mpc(lnr, mp.pi), "L": mp.mpc(0, y)})
res = max(abs(mp.exp(r["u"] * r["L"]) - r["L"]) for r in reell)
lauf.pruefe("reelle Paare: e^{π/2+2πk} teilen L = i; −exp((π/2+2πk)π/W(π)) teilen L = i·W(π)/π", res < mp.mpf(10) ** -30, mp.nstr(res, 3))

# gemischt: b = e^β > 0 und d = −s < 0 mit gemeinsamem L ⇔ Re(1/L) = −1/(2m), dann ln s = β − Re(2πim/L)
gemischt = []
for m in (1, 2, -1, -2):
    for kz in range(-2, 3):
        betas = np.linspace(-6, 6, 4801)
        betas = betas[np.abs(betas) > 1e-6]
        with np.errstate(all="ignore"):
            Lk = -lambertw(-betas, kz) / betas
            F = (1 / Lk).real + 1 / (2 * m)
        for i in np.where(np.sign(F[:-1]) * np.sign(F[1:]) < 0)[0]:
            if not (np.isfinite(F[i]) and np.isfinite(F[i + 1])) or abs(F[i] - F[i + 1]) > 0.5:
                continue
            fb = lambda be: (1 / (-mp.lambertw(-be, kz) / be)).real + mp.mpf(1) / (2 * m)
            try:
                be = mp.findroot(fb, (mp.mpf(betas[i]), mp.mpf(betas[i + 1])), solver="anderson")
            except (ValueError, ZeroDivisionError):
                continue
            L = -mp.lambertw(-be, kz) / be
            lns = be - (ZWEI_PI_I * m / L).real
            ud = mp.mpc(lns, mp.pi)
            if abs(mp.exp(ud * L) - L) < mp.mpf(10) ** -25 and abs(mp.exp(be * L) - L) < mp.mpf(10) ** -25:
                gemischt.append({"art": f"gemischt m={m}, Zweig {kz}", "b": mp.e ** be, "u": be, "d": -mp.e ** lns, "ud": ud, "L": L})
lauf.pruefe("gemischte Paare (positive und negative reelle Basis mit gleichem L) gefunden und nachgerechnet", len(gemischt) >= 2,
            f"{len(gemischt)} Paare")

kategorien = []
for r in reell:
    q = periode(complex(r["u"]))
    kategorien.append({"basis": mp.nstr(r["b"], 8), "L": mp.nstr(r["L"], 8), "lam": mp.nstr(abs(r["L"] * r["u"]), 5), "turm": name_periode(q)})
for g in gemischt:
    kategorien.append({"basis": mp.nstr(g["b"], 8), "partner": mp.nstr(g["d"], 8), "L": mp.nstr(g["L"], 8),
                       "turm_b": name_periode(periode(complex(g["u"]))), "turm_d": name_periode(periode(complex(g["ud"])))})
for z in kategorien:
    print("  ", z)

# ------------------------------------------------------------------ (C) Basiswechsel am gemeinsamen Fixpunkt
def basiswechsel(ub, ud, L, xs, normiere=True):
    B, D = AmFixpunkt(ub, L), AmFixpunkt(ud, L)
    w = mp.mpc(1)
    if normiere:
        try:
            w = D.sigma(mp.mpc(1))              # T_d(0) = 1
        except ValueError:
            normiere = False
    G, rest = [], []
    for x in xs:
        v = mp.exp(x * D.log_lam) * w
        T = D.Psi(v)
        hv = B.sigma(T)
        G.append((x * D.log_lam + mp.log(w) + mp.log(hv / v)) / B.log_lam)
        rest.append(mp.log(hv / v) / B.log_lam)
    C1, C2 = mp.log(w) / B.log_lam, D.log_lam / B.log_lam
    a1 = D.u / (2 * (D.lam - 1)) - B.u / (2 * (B.lam - 1))
    pred = [a1 * mp.exp(x * D.log_lam) * w / B.log_lam for x in xs]
    return {"B": B, "D": D, "G": G, "rest": rest, "pred": pred, "C1": C1, "C2": C2, "normiert": normiere, "T0": D.Psi(w)}


xs = [mp.mpf(x) / 4 for x in range(-48, -3)]          # x = −12 … −1
paare = {"e^{π/2} ← e^{5π/2} (L = i)": (mp.pi / 2, 5 * mp.pi / 2, mp.mpc(0, 1)),
         "e^{π/2} ← e^{−3π/2} (L = i, Partner links)": (mp.pi / 2, -3 * mp.pi / 2, mp.mpc(0, 1)),
         "√2 ← −√2 (L = 2)": (mp.log(2) / 2, mp.mpc(mp.log(2) / 2, mp.pi), mp.mpf(2)),
         "e ← e^{1−2πi/L_e} (L = L_e)": (mp.mpf(1), fam[2], Le)}
if gemischt:
    g0 = min(gemischt, key=lambda g: abs(g["b"] - 3))
    paare[f"{mp.nstr(g0['b'], 5)} ← {mp.nstr(g0['d'], 5)} (gemischt)"] = (g0["u"], g0["ud"], g0["L"])
ergebnisse = {}
for name, (ub, ud, L) in paare.items():
    try:
        r = basiswechsel(ub, ud, L, xs)
    except ValueError as fehler:
        print("  übersprungen:", name, fehler)
        continue
    ergebnisse[name] = r
    D = r["D"]
    # (1) Linearität: Rest fällt wie |λ_d|^x, (2) erste Korrektur stimmt
    i1, i2 = 4, 20                                    # x = −11, −7
    rate = (mp.log(abs(r["rest"][i2])) - mp.log(abs(r["rest"][i1]))) / (xs[i2] - xs[i1])
    quot = max(abs(r["rest"][i] / r["pred"][i] - 1) for i in range(0, 8))
    lauf.pruefe(f"{name}: α_b(T_d(x)) − (C1 + C2 x) fällt wie |λ_d|^x (Rate {mp.nstr(rate, 6)} vs ln|λ_d| = {mp.nstr(mp.log(abs(D.lam)), 6)})",
                abs(rate - mp.log(abs(D.lam))) < 0.02, f"Rest bei x = −12: {mp.nstr(abs(r['rest'][0]), 3)}")
    lauf.pruefe(f"{name}: erste Korrektur a1·λ_d^x·w/log λ_b stimmt (x ≤ −10)", quot < 0.01, mp.nstr(quot, 3))
    # (3) lineare Ausgleichsgerade auf [−12, −8] trifft C2 = log λ_d/log λ_b
    n = 17
    X = mp.matrix([[1, xs[i]] for i in range(n)])
    Y = mp.matrix([r["G"][i] for i in range(n)])
    koef = mp.lu_solve(X.T * X, X.T * Y)
    lauf.pruefe(f"{name}: Ausgleichsgerade auf [−12, −8] hat Steigung log λ_d/log λ_b", abs(koef[1] - r["C2"]) < mp.mpf(10) ** -6,
                f"C2 = {mp.nstr(r['C2'], 10)}, Fit {mp.nstr(koef[1], 10)}")

# Gegenprobe: verschiedene Fixpunkte (b = e^{π/2} mit L = i, d = e mit L_e)
B = AmFixpunkt(mp.pi / 2, mp.mpc(0, 1))
D = AmFixpunkt(mp.mpf(1), Le)
steig = []
for x in (mp.mpf(-20), mp.mpf(-40)):
    werte = [mp.log(B.sigma(D.Psi(mp.exp((x + s) * D.log_lam)))) / B.log_lam for s in (0, mp.mpf("0.01"))]
    steig.append(abs((werte[1] - werte[0]) / mp.mpf("0.01")))
erwartet = abs(D.lam) ** -20
lauf.pruefe("Gegenprobe verschiedene Fixpunkte (e^{π/2} gegen e): Steigung fällt wie |λ_e|^x gegen 0 — α_b(T_d(x)) wird konstant, nicht linear",
            steig[1] < 1e-4 and abs(steig[1] / steig[0] / erwartet - 1) < 0.05, [mp.nstr(s, 3) for s in steig])

# Exakte Linearität nur für b = d: sonst wäre f_d eine gebrochene Iterierte von f_b; Kommutator ≠ 0 als Zeuge
zb, zd = mp.pi / 2, 5 * mp.pi / 2
z0 = mp.mpc("0.3", "0.9")
komm = abs(mp.exp(zd * mp.exp(zb * z0)) - mp.exp(zb * mp.exp(zd * z0)))
lauf.pruefe("exakt linear unmöglich für b ≠ d: f_b und f_d vertauschen nicht (Zeuge bei z = 0.3 + 0.9i)", komm > 1e-3, mp.nstr(komm, 4))

lauf.daten = {"kategorien": kategorien,
              "paare": {k: {"C1": mp.nstr(v["C1"], 15), "C2": mp.nstr(v["C2"], 15), "T_d(0)": mp.nstr(v["T0"], 12), "normiert": v["normiert"],
                            "rest_x-12": mp.nstr(abs(v["rest"][0]), 4), "rest_x-1": mp.nstr(abs(v["rest"][-1]), 4)} for k, v in ergebnisse.items()},
              "gegenprobe_steigung": [mp.nstr(s, 4) for s in steig]}

# ------------------------------------------------------------------ Bild
fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
th = np.linspace(0, 2 * np.pi, 800)
tt = np.exp(1j * th)
rand = np.exp(tt * np.exp(-tt))
ax[0].plot(rand.real, rand.imag, color="#b5523b", lw=1.1, label="Shell–Thron-Rand")
for r in reell:
    b = complex(r["b"])
    if abs(b) < 6:
        ax[0].plot([b.real], [b.imag], "o", ms=5, color="#3d6fb6" if b.real > 0 else "#2f8f83")
        ax[0].annotate(mp.nstr(r["b"], 4), (b.real, b.imag), fontsize=6, xytext=(3, 4), textcoords="offset points")
for g in gemischt:
    for b in (complex(g["b"]), complex(g["d"])):
        if abs(b) < 6:
            ax[0].plot([b.real], [b.imag], "s", ms=4, mfc="none", color="#7a4fb0")
famb = [complex(mp.exp(u)) for u in fam]
ax[0].plot([b.real for b in famb if abs(b) < 6], [b.imag for b in famb if abs(b) < 6], "^", ms=5, color="#c08a2b", label="Familie zu L_e")
ax[0].set_xlim(-4, 6); ax[0].set_ylim(-3, 3); ax[0].set_aspect("equal")
ax[0].set_title("Basen mit gemeinsamem Fixpunkt"); ax[0].legend(fontsize=7)
for name, r in ergebnisse.items():
    ax[1].plot([float(x) for x in xs], [float(mp.re(gg)) for gg in r["G"]], label=f"Re α_b(T_d(x)): {name}")
    ax[2].semilogy([float(x) for x in xs], [float(abs(q)) + 1e-40 for q in r["rest"]], label=name)
    ax[2].semilogy([float(x) for x in xs], [float(abs(q)) for q in r["pred"]], ":", color="#888")
ax[1].set_xlabel("x"); ax[1].set_title("Basiswechsel am gemeinsamen L: eine Gerade"); ax[1].legend(fontsize=6)
ax[2].set_xlabel("x"); ax[2].set_title("|Abweichung von C1 + C2·x| (gepunktet: erste Korrektur)"); ax[2].legend(fontsize=6)
lauf.bild(fig)
lauf.zusammenfassung = ("Basen mit gemeinsamem L bilden Familien u_k = (Log L + 2πik)/L. Reell: e^{π/2+2πk} (L = i), "
                        f"−exp((π/2+2πk)π/W(π)) (L = i·W(π)/π) und {len(gemischt)} gemischte Paare. Am gemeinsamen L ist α_b(T_d(x)) = C1 + C2 x "
                        "+ O(λ_d^x), C2 = log λ_d/log λ_b — keine bleibende Welle; mit verschiedenen L wird α_b(T_d) konstant.")
lauf.schreibe()
