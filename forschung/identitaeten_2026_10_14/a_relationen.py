"""Etappe 1, Teil A4: Suche nach ganzzahligen Relationen (LLL über python-flint, Kontrolle mit mpmath.pslq/findpoly/identify).

Für jede Zielkonstante x (Genauigkeit P_ok aus dem Abgleich zweier Rechenstufen):
  lin   : x gegen die klassische Basis (22 Zahlen)                        → Schranke für die Norm jeder Relation
  klein : x gegen kleine Teilbasen (höhere Koeffizienten möglich)
  log   : ln|x| gegen Logarithmen der Basis (Produkte von Potenzen)
  poly  : 1, x, …, x^6 (algebraisch vom Grad ≤ 6)
Suche mit P_find Stellen, Treffer zählt nur, wenn der Rest bei P_ok ≥ P_find + 15 Stellen klein bleibt.
Schranke: nach LLL ist jeder Gittervektor ≥ min‖b_i*‖; eine echte Relation m ergibt einen Vektor der Länge ≤ ‖m‖·√(1+n/4),
also ‖m‖₂ ≥ min‖b_i*‖/√(1+n/4).
Aufruf: python a_relationen.py   → daten/relationen.json, Kurzfassung auf stdout
"""
import json
import math
import sys
import time
from pathlib import Path

import flint
import mpmath as mp

HIER = Path(__file__).resolve().parent
D = HIER / "daten"
mp.mp.dps = 260
t0 = time.time()


def log(*a):
    print(f"[{time.time()-t0:6.1f}s]", *a, flush=True)


# ------------------------------------------------------------------ Zielkonstanten und ihre Genauigkeit
def lade(dps):
    p = D / f"konstanten_{dps}.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else None


stufen = {d: lade(d) for d in (120, 160, 200)}
stufen = {d: v for d, v in stufen.items() if v}
hoch = max(stufen)
tief = sorted(stufen)[-2] if len(stufen) > 1 else None
ziele = {}
for k, v in stufen[hoch].items():
    if k.startswith("_"):
        continue
    x = mp.mpf(v)
    if tief and k in stufen[tief]:
        d = abs(x - mp.mpf(stufen[tief][k]))
        p_tief = int(-mp.log10(d / max(abs(x), 1))) if d > 0 else tief
    else:
        p_tief = None
    ziele[k] = {"x": x, "p_abgleich": p_tief}

# Genauigkeit der höchsten Stufe: Fork liefert ≥ 0.95·dps (gp-tetration METHODS); mit 200-Stufe gemessen, sonst konservativ
P_MOTOR = {120: 108, 160: 140, 200: 160}[hoch]   # 200 gegen 160: Abweichung ≤ 4e−159 bei allen Grössen (gemessen) if hoch in (120, 160, 200) else 100
fatou = json.loads((D / "fatou.json").read_text(encoding="utf-8"))
for k in ("xexp.alpha(1)", "xexp.halb(1)", "expm1.alpha(1)", "expm1.halb(1)"):
    ziele["fatou." + k] = {"x": mp.mpf(fatou[k]), "p_abgleich": 180}
th = D / f"theta_{hoch}.json"
if th.is_file():
    tj = json.loads(th.read_text(encoding="utf-8"))
    tt = json.loads((D / "theta_160.json").read_text(encoding="utf-8"))
    for n in range(9):      # Genauigkeit aus dem Vergleich der Motorstufen 160/200 statt aus der inneren Konsistenz
        d = abs(mp.mpc(*map(mp.mpf, tj[f"a_{n}"])) - mp.mpc(*map(mp.mpf, tt[f"a_{n}"])))
        if n == 0:
            d = min(abs(d + k * 2j * mp.pi / mp.conj(-mp.lambertw(-1))) for k in range(-3, 4))
        tj[f"abgleich_a_{n}"] = {"N": mp.nstr(d, 5), "y": "0"}
        print("theta a_%d: |a(200) - a(160)| =" % n, mp.nstr(d, 3))
    Lp = mp.conj(-mp.lambertw(-1))
    per = 2j * mp.pi / Lp
    a0 = mp.mpc(*map(mp.mpf, tj["a_0"]))
    a0 = min((a0 + k * per for k in range(-200, 201)), key=abs)       # Vertreter modulo 2πi/L mit kleinstem Betrag
    p0 = int(-mp.log10(max(mp.mpf(tj["abgleich_a_0"]["N"]), mp.mpf(tj["abgleich_a_0"]["y"]), mp.mpf(10) ** -200)))
    ziele["theta.Re a_0(red)"] = {"x": mp.re(a0), "p_abgleich": p0}
    ziele["theta.Im a_0(red)"] = {"x": mp.im(a0), "p_abgleich": p0}
    for n in (1, 2, 3):
        re_, im_ = map(mp.mpf, tj[f"a_{n}"])
        p = int(-mp.log10(max(mp.mpf(tj[f"abgleich_a_{n}"]["N"]), mp.mpf(tj[f"abgleich_a_{n}"]["y"]), mp.mpf(10) ** -200)))
        ziele[f"theta.Re a_{n}"] = {"x": re_, "p_abgleich": p}
        ziele[f"theta.Im a_{n}"] = {"x": im_, "p_abgleich": p}
        ziele[f"theta.|a_{n}|"] = {"x": abs(mp.mpc(re_, im_)), "p_abgleich": p}
        ziele[f"theta.arg a_{n}"] = {"x": mp.arg(mp.mpc(re_, im_)), "p_abgleich": p}

# redundante/triviale Einträge weg
for k in list(ziele):
    if k in ("e.taylor0_a0",) or "GL110" in k or k.endswith("_k5") or "F(1)_" in k or "F(2)_" in k:
        ziele.pop(k)
for k, z in ziele.items():
    pa = z["p_abgleich"]
    z["P_ok"] = min(P_MOTOR, pa if pa else P_MOTOR) if not k.startswith("fatou.") else 175
    if k.startswith("theta."):
        z["P_ok"] = min(z["P_ok"], pa)

# ------------------------------------------------------------------ Basen
L = mp.conj(-mp.lambertw(-1))      # Im L > 0
L2 = -mp.lambertw(-mp.log(2)) / mp.log(2)          # Fixpunkt von 2^z
K = {
    "1": mp.mpf(1), "pi": mp.pi, "e": mp.e, "gamma": mp.euler, "ln2": mp.log(2), "ln3": mp.log(3), "lnpi": mp.log(mp.pi),
    "zeta3": mp.zeta(3), "G": mp.catalan, "Gamma(1/4)": mp.gamma(mp.mpf(1) / 4), "Gamma(1/3)": mp.gamma(mp.mpf(1) / 3),
    "Omega": mp.lambertw(1).real, "ReL": mp.re(L), "ImL": mp.im(L), "|L|": abs(L), "lngamma": mp.log(mp.euler),
    "lnzeta3": mp.log(mp.zeta(3)), "lnG": mp.log(mp.catalan), "lnGamma(1/4)": mp.log(mp.gamma(mp.mpf(1) / 4)),
    "lnGamma(1/3)": mp.log(mp.gamma(mp.mpf(1) / 3)), "e^e": mp.e ** mp.e, "pi^2": mp.pi ** 2,
}
K2 = {"ReL2": mp.re(L2), "ImL2": mp.im(L2), "|L2|": abs(L2), "lnln2": mp.log(mp.log(2)), "1/ln2": 1 / mp.log(2)}
LOGB = {   # Logarithmen für Produktrelationen x = ∏ b^m
    "ln2": mp.log(2), "ln3": mp.log(3), "ln5": mp.log(5), "lnpi": mp.log(mp.pi), "1(=ln e)": mp.mpf(1), "lngamma": mp.log(mp.euler),
    "lnzeta3": mp.log(mp.zeta(3)), "lnG": mp.log(mp.catalan), "lnGamma(1/4)": mp.log(mp.gamma(mp.mpf(1) / 4)),
    "lnGamma(1/3)": mp.log(mp.gamma(mp.mpf(1) / 3)), "lnOmega(=-Omega)": -mp.lambertw(1).real, "ln|L|(=ReL)": mp.re(L),
    "lnImL": mp.log(mp.im(L)), "e(=ln e^e)": mp.e,
}
KLEIN = {
    "elementar": ["1", "pi", "e", "gamma", "ln2"],
    "fixpunkt": ["1", "pi", "ReL", "ImL", "|L|"],
    "omega": ["1", "Omega", "e", "pi", "gamma"],
    "zeta_catalan": ["1", "pi", "zeta3", "G", "pi^2"],
}


# ------------------------------------------------------------------ LLL-Relationssuche
def lll_relation(xs, P_find):
    """xs: Liste mpf. Rückgabe (kürzester Kandidat m, Rest, Schranke für ‖m‖₂)."""
    n = len(xs)
    C = mp.mpf(10) ** P_find
    rows = []
    for i, x in enumerate(xs):
        r = [0] * n + [int(mp.nint(C * x))]
        r[i] = 1
        rows.append(r)
    M = flint.fmpz_mat(rows).lll()
    B = [[int(M[i, j]) for j in range(n + 1)] for i in range(n)]
    # Gram–Schmidt-Normen
    with mp.workdps(P_find + 60):
        V = [[mp.mpf(v) for v in row] for row in B]
        Bs = []
        norms = []
        for v in V:
            w = list(v)
            for u, uu in zip(Bs, norms):
                c = mp.fsum(a * b for a, b in zip(v, u)) / uu
                w = [a - c * b for a, b in zip(w, u)]
            Bs.append(w)
            norms.append(mp.fsum(a * a for a in w))
        schranke = mp.sqrt(min(norms)) / mp.sqrt(1 + mp.mpf(n) / 4)
    m = B[0][:n]
    return m, schranke


def pruefe_relation(m, xs_voll, P_ok):
    with mp.workdps(P_ok + 20):
        r = abs(mp.fsum(mi * x for mi, x in zip(m, xs_voll)))
        nm = mp.sqrt(sum(mi * mi for mi in m))
        return r, nm


ERG = {}


def suche(name, ziel, basis_namen, basis_werte, art):
    z = ziele[ziel]
    P_ok = z["P_ok"]
    P_find = P_ok - 20
    xs = [z["x"]] + basis_werte
    m, schranke = lll_relation(xs, P_find)
    r, nm = pruefe_relation(m, xs, P_ok)
    rel = (r < mp.mpf(10) ** (-(P_find + 15)) * max(nm, 1)) and m[0] != 0
    eintrag = {"art": art, "basis": basis_namen, "P_find": P_find, "P_ok": P_ok,
               "schranke_norm": mp.nstr(schranke, 4), "log10_schranke": float(mp.log10(schranke)),
               "kandidat_norm": mp.nstr(nm, 4), "kandidat_rest": mp.nstr(r, 3), "relation": rel}
    if rel:
        eintrag["koeffizienten"] = m
        log("!!! RELATION", ziel, art, list(zip(["x"] + basis_namen, m)), "Rest", mp.nstr(r, 3))
    ERG.setdefault(ziel, []).append(eintrag)
    return eintrag


def basis_fuer(ziel):
    b = dict(K)
    if ziel.startswith("2."):
        b.update(K2)
    return b


# Positivkontrollen: bekannte Relationen müssen gefunden werden
with mp.workdps(300):
    ziele["kontrolle.lin 3pi-2gamma+5ln2/7+ImL/11"] = {"x": 3 * mp.pi - 2 * mp.euler + 5 * mp.log(2) / 7 + mp.im(L) / 11, "P_ok": 160, "p_abgleich": None}
    ziele["kontrolle.log Gamma(1/4)^2 e^3/pi^(3/2)"] = {"x": mp.gamma(mp.mpf(1) / 4) ** 2 * mp.e ** 3 / mp.pi ** 1.5, "P_ok": 160, "p_abgleich": None}
    ziele["kontrolle.poly x^5-x-1"] = {"x": mp.findroot(lambda t: t ** 5 - t - 1, 1.17), "P_ok": 160, "p_abgleich": None}
for ziel in list(ziele):
    b = basis_fuer(ziel)
    e = suche(ziel, ziel, list(b), list(b.values()), "lin")
    for kn, namen in KLEIN.items():
        suche(ziel, ziel, namen, [b[q] for q in namen], "klein:" + kn)
    x = ziele[ziel]["x"]
    if x != 0:
        lb = dict(LOGB)
        if ziel.startswith("2."):
            lb["lnln2"] = mp.log(mp.log(2))
        ziele["_ln"] = {"x": mp.log(abs(x)), "P_ok": ziele[ziel]["P_ok"]}
        e2 = suche(ziel, "_ln", list(lb), list(lb.values()), "log")
        ERG[ziel] = ERG.get(ziel, []) + ERG.pop("_ln")
        ziele.pop("_ln")
    # algebraisch Grad ≤ 6
    ziele["_poly"] = {"x": x, "P_ok": ziele[ziel]["P_ok"]}
    for deg in (2, 4, 6):
        pw = [x ** k for k in range(2, deg + 1)]
        e3 = suche(ziel, "_poly", ["1"] + [f"x^{k}" for k in range(2, deg + 1)], [mp.mpf(1)] + pw, f"poly{deg}")
    ERG[ziel] = ERG.get(ziel, []) + ERG.pop("_poly")
    ziele.pop("_poly")
    kurz = {e["art"]: round(e["log10_schranke"], 1) for e in ERG[ziel]}
    log(ziel, "P_ok", ziele[ziel]["P_ok"], "log10-Schranken", kurz)

# ------------------------------------------------------------------ gemeinsame Relation unter den Kneser-Zahlen (Basis e)
gemeinsam = [k for k in ziele if k.startswith("e.") and ziele[k]["P_ok"] >= 100][:18]
ziele["_joint"] = {"x": mp.mpf(1), "P_ok": min(ziele[k]["P_ok"] for k in gemeinsam)}
werte = [ziele[k]["x"] for k in gemeinsam] + [mp.pi, mp.re(L), mp.im(L)]
e = suche("joint", "_joint", gemeinsam + ["pi", "ReL", "ImL"], werte, "gemeinsam")
ERG["gemeinsam(1 + Kneser-e-Zahlen)"] = ERG.pop("_joint")
log("gemeinsam", e["log10_schranke"], e["relation"])

# ------------------------------------------------------------------ Gegenprobe mit mpmath.identify bei 15 Stellen → bei voller Genauigkeit prüfen
near = {}
konst = {"pi": mp.pi, "e": mp.e, "euler": mp.euler, "log(2)": mp.log(2), "catalan": mp.catalan, "zeta(3)": mp.zeta(3)}
for ziel, z in [(k, v) for k, v in ziele.items() if not k.startswith("_")]:
    with mp.workdps(15):
        try:
            cands = mp.identify(z["x"], ["pi", "e", "euler", "log(2)"], tol=mp.mpf(10) ** -12)
        except Exception:
            cands = None
    if cands:
        with mp.workdps(z["P_ok"] + 10):
            loc = dict(konst)
            loc.update({"pi": mp.pi, "e": mp.e, "euler": mp.euler, "catalan": mp.catalan})
            try:
                wert = eval(cands.replace("^", "**"), {"__builtins__": {}}, {**{k: v for k, v in loc.items()}, "log": mp.log, "sqrt": mp.sqrt, "zeta": mp.zeta, "exp": mp.exp})
                d = abs(wert - z["x"])
                near[ziel] = {"identify@15": cands, "abweichung_voll": mp.nstr(d, 3)}
            except Exception as ex:
                near[ziel] = {"identify@15": cands, "fehler": str(ex)}
# wörtliche Beinahe-Treffer: x ≈ (p/q)·b mit kleinen p, q
beinahe = {}
for ziel, z in [(k, v) for k, v in ziele.items() if not k.startswith("_")]:
    best = []
    for bn, bv in K.items():
        for q in range(1, 13):
            for p in range(-36, 37):
                if p == 0:
                    continue
                c = mp.mpf(p) / q * bv
                rel = abs(z["x"] - c) / abs(z["x"])
                if rel < mp.mpf("2e-4"):
                    best.append((float(rel), f"{p}/{q}·{bn}" if bn != "1" else f"{p}/{q}"))
    best.sort()
    if best:
        beinahe[ziel] = [(f"{r:.2e}", s) for r, s in best[:3]]
log("identify-Kandidaten:", near)
log("Beinahe-Treffer:", beinahe)

out = {"P_motor_annahme": P_MOTOR, "stufen": list(stufen),
       "ziele": {k: {"wert": mp.nstr(z["x"], 60), "P_ok": z["P_ok"], "p_abgleich": z.get("p_abgleich")} for k, z in ziele.items() if not k.startswith("_")},
       "suche": ERG, "identify": near, "beinahe": beinahe}
(D / "relationen.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=str), encoding="utf-8")
log("fertig")
