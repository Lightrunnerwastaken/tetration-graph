"""Ideen-Landkarte: Kopplung, Lage, Haupt-/Hilfsaussagen, Füllstand und Herkunft — reine Rechnung, deterministisch.

Ergebnis von berechne() geht als `ideen` in window.GRAPH; text_md() ist dieselbe Struktur als Text für Agenten.
"""
import math
from collections import defaultdict

STARK = {"stuetzt_sich_auf", "verallgemeinert", "ersetzt", "korrigiert", "widerlegt"}
MITTEL = {"stuetzt", "motiviert"}
GRUPPEN = {
    "bewiesen": {"definition", "bewiesen", "bewiesen_literatur", "entlastet"},
    "bedingt": {"bedingt", "skizze"},
    "numerisch": {"numerisch"},
    "vermutung": {"vermutung", "annahme"},
    "offen": {"offen"},
    "verworfen": {"widerlegt", "gescheitert", "ueberholt"},
}
GRUPPE_VON = {s: g for g, ss in GRUPPEN.items() for s in ss}
GRUPPE_NAME = {"bewiesen": "bewiesen", "bedingt": "bedingt", "numerisch": "numerisch", "vermutung": "Vermutung",
               "offen": "offen", "verworfen": "verworfen"}
HAUPT_ARTEN = {"satz", "vermutung", "frage", "konstruktion", "definition"}
GOLD = 2.399963229728653
OHNE = "_ohne"


def _gewicht(typ):
    return 1.0 if typ in STARK else 0.5 if typ in MITTEL else 0.25


def kopplung(g, unter_von):
    """w[(U, V)] mit U < V: Kanten zwischen den Aussagen zweier Unterideen plus Anteil gemeinsamer Objekte."""
    w = defaultdict(float)
    for k in g.kanten:
        u, v = unter_von.get(k["von"]), unter_von.get(k["nach"])
        if u and v and u != v:
            w[tuple(sorted((u, v)))] += _gewicht(k["typ"])
    obj, n = defaultdict(set), defaultdict(int)
    for i, u in unter_von.items():
        if u:
            n[u] += 1
            obj[u] |= set(g.knoten[i].get("objekte") or [])
    ids = sorted(n)
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            U, V = ids[a], ids[b]
            if obj[U] and obj[V]:
                j = len(obj[U] & obj[V]) / len(obj[U] | obj[V])
                if j > 0:
                    w[(U, V)] += 0.5 * j * min(n[U], n[V])
    return dict(w)


def _gruppe(a):
    return GRUPPE_VON.get(a.get("_status_eff") or a.get("status"), "offen")


def _herkunft(a, quellen):
    arten = {(quellen.get(q.get("ref")) or {}).get("art") for q in (a.get("quellen") or []) if isinstance(q, dict)}
    return "literatur" if "literatur" in arten else "forum" if arten & {"forum", "notiz"} else "eigen"


def _zaehle(werte):
    z = defaultdict(int)
    for v in werte:
        z[v] += 1
    return dict(sorted(z.items()))


def _trennen(P, R, luecke, schritte):
    """Überlappungen auflösen (Paare auf Abstand r_i + r_j + luecke schieben)."""
    for _ in range(schritte):
        bewegt = False
        for i in range(len(P)):
            for j in range(i + 1, len(P)):
                dx, dy = P[j][0] - P[i][0], P[j][1] - P[i][1]
                d = math.hypot(dx, dy) or 1e-9
                soll = R[i] + R[j] + luecke
                if d < soll:
                    s = (soll - d) / 2 / d
                    P[i][0] -= dx * s; P[i][1] -= dy * s
                    P[j][0] += dx * s; P[j][1] += dy * s
                    bewegt = True
        if not bewegt:
            break


def _layout_haupt(ids, R, w):
    n = len(ids)
    P = [[math.cos(2 * math.pi * i / n) * 400 * math.sqrt(n), math.sin(2 * math.pi * i / n) * 400 * math.sqrt(n)] for i in range(n)]
    for t in range(600):
        lr = 0.5 * (1 - t / 600) + 0.02
        F = [[0.0, 0.0] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                dx, dy = P[j][0] - P[i][0], P[j][1] - P[i][1]
                d = math.hypot(dx, dy) or 1e-9
                wij = w.get((ids[i], ids[j]), 0) + w.get((ids[j], ids[i]), 0)
                f = 0.0
                if wij > 0:
                    f += 0.08 * (d - (R[i] + R[j] + 60 / (1 + wij / 5))) * min(1.0, 0.25 + wij / 20)
                if d < R[i] + R[j] + 40:
                    f -= 0.5 * (R[i] + R[j] + 40 - d)
                F[i][0] += f * dx / d; F[i][1] += f * dy / d
                F[j][0] -= f * dx / d; F[j][1] -= f * dy / d
        for i in range(n):
            F[i][0] -= 0.002 * P[i][0]; F[i][1] -= 0.002 * P[i][1]
            P[i][0] += lr * F[i][0]; P[i][1] += lr * F[i][1]
    _trennen(P, R, 40, 400)
    cx, cy = sum(p[0] for p in P) / n, sum(p[1] for p in P) / n
    return [[p[0] - cx, p[1] - cy] for p in P]


def _layout_unter(ids, r, Rh, w_innen, zug):
    """ids in einer Hauptidee mit Radius Rh; zug[i] = (dx, dy, gewicht) Richtung gekoppelter fremder Hauptideen."""
    n = len(ids)
    P = [[(Rh * 0.55) * math.sqrt((i + 0.5) / n) * math.cos(i * GOLD), (Rh * 0.55) * math.sqrt((i + 0.5) / n) * math.sin(i * GOLD)]
         for i in range(n)]
    for t in range(300):
        lr = 0.5 * (1 - t / 300) + 0.02
        F = [[0.0, 0.0] for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                dx, dy = P[j][0] - P[i][0], P[j][1] - P[i][1]
                d = math.hypot(dx, dy) or 1e-9
                soll = r[i] + r[j] + 8
                wij = w_innen.get(tuple(sorted((ids[i], ids[j]))), 0)
                f = (0.05 * min(1.0, wij / 4) * (d - soll)) if wij > 0 else 0.0
                if d < soll:
                    f -= 0.5 * (soll - d)
                F[i][0] += f * dx / d; F[i][1] += f * dy / d
                F[j][0] -= f * dx / d; F[j][1] -= f * dy / d
        for i in range(n):
            zx, zy, zw = zug[i]
            if zw > 0:
                ziel = Rh - r[i] - 10
                F[i][0] += 0.3 * min(1.0, zw / 4) * (zx * ziel - P[i][0]) * 0.05
                F[i][1] += 0.3 * min(1.0, zw / 4) * (zy * ziel - P[i][1]) * 0.05
            P[i][0] += lr * F[i][0]; P[i][1] += lr * F[i][1]
            m = Rh - r[i] - 10
            d = math.hypot(*P[i])
            if d > m > 0:
                P[i][0] *= m / d; P[i][1] *= m / d
    for _ in range(200):                     # Überlappung lösen, ohne die Blase zu verlassen
        _trennen(P, r, 6, 1)
        for i in range(n):
            m = max(0.0, Rh - r[i] - 10)
            d = math.hypot(*P[i])
            if d > m:
                P[i][0] *= m / d if d else 0; P[i][1] *= m / d if d else 0
    return P


def berechne(g, roh):
    quellen = roh.get("quellen") or {}
    ideen = [dict(i) for i in (roh.get("ideen") or [])]
    unter_von = {k: a.get("_idee") for k, a in g.knoten.items()}
    if any(v is None for v in unter_von.values()):          # Sammelidee, damit keine Aussage verloren geht
        ideen.append({"id": OHNE + "_h", "name": "Ohne Zuordnung", "kurz": "Aussagen, die noch keiner Idee zugeordnet sind."})
        ideen.append({"id": OHNE, "eltern": OHNE + "_h", "name": "Ohne Zuordnung", "kurz": "noch keiner Idee zugeordnet"})
        unter_von = {k: (v or OHNE) for k, v in unter_von.items()}
    haupt = [i for i in ideen if not i.get("eltern")]
    unter = [i for i in ideen if i.get("eltern")]
    mitglieder = defaultdict(list)
    for k in sorted(g.knoten):
        mitglieder[unter_von[k]].append(k)
    eltern = {u["id"]: u["eltern"] for u in unter}
    w = kopplung(g, unter_von)
    wh = defaultdict(float)
    for (a, b), x in w.items():
        if eltern.get(a) != eltern.get(b):
            wh[tuple(sorted((eltern[a], eltern[b])))] += x

    # --- Aussagen: Haupt/Hilfs, Rang, Radius
    aussagen = {}
    rang_in = {}
    for u in unter:
        ms = mitglieder[u["id"]]
        mset = set(ms)
        def rang(i):
            a = g.knoten[i]
            tr = len([j for j in a.get("_traegt") or [] if j in mset])
            vs = len([j for j in a.get("_abhaengig") or [] if j in mset])
            return (tr - vs, a.get("_tragweite", 0), i)
        rang_in[u["id"]] = sorted(ms, key=lambda i: (-rang(i)[0], -rang(i)[1], i))
        kand = [i for i in ms if g.knoten[i].get("art") in HAUPT_ARTEN and
                (g.knoten[i].get("_tragweite", 0) >= 1 or (g.knoten[i].get("bewertung") or {}).get("nutzen_stufe") == "hoch")]
        kand.sort(key=lambda i: (-g.knoten[i].get("_tragweite", 0),
                                 (g.knoten[i].get("bewertung") or {}).get("nutzen_stufe") != "hoch", i))
        if not kand and ms:
            kand = sorted(ms, key=lambda i: (-g.knoten[i].get("_tragweite", 0), i))
        wichtig = set(kand[:max(1, math.ceil(len(ms) / 3))]) if ms else set()
        for i in ms:
            r = 5 + 2 * math.sqrt(min(g.knoten[i].get("_tragweite", 0), 25))
            aussagen[i] = {"haupt": i in wichtig, "r": round(r * (1.4 if i in wichtig else 1), 1)}

    # --- Radien
    ru = {u["id"]: 26 + 15 * math.sqrt(len(mitglieder[u["id"]])) for u in unter}
    kinder = defaultdict(list)
    for u in unter:
        kinder[u["eltern"]].append(u["id"])
    rh = {h["id"]: max(120.0, 1.35 * math.sqrt(sum(ru[c] ** 2 for c in kinder[h["id"]]) or 1)) + 20 for h in haupt}

    # --- Hauptideen
    hid = sorted((h["id"] for h in haupt), key=lambda i: (-sum(x for k, x in wh.items() if i in k), i))
    PH = _layout_haupt(hid, [rh[i] for i in hid], dict(wh))
    pos_h = {i: PH[n] for n, i in enumerate(hid)}

    # --- Unterideen
    pos_u = {}
    for h in hid:
        ids = sorted(kinder[h])
        if not ids:
            continue
        w_innen = {k: x for k, x in w.items() if eltern.get(k[0]) == h and eltern.get(k[1]) == h}
        zug = []
        for u in ids:
            zx = zy = zw = 0.0
            for (a, b), x in w.items():
                if u not in (a, b):
                    continue
                anderes = b if a == u else a
                eh = eltern.get(anderes)
                if eh == h:
                    continue
                dx, dy = pos_h[eh][0] - pos_h[h][0], pos_h[eh][1] - pos_h[h][1]
                d = math.hypot(dx, dy) or 1
                zx += x * dx / d; zy += x * dy / d; zw += x
            lz = math.hypot(zx, zy)
            zug.append((zx / lz, zy / lz, zw) if lz > 1e-9 else (0.0, 0.0, 0.0))
        P = _layout_unter(ids, [ru[u] for u in ids], rh[h], w_innen, zug)
        for n, u in enumerate(ids):
            pos_u[u] = [pos_h[h][0] + P[n][0], pos_h[h][1] + P[n][1]]

    # --- Aussagen in der Unteridee: tragende innen (Spirale)
    for u in unter:
        ms = rang_in[u["id"]]
        n = len(ms)
        cx, cy = pos_u[u["id"]]
        for k, i in enumerate(ms):
            rho = (ru[u["id"]] - 12) * math.sqrt((k + 0.5) / n) if n > 1 else 0
            aussagen[i]["x"] = round(cx + rho * math.cos(k * GOLD), 1)
            aussagen[i]["y"] = round(cy + rho * math.sin(k * GOLD), 1)
            aussagen[i]["idee"] = u["id"]

    def stand(ids):
        fuell = _zaehle(_gruppe(g.knoten[i]) for i in ids)
        herk = _zaehle(_herkunft(g.knoten[i], quellen) for i in ids)
        return fuell, herk

    def offen(ids):
        o = [i for i in ids if _gruppe(g.knoten[i]) in ("vermutung", "offen")
             or (g.knoten[i].get("_status_eff") or g.knoten[i].get("status")) == "bedingt"]
        return sorted(o, key=lambda i: (-g.knoten[i].get("_tragweite", 0), i))

    aus_h, aus_u = [], []
    for h in haupt:
        ids = [i for c in kinder[h["id"]] for i in mitglieder[c]]
        f, hk = stand(ids)
        aus_h.append({"id": h["id"], "name": h["name"], "kurz": h.get("kurz", ""), "x": round(pos_h[h["id"]][0], 1),
                      "y": round(pos_h[h["id"]][1], 1), "r": round(rh[h["id"]], 1), "n": len(ids), "fuell": f, "herkunft": hk,
                      "kinder": sorted(kinder[h["id"]])})
    for u in unter:
        ids = mitglieder[u["id"]]
        f, hk = stand(ids)
        aus_u.append({"id": u["id"], "eltern": u["eltern"], "name": u["name"], "kurz": u.get("kurz", ""),
                      "x": round(pos_u[u["id"]][0], 1), "y": round(pos_u[u["id"]][1], 1), "r": round(ru[u["id"]], 1),
                      "n": len(ids), "fuell": f, "herkunft": hk, "offen": offen(ids), "mitglieder": rang_in[u["id"]]})
    return {
        "haupt": aus_h,
        "unter": aus_u,
        "kopplung_haupt": [[a, b, round(x, 2)] for (a, b), x in sorted(wh.items())],
        "kopplung_unter": [[a, b, round(x, 2)] for (a, b), x in sorted(w.items())],
        "aussagen": aussagen,
    }


def _stand_text(f, hk):
    teile = [f"{GRUPPE_NAME[k]} {f[k]}" for k in GRUPPEN if f.get(k)]
    herk = " / ".join(f"{n} {hk.get(k, 0)}" for k, n in (("literatur", "Literatur"), ("forum", "Forum"), ("eigen", "eigen")))
    return ", ".join(teile) + " · " + herk


def text_md(d):
    name = {i["id"]: i["name"] for i in d["haupt"] + d["unter"]}
    z = ["# Landkarte der Tetration", "",
         "Hauptideen und Unterideen mit Füllstand (Status der Aussagen), Herkunft (Literatur / Forum / nur eigene Arbeit) "
         "und den offenen Aussagen (Vermutungen, Fragen, Annahmen, bedingte Sätze; nach Tragweite). Erzeugt von baue.py.", ""]
    for h in sorted(d["haupt"], key=lambda h: -h["n"]):
        z += [f"## {h['name']} ({h['n']} Aussagen)", "", h["kurz"], "", f"Stand: {_stand_text(h['fuell'], h['herkunft'])}", ""]
        for u in sorted((u for u in d["unter"] if u["eltern"] == h["id"]), key=lambda u: -u["n"]):
            o = ", ".join(u["offen"][:8]) + (" …" if len(u["offen"]) > 8 else "")
            z.append(f"- **{u['name']}** ({u['n']}): {_stand_text(u['fuell'], u['herkunft'])}" + (f" · offen: {o}" if o else ""))
        z.append("")
    z += ["## Stärkste Kopplungen zwischen Hauptideen", ""]
    for a, b, x in sorted(d["kopplung_haupt"], key=lambda t: -t[2])[:10]:
        z.append(f"- {name[a]} — {name[b]}: {x}")
    eltern = {u["id"]: u["eltern"] for u in d["unter"]}
    z += ["", "## Stärkste Kopplungen über Ideengrenzen", ""]
    quer = [t for t in d["kopplung_unter"] if eltern.get(t[0]) != eltern.get(t[1])]
    for a, b, x in sorted(quer, key=lambda t: -t[2])[:15]:
        z.append(f"- {name[a]} ({name[eltern[a]]}) — {name[b]} ({name[eltern[b]]}): {x}")
    return "\n".join(z) + "\n"
