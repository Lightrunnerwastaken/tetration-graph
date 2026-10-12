"""Englische Fassung der Strukturtexte: Gebiete, Stränge, Ideen, Lehrbuch-Kapitel, Objekte, Experimente.

Gleiche Idee wie wissen/uebersetzung.py (Aussagen), nur andere Gegenstände. Die Übersetzungen stehen in
`daten/uebersetzung_struktur/<art>.yaml`, ein Eintrag je Gegenstand:

    - art: strang
      id: gr.werte
      quelle_hash: 3fa9c1d20b7e        # Hash der übersetzten deutschen Texte (Pfad -> Text)
      herkunft: {modell: sonnet, datum: '2026-10-13'}
      geprueft: nein                   # nein | zweitpruefung | mensch
      en: {name: Values and magnitudes, verlauf: ...}

Übersetzbare Felder je Art (`pfade`): gebiet name, kurz · strang name, verlauf · idee name, kurz · objekt name ·
experiment titel, frage, ergebnis.zusammenfassung, ergebnis.pruefungen[i].name · quelle notiz (und titel bei eigenen Arbeiten) ·
kapitel titel, untertitel und je Teil `abschnitt`, `text`, `unterschrift` (Liste `teile` in deutscher Länge, nicht Übersetztes bleibt `{}`).

Befehle (aus dem Wurzelverzeichnis):
    python -m wissen.uebersetzung struktur-stand
    python -m wissen.uebersetzung struktur-export --art strang -o stapel.json
    python -m wissen.uebersetzung struktur-einlesen uebersetzt.json --modell sonnet
"""
import datetime
import json
from collections import Counter, defaultdict
from pathlib import Path

from . import schema as S
from .uebersetzung import flach, hole, pruefe_text, quell_hash

EINFACH = {
    "gebiet": ("name", "kurz"),
    "strang": ("name", "verlauf"),
    "idee": ("name", "kurz"),
    "objekt": ("name",),
}
ARTEN = [*EINFACH, "experiment", "quelle", "kapitel"]
# Titel von Quellen werden nur übersetzt, wenn es eigene oder deutschsprachig beschriebene Einträge sind (Titel fremder Werke bleiben)
QUELLE_TITEL_ARTEN = {"eigene", "notiz", "software"}
QUELLE_TITEL_IDS = {"tetrationforum"}
TEIL_FELDER = ("abschnitt", "text", "unterschrift")


def gegenstaende(roh: dict) -> dict:
    """{art: {id: deutsches Wörterbuch}}; die Wörterbücher sind die Objekte aus roh (Übersetzung wird dort angehängt)."""
    return {
        "gebiet": {x["id"]: x for x in roh.get("gebiete") or []},
        "strang": {x["id"]: x for x in roh.get("straenge") or []},
        "idee": {x["id"]: x for x in roh.get("ideen") or []},
        "objekt": {k: v for k, v in (roh.get("objekte") or {}).items() if isinstance(v, dict)},
        "experiment": {x["id"]: x for x in roh.get("experimente") or []},
        "quelle": {k: v for k, v in (roh.get("quellen") or {}).items() if isinstance(v, dict)},
        "kapitel": {k["id"]: k for k in (roh.get("lehrbuch") or {}).get("kapitel") or []},
    }


def deutsch_ansicht(art: str, item: dict) -> dict:
    """Experimente tragen ihr Ergebnis unter `_ergebnis` (aus experimente/ergebnisse/*.json); für die Übersetzung heisst es `ergebnis`."""
    if art == "experiment":
        return {**item, "ergebnis": item.get("_ergebnis") or {}}
    return item


def pfade(art: str, item: dict, ident: str = "") -> list:
    """Übersetzbare Textpfade eines deutschen Gegenstands (nur was es gibt und nicht leer ist)."""
    item = deutsch_ansicht(art, item)

    def text(p):
        v = hole(item, p)
        return isinstance(v, str) and bool(v.strip())
    if art == "experiment":
        p = [f for f in ("titel", "frage", "ergebnis.zusammenfassung") if text(f)]
        return p + [f"ergebnis.pruefungen[{i}].name" for i, x in enumerate(item["ergebnis"].get("pruefungen") or [])
                    if isinstance(x, dict) and isinstance(x.get("name"), str) and x["name"].strip()]
    if art == "quelle":
        p = ["notiz"] if text("notiz") else []
        if text("titel") and (item.get("art") in QUELLE_TITEL_ARTEN or ident in QUELLE_TITEL_IDS):
            p.append("titel")
        return p
    if art in EINFACH:
        return [p for p in EINFACH[art] if text(p)]
    p = [f for f in ("titel", "untertitel") if text(f)]
    for i, t in enumerate(item.get("teile") or []):
        p += [f"teile[{i}].{f}" for f in TEIL_FELDER if isinstance(t, dict) and isinstance(t.get(f), str) and t[f].strip()]
    return p


def vorlage(art: str, item: dict, ident: str = "") -> dict:
    """Die deutsche Struktur, die der Übersetzer ausfüllt (Listen in voller Länge)."""
    if art == "experiment":
        from .uebersetzung import _bauen
        d = deutsch_ansicht(art, item)
        out = _bauen({p: hole(d, p) for p in pfade(art, item)}, None)
        pr = d["ergebnis"].get("pruefungen") or []
        if "ergebnis" in out and "pruefungen" in out["ergebnis"]:
            out["ergebnis"]["pruefungen"] += [{} for _ in range(len(pr) - len(out["ergebnis"]["pruefungen"]))]
        return out
    if art != "kapitel":
        return {p: hole(item, p) for p in pfade(art, item, ident)}
    out = {f: item[f] for f in ("titel", "untertitel") if f in pfade(art, item)}
    out["teile"] = [{f: t[f] for f in TEIL_FELDER if isinstance(t, dict) and f in t and isinstance(t[f], str) and t[f].strip()}
                    for t in item.get("teile") or []]
    return out


def pruefe_form(art: str, item: dict, en: dict, ident: str = "") -> list:
    if not isinstance(en, dict):
        return [("FEHLER", "uebersetzung-struktur", "'en' muss ein Wörterbuch sein")]
    b = []
    erlaubt = set(pfade(art, item, ident))
    if art == "kapitel" and "teile" in en:
        if not isinstance(en["teile"], list) or len(en["teile"]) != len(item.get("teile") or []):
            return [("FEHLER", "uebersetzung-struktur", f"'teile' muss {len(item.get('teile') or [])} Einträge haben")]
        for i, t in enumerate(en["teile"]):
            if not isinstance(t, dict):
                b.append(("FEHLER", "uebersetzung-struktur", f"teile[{i}] muss ein Wörterbuch sein"))
    for p, v in flach(en).items():
        if p not in erlaubt:
            b.append(("FEHLER", "uebersetzung-struktur", f"'{p}' ist nicht übersetzbar oder im Deutschen leer"))
        elif not v.strip():
            b.append(("FEHLER", "uebersetzung-struktur", f"'{p}' ist leer"))
    return b


def pruefe_eintrag(art: str, item: dict, e: dict) -> list:
    ident = e.get("id", "")
    de = deutsch_ansicht(art, item)
    b = []
    gp = e.get("geprueft", "nein")
    if gp not in S.UEBERSETZ_GEPRUEFT:
        b.append(("FEHLER", "uebersetzung-struktur", f"geprueft='{gp}' unbekannt (erlaubt: {sorted(S.UEBERSETZ_GEPRUEFT)})"))
    en = e.get("en")
    if en is None:
        return b + [("FEHLER", "uebersetzung-struktur", "Eintrag ohne 'en'")]
    form = pruefe_form(art, item, en, ident)
    b += form
    if not form:
        for p, txt in flach(en).items():
            b += [(s, r, f"{p}: {t}") for s, r, t in pruefe_text(hole(de, p), txt)]
        aktuell = quell_hash({p: hole(de, p) for p in flach(en)})
        if not e.get("quelle_hash"):
            b.append(("WARNUNG", "uebersetzung-veraltet", "ohne Quell-Hash (mit `struktur-einlesen` stempeln)"))
        elif e["quelle_hash"] != aktuell:
            b.append(("WARNUNG", "uebersetzung-veraltet", "deutscher Text hat sich seit der Übersetzung geändert"))
    return b


def einbinden(g):
    """Hängt die Übersetzungen als item['en'] an die Gegenstände (roh-Objekte) und prüft sie."""
    alle = gegenstaende(g.roh)
    gesehen = set()
    for e in g.roh.get("uebersetzung_struktur") or []:
        art, i = e.get("art"), e.get("id")
        if art not in alle:
            g.befund("FEHLER", "uebersetzung-struktur", f"{art}:{i}", f"unbekannte Art (erlaubt: {', '.join(ARTEN)})")
            continue
        if i not in alle[art]:
            g.befund("FEHLER", "uebersetzung-struktur", f"{art}:{i}", "Übersetzung für unbekannten Gegenstand")
            continue
        if (art, i) in gesehen:
            g.befund("FEHLER", "uebersetzung-struktur", f"{art}:{i}", "mehrere Einträge für denselben Gegenstand")
            continue
        gesehen.add((art, i))
        item = alle[art][i]
        b = pruefe_eintrag(art, item, e)
        for s, r, t in b:
            g.befund(s, r if r.startswith("uebersetzung") else "uebersetzung-struktur", f"{art}:{i}", t)
        if not any(s == "FEHLER" for s, _, _ in b):
            item["en"] = e["en"]
            item["_en_stand"] = "veraltet" if any(r == "uebersetzung-veraltet" for _, r, _ in b) else "aktuell"
    if gesehen:
        for art in ARTEN:
            n = len(alle[art])
            mit = sum(1 for (a, _) in gesehen if a == art)
            if mit < n:
                g.befund("HINWEIS", "uebersetzung-fortschritt", f"{art}", f"{mit} von {n} {art}-Gegenständen übersetzt")
    return g


# ------------------------------------------------------------------ Befehle
def lade_roh(wurzel: Path) -> dict:
    """roh nach vollem Lauf des Prüfprogramms (Experimente tragen dann `_ergebnis`)."""
    from .graph import Graph
    from .laden import lade
    roh = lade(wurzel)
    Graph(roh, wurzel).alles()
    return roh


def befehl_stand(wurzel: Path):
    roh = lade_roh(wurzel)
    alle = gegenstaende(roh)
    vorh = {(e["art"], e["id"]): e for e in roh.get("uebersetzung_struktur") or []}
    print(f"{'Art':11} {'gesamt':>6} {'übersetzt':>9} {'veraltet':>8} {'gegengelesen':>12}")
    for art in ARTEN:
        n = len(alle[art])
        mit = [(art, i) for i in alle[art] if (art, i) in vorh]
        alt = sum(1 for (a, i) in mit
                  if any(r == "uebersetzung-veraltet" for _, r, _ in pruefe_eintrag(art, alle[art][i], vorh[(a, i)])))
        gel = sum(1 for k in mit if vorh[k].get("geprueft", "nein") != "nein")
        print(f"{art:11} {n:>6} {len(mit):>9} {alt:>8} {gel:>12}")


def befehl_export(wurzel: Path, art: str, maximum: int, nur_fehlend: bool, ausgabe):
    roh = lade_roh(wurzel)
    alle = gegenstaende(roh)[art]
    vorh = {(e["art"], e["id"]) for e in roh.get("uebersetzung_struktur") or []}
    stapel = []
    for i, item in alle.items():
        if not pfade(art, item, i) or (nur_fehlend and (art, i) in vorh):
            continue
        stapel.append({"art": art, "id": i, "de": vorlage(art, item, i)})
    if maximum:
        stapel = stapel[:maximum]
    text = json.dumps(stapel, ensure_ascii=False, indent=1)
    if ausgabe:
        Path(ausgabe).write_text(text, encoding="utf-8")
        print(f"{len(stapel)} Einträge → {ausgabe}")
    else:
        print(text)


def befehl_einlesen(wurzel: Path, datei: str, modell: str, geprueft: str, trocken: bool) -> int:
    import yaml

    from .uebersetzung import _yaml_schreiben
    roh = lade_roh(wurzel)
    alle = gegenstaende(roh)
    eingabe = json.loads(Path(datei).read_text(encoding="utf-8"))
    fehler = 0
    neu = defaultdict(list)
    for z in eingabe:
        art, i = z.get("art"), z.get("id")
        if art not in alle or i not in alle[art]:
            print(f"FEHLER {art}:{i}: unbekannter Gegenstand")
            fehler += 1
            continue
        item = alle[art][i]
        e = {"art": art, "id": i, "quelle_hash": quell_hash({p: hole(deutsch_ansicht(art, item), p) for p in flach(z["en"])}),
             "herkunft": {"modell": modell, "datum": datetime.date.today().isoformat()}, "geprueft": geprueft, "en": z["en"]}
        b = pruefe_eintrag(art, item, e)
        for s, r, t in b:
            if s != "HINWEIS":
                print(f"{s:8} {art}:{i}: {r}: {t}")
        if any(s == "FEHLER" for s, _, _ in b):
            fehler += 1
            continue
        neu[art].append(e)
    if fehler:
        print(f"{fehler} Einträge mit FEHLER nicht übernommen.")
    if trocken:
        print(f"Trockenlauf: {sum(len(v) for v in neu.values())} Einträge wären geschrieben worden.")
        return 1 if fehler else 0
    ordner = wurzel / "daten" / "uebersetzung_struktur"
    for art, eintraege in neu.items():
        ziel = ordner / f"{art}.yaml"
        bestand = {e["id"]: e for e in (yaml.safe_load(ziel.read_text(encoding="utf-8")) or [])} if ziel.exists() else {}
        for e in eintraege:
            bestand[e["id"]] = e
        _yaml_schreiben(ziel, list(bestand.values()))
        print(f"{ziel.relative_to(wurzel)}: {len(eintraege)} Einträge geschrieben")
    return 1 if fehler else 0
