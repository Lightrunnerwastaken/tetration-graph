"""Englische Fassung: Aufbau der Übersetzungsdateien, Prüfregeln, Stand, Export und Einlesen.

Die Übersetzungen stehen NICHT in den deutschen Aussagedateien, sondern in `daten/uebersetzung/<gebiet>.yaml` (maschinell geführt):

    - id: gr.satz.kettenregel
      quelle_hash: 3fa9c1d20b7e        # Hash der deutschen Texte, die hier übersetzt sind (Pfad -> Text)
      theorie_hash: 91ce44aa07d2       # nur wenn en_theorie vorhanden
      herkunft: {modell: sonnet, datum: '2026-10-13'}
      geprueft: nein                   # nein | zweitpruefung | mensch – oder je Feldgruppe {kern: …, beweis: …, intern: …, theorie: …}
      en:
        titel: Chain rule of tetration
        aussage: ...
        beweis: {text: ...}            # Unterstruktur wie im Deutschen, nur die Textfelder
        formen: [{notiz: ...}, {}]     # gleiche Länge wie deutsch; `tex` wird nie übersetzt
      en_theorie: {worum: ..., genau: ..., warum: ..., beispiel: ..., einordnung: ...}

Ändert sich ein übersetzter deutscher Text, passt der Hash nicht mehr und die Prüfung meldet „uebersetzung-veraltet“.

Befehle (aus dem Wurzelverzeichnis):
    python -m wissen.uebersetzung stand
    python -m wissen.uebersetzung export --gebiet fp --felder kern --max 20 -o stapel.json
    python -m wissen.uebersetzung export --gebiet fp --felder theorie -o stapel.json
    python -m wissen.uebersetzung einlesen uebersetzt.json --modell sonnet [--ersetzen]
    python -m wissen.uebersetzung abgleich sonnet.json haiku.json      # wo zwei unabhängige Übersetzungen auseinanderliegen
    python -m wissen.uebersetzung struktur-stand | struktur-export --art strang | struktur-einlesen datei.json --modell sonnet
        (Gebiete, Stränge, Ideen, Lehrbuch-Kapitel, Objekte, Experimente: siehe wissen/uebersetzung_struktur.py)
"""
import argparse
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from . import schema as S

# ------------------------------------------------------------------ Muster für die maschinellen Prüfungen
MATHE = re.compile(r"\$\$.+?\$\$|\$.+?\$", re.S)
VERWEIS = re.compile(r"\[\[.+?\]\]|!\[\[.+?\]\]")
KENNUNG = re.compile(r"\b(?:gr|fg|it|fp|kl|kn|bw|ar|br|exp)\.[a-z0-9_]+(?:\.[A-Za-z0-9_]+)*")
ZAHL = re.compile(r"(?<![\w.])\d+(?:\.\d+)?")
# (deutsches Signalwort, was im Englischen vorkommen muss)
SIGNAL = [
    (r"\bnicht\b|\bkein", r"\bnot\b|\bno\b|\bnever\b|\bneither\b|n't|\bnon-?|\bwithout\b|\bnor\b|cannot|\bimpossible\b|\bnone\b|\bnowhere\b|\bun[a-z]+"),
    (r"\bnur\b", r"\bonly\b|\bmerely\b|\bsolely\b|\balone\b"),
    (r"\bgenau\b", r"\bexactly\b|\bprecisely\b|\bjust\b|\bif and only if\b|\bexact\b|\bsharp\b|\bunique|\baccura|\bto within\b"),
    (r"\bjede[rsmn]?\b|\balle\b", r"\bevery\b|\bany\b|\ball\b|\beach\b"),
    (r"\bvermutlich\b|\bwohl\b(?! aber)|\bvielleicht\b", r"\bpresumably\b|\bprobably\b|\blikely\b|\bperhaps\b|\bsuppos|\bmay\b|\bmight\b|\bapparently\b|\bexpect"),
    (r"\bbedingt", r"\bconditional"),
    (r"\bVermutung", r"\bconjecture"),
    (r"\bwiderlegt", r"\brefut|\bdisprov|\bfalsif"),
]
BRITISCH = re.compile(r"neighbour|behaviour|colour|analys(?:e|ing)\b|normalis|favour|centre\b|\bgrey\b|organis|recognis|honour|"
                      r"generalis|characteris|minimis|maximis|parametris|neighbouring", re.I)
DEUTSCHE_WOERTER = re.compile(r" (?:und|der|die|das|ist|wird|mit|für|nicht|auch|oder|wenn|dann|aber) ")
UMLAUT = re.compile(r"\w*[äöüÄÖÜß]\w*")


# ------------------------------------------------------------------ Pfade in verschachtelten Feldern
def flach(o, pfad=""):
    """Alle Textblätter als {pfad: text}. Listen: pfad[i]."""
    out = {}
    if isinstance(o, str):
        out[pfad] = o
    elif isinstance(o, dict):
        for k, v in o.items():
            out.update(flach(v, f"{pfad}.{k}" if pfad else str(k)))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            out.update(flach(v, f"{pfad}[{i}]"))
    return out


_TEIL = re.compile(r"([^.\[\]]+)|\[(\d+)\]")


def hole(o, pfad):
    """Wert zu 'beweis.text' / 'formen[0].notiz'; None wenn es ihn nicht gibt."""
    for name, idx in _TEIL.findall(pfad):
        if idx != "":
            if not isinstance(o, list) or int(idx) >= len(o):
                return None
            o = o[int(idx)]
        else:
            if not isinstance(o, dict) or name not in o:
                return None
            o = o[name]
    return o


FELDGRUPPEN = {
    "kern": ["titel", "aussage", "bereich", "einfach", "formen"],
    "beweis": ["beweis", "gegenbeispielsuche", "naechste_frage", "bedingungen_text"],
    "intern": ["pruefung", "bewertung", "voraussetzungen"],
    "zusatz": ["versionen", "nachgeprueft", "quellen"],
}
GRUPPE_VON_FELD = {f: g for g, fs in FELDGRUPPEN.items() for f in fs}
RANG = {"nein": 0, "zweitpruefung": 1, "mensch": 2}


def gruppen_im_eintrag(e: dict) -> set:
    """Feldgruppen, die ein Übersetzungseintrag enthält."""
    gr = {GRUPPE_VON_FELD[f] for f in (e.get("en") or {}) if f in GRUPPE_VON_FELD}
    return gr | ({"theorie"} if e.get("en_theorie") else set())


def geprueft_je_gruppe(e: dict) -> dict:
    """Prüfstand je enthaltener Feldgruppe; ein einzelner Wert gilt für alle Gruppen des Eintrags."""
    g = e.get("geprueft", "nein")
    gr = gruppen_im_eintrag(e)
    if isinstance(g, dict):
        return {k: g.get(k, "nein") for k in sorted(gr)}
    return {k: g for k in sorted(gr)}


def schlechtester_stand(e: dict) -> str:
    werte = list(geprueft_je_gruppe(e).values()) or [e.get("geprueft", "nein") if isinstance(e.get("geprueft", "nein"), str) else "nein"]
    return min(werte, key=lambda w: RANG.get(w, 0))


def quell_hash(texte: dict) -> str:
    """Hash der deutschen Texte (Pfad -> Text), unabhängig von der Reihenfolge."""
    return hashlib.sha256(json.dumps(sorted(texte.items()), ensure_ascii=False).encode("utf-8")).hexdigest()[:12]


def deutsch(a: dict, pfade) -> dict:
    return {p: hole(a, p) for p in pfade}


def deutsch_theorie(a: dict, teile) -> dict:
    return {t: (a.get("theorie") or {}).get(t) for t in teile}


# ------------------------------------------------------------------ Form des en-Blocks
def pruefe_form(a: dict, en: dict) -> list:
    """Befunde (stufe, regel, text) zur Gestalt eines en-Blocks gegen die deutsche Aussage."""
    b = []
    if not isinstance(en, dict):
        return [("FEHLER", "uebersetzung", "'en' muss ein Feld-Wörterbuch sein")]
    for f, v in en.items():
        if f not in S.UEBERSETZBAR:
            b.append(("FEHLER", "uebersetzung", f"Feld '{f}' ist nicht übersetzbar"))
        elif f in S.UEBERSETZ_TEXT:
            if not isinstance(v, str) or not v.strip():
                b.append(("FEHLER", "uebersetzung", f"'{f}' muss ein nichtleerer Text sein"))
        elif f in S.UEBERSETZ_LISTEN:
            de = a.get(f) or []
            if not isinstance(v, list):
                b.append(("FEHLER", "uebersetzung", f"'{f}' muss eine Liste sein"))
            elif len(v) != len(de):
                b.append(("FEHLER", "uebersetzung", f"'{f}': {len(v)} Einträge, deutsch {len(de)}"))
            else:
                erlaubt = S.UEBERSETZ_BAU[f]
                for i, e in enumerate(v):
                    if not isinstance(e, dict):
                        b.append(("FEHLER", "uebersetzung", f"'{f}[{i}]' muss ein Wörterbuch sein"))
                        continue
                    for k in e:
                        if k not in erlaubt:
                            b.append(("FEHLER", "uebersetzung", f"'{f}[{i}].{k}' ist nicht übersetzbar (erlaubt: {sorted(erlaubt)})"))
        else:
            if not isinstance(v, dict):
                b.append(("FEHLER", "uebersetzung", f"'{f}' muss ein Wörterbuch sein"))
                continue
            for k in v:
                if k not in S.UEBERSETZ_BAU[f]:
                    b.append(("FEHLER", "uebersetzung", f"'{f}.{k}' ist nicht übersetzbar (erlaubt: {sorted(S.UEBERSETZ_BAU[f])})"))
            if f == "beweis" and "luecken" in v:
                de = (a.get("beweis") or {}).get("luecken") or []
                if not isinstance(v["luecken"], list) or len(v["luecken"]) != len(de):
                    b.append(("FEHLER", "uebersetzung", f"'beweis.luecken' muss eine Liste mit {len(de)} Einträgen sein"))
    for p, v in flach(en).items():
        if hole(a, p) is None:
            b.append(("FEHLER", "uebersetzung", f"'{p}' hat im Deutschen keine Entsprechung"))
        elif not v.strip():
            b.append(("FEHLER", "uebersetzung", f"'{p}' ist leer"))
    return b


def pruefe_form_theorie(a: dict, en: dict) -> list:
    b = []
    if not isinstance(en, dict):
        return [("FEHLER", "uebersetzung", "'en_theorie' muss ein Wörterbuch sein")]
    for t, v in en.items():
        if t not in S.UEBERSETZ_THEORIE:
            b.append(("FEHLER", "uebersetzung", f"Theorie-Teil '{t}' ist nicht übersetzbar"))
        elif not isinstance(v, str) or not v.strip():
            b.append(("FEHLER", "uebersetzung", f"Theorie-Teil '{t}' muss ein nichtleerer Text sein"))
        elif not isinstance((a.get("theorie") or {}).get(t), str):
            b.append(("FEHLER", "uebersetzung", f"Theorie-Teil '{t}' hat im Deutschen keine Entsprechung"))
    return b


# ------------------------------------------------------------------ Invarianten Text gegen Text
def _mathe(s: str) -> Counter:
    """Formeln als Mehrmenge; Prosa in \\text{…} darf übersetzt werden und wird eingeebnet."""
    return Counter(re.sub(r"\\text\{[^}]*\}", r"\\text{}", m) for m in MATHE.findall(s))


def _verweis_ids(s: str) -> list:
    """[[id|Linktext]] -> [[id]]: nur das Ziel muss gleich bleiben, der angezeigte Text wird übersetzt."""
    return [re.sub(r"\|[^\]]*\]\]$", "]]", v) for v in VERWEIS.findall(s)]


def _dezimalpunkt(s: str) -> str:
    """Dezimalkomma -> Punkt, ausser in Intervallen und Mengen wie (0,1), [0,1] oder {2,10}."""
    def ersetze(m):
        return m.group(0) if m.group(1) and m.group(4) else f"{m.group(1) or ''}{m.group(2)}.{m.group(3)}{m.group(4) or ''}"
    return re.sub(r"([\[({])?[−-]?\b(\d+),(\d+)\b(?!\.\d)([)\]}])?", ersetze, s)


def pruefe_text(de: str, en: str) -> list:
    """Maschinelle Prüfung eines Textpaares. FEHLER = Inhalt kann verfälscht sein, WARNUNG = bitte ansehen."""
    b = []
    qm, um = _mathe(de), _mathe(en)
    if qm != um:
        b.append(("FEHLER", "uebersetzung-invariante",
                  f"Formeln verändert: fehlt {list((qm - um).elements())[:2]}, neu {list((um - qm).elements())[:2]}"))
    if Counter(_verweis_ids(de)) != Counter(_verweis_ids(en)):          # Linktext hinter | darf übersetzt werden
        b.append(("FEHLER", "uebersetzung-invariante", "Verweise [[…]] verändert"))
    if Counter(KENNUNG.findall(de)) != Counter(KENNUNG.findall(en)):
        b.append(("FEHLER", "uebersetzung-invariante", "Kennungen (z. B. exp.…, kn.…) verändert"))
    # deutsches Dezimalkomma (1,5) darf im Englischen zum Punkt werden und soll es auch
    zq = Counter(ZAHL.findall(_dezimalpunkt(MATHE.sub(" ", de))))
    zu = Counter(ZAHL.findall(MATHE.sub(" ", en)))
    if zq != zu:
        b.append(("WARNUNG", "uebersetzung-zahlen", f"Zahlen ausserhalb von Formeln verändert: fehlt {dict(zq - zu)}, neu {dict(zu - zq)}"))
    if de.count("\n\n") != en.count("\n\n"):
        b.append(("WARNUNG", "uebersetzung-absaetze", f"Absätze: {de.count(chr(10) * 2)} deutsch, {en.count(chr(10) * 2)} englisch"))
    de_s, en_s = KENNUNG.sub(" ", de), KENNUNG.sub(" ", en)       # Kennungen wie ar.satz.keine_arith_info enthalten Signalwörter
    for muster_de, muster_en in SIGNAL:
        if re.search(muster_de, de_s, re.I) and not re.search(muster_en, en_s, re.I):
            b.append(("WARNUNG", "uebersetzung-wortwahl", f"Signalwort '{muster_de}' im Deutschen, im Englischen nichts Entsprechendes"))
    if BRITISCH.search(en):
        b.append(("WARNUNG", "uebersetzung-wortwahl", f"britische Schreibweise: {BRITISCH.search(en).group(0)}"))
    ohne = MATHE.sub(" ", en).replace("van der", "van_der")
    uml = [w for w in UMLAUT.findall(ohne) if w not in de]
    if uml or DEUTSCHE_WOERTER.search(ohne):
        b.append(("WARNUNG", "uebersetzung-wortwahl", "Deutsch im englischen Text" + (f": {uml[:3]}" if uml else "")))
    return b


def pruefe_eintrag(a: dict, eintrag: dict) -> list:
    """Alle Befunde (stufe, regel, text) zu einem Übersetzungseintrag; a = deutsche Aussage (mit Theorie-Block)."""
    b = []
    en, th = eintrag.get("en"), eintrag.get("en_theorie")
    if en is None and th is None:
        return [("FEHLER", "uebersetzung", "Eintrag ohne 'en' und ohne 'en_theorie'")]
    gp = eintrag.get("geprueft", "nein")
    erlaubt_gruppen = set(FELDGRUPPEN) | {"theorie"}
    if isinstance(gp, dict):
        schlecht = [f"{k}={v}" for k, v in gp.items() if v not in S.UEBERSETZ_GEPRUEFT or k not in erlaubt_gruppen]
    else:
        schlecht = [] if gp in S.UEBERSETZ_GEPRUEFT else [str(gp)]
    if schlecht:
        b.append(("FEHLER", "uebersetzung", f"geprueft unbekannt: {', '.join(schlecht)} (Werte: {sorted(S.UEBERSETZ_GEPRUEFT)}, Gruppen: {sorted(erlaubt_gruppen)})"))
    if en is not None:
        form = pruefe_form(a, en)
        b += form
        if not form:
            for p, txt in flach(en).items():
                b += [(s, r, f"{p}: {t}") for s, r, t in pruefe_text(hole(a, p), txt)]
            b += _stand(eintrag.get("quelle_hash"), quell_hash(deutsch(a, flach(en))), "Aussage")
    if th is not None:
        form = pruefe_form_theorie(a, th)
        b += form
        if not form:
            for t, txt in th.items():
                b += [(s, r, f"theorie.{t}: {x}") for s, r, x in pruefe_text(a["theorie"][t], txt)]
            b += _stand(eintrag.get("theorie_hash"), quell_hash(deutsch_theorie(a, th)), "Theorie")
    return b


def _stand(gespeichert, aktuell, was):
    if not gespeichert:
        return [("WARNUNG", "uebersetzung-veraltet", f"{was}: ohne Quell-Hash (mit `einlesen` stempeln)")]
    if gespeichert != aktuell:
        return [("WARNUNG", "uebersetzung-veraltet", f"{was}: deutscher Text hat sich seit der Übersetzung geändert")]
    return []


def stand_eines(a: dict, eintrag: dict) -> str:
    """'aktuell' | 'veraltet' | 'ungestempelt' | 'fehlerhaft' für die Ansicht (eine Zeile pro Aussage)."""
    b = pruefe_eintrag(a, eintrag)
    if any(s == "FEHLER" for s, _, _ in b):
        return "fehlerhaft"
    if any(r == "uebersetzung-veraltet" and "ohne Quell-Hash" in t for _, r, t in b):
        return "ungestempelt"
    if any(r == "uebersetzung-veraltet" for _, r, _ in b):
        return "veraltet"
    return "aktuell"


# ------------------------------------------------------------------ Einbindung in den Graphen
def einbinden(g):
    """Hängt die Übersetzungseinträge an die Aussagen (a['en'], a['en_theorie'], a['_en_stand']) und prüft sie."""
    gesehen = {}
    for e in g.roh.get("uebersetzung") or []:
        i = e.get("id")
        if i not in g.knoten:
            g.befund("FEHLER", "uebersetzung", str(i), "Übersetzung für unbekannte Aussage")
            continue
        if i in gesehen:
            g.befund("FEHLER", "uebersetzung", i, "mehrere Übersetzungseinträge für dieselbe Aussage")
            continue
        a = g.knoten[i]
        if a.get("en") is not None and e.get("en") is not None:
            g.befund("FEHLER", "uebersetzung", i, "'en' steht sowohl in der Aussagedatei als auch in daten/uebersetzung/")
            continue
        gesehen[i] = e
        for s, r, t in pruefe_eintrag(a, e):
            g.befund(s, r, i, t)
        if e.get("en") is not None:
            a["en"] = e["en"]
        if e.get("en_theorie") is not None:
            a["en_theorie"] = e["en_theorie"]
        a["_en_stand"] = stand_eines(a, e)
        a["_en_geprueft"] = schlechtester_stand(e)
    # Fortschritt je Gebiet (nur wenn schon übersetzt wird)
    if gesehen:
        for geb in g.roh.get("gebiete") or []:
            ids = [i for i, a in g.knoten.items() if a.get("gebiet") == geb["id"]]
            if not ids:
                continue
            mit = [i for i in ids if i in gesehen and g.knoten[i].get("en")]
            alt = sum(1 for i in mit if g.knoten[i]["_en_stand"] != "aktuell")
            th = sum(1 for i in ids if i in gesehen and gesehen[i].get("en_theorie"))
            geprueft = sum(1 for i in mit if schlechtester_stand(gesehen[i]) != "nein")
            if len(mit) < len(ids) or alt or th < len(ids):
                g.befund("HINWEIS", "uebersetzung-fortschritt", geb["id"],
                         f"{len(mit)} von {len(ids)} Aussagen übersetzt ({alt} nicht aktuell, {geprueft} gegengelesen), "
                         f"{th} von {len(ids)} Theorie-Blöcken")
    return g


# ------------------------------------------------------------------ Werkzeuge (Befehlszeile)
def _blatt_pfade(a, feld):
    """Übersetzbare Pfade einer Aussage für ein Feld (nur was es im Deutschen gibt)."""
    v = a.get(feld)
    if v is None:
        return []
    if feld in S.UEBERSETZ_TEXT:
        return [feld] if isinstance(v, str) and v.strip() else []
    if feld in S.UEBERSETZ_LISTEN:
        return [f"{feld}[{i}].{k}" for i, e in enumerate(v) if isinstance(e, dict)
                for k in sorted(S.UEBERSETZ_BAU[feld]) if isinstance(e.get(k), str) and e[k].strip()]
    pf = [f"{feld}.{k}" for k in sorted(S.UEBERSETZ_BAU[feld]) if isinstance(v, dict) and isinstance(v.get(k), str) and v[k].strip()]
    if feld == "beweis" and isinstance(v, dict):
        pf += [f"beweis.luecken[{i}]" for i, x in enumerate(v.get("luecken") or []) if isinstance(x, str) and x.strip()]
    return pf


def _bauen(pfade_texte, vorlage_aussage):
    """Baut aus {pfad: text} die verschachtelte en-Struktur (Listen in deutscher Länge)."""
    out = {}
    for p, t in pfade_texte.items():
        teile = [(n, i) for n, i in _TEIL.findall(p)]
        cur = out
        for j, (name, idx) in enumerate(teile):
            letzte = j == len(teile) - 1
            if idx != "":
                lst = cur
                while len(lst) <= int(idx):
                    lst.append({})
                if letzte:
                    lst[int(idx)] = t
                else:
                    cur = lst[int(idx)]
            else:
                nxt = teile[j + 1] if not letzte else None
                if letzte:
                    cur[name] = t
                else:
                    cur = cur.setdefault(name, [] if nxt[1] != "" else {})
    return out


def _auffuellen(out: dict, a: dict) -> dict:
    """Listen im Export auf die deutsche Länge bringen (fehlende Einträge `{}`), damit der Übersetzer die Ausrichtung behält."""
    for f in S.UEBERSETZ_LISTEN:
        if f in out and isinstance(a.get(f), list):
            out[f] += [{} for _ in range(len(a[f]) - len(out[f]))]
    luecken = (a.get("beweis") or {}).get("luecken") if isinstance(a.get("beweis"), dict) else None
    if "beweis" in out and "luecken" in out["beweis"] and luecken:
        out["beweis"]["luecken"] += [""] * (len(luecken) - len(out["beweis"]["luecken"]))
    return out


def lade_graph(wurzel: Path):
    from .graph import Graph
    from .laden import lade
    roh = lade(wurzel)
    return Graph(roh, wurzel).alles()


def befehl_stand(wurzel):
    """Fortschritt je Gebiet; je Feldgruppe 'übersetzt/gegengelesen' (gegengelesen = zweitpruefung oder mensch)."""
    g = lade_graph(wurzel)
    gruppen = [*FELDGRUPPEN, "theorie"]
    zeilen = defaultdict(Counter)
    eintraege = {e["id"]: e for e in g.roh.get("uebersetzung") or []}
    for i, a in g.knoten.items():
        z = zeilen[a.get("gebiet")]
        z["gesamt"] += 1
        e = eintraege.get(i)
        if not e:
            continue
        if a.get("en"):
            z["veraltet" if a["_en_stand"] != "aktuell" else "aktuell"] += 1
        for gr, wert in geprueft_je_gruppe(e).items():
            z[gr + ":uebersetzt"] += 1
            z[gr + ":geprueft"] += RANG.get(wert, 0) > 0
    kopf = f"{'Gebiet':8} {'gesamt':>6} {'veraltet':>8}  " + "  ".join(f"{gr:>9}" for gr in gruppen) + "     (je Gruppe: übersetzt/gegengelesen)"
    print(kopf)
    summe = Counter()
    for geb, z in [*sorted(zeilen.items()), ("Summe", None)]:
        if z is None:
            z = summe
        else:
            summe.update(z)
        print(f"{geb:8} {z['gesamt']:>6} {z['veraltet']:>8}  " +
              "  ".join(f"{str(z[gr + ':uebersetzt']) + '/' + str(z[gr + ':geprueft']):>9}" for gr in gruppen))


def befehl_export(wurzel, gebiet, felder, maximum, nur_fehlend, ausgabe):
    g = lade_graph(wurzel)
    stapel = []
    for i, a in g.knoten.items():
        if gebiet and a.get("gebiet") != gebiet:
            continue
        if felder == "theorie":
            teile = {t: (a.get("theorie") or {}).get(t) for t in sorted(S.UEBERSETZ_THEORIE)}
            teile = {t: v for t, v in teile.items() if isinstance(v, str) and v.strip()}
            if not teile or (nur_fehlend and a.get("en_theorie")):
                continue
            stapel.append({"id": i, "art": a.get("art"), "status": a.get("status"), "de_theorie": teile})
            continue
        pfade = []
        for f in FELDGRUPPEN[felder]:
            pfade += _blatt_pfade(a, f)
        if not pfade:
            continue
        if nur_fehlend and a.get("en") and not (set(pfade) - set(flach(a["en"]))):
            continue
        texte = {p: hole(a, p) for p in pfade}
        if nur_fehlend and a.get("en"):
            texte = {p: t for p, t in texte.items() if p not in flach(a["en"])}
        stapel.append({"id": i, "art": a.get("art"), "status": a.get("status"), "de": _auffuellen(_bauen(texte, a), a)})
    if maximum:
        stapel = stapel[:maximum]
    text = json.dumps(stapel, ensure_ascii=False, indent=1)
    if ausgabe:
        Path(ausgabe).write_text(text, encoding="utf-8")
        print(f"{len(stapel)} Einträge → {ausgabe}")
    else:
        print(text)


def _yaml_schreiben(pfad: Path, eintraege: list):
    import yaml

    class Dumper(yaml.SafeDumper):
        pass

    def text(dumper, daten):
        return dumper.represent_scalar("tag:yaml.org,2002:str", daten, style="|" if "\n" in daten else None)

    Dumper.add_representer(str, text)
    kopf = ("# Maschinell geführt (python -m wissen.uebersetzung einlesen …). Englische Fassung je Aussage; quelle_hash = Hash der\n"
            "# übersetzten deutschen Texte. Nicht von Hand ändern, ausser Fehler in einzelnen Texten – dann `geprueft` anpassen.\n")
    pfad.parent.mkdir(parents=True, exist_ok=True)
    pfad.write_text(kopf + yaml.dump(sorted(eintraege, key=lambda e: e["id"]), Dumper=Dumper, allow_unicode=True,
                                     sort_keys=False, width=140, default_flow_style=False), encoding="utf-8",
                    newline="\n")                       # auch unter Windows LF, sonst wird jede Datei ganz neu eingecheckt


def befehl_einlesen(wurzel, datei, modell, ersetzen, trocken, geprueft):
    import datetime

    import yaml
    g = lade_graph(wurzel)
    eingabe = json.loads(Path(datei).read_text(encoding="utf-8"))
    vorhanden = {e["id"]: e for e in g.roh.get("uebersetzung") or []}
    herkunft_datei = {}
    daten = wurzel / "daten" / "uebersetzung"
    for pfad in sorted(daten.glob("*.yaml")) if daten.is_dir() else []:
        for e in yaml.safe_load(pfad.read_text(encoding="utf-8")) or []:
            herkunft_datei[e["id"]] = pfad
    neu = defaultdict(list)       # Datei -> Einträge
    fehler = 0
    for z in eingabe:
        i = z["id"]
        if i not in g.knoten:
            print(f"FEHLER {i}: unbekannte Aussage")
            fehler += 1
            continue
        a = g.knoten[i]
        alt = vorhanden.get(i)
        e = {"id": i}
        if alt and not ersetzen:
            for k, v in alt.items():
                e[k] = json.loads(json.dumps(v))
        if z.get("en") is not None:
            if "en" in e and not ersetzen:
                if any(r == "uebersetzung-veraltet" for _, r, _ in pruefe_eintrag(a, alt)):
                    print(f"FEHLER {i}: vorhandene Übersetzung ist veraltet – mit --ersetzen neu schreiben")
                    fehler += 1
                    continue
                e["en"] = _verschmelzen(e["en"], z["en"])
            else:
                e["en"] = z["en"]
        if z.get("en_theorie") is not None:
            e["en_theorie"] = {**(e.get("en_theorie") or {}), **z["en_theorie"]} if not ersetzen else z["en_theorie"]
        if "en" in e:
            e["quelle_hash"] = quell_hash(deutsch(a, flach(e["en"])))
        if "en_theorie" in e:
            e["theorie_hash"] = quell_hash(deutsch_theorie(a, e["en_theorie"]))
        e["herkunft"] = {"modell": modell, "datum": datetime.date.today().isoformat()}
        neue_gruppen = gruppen_im_eintrag({"en": z.get("en"), "en_theorie": z.get("en_theorie")})
        alte = geprueft_je_gruppe(e) if (alt and not ersetzen) else {}
        stand = {**alte, **{g: geprueft for g in neue_gruppen}}
        e["geprueft"] = next(iter(stand.values())) if len(set(stand.values())) == 1 else dict(sorted(stand.items()))
        b = pruefe_eintrag(a, e)
        for s, r, t in b:
            if s != "HINWEIS":
                print(f"{s:8} {i}: {r}: {t}")
        if any(s == "FEHLER" for s, _, _ in b):
            fehler += 1
            continue
        ziel = herkunft_datei.get(i) or daten / f"{a.get('gebiet')}.yaml"
        reihenfolge = ["id", "quelle_hash", "theorie_hash", "herkunft", "geprueft", "en", "en_theorie"]
        neu[ziel].append({k: e[k] for k in reihenfolge if k in e})
    if fehler:
        print(f"{fehler} Einträge mit FEHLER nicht übernommen.")
    if trocken:
        print(f"Trockenlauf: {sum(len(v) for v in neu.values())} Einträge wären geschrieben worden.")
        return 1 if fehler else 0
    for ziel, eintraege in neu.items():
        bestand = {e["id"]: e for e in (yaml.safe_load(ziel.read_text(encoding="utf-8")) or [])} if ziel.exists() else {}
        for e in eintraege:
            bestand[e["id"]] = e
        _yaml_schreiben(ziel, list(bestand.values()))
        print(f"{ziel.relative_to(wurzel)}: {len(eintraege)} Einträge geschrieben")
    return 1 if fehler else 0


def _verschmelzen(alt, neu):
    """Neue Felder zu einem vorhandenen en-Block; Listen werden elementweise ergänzt."""
    out = json.loads(json.dumps(alt))
    for k, v in neu.items():
        if isinstance(v, list) and isinstance(out.get(k), list):
            for i, x in enumerate(v):
                if i < len(out[k]) and isinstance(x, dict) and isinstance(out[k][i], dict):
                    out[k][i] = {**out[k][i], **x}
                elif i >= len(out[k]):
                    out[k].append(x)
                elif x:
                    out[k][i] = x
        elif isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = {**out[k], **v}
        else:
            out[k] = v
    return out


def _woerter(s: str) -> list:
    """Englische Wörter ohne Formeln, kleingeschrieben (für den Abgleich zweier Übersetzungen)."""
    return re.findall(r"[a-zA-Z\u00c0-\u024f][\w'-]*", MATHE.sub(" ", s).lower())


def abgleich(a, b, schwelle=0.85):
    """Vergleicht zwei Übersetzungsstapel [{id, en | en_theorie}]; liefert [(id, pfad, aehnlichkeit, nur_a, nur_b)] unter der Schwelle."""
    import difflib
    zu_b = {e["id"]: e for e in b}
    ausgabe = []
    for e in a:
        f = zu_b.get(e["id"])
        if f is None:
            ausgabe.append((e["id"], "(fehlt in B)", 0.0, [], []))
            continue
        for schluessel in ("en", "en_theorie"):
            if schluessel not in e:
                continue
            la, lb = flach(e[schluessel]), flach(f.get(schluessel) or {})
            for pfad, ta in la.items():
                tb = lb.get(pfad)
                if tb is None:
                    ausgabe.append((e["id"], f"{schluessel}.{pfad}", 0.0, [], ["(fehlt in B)"]))
                    continue
                wa, wb = _woerter(ta), _woerter(tb)
                sm = difflib.SequenceMatcher(None, wa, wb)
                if sm.ratio() < schwelle:
                    nur_a = [w for op, i1, i2, _, _ in sm.get_opcodes() if op != "equal" for w in wa[i1:i2]]
                    nur_b = [w for op, _, _, j1, j2 in sm.get_opcodes() if op != "equal" for w in wb[j1:j2]]
                    ausgabe.append((e["id"], f"{schluessel}.{pfad}", round(sm.ratio(), 2), nur_a, nur_b))
    return sorted(ausgabe, key=lambda x: x[2])


def befehl_abgleich(datei_a, datei_b, schwelle):
    a = json.loads(Path(datei_a).read_text(encoding="utf-8"))
    b = json.loads(Path(datei_b).read_text(encoding="utf-8"))
    abw = abgleich(a, b, schwelle)
    print(f"{len(abw)} Texte unter Ähnlichkeit {schwelle}")
    for i, pfad, r, na, nb in abw:
        print(f"- {i} · {pfad} · {r}\n    A: {' '.join(na)[:200]}\n    B: {' '.join(nb)[:200]}")


def main(argv=None):
    wurzel = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(prog="python -m wissen.uebersetzung")
    sub = ap.add_subparsers(dest="befehl", required=True)
    sub.add_parser("stand", help="Fortschritt je Gebiet")
    e = sub.add_parser("export", help="deutsche Texte als Stapel für Übersetzer")
    e.add_argument("--gebiet")
    e.add_argument("--felder", choices=[*FELDGRUPPEN, "theorie"], default="kern")
    e.add_argument("--max", type=int, default=0)
    e.add_argument("--nur-fehlend", action="store_true")
    e.add_argument("-o", "--ausgabe")
    i = sub.add_parser("einlesen", help="übersetzten Stapel prüfen, stempeln und in daten/uebersetzung/ schreiben")
    i.add_argument("datei")
    i.add_argument("--modell", required=True)
    i.add_argument("--geprueft", choices=sorted(S.UEBERSETZ_GEPRUEFT), default="nein")
    i.add_argument("--ersetzen", action="store_true", help="vorhandenen Eintrag ganz ersetzen (nötig bei veralteten)")
    i.add_argument("--trocken", action="store_true")
    v = sub.add_parser("abgleich", help="zwei unabhängige Übersetzungen vergleichen, Abweichungen zeigen")
    v.add_argument("datei_a")
    v.add_argument("datei_b")
    v.add_argument("--schwelle", type=float, default=0.85)
    sub.add_parser("struktur-stand", help="Fortschritt der Strukturtexte")
    se = sub.add_parser("struktur-export", help="deutsche Strukturtexte als Stapel für Übersetzer")
    se.add_argument("--art", required=True, choices=["gebiet", "strang", "idee", "objekt", "experiment", "quelle", "kapitel"])
    se.add_argument("--max", type=int, default=0)
    se.add_argument("--nur-fehlend", action="store_true")
    se.add_argument("-o", "--ausgabe")
    si = sub.add_parser("struktur-einlesen", help="übersetzten Stapel prüfen, stempeln und in daten/uebersetzung_struktur/ schreiben")
    si.add_argument("datei")
    si.add_argument("--modell", required=True)
    si.add_argument("--geprueft", choices=sorted(S.UEBERSETZ_GEPRUEFT), default="nein")
    si.add_argument("--trocken", action="store_true")
    a = ap.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except AttributeError:
        pass
    if a.befehl.startswith("struktur-"):
        from . import uebersetzung_struktur as U2
        if a.befehl == "struktur-stand":
            U2.befehl_stand(wurzel)
        elif a.befehl == "struktur-export":
            U2.befehl_export(wurzel, a.art, a.max, a.nur_fehlend, a.ausgabe)
        else:
            return U2.befehl_einlesen(wurzel, a.datei, a.modell, a.geprueft, a.trocken)
        return 0
    if a.befehl == "stand":
        befehl_stand(wurzel)
    elif a.befehl == "abgleich":
        befehl_abgleich(a.datei_a, a.datei_b, a.schwelle)
    elif a.befehl == "export":
        befehl_export(wurzel, a.gebiet, a.felder, a.max, a.nur_fehlend, a.ausgabe)
    else:
        return befehl_einlesen(wurzel, a.datei, a.modell, a.ersetzen, a.trocken, a.geprueft)
    return 0


if __name__ == "__main__":
    sys.exit(main())
