"""Prüfen und Ableiten.

Eingabe: das Ergebnis von laden.lade(). Ausgabe: ein Graph-Dict mit abgeleiteten
Feldern (effektiver Status, Neu-prüfen-Markierungen, Tragweite, Verwandte) und
eine Liste von Befunden. Die Regeln stehen jeweils als ein Satz im Kommentar.
"""
import json
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path

from . import schema as S
from .laden import datei_hash  # noqa: F401

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experimente"))
from hashwert import experiment_hash  # noqa: E402


@dataclass
class Befund:
    stufe: str      # FEHLER | WARNUNG | HINWEIS
    regel: str
    id: str
    text: str


def _liste(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def _dep_ids(a):
    """stuetzt_sich_auf darf Strings oder {id, wofuer} enthalten."""
    out = []
    for d in _liste(a.get("stuetzt_sich_auf")):
        out.append(d["id"] if isinstance(d, dict) else d)
    return out


class Graph:
    def __init__(self, roh: dict, wurzel: Path | None = None):
        self.roh = roh
        self.wurzel = wurzel
        self.befunde: list[Befund] = []
        self.knoten = {}
        self.kanten = []          # dicts: von, nach, typ, notiz
        self.gebiete = {g["id"]: g for g in roh["gebiete"]}
        self.objekte = roh["objekte"]
        self.quellen = roh["quellen"]
        self.experimente = {e["id"]: e for e in roh["experimente"]}
        self.straenge = roh.get("straenge") or []
        self.ideen = roh.get("ideen") or []

    # ------------------------------------------------------------------ util
    def befund(self, stufe, regel, id_, text):
        self.befunde.append(Befund(stufe, regel, id_, text))

    # ------------------------------------------------------------- aufbauen
    def aufbauen(self):
        for a in self.roh["aussagen"]:
            i = a.get("id")
            if not i:
                self.befund("FEHLER", "pflichtfeld", "?", f"Aussage ohne id in {a.get('_datei')}")
                continue
            if i in self.knoten:
                self.befund("FEHLER", "doppelt", i, "id kommt mehrfach vor")
                continue
            self.knoten[i] = a
        for i, a in self.knoten.items():
            for d in _dep_ids(a):
                self.kanten.append({"von": i, "nach": d, "typ": "stuetzt_sich_auf"})
            for k in _liste(a.get("kanten")):
                self.kanten.append({"von": i, "nach": k.get("ziel"), "typ": k.get("typ"),
                                    "notiz": k.get("notiz", "")})
        return self

    # --------------------------------------------------------------- prüfen
    def pruefen_schema(self):
        for i, a in self.knoten.items():
            for f in S.PFLICHTFELDER:
                if not a.get(f):
                    self.befund("FEHLER", "pflichtfeld", i, f"Feld '{f}' fehlt")
            art, st, eb = a.get("art"), a.get("status"), a.get("ebene")
            if art not in S.ARTEN:
                self.befund("FEHLER", "vokabular", i, f"unbekannte Art '{art}'")
            if st not in S.STATUS:
                self.befund("FEHLER", "vokabular", i, f"unbekannter Status '{st}'")
            if eb not in S.EBENEN:
                self.befund("FEHLER", "vokabular", i, f"unbekannte Ebene '{eb}'")
            if a.get("gebiet") not in self.gebiete:
                self.befund("FEHLER", "verweis", i, f"unbekanntes Gebiet '{a.get('gebiet')}'")
            elif not i.startswith(a["gebiet"] + "."):
                self.befund("WARNUNG", "id-praefix", i, f"id beginnt nicht mit '{a['gebiet']}.'")
            for o in _liste(a.get("objekte")):
                if o not in self.objekte:
                    self.befund("FEHLER", "verweis", i, f"unbekanntes Objekt '{o}'")
            for q in _liste(a.get("quellen")):
                ref = q.get("ref") if isinstance(q, dict) else q
                if ref not in self.quellen:
                    self.befund("FEHLER", "verweis", i, f"unbekannte Quelle '{ref}'")
                elif isinstance(q, dict) and q.get("geprueft") and q["geprueft"] not in S.QUELLEN_PRUEFUNG:
                    self.befund("FEHLER", "vokabular", i, f"Quellenprüfung '{q['geprueft']}' unbekannt")
            for e in _liste(a.get("experimente")):
                if e not in self.experimente:
                    self.befund("FEHLER", "verweis", i, f"unbekanntes Experiment '{e}'")
            for v in _liste(a.get("voraussetzungen")):
                if not isinstance(v, dict) or not v.get("text"):
                    self.befund("FEHLER", "voraussetzung", i, "Voraussetzung ohne Text")
                    continue
                b = v.get("benutzt")
                if b not in S.BENUTZT:
                    self.befund("FEHLER", "voraussetzung", i, f"benutzt='{b}' (erlaubt: ja/nein/unklar)")
                elif b == "nein":
                    # Regel: eine im Beweis unbenutzte Voraussetzung ist ein Kandidat zum Streichen.
                    self.befund("HINWEIS", "abschwaechbar", i, f"Voraussetzung unbenutzt: {v['text']}")
                elif b == "unklar":
                    self.befund("HINWEIS", "voraussetzung-offen", i, f"Benutzung unklar: {v['text']}")
                if v.get("ref") and v["ref"] not in self.knoten:
                    self.befund("FEHLER", "verweis", i, f"Voraussetzung verweist auf unbekannte Aussage '{v['ref']}'")
            for b in _liste(a.get("bedingungen")):
                if b not in self.knoten:
                    self.befund("FEHLER", "verweis", i, f"Bedingung '{b}' ist keine Aussage")
            if a.get("niveau") and a["niveau"] not in S.NIVEAU:
                self.befund("FEHLER", "vokabular", i, f"unbekanntes Niveau '{a['niveau']}'")
            for f in _liste(a.get("formen")):
                if not (isinstance(f, str) or (isinstance(f, dict) and f.get("tex"))):
                    self.befund("FEHLER", "formen", i, "Form braucht 'tex' (und optional 'notiz')")
            if a.get("en") is not None:
                if not isinstance(a["en"], dict):
                    self.befund("FEHLER", "uebersetzung", i, "'en' muss ein Feld-Wörterbuch sein")
                else:
                    for f in a["en"]:
                        if f not in S.UEBERSETZBAR:
                            self.befund("WARNUNG", "uebersetzung", i, f"Feld '{f}' ist nicht als übersetzbar vorgesehen")
        for k in self.kanten:
            if k["typ"] not in S.KANTEN:
                self.befund("FEHLER", "vokabular", k["von"], f"unbekannter Kantentyp '{k['typ']}'")
            if k["nach"] not in self.knoten:
                self.befund("FEHLER", "verweis", k["von"], f"Kante {k['typ']} → unbekanntes Ziel '{k['nach']}'")
            if k["nach"] == k["von"]:
                self.befund("FEHLER", "schleife", k["von"], "Kante auf sich selbst")
        for e in self.experimente.values():
            for s in _liste(e.get("stuetzt")):
                if s not in self.knoten:
                    self.befund("FEHLER", "verweis", e["id"], f"Experiment stützt unbekannte Aussage '{s}'")
        return self

    def pruefen_inhalt(self):
        for i, a in self.knoten.items():
            art, st = a.get("art"), a.get("status")
            if art in S.PRUEFPFLICHT:
                p = a.get("pruefung") or {}
                fehlend = [S.PRUEFFRAGEN[k].split(" ", 1)[0] for k in S.PRUEFFRAGEN if not p.get(k)]
                if fehlend:
                    self.befund("WARNUNG", "prueffragen", i, "Prüffragen unbeantwortet: " + ", ".join(fehlend))
                bw = a.get("bewertung") or {}
                for f in ("richtigkeit", "neuheit", "nutzen"):
                    if not bw.get(f):
                        self.befund("WARNUNG", "bewertung", i, f"Bewertung '{f}' fehlt")
                if bw.get("neuheit") and bw["neuheit"] not in S.NEUHEIT:
                    self.befund("FEHLER", "vokabular", i, f"Neuheit '{bw['neuheit']}' unbekannt")
                if bw.get("nutzen_stufe") and bw["nutzen_stufe"] not in S.NUTZEN:
                    self.befund("FEHLER", "vokabular", i, f"Nutzenstufe '{bw['nutzen_stufe']}' unbekannt")
                # Regel: wer "vermutlich neu" sagt, muss die Recherche nennen.
                if bw.get("neuheit") == "vermutlich_neu" and not bw.get("recherche"):
                    self.befund("WARNUNG", "neuheit", i, "'vermutlich_neu' ohne dokumentierte Recherche")
            # Regel: Einstiegsaussagen erklären sich zusätzlich in einfachen Worten.
            if a.get("niveau") == "einstieg" and art != "definition" and not a.get("einfach"):
                self.befund("HINWEIS", "einfach", i, "Einstiegsaussage ohne Erklärung in einfachen Worten")
            # Regel: Vermutungen brauchen eine dokumentierte Gegenbeispielsuche.
            if art == "vermutung" and not a.get("gegenbeispielsuche"):
                self.befund("WARNUNG", "gegenbeispiel", i, "Vermutung ohne dokumentierte Gegenbeispielsuche")
            # Regel: Beobachtungen brauchen Experiment oder Quelle.
            if art == "beobachtung" and not (a.get("experimente") or a.get("quellen")):
                self.befund("WARNUNG", "beleg", i, "Beobachtung ohne Experiment und ohne Quelle")
            # Regel: Skizzen nennen ihre Lücken.
            if st == "skizze" and not (a.get("beweis") or {}).get("luecken"):
                self.befund("WARNUNG", "luecken", i, "Status 'skizze' ohne benannte Beweislücken")
            # Regel: bedingte Beweise nennen ihre Bedingungen.
            if st == "bedingt" and not (a.get("bedingungen") or a.get("bedingungen_text")):
                self.befund("WARNUNG", "bedingungen", i, "Status 'bedingt' ohne Bedingungen")
            if st in S.BEWIESEN and st != "definition" and art in S.BEHAUPTUNGEN and not a.get("beweis"):
                self.befund("WARNUNG", "beweis", i, "als bewiesen markiert, aber ohne Beweisangabe")
            if st == "bewiesen_literatur" and not a.get("quellen"):
                self.befund("FEHLER", "quelle", i, "'bewiesen_literatur' ohne Quelle")
        return self


    def pruefen_straenge(self):
        zuordnung = {}
        for s in self.straenge:
            sid = s.get("id", "?")
            if s.get("gebiet") not in self.gebiete:
                self.befund("FEHLER", "strang", sid, f"unbekanntes Gebiet '{s.get('gebiet')}'")
            for f in ("name", "verlauf"):
                if not s.get(f):
                    self.befund("WARNUNG", "strang", sid, f"Feld '{f}' fehlt")
            for m in _liste(s.get("mitglieder")):
                if m not in self.knoten:
                    self.befund("FEHLER", "strang", sid, f"Mitglied '{m}' ist keine Aussage")
                    continue
                if self.knoten[m].get("gebiet") != s.get("gebiet"):
                    self.befund("WARNUNG", "strang", sid, f"Mitglied {m} gehört zu einem anderen Gebiet")
                if m in zuordnung:
                    self.befund("FEHLER", "strang", m, f"in zwei Strängen: {zuordnung[m]} und {sid}")
                zuordnung[m] = sid
        if self.straenge:
            for i, a in self.knoten.items():
                a["_strang"] = zuordnung.get(i)
                # Regel: auf der Karte soll jede Aussage einem Verlauf angehören.
                if i not in zuordnung:
                    self.befund("HINWEIS", "ohne-strang", i, "gehört keinem Forschungsstrang an")
        return self

    def pruefen_ideen(self):
        """Ideen-Landkarte: Hauptideen (ohne eltern) und Unterideen; Zuordnung über Stränge, Ausnahmen über 'aussagen'."""
        if not self.ideen:
            return self
        haupt = {i["id"] for i in self.ideen if not i.get("eltern")}
        strang_ids = {s.get("id") for s in self.straenge}
        von_strang, von_aussage = {}, {}
        for i in self.ideen:
            iid = i.get("id", "?")
            for f in ("name", "kurz"):
                if not i.get(f):
                    self.befund("WARNUNG", "idee", iid, f"Feld '{f}' fehlt")
            if not i.get("eltern"):
                continue
            if i["eltern"] not in haupt:
                self.befund("FEHLER", "idee", iid, f"Eltern '{i['eltern']}' ist keine Hauptidee")
            for s in _liste(i.get("straenge")):
                if s not in strang_ids:
                    self.befund("FEHLER", "idee", iid, f"unbekannter Strang '{s}'")
                elif s in von_strang:
                    self.befund("FEHLER", "idee-doppelt", s, f"Strang in {von_strang[s]} und {iid}")
                else:
                    von_strang[s] = iid
            for a in _liste(i.get("aussagen")):
                if a not in self.knoten:
                    self.befund("FEHLER", "idee", iid, f"unbekannte Aussage '{a}'")
                elif a in von_aussage:
                    self.befund("FEHLER", "idee-doppelt", a, f"Aussage in {von_aussage[a]} und {iid}")
                else:
                    von_aussage[a] = iid
        for k, a in self.knoten.items():
            a["_idee"] = von_aussage.get(k) or von_strang.get(a.get("_strang"))
            if not a["_idee"]:
                self.befund("HINWEIS", "ohne-idee", k, "gehört keiner Idee der Landkarte an")
        return self

    # ------------------------------------------------------------ ableiten
    def _ordnungskanten(self):
        """Kanten, die die Auswertungsreihenfolge bestimmen: (vorher, nachher)."""
        paare = []
        for k in self.kanten:
            if k["nach"] not in self.knoten:
                continue
            if S.KANTEN.get(k["typ"], {}).get("beweis"):
                paare.append((k["nach"], k["von"]))      # Abhängigkeit zuerst
            elif k["typ"] == "entlastet":
                paare.append((k["von"], k["nach"]))      # Beweis vor Annahme
        return paare

    def topologisch(self):
        nach = defaultdict(list)
        grad = {i: 0 for i in self.knoten}
        for vor, nachher in self._ordnungskanten():
            nach[vor].append(nachher)
            grad[nachher] += 1
        reihe = [i for i, g in grad.items() if g == 0]
        out = []
        while reihe:
            i = reihe.pop()
            out.append(i)
            for j in nach[i]:
                grad[j] -= 1
                if grad[j] == 0:
                    reihe.append(j)
        if len(out) < len(self.knoten):
            rest = [i for i in self.knoten if i not in set(out)]
            # Regel: Beweisketten dürfen keine Kreise enthalten.
            self.befund("FEHLER", "kreis", rest[0], "Kreis in Beweisabhängigkeiten: " + ", ".join(sorted(rest)[:8]))
            out += rest
        return out

    def ableiten(self):
        deps = {i: [d for d in _dep_ids(a) if d in self.knoten] for i, a in self.knoten.items()}
        for k in self.kanten:
            if k["typ"] == "folgt_unmittelbar_aus" and k["nach"] in self.knoten:
                deps[k["von"]].append(k["nach"])
        eingehend = defaultdict(list)
        for k in self.kanten:
            eingehend[k["nach"]].append(k)
        eff = {}
        wurzel = {}
        auto_bed = defaultdict(list)
        for i in self.topologisch():
            a = self.knoten[i]
            st = a.get("status")
            art = a.get("art")
            if st == "annahme":
                entl = [k["von"] for k in eingehend[i] if k["typ"] == "entlastet"]
                # Regel: eine Annahme gilt als entlastet, sobald ein bewiesener Satz sie entlastet.
                if any(eff.get(e) in S.BEWIESEN for e in entl):
                    st = "entlastet"
                elif entl:
                    a["_teilweise_entlastet_durch"] = entl
            if st in S.BEWIESEN and art in S.BEHAUPTUNGEN:
                schwach = [d for d in deps[i] if S.STATUS.get(eff.get(d), {}).get("rang", 5) < 5
                           and eff.get(d) not in S.TOT]
                if schwach:
                    # Regel: "bewiesen" auf schwächerem Fundament wird zu "bedingt".
                    st = "bedingt"
                    auto_bed[i] += schwach
                    self.befund("WARNUNG", "status-zu-stark", i,
                                "als bewiesen markiert, stützt sich aber auf " +
                                ", ".join(f"{d} ({eff.get(d)})" for d in schwach))
            eff[i] = st
            # Wurzelbedingungen: die unbewiesenen Enden der Beweiskette (bedingte Zwischenglieder zählen nicht)
            wurzeln = set()
            for d in deps[i]:
                sd = eff.get(d)
                if sd in S.TOT or sd is None:
                    continue
                if sd == "bedingt":
                    wurzeln |= wurzel.get(d, set())
                elif S.STATUS.get(sd, {}).get("rang", 5) < 5:
                    wurzeln.add(d)
            if st == "bedingt":
                for d in _liste(a.get("bedingungen")):
                    sd = eff.get(d)
                    if sd == "bedingt":
                        wurzeln |= wurzel.get(d, set())
                    elif sd not in S.BEWIESEN:
                        wurzeln.add(d)
                genannt = set()
                for d in _liste(a.get("bedingungen")):
                    genannt |= wurzel.get(d, set()) | {d}
                for d in sorted(wurzeln - genannt):
                    # Regel: jede offene Grundannahme der Kette muss als Bedingung sichtbar sein.
                    self.befund("HINWEIS", "bedingung-implizit", i,
                                f"Grundannahme {d} ({eff.get(d)}) der Beweiskette ist nicht als Bedingung genannt")
            wurzel[i] = wurzeln
        # widerlegt / ersetzt / korrigiert
        for k in self.kanten:
            ziel = k["nach"]
            if ziel not in self.knoten:
                continue
            if k["typ"] == "widerlegt" and self.knoten[ziel].get("status") not in S.TOT | {"gescheitert"}:
                self.befund("WARNUNG", "widerspruch", ziel,
                            f"wird von {k['von']} widerlegt, Status ist aber '{self.knoten[ziel].get('status')}'")
            schon = {n.get("gegen") for n in _liste(self.knoten[ziel].get("nachgeprueft"))}
            if (k["typ"] in ("ersetzt", "korrigiert") and self.knoten[ziel].get("status") not in S.TOT
                    and k["von"] not in schon):
                self.befund("HINWEIS", "fassung", ziel,
                            f"wird von {k['von']} {k['typ']}; Status '{self.knoten[ziel].get('status')}' prüfen")
        # Neu prüfen: Abhängigkeit (transitiv) von Totem, sofern nicht nachgeprüft
        nachgeprueft = {i: {n.get("gegen") for n in _liste(a.get("nachgeprueft"))} for i, a in self.knoten.items()}
        tot = {i for i, s in eff.items() if s in S.TOT}
        for k in self.kanten:
            if k["typ"] == "widerlegt" and k["nach"] in self.knoten:
                tot.add(k["nach"])
        neu_pruefen = {}
        for i in self.topologisch():
            gruende = []
            for d in deps[i]:
                if d in tot and d not in nachgeprueft[i]:
                    gruende.append([d])
                elif d in neu_pruefen:
                    for pfad in neu_pruefen[d]:
                        if pfad[-1] not in nachgeprueft[i]:
                            gruende.append([d] + pfad)
            if gruende and i not in tot:
                neu_pruefen[i] = gruende[:5]
                self.befund("WARNUNG", "neu-pruefen", i,
                            "hängt an widerlegter/überholter Aussage: " +
                            "; ".join(" → ".join(p) for p in gruende[:3]))
        # Tragweite: wie viele Aussagen hängen (transitiv) daran
        rueck = defaultdict(set)
        for i, ds in deps.items():
            for d in ds:
                rueck[d].add(i)
        tragweite = {}
        for i in self.knoten:
            seen, stack = set(), [i]
            while stack:
                for j in rueck[stack.pop()]:
                    if j not in seen:
                        seen.add(j)
                        stack.append(j)
            tragweite[i] = len(seen)
        for i, a in self.knoten.items():
            # Standard-Niveau, wenn nicht eingetragen: Grundlagen = Einstieg, Forschungsarten und Basiswechsel = Forschung
            if a.get("niveau"):
                a["_niveau"] = a["niveau"]
            elif a.get("gebiet") == "gr":
                a["_niveau"] = "einstieg"
            elif a.get("gebiet") == "bw" or a.get("art") in {"vermutung", "frage", "beobachtung", "gescheitert", "korrektur", "gegenbeispiel", "annahme"}:
                a["_niveau"] = "forschung"
            else:
                a["_niveau"] = "vertiefung"
            a["_status_eff"] = eff.get(i, a.get("status"))
            a["_bedingungen_auto"] = auto_bed.get(i, [])
            a["_wurzelbedingungen"] = sorted(wurzel.get(i, set()))
            a["_neu_pruefen"] = neu_pruefen.get(i, [])
            a["_tragweite"] = tragweite[i]
            a["_abhaengig"] = sorted(set(deps[i]))
            a["_traegt"] = sorted(rueck[i])
        self.eff = eff
        return self

    def verwandte(self, n=8):
        nachbarn = defaultdict(lambda: defaultdict(float))
        gruende = defaultdict(lambda: defaultdict(list))
        for k in self.kanten:
            if k["nach"] not in self.knoten or k["typ"] not in S.KANTEN:
                continue
            w = S.KANTEN[k["typ"]]["gewicht"]
            u, v = k["von"], k["nach"]
            nachbarn[u][v] += w
            nachbarn[v][u] += w
            gruende[u][v].append(S.KANTEN[k["typ"]]["label"])
            gruende[v][u].append(S.KANTEN[k["typ"]]["label"] + " (rückwärts)")
        for e in self.experimente.values():
            st = [s for s in _liste(e.get("stuetzt")) if s in self.knoten]
            for x in st:
                for y in st:
                    if x < y:
                        nachbarn[x][y] += 1.0
                        nachbarn[y][x] += 1.0
                        gruende[x][y].append(f"gemeinsames Experiment {e['id']}")
                        gruende[y][x].append(f"gemeinsames Experiment {e['id']}")
        obj = {i: set(_liste(a.get("objekte"))) for i, a in self.knoten.items()}
        ids = list(self.knoten)
        for u in ids:
            score = defaultdict(float)
            why = defaultdict(list)
            for v, w in nachbarn[u].items():
                score[v] += w
                why[v] += gruende[u][v]
            for v, w in list(nachbarn[u].items()):
                for x, w2 in nachbarn[v].items():
                    if x != u:
                        score[x] += 0.35 * min(w, w2)
                        if f"über {v}" not in why[x]:
                            why[x].append(f"über {v}")
            for v in ids:
                if v == u or not obj[u] or not obj[v]:
                    continue
                j = len(obj[u] & obj[v]) / len(obj[u] | obj[v])
                if j > 0:
                    score[v] += 2.5 * j
                    if j >= 0.34:
                        why[v].append("gemeinsame Objekte: " + ", ".join(sorted(obj[u] & obj[v])))
            liste = sorted((s, v) for v, s in score.items() if s > 0.6)[::-1][:n]
            self.knoten[u]["_verwandt"] = [{"id": v, "score": round(s, 2), "warum": why[v][:3]} for s, v in liste]
        return self

    # ---------------------------------------------------------- Experimente
    def experimente_einlesen(self):
        for e in self.experimente.values():
            e["_zustand"] = "fehlt"
            if not self.wurzel:
                continue
            skript = self.wurzel / e.get("skript", "")
            erg = self.wurzel / e.get("ergebnis", "")
            if not skript.is_file():
                self.befund("FEHLER", "experiment", e["id"], f"Skript fehlt: {e.get('skript')}")
                continue
            if not erg.is_file():
                self.befund("WARNUNG", "experiment", e["id"], "noch nicht gelaufen (kein Ergebnis)")
                continue
            daten = json.loads(erg.read_text(encoding="utf-8"))
            e["_ergebnis"] = daten
            # Regel: ein Ergebnis gilt nur, wenn es vom aktuellen Skript stammt.
            if daten.get("skript_hash") != experiment_hash(skript):
                e["_zustand"] = "veraltet"
                self.befund("WARNUNG", "experiment", e["id"], "Skript seit dem letzten Lauf geändert (Ergebnis veraltet)")
            else:
                e["_zustand"] = "aktuell"
            for f in daten.get("fehlgeschlagen") or []:
                self.befund("WARNUNG", "experiment-pruefung", e["id"], f"Prüfung fehlgeschlagen: {f}")
            if e.get("bild"):
                b = self.wurzel / e["bild"]
                if b.is_file():
                    e["_svg"] = b.read_text(encoding="utf-8")
                else:
                    self.befund("WARNUNG", "experiment", e["id"], f"Bild fehlt: {e['bild']}")
        return self

    # ------------------------------------------------------------- gesamt
    def alles(self):
        from . import lehrbuch, uebersetzung
        g = (self.aufbauen().pruefen_schema().pruefen_inhalt().pruefen_straenge().pruefen_ideen()
             .ableiten().verwandte().experimente_einlesen())
        g = lehrbuch.pruefen(g)     # Theorie-Blöcke sind schon angehängt (laden.py)
        return uebersetzung.einbinden(g)

    def statistik(self):
        z = defaultdict(int)
        for a in self.knoten.values():
            z["niveau:" + a.get("_niveau", "?")] += 1
            z["status:" + a["_status_eff"]] += 1
            z["art:" + a.get("art", "?")] += 1
            if a.get("_en_stand"):
                z["uebersetzung:" + a["_en_stand"]] += 1
            z["ebene:" + a.get("ebene", "?")] += 1
        qn = 0
        for a in self.knoten.values():
            for q in _liste(a.get("quellen")):
                ref = q.get("ref") if isinstance(q, dict) else q
                g = (q.get("geprueft") if isinstance(q, dict) else None) or self.quellen.get(ref, {}).get("geprueft")
                if g == "nein":
                    qn += 1
        z["zitate_ungeprueft"] = qn
        for b in self.befunde:
            z["befund:" + b.stufe] += 1
        return dict(sorted(z.items()))

    def befunde_dict(self):
        ordnung = {"FEHLER": 0, "WARNUNG": 1, "HINWEIS": 2}
        return [asdict(b) for b in sorted(self.befunde, key=lambda b: (ordnung[b.stufe], b.regel, b.id))]
