# Ideen-Landkarte Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eine neue Startansicht „Ideen“, die die Tetration als Netz benannter Haupt- und Unterideen zeigt (Lage aus der Kopplung berechnet, Füllstand und Herkunft als Ringe, aufklappbar bis zur Aussage), dazu eine Textfassung für Agenten.

**Architecture:** Daten in `daten/ideen.yaml` (Lader wie bei den Strängen). `wissen/graph.py` prüft die Zuordnung und setzt `_idee`. Das neue reine Rechenmodul `wissen/landkarte.py` berechnet Kopplung, Lage, Haupt/Hilfs-Einstufung, Füllstand, Herkunft und den Text `ausgabe/landkarte.md`; `baue.py` legt das Ergebnis als `ideen` in `window.GRAPH`. Die Ansicht ist ein eigenes Canvas-Modul `ansicht/js/02b_ideen.js` neben der bestehenden Karte.

**Tech Stack:** Python 3.13 (Miniconda, PyYAML, pytest), Canvas 2D in einer Offline-HTML-Datei.

**Spec:** `docs/superpowers/specs/2026-10-12-ideen-landkarte-design.md`

## Global Constraints

- Bauen und testen nur mit `C:/Users/janis/miniconda3/python.exe` (plain `python` hat kein PyYAML).
- Backslashes in YAML/Skripten nur mit Write/Edit schreiben, nie per Bash-Heredoc.
- `experimente/gemeinsam.py` nicht anfassen (sonst werden alle Experimente „veraltet“).
- Lage deterministisch: gleiche Daten ⇒ gleiche Koordinaten (kein `random` ohne festen Startwert).
- Sprache der Oberfläche Deutsch, einfache Wörter; Karte immer dunkel (Janis' Wahl), Farben wie `KF` in `02_karte.js`.
- Die alte Karte (Tab „Karte“) bleibt unverändert erreichbar.

## Review Focus

- Aussage, deren Strang keiner Unteridee zugeordnet ist → Hinweis „ohne-idee“, Ansicht legt sie in eine Sammel-Unteridee „Ohne Zuordnung“ statt sie zu verlieren (Test in Task 2).
- Unteridee ohne Aussagen (leerer Strang) → wird gezeichnet, Radius > 0, kein Absturz durch Division (Test in Task 2).
- Sehr ungleiche Grössen (Basiswechsel ~60 Aussagen, Reell 1 Aussage) → keine Überlappung, kleine Kreise bleiben klickbar (Test in Task 2: Überlappung; Mindestradius).
- Suche bzw. `#a=<id>` im Link führt in der Ideen-Ansicht direkt zur Aussage (Pfad öffnet sich) (Prüfung im Browser, Task 4).
- Fenstergrösse ändern oder Detailtafel öffnen, während eine Idee offen ist → Kamera passt sich an, nichts verschwindet (Browser, Task 4).

---

### Task 1: Ideen-Daten, Lader und Prüfregeln

**Files:**
- Create: `daten/ideen.yaml`
- Modify: `wissen/laden.py` (Eintrag `"ideen": _liste(daten, "ideen", h)` in `lade`)
- Modify: `wissen/graph.py` (neue Methode `pruefen_ideen`, in `alles()` nach `pruefen_straenge`)
- Test: `tests/test_ideen.py`

**Interfaces:**
- Produces: `roh["ideen"]`: Liste von Dicts `{id, name, kurz, eltern?, straenge?, aussagen?}`; nach `Graph.alles()` trägt jede Aussage `a["_idee"]` (id einer Unteridee oder `None`); `Graph.ideen` = Liste der gültigen Ideen.

- [ ] **Step 1: Test schreiben** — `tests/test_ideen.py` mit den Hilfen aus `tests/test_graph.py` (`knoten`, `regeln`; `graph` um Parameter `ideen` erweitert):

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tests.test_graph import knoten, regeln  # noqa: E402
from wissen.graph import Graph  # noqa: E402


def graph(*aussagen, straenge=None, ideen=None):
    roh = {"gebiete": [{"id": "g", "name": "G"}], "objekte": {"o": {"name": "o"}},
           "quellen": {"q": {"geprueft": "nein"}}, "experimente": [], "aussagen": list(aussagen),
           "straenge": straenge or [], "ideen": ideen or []}
    return Graph(roh).alles()


STR = [{"id": "g.s1", "gebiet": "g", "name": "S1", "verlauf": "x", "mitglieder": ["g.a", "g.b"]},
       {"id": "g.s2", "gebiet": "g", "name": "S2", "verlauf": "x", "mitglieder": ["g.c"]}]
IDEEN = [{"id": "h", "name": "H", "kurz": "x"},
         {"id": "h.u", "eltern": "h", "name": "U", "kurz": "x", "straenge": ["g.s1"]},
         {"id": "h.v", "eltern": "h", "name": "V", "kurz": "x", "straenge": ["g.s2"], "aussagen": ["g.b"]}]


def test_zuordnung_ueber_strang_und_ausnahme():
    g = graph(knoten("a"), knoten("b"), knoten("c"), straenge=STR, ideen=IDEEN)
    assert g.knoten["g.a"]["_idee"] == "h.u"
    assert g.knoten["g.b"]["_idee"] == "h.v"      # Ausnahme schlägt Strang
    assert g.knoten["g.c"]["_idee"] == "h.v"
    assert not {r for r, _ in regeln(g) if r.startswith("idee") or r == "ohne-idee"}


def test_ohne_idee_ist_hinweis():
    g = graph(knoten("a"), knoten("b"), knoten("c"), straenge=STR, ideen=IDEEN[:2])
    assert ("ohne-idee", "g.c") in regeln(g)
    assert g.knoten["g.c"]["_idee"] is None


def test_fehlerhafte_ideen():
    kaputt = IDEEN + [{"id": "x.u", "eltern": "fehlt", "name": "X", "kurz": "x", "straenge": ["g.s1", "g.nix"],
                       "aussagen": ["g.nix"]}]
    g = graph(knoten("a"), knoten("b"), knoten("c"), straenge=STR, ideen=kaputt)
    r = regeln(g)
    assert ("idee", "x.u") in r                       # Eltern unbekannt / Strang unbekannt / Aussage unbekannt
    assert any(rr == "idee-doppelt" for rr, _ in r)   # g.s1 in zwei Unterideen


def test_ohne_ideendatei_keine_hinweise():
    g = graph(knoten("a"), straenge=[{"id": "g.s", "gebiet": "g", "name": "S", "verlauf": "x", "mitglieder": ["g.a"]}])
    assert not {r for r, _ in regeln(g) if r in ("ohne-idee", "idee", "idee-doppelt")}
```

- [ ] **Step 2: Test laufen lassen, muss scheitern** — `C:/Users/janis/miniconda3/python.exe -m pytest -q tests/test_ideen.py` → FAIL (`_idee` fehlt).

- [ ] **Step 3: Umsetzen.** In `Graph.__init__` `self.ideen = roh.get("ideen") or []`. Methode:

```python
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
```

In `alles()`: `.pruefen_straenge().pruefen_ideen().ableiten()`. In `wissen/laden.py`: `"ideen": _liste(daten, "ideen", h),` neben `"straenge"`. In `baue.py` `g_export`: `"_idee"` in die Liste der exportierten `_`-Felder.

- [ ] **Step 4: `daten/ideen.yaml` anlegen** — 8 Hauptideen, Unterideen mit den Strängen (vollständige Zuordnung aller 48 Stränge):

| Unteridee | Stränge | Ausnahmen (`aussagen`) |
|---|---|---|
| gleichung.begriffe | fg.begriffe, fg.rand | |
| gleichung.freiheit | fg.freiheit | |
| gleichung.festfrei | gr.identitaeten | |
| gleichung.werte | gr.werte, gr.kurve, gr.ableitung, gr.basis | |
| gleichung.integrale | gr.integrale, gr.turmintegrale | |
| gleichung.slog | gr.slog | |
| gleichung.superwurzeln | gr.umkehrungen | |
| dynamik.fixpunkte | fp.fixpunkte | |
| dynamik.koenigs | it.hyperbolisch | |
| dynamik.parabolisch | it.parabolisch | |
| dynamik.indifferent | it.indifferent | |
| turm.konvergenz | fp.turm | |
| turm.periode2 | fp.periode2 | |
| turm.regionen | fp.parameterraum | |
| turm.spiralen | fp.spiralen_mittel | |
| fortsetzung.regulaer | kl.regulaer, kl.zwei_abel, kl.bruch | |
| fortsetzung.eta | kl.eta, br.engpass | |
| fortsetzung.kneser | kn.konstruktion, kn.rand | |
| fortsetzung.andere | kn.andere | |
| fortsetzung.eindeutigkeit | kn.eindeutigkeit | |
| fortsetzung.reell | br.reell | |
| fortsetzung.halbexp | kn.halbe_exp | |
| basiswechsel.gleiches_l | kl.gleiches_L | |
| basiswechsel.abschaelen | — | bw.def.F_bd, bw.prop.gruppoid, bw.prop.ganzzahl, bw.ann.hyp22, bw.ann.hyp24, bw.lem.wachstumsluecke, bw.lem.stoppregel, bw.satz.schwanzkonvergenz |
| basiswechsel.phase | bw.phasengesetz | |
| basiswechsel.kreis | bw.kreis | |
| basiswechsel.moden | bw.moden | |
| basiswechsel.antwort | bw.antwort, bw.rh | |
| basiswechsel.rekonstruktion | bw.rekonstruktion | |
| stufen.pentation | kn.hoehere_stufen | |
| stufen.etas | kn.etas | |
| stufen.hyperop | — | br.def.hyperoperation, br.satz.ackermann, br.satz.ordinal_tetration, br.satz.hyperop_varianten |
| zahlen.transzendenz | ar.klassisch, ar.gleichung, ar.zaehlen | |
| zahlen.ks | ar.ks | |
| zahlen.kongruenzen | ar.kongruenzen, br.diskret | |
| bruecken.logik | br.logik | |
| bruecken.fourier | br.kreis | |

(Abweichung vom Entwurf: Ausnahmen stehen als Liste `aussagen:` in der Unteridee statt als Feld `idee:` in jeder Aussage — gleiche Wirkung, aber die Aussagen-Dateien bleiben unberührt.) Jede Idee bekommt `name` und einen Satz `kurz` in einfacher Sprache.

- [ ] **Step 5: Tests und Prüfung** — `pytest -q tests` → alle grün; `baue.py --nur-pruefen` → 0 Fehler, kein „ohne-idee“.

- [ ] **Step 6: Commit** `git add daten/ideen.yaml wissen tests baue.py; git commit -m "Ideen-Landkarte: Ideen-Daten und Prüfregeln"`

### Task 2: Rechenmodul `wissen/landkarte.py` und Export

**Files:**
- Create: `wissen/landkarte.py`
- Modify: `baue.py` (`g_export` bekommt `"ideen": landkarte.berechne(g, roh)`; nach dem Bau `ausgabe/landkarte.md` schreiben)
- Test: `tests/test_landkarte.py`

**Interfaces:**
- Consumes: `g.knoten` (mit `_idee`, `_status_eff`, `_tragweite`, `objekte`, `quellen`, `art`, `bewertung`), `g.kanten` (`von`, `nach`, `typ`), `roh["ideen"]`, `roh["quellen"]`.
- Produces: `berechne(g, roh) -> dict` mit
  `{"haupt": [{id,name,kurz,x,y,r,n,fuell:{gruppe:anzahl},herkunft:{literatur,forum,eigen}}],
    "unter": [{id,eltern,name,kurz,x,y,r,n,fuell,herkunft,offen:[ids]}],
    "kopplung_haupt": [[a,b,w]], "kopplung_unter": [[a,b,w]],
    "aussagen": {id: {"x","y","r","haupt": bool}}}` — Koordinaten absolut in Weltkoordinaten;
  `kopplung(g, unter_von) -> dict[(a,b)] -> float`; `text_md(daten) -> str`.

- [ ] **Step 1: Test schreiben** (`tests/test_landkarte.py`): Mini-Graph mit zwei Hauptideen, drei Unterideen.

```python
from wissen import landkarte as L
from tests.test_ideen import graph, knoten


def mini():
    ideen = [{"id": "h1", "name": "H1", "kurz": "x"}, {"id": "h2", "name": "H2", "kurz": "x"},
             {"id": "h1.a", "eltern": "h1", "name": "A", "kurz": "x", "aussagen": ["g.a1", "g.a2"]},
             {"id": "h1.b", "eltern": "h1", "name": "B", "kurz": "x", "aussagen": ["g.b1"]},
             {"id": "h2.c", "eltern": "h2", "name": "C", "kurz": "x", "aussagen": ["g.c1"]},
             {"id": "h2.leer", "eltern": "h2", "name": "Leer", "kurz": "x"}]
    a1 = knoten("a1", objekte=["o"]); a2 = knoten("a2", stuetzt_sich_auf=["g.a1"])
    b1 = knoten("b1", stuetzt_sich_auf=["g.a1"], objekte=["o"])
    c1 = knoten("c1", art="vermutung", status="vermutung", kanten=[{"typ": "verwandt", "ziel": "g.a2"}])
    roh_ideen = ideen
    return graph(a1, a2, b1, c1, ideen=roh_ideen), roh_ideen


def test_kopplung_von_hand():
    g, ideen = mini()
    w = L.kopplung(g, {k: a["_idee"] for k, a in g.knoten.items()})
    # a↔b: eine stuetzt_sich_auf-Kante (1.0) + Objekte: Jaccard({o},{o}) = 1 · min(2,1) · 0.5 = 0.5
    assert abs(w[("h1.a", "h1.b")] - 1.5) < 1e-9
    assert abs(w[("h1.a", "h2.c")] - 0.25) < 1e-9          # verwandt
    assert ("h1.b", "h2.c") not in w


def test_lage_deterministisch_ohne_ueberlappung():
    g, ideen = mini()
    d1 = L.berechne(g, {"ideen": ideen, "quellen": {"q": {"art": "literatur"}}})
    d2 = L.berechne(g, {"ideen": ideen, "quellen": {"q": {"art": "literatur"}}})
    assert d1 == d2
    H = d1["haupt"]
    for i in range(len(H)):
        for j in range(i + 1, len(H)):
            assert ((H[i]["x"] - H[j]["x"]) ** 2 + (H[i]["y"] - H[j]["y"]) ** 2) ** 0.5 >= H[i]["r"] + H[j]["r"] - 1e-6
    hp = {h["id"]: h for h in H}
    for u in d1["unter"]:
        p = hp[u["eltern"]]
        assert ((u["x"] - p["x"]) ** 2 + (u["y"] - p["y"]) ** 2) ** 0.5 + u["r"] <= p["r"] + 1e-6
        assert u["r"] > 0
    leer = next(u for u in d1["unter"] if u["id"] == "h2.leer")
    assert leer["n"] == 0


def test_hauptaussagen_hoechstens_ein_drittel_mindestens_eine():
    g, ideen = mini()
    d = L.berechne(g, {"ideen": ideen, "quellen": {}})
    in_a = [d["aussagen"][i]["haupt"] for i in ("g.a1", "g.a2")]
    assert sum(in_a) == 1 and d["aussagen"]["g.a1"]["haupt"]   # a1 trägt zwei Aussagen


def test_fuellstand_und_offen():
    g, ideen = mini()
    d = L.berechne(g, {"ideen": ideen, "quellen": {}})
    c = next(u for u in d["unter"] if u["id"] == "h2.c")
    assert c["fuell"] == {"vermutung": 1} and c["offen"] == ["g.c1"]
    assert "Leer" in L.text_md(d) and "g.c1" in L.text_md(d)


def test_ohne_idee_landet_in_sammelidee():
    g, ideen = mini()
    g.knoten["g.c1"]["_idee"] = None
    d = L.berechne(g, {"ideen": ideen, "quellen": {}})
    assert d["aussagen"]["g.c1"]
    assert any(u["id"] == "_ohne" for u in d["unter"])
```

- [ ] **Step 2: scheitern lassen** — `pytest -q tests/test_landkarte.py` → FAIL (Modul fehlt).

- [ ] **Step 3: Umsetzen** (`wissen/landkarte.py`). Regeln:
  - Kantengewichte: `{stuetzt_sich_auf, verallgemeinert, ersetzt, korrigiert, widerlegt}: 1.0`, `{stuetzt, motiviert}: 0.5`, sonst `0.25`; Kanten aus `g.kanten` plus `stuetzt_sich_auf`-Listen (wie `_abhaengig`). Zwischen Aussagen derselben Unteridee zählt nichts.
  - Objekt-Anteil: `0.5 · |O_U ∩ O_V| / |O_U ∪ O_V| · min(n_U, n_V)`, nur wenn > 0 und beide ≥ 1 Aussage.
  - Statusgruppen: `bewiesen`={definition, bewiesen, bewiesen_literatur, entlastet}, `bedingt`={bedingt, skizze}, `numerisch`, `vermutung`={vermutung, annahme}, `offen`={offen}, `verworfen`={widerlegt, gescheitert, ueberholt}.
  - Herkunft je Aussage: `literatur` wenn eine Quelle `art: literatur`, sonst `forum` wenn eine `forum`/`notiz`, sonst `eigen`.
  - `offen`: Aussagen mit Gruppe vermutung/offen oder Status bedingt, sortiert nach `_tragweite` absteigend.
  - Radien: Aussage `r = 5 + 2·√min(tragweite, 25)` (Haupt ×1.4); Unteridee `r_u = 26 + 15·√n`; Hauptidee `r_h = max(120, 1.35·√Σ r_u²) + 20`.
  - Unterideen in der Hauptidee: Startlage Goldener-Winkel-Spirale (Reihenfolge nach id), 300 Schritte: Abstossung bei Überlappung (Paar trennen auf `r_i+r_j+8`), Anziehung zu Geschwistern ∝ w (Ziel `r_i+r_j+8`), Zug zum Rand Richtung gekoppelter fremder Hauptideen (Gewicht 0.3·w, erst nach dem Hauptideen-Layout — daher zweistufig: Hauptideen zuerst mit `r_h`), danach Projektion in `|p| ≤ r_h − r_u − 10`; Schrittweite abklingend 0.5 → 0.02.
  - Hauptideen: Startlage auf Kreis (Reihenfolge nach Gesamtkopplung absteigend), 600 Schritte: Paare mit w > 0 ziehen auf `r_i+r_j+60/(1+w/5)`, alle Paare stossen ab, wenn `d < r_i+r_j+40`; schwacher Zug zur Mitte (0.002); Schluss: 200 Schritte reine Überlappungsauflösung.
  - Aussagen in der Unteridee: Rang = (#Abhängige in U − #Voraussetzungen in U, tragweite) absteigend; Platz i auf Spirale `ρ = (r_u − 12)·√((i+0.5)/n)`, Winkel `i·2.39996`.
  - Haupt-Aussagen: Kandidaten = Art in {satz, vermutung, frage, konstruktion, definition} mit `tragweite ≥ 1` oder `nutzen_stufe == "hoch"`; sortiert nach (tragweite, nutzen hoch) absteigend; davon die ersten `max(1, ceil(n/3))`; gibt es keinen Kandidaten, die Aussage mit grösster Tragweite.
  - Aussagen ohne `_idee` → Unteridee `_ohne` („Ohne Zuordnung“) unter einer Hauptidee `_ohne` (nur wenn nötig).
  - `text_md(d)`: Überschrift, je Hauptidee (nach Grösse) Zeile mit Füllstand und Herkunft, darunter je Unteridee `- **Name** (n): bewiesen 5, bedingt 2, … · Literatur 3 / Forum 1 / eigen 3 · offen: id1, id2, id3`; am Ende „Stärkste Kopplungen zwischen Hauptideen“ (Top 10) und „Stärkste Kopplungen über Ideengrenzen“ (Top 15 Unterideen-Paare verschiedener Hauptideen).
  - Alle Zahlen auf 1 Nachkommastelle runden (stabile JSON-Ausgabe).

- [ ] **Step 4: grün** — `pytest -q tests` → alles grün.

- [ ] **Step 5: in den Bau** — `baue.py`: `from wissen import landkarte`; in `g_export` `"ideen": landkarte.berechne(g, roh)`; in `main` nach dem HTML: `(aus / "landkarte.md").write_text(landkarte.text_md(daten["ideen"]), encoding="utf-8")`. `baue.py` laufen lassen, `ausgabe/landkarte.md` lesen und auf Plausibilität prüfen (Basiswechsel hat die meisten Aussagen; Welle/Phase gekoppelt mit Fortsetzung).

- [ ] **Step 6: Commit** `git commit -m "Ideen-Landkarte: Kopplung, Lage, Füllstand, Textfassung"`

### Task 3: Ansicht `02b_ideen.js`

**Files:**
- Create: `ansicht/js/02b_ideen.js`
- Modify: `ansicht/vorlage.html` (Tab „Ideen“ zuerst; `<canvas id="ideen-flaeche" class="karte" hidden>`)
- Modify: `ansicht/js/06_steuerung.js` (Ansicht `ideen`: Start, `setzeAnsicht`, `waehle`, `pfadZeigen`, `filterLeiste`, Zoomknöpfe, Esc, Hash)
- Modify: `ansicht/stil.css` (nur was das Canvas und die Pfadleiste brauchen)

**Interfaces:**
- Consumes: `G.ideen` (Task 2), `K`, `zustand`, `zeigeKnoten`, `tafelSchliessen`, `waehle`, `esc`, `reinText`, `kurz`, `BEWEGUNG`.
- Produces: globales `ideen` mit `zeige()`, `verberge()`, `zuKnoten(id)`, `oeffne(id)` (Haupt- oder Unteridee), `hoch()`, `zoom(f)`, `pfadHtml()`, `anstossen()`.

Verhalten (vollständig im Modul umzusetzen):
- Zustand `offen = { haupt: null|id, unter: null|id }`; Kamera `{x,y,k}` mit Flug wie `fliege()` in `02_karte.js` (Easing kubisch, 820 ms, ohne Bewegung sofort).
- Stufe 0: alle Hauptideen; Verbindungen `kopplung_haupt` als Linien, Breite `1 + 2·√w` px, Deckkraft 0.25–0.6. Kreis: Füllung dunkel getönt, **Füllstand-Ring** aussen (Bogenstücke in Statusfarben `KF.status`, Gruppenfarbe = Farbe des ersten Status der Gruppe; `offen` als gestrichelter Bogen), **Herkunfts-Ring** innen (2 px: Literatur `#93a6d4`, Forum `#e3ab3d`, eigen `#7fe3c4`), Name (600 18px, Serifenschrift wie Karte), darunter „n Aussagen“.
- Stufe 1 (Hauptidee offen): Kamera passt die Hauptidee ein; ihre Unterideen als Kreise mit beiden Ringen und Namen; Linien `kopplung_unter` zwischen Geschwistern; Kopplungen zu Unterideen anderer Hauptideen als Linie zum **Randchip** der fremden Hauptidee. Randchips: andere Hauptideen als Pille am Bildschirmrand in Richtung ihres Zentrums (Name, klickbar → öffnet sie).
- Stufe 2 (Unteridee offen): Kamera passt die Unteridee ein; Aussagen als Punkte: Haupt-Aussagen gross, hell, mit Titel (bis 40 Zeichen), Hilfsaussagen klein, Deckkraft 0.55, Titel nur bei Hover; `stuetzt_sich_auf`-Linien innerhalb der Unteridee dünn; Geschwister-Unterideen gedämpft sichtbar und klickbar.
- Klick auf leere Fläche: nichts; Klick auf Aussage: `waehle(id)` (öffnet Tafel), Auswahl hell umrandet.
- Hover: Tipp mit Name, `kurz` und Füllstand in Worten („12 bewiesen, 3 bedingt, 4 Vermutungen, 2 offen“).
- Mausrad: Zoom um den Zeiger, begrenzt auf 0.6×–4× der Einpass-Grösse der offenen Stufe; Ziehen verschiebt.
- `hoch()`: Stufe 2 → 1 → 0; Esc ruft `hoch()`, wenn die Ideen-Ansicht aktiv ist und keine Eingabe fokussiert ist.
- `zuKnoten(id)`: öffnet Haupt- und Unteridee der Aussage und markiert sie.
- `pfadHtml()`: `Ideen › Haupt › Unter › Aussage` als Knöpfe mit `data-idee="…"` (Ideen = `data-idee="_alle"`); Klick-Behandlung in `06_steuerung.js`.
- `ResizeObserver` auf das Elternelement, Gerätepixel ≤ 2 wie in der Karte.

Steuerung:
- Startansicht ohne Hash: `ideen`; `#v=karte` öffnet weiter die alte Karte; gültige Ansichten `["ideen","karte","theorie","forschung"]`.
- `setzeAnsicht`: in jedem Zweig `ideen.verberge()`, für `ideen` alle anderen verbergen und `ideen.zeige()`.
- `waehle(id)`: in Ansicht `ideen` → `zeigeKnoten(id); ideen.zuKnoten(id)`.
- `filterLeiste`: in `ideen` leer.
- Zoomknöpfe: `ideen.zoom(1.4|1/1.4)`, Einpassen → `ideen.oeffne(null)`.

- [ ] **Step 1:** Modul schreiben, Vorlage und Steuerung anpassen.
- [ ] **Step 2:** `baue.py` laufen lassen (0 Fehler); `node --check` auf das zusammengesetzte JS (aus `ausgabe/tetration-graph.html` den App-Teil extrahieren) — keine Syntaxfehler.
- [ ] **Step 3: Commit** `git commit -m "Ideen-Landkarte: neue Startansicht"`

### Task 4: Stücke einer Vermutung in der Tafel, Prüfung im Browser

**Files:**
- Modify: `ansicht/js/04_tafel.js` (in `zeigeKnoten`: Abschnitt „Woraus es besteht“)
- Modify: `ansicht/stil.css` (Klassen `.stuecke`, `.stueck-voll`, `.stueck-halb`, `.stueck-fehlt`)

- [ ] **Step 1:** In `zeigeKnoten(id)`: wenn `k._status_eff` in {vermutung, bedingt, skizze}, Baum aus `_abhaengig` rekursiv (Tiefe ≤ 3, jede Aussage nur einmal) als verschachtelte Liste; Klasse nach Status der Stücke: bewiesen-Gruppe → `stueck-voll`, bedingt/skizze/numerisch → `stueck-halb` (Deckkraft 0.6), vermutung/annahme/offen → `stueck-fehlt` (hohler Punkt, Zusatz „hier fehlt es“); Einträge sind `data-id`-Knöpfe. Überschrift „Woraus es besteht“ und Zählzeile („7 Stücke: 4 sicher, 2 bedingt, 1 fehlt“).
- [ ] **Step 2: Browser** — `.claude/launch.json` Eintrag „tetration-graph“ (Port 5196) mit `preview_start`; prüfen:
  1. Startansicht zeigt 8 Hauptideen mit Ringen, keine Konsolenfehler.
  2. Klick Hauptidee → Unterideen + Randchips; Klick Unteridee → Aussagen, Haupt/Hilfs unterscheidbar; Klick Aussage → Tafel.
  3. Pfadleiste und Esc führen zurück; Tab „Karte“ zeigt die alte Karte.
  4. `#a=kn.verm.polynom_halbexp_schranke` (oder eine andere bedingte Vermutung) öffnet direkt den Pfad, Tafel zeigt „Woraus es besteht“ mit fehlenden Stücken.
  5. Fenster schmal ziehen (Tablet-Grösse) → nichts verschwindet.
  Bildschirmfotos von Stufe 0, 1, 2 als Beleg.
- [ ] **Step 3:** Fehler beheben, erneut prüfen.
- [ ] **Step 4: Commit und Push** `git commit -m "Ideen-Landkarte: Stücke einer Vermutung, geprüft"; git push origin main`
