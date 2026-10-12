"""Prüft die Regeln der englischen Fassung (wissen/uebersetzung.py)."""
import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from test_graph import knoten, regeln  # noqa: E402
from wissen import uebersetzung as U  # noqa: E402
from wissen.graph import Graph  # noqa: E402


def graph(aussagen, uebersetzung):
    roh = {"gebiete": [{"id": "g", "name": "G"}], "objekte": {"o": {"name": "o"}},
           "quellen": {"q": {"geprueft": "nein"}}, "experimente": [], "aussagen": copy.deepcopy(aussagen),
           "straenge": [], "uebersetzung": uebersetzung}
    return Graph(roh).alles()


def eintrag(a, en, **kw):
    e = {"id": a["id"], "en": en, "geprueft": "nein", **kw}
    e["quelle_hash"] = U.quell_hash(U.deutsch(a, U.flach(en)))
    return e


def stufen(befunde, regel=None):
    return {(s, r) for s, r, _ in befunde if regel is None or r == regel}


# ------------------------------------------------------------------ Textpaare
def test_formeln_muessen_gleich_bleiben():
    assert stufen(U.pruefe_text("Sei $x^2$.", "Let $x^3$.")) == {("FEHLER", "uebersetzung-invariante")}
    assert not U.pruefe_text("Sei $x^2$ und $$y$$.", "Let $x^2$ and $$y$$.")


def test_prosa_in_text_darf_uebersetzt_werden():
    assert not U.pruefe_text(r"$$a\qquad\text{mit}\quad b$$", r"$$a\qquad\text{with}\quad b$$")
    assert stufen(U.pruefe_text(r"$$a\text{mit}b$$", r"$$a\text{with}c$$")) == {("FEHLER", "uebersetzung-invariante")}


def test_verweise_und_kennungen():
    assert stufen(U.pruefe_text("Siehe [[kn.konstr.kneser]].", "See [[kn.konstr.kneser2]].")) >= {("FEHLER", "uebersetzung-invariante")}
    assert stufen(U.pruefe_text("Siehe exp.funde.", "See exp.funden.")) >= {("FEHLER", "uebersetzung-invariante")}
    assert not U.pruefe_text("Siehe [[kn.konstr.kneser]] und exp.funde.", "See [[kn.konstr.kneser]] and exp.funde.")


def test_zahlen_ausserhalb_von_formeln_warnen():
    assert stufen(U.pruefe_text("Zeit-1-Abbildung", "time-one map")) == {("WARNUNG", "uebersetzung-zahlen")}
    assert not U.pruefe_text("Zeit-1-Abbildung", "time-1 map")
    assert not U.pruefe_text("etwa 1,5 mal", "about 1.5 times")                   # Dezimalkomma -> Punkt ist gewollt
    assert stufen(U.pruefe_text("etwa 1,5 mal", "about 1,5 times")) == {("WARNUNG", "uebersetzung-zahlen")}


def test_signalwoerter_britisch_und_deutsch():
    assert ("WARNUNG", "uebersetzung-wortwahl") in stufen(U.pruefe_text("Das gilt nur für reelle b.", "This holds for real b."))
    assert not stufen(U.pruefe_text("Das gilt nur für reelle b.", "This holds only for real b."), "uebersetzung-wortwahl")
    assert ("WARNUNG", "uebersetzung-wortwahl") in stufen(U.pruefe_text("Nachbarn", "neighbouring points"))
    assert ("WARNUNG", "uebersetzung-wortwahl") in stufen(U.pruefe_text("Beispiel", "Das ist ein Beispiel"))
    assert not stufen(U.pruefe_text("nach Müller–Schleicher", "after Müller–Schleicher"), "uebersetzung-wortwahl")
    assert not stufen(U.pruefe_text("Verallgemeinert ar.satz.keine_arith_info", "Generalizes ar.satz.keine_arith_info"), "uebersetzung-wortwahl")
    assert not stufen(U.pruefe_text("Sie sind nicht reell, wohl aber periodisch.", "They are not real, but they are periodic."), "uebersetzung-wortwahl")
    assert ("WARNUNG", "uebersetzung-wortwahl") in stufen(U.pruefe_text("Das gilt wohl.", "This holds."))
    assert not stufen(U.pruefe_text("Hardy und van der Hoeven", "Hardy and van der Hoeven"), "uebersetzung-wortwahl")


def test_absaetze():
    assert ("WARNUNG", "uebersetzung-absaetze") in stufen(U.pruefe_text("a\n\nb", "a b"))


# ------------------------------------------------------------------ Gestalt des en-Blocks
def test_form_des_en_blocks():
    a = knoten("a", formen=[{"tex": "x", "notiz": "n"}, {"tex": "y"}], einfach="e")
    ok = {"titel": "A", "formen": [{"notiz": "n"}, {}], "beweis": {"text": "p"}}
    assert not U.pruefe_form(a, ok)
    assert U.pruefe_form(a, {"nonsense": "x"})
    assert U.pruefe_form(a, {"formen": [{"notiz": "n"}]})                   # falsche Länge
    assert U.pruefe_form(a, {"formen": [{"tex": "z"}, {}]})                 # tex nie übersetzen
    assert U.pruefe_form(a, {"formen": [{}, {"notiz": "n2"}]})              # im Deutschen keine Notiz
    assert U.pruefe_form(a, {"naechste_frage": "x"})                        # im Deutschen nicht vorhanden
    assert U.pruefe_form(a, {"titel": " "})                                 # leer
    assert U.pruefe_form(a, {"beweis": {"art": "vollstaendig"}})            # art ist kein Text


# ------------------------------------------------------------------ Hash und Veralten im Graphen
def test_aktuell_und_veraltet():
    a = knoten("a", titel="Satz")
    e = eintrag(a, {"titel": "Theorem"})
    g = graph([dict(a)], [e])
    assert g.knoten["g.a"]["_en_stand"] == "aktuell"
    assert g.knoten["g.a"]["en"]["titel"] == "Theorem"
    assert not [b for b in g.befunde if b.regel.startswith("uebersetzung") and b.stufe != "HINWEIS"]
    a2 = dict(a, titel="Satz geändert")
    g = graph([a2], [e])
    assert g.knoten["g.a"]["_en_stand"] == "veraltet"
    assert ("uebersetzung-veraltet", "g.a") in regeln(g)


def test_nicht_uebersetzte_felder_machen_nicht_veraltet():
    a = knoten("a", titel="Satz", aussage="alt")
    e = eintrag(a, {"titel": "Theorem"})
    g = graph([dict(a, aussage="neu")], [e])     # aussage ist nicht übersetzt, also nichts veraltet
    assert g.knoten["g.a"]["_en_stand"] == "aktuell"


def test_ohne_hash_ungestempelt():
    a = knoten("a", titel="Satz")
    e = eintrag(a, {"titel": "Theorem"})
    del e["quelle_hash"]
    g = graph([a], [e])
    assert g.knoten["g.a"]["_en_stand"] == "ungestempelt"


def test_fehlerhafte_eintraege():
    a, b = knoten("a"), knoten("b")
    assert ("uebersetzung", "g.fehlt") in regeln(graph([a], [{"id": "g.fehlt", "en": {"titel": "x"}}]))
    e = eintrag(a, {"titel": "T"})
    assert ("uebersetzung", "g.a") in regeln(graph([a], [e, dict(e)]))                       # doppelt
    assert ("uebersetzung", "g.a") in regeln(graph([dict(a, en={"titel": "T"})], [e]))       # inline und Datei
    assert ("uebersetzung", "g.a") in regeln(graph([a], [{"id": "g.a"}]))                    # leer
    assert ("uebersetzung", "g.a") in regeln(graph([a], [dict(e, geprueft="vielleicht")]))
    kaputt = eintrag(a, {"aussage": "$x^2$"})
    kaputt["en"]["aussage"] = "$x^3$"
    g = graph([dict(a, aussage="$x^2$")], [kaputt])
    assert ("uebersetzung-invariante", "g.a") in regeln(g)
    assert any(x.stufe == "FEHLER" for x in g.befunde if x.regel == "uebersetzung-invariante")


def test_fortschritt_nur_wenn_uebersetzt_wird():
    a, b = knoten("a"), knoten("b")
    assert not [x for x in graph([a, b], []).befunde if x.regel == "uebersetzung-fortschritt"]
    g = graph([a, b], [eintrag(a, {"titel": "T"})])
    f = [x for x in g.befunde if x.regel == "uebersetzung-fortschritt"]
    assert len(f) == 1 and "1 von 2 Aussagen übersetzt" in f[0].text


def test_theorie_uebersetzung():
    a = knoten("a", theorie={"worum": "Hallo", "genau": "Welt"})
    e = {"id": "g.a", "en_theorie": {"worum": "Hello"}, "theorie_hash": U.quell_hash(U.deutsch_theorie(a, ["worum"]))}
    g = graph([a], [e])
    assert g.knoten["g.a"]["en_theorie"]["worum"] == "Hello"
    assert g.knoten["g.a"]["_en_stand"] == "aktuell"
    a2 = knoten("a", theorie={"worum": "Hallo geändert", "genau": "Welt"})
    assert ("uebersetzung-veraltet", "g.a") in regeln(graph([a2], [e]))
    assert ("uebersetzung", "g.a") in regeln(graph([a], [dict(e, en_theorie={"bilder": "x"})]))


# ------------------------------------------------------------------ Hilfsfunktionen für Export/Einlesen
def test_pfade_und_bauen_sind_umkehrbar():
    en = {"titel": "T", "beweis": {"text": "B"}, "formen": [{"notiz": "n1"}, {}, {"notiz": "n3"}], "pruefung": {"f5_neue_frage": "F"}}
    assert U._bauen(U.flach(en), None) == en
    a = {"formen": [{"tex": "x", "notiz": "n1"}, {"tex": "y"}, {"tex": "z", "notiz": "n3"}], "titel": "t"}
    assert U._blatt_pfade(a, "formen") == ["formen[0].notiz", "formen[2].notiz"]
    assert U.hole(a, "formen[2].notiz") == "n3" and U.hole(a, "formen[1].notiz") is None and U.hole(a, "x.y") is None


def test_hash_reihenfolge_und_aenderung():
    assert U.quell_hash({"a": "1", "b": "2"}) == U.quell_hash({"b": "2", "a": "1"})
    assert U.quell_hash({"a": "1"}) != U.quell_hash({"a": "2"})


def test_verschmelzen():
    alt = {"titel": "T", "formen": [{"notiz": "n1"}, {}], "beweis": {"text": "B"}}
    neu = {"aussage": "A", "formen": [{}, {"notiz": "n2"}], "beweis": {"text": "B2"}}
    m = U._verschmelzen(alt, neu)
    assert m["titel"] == "T" and m["aussage"] == "A" and m["beweis"]["text"] == "B2"
    assert m["formen"] == [{"notiz": "n1"}, {"notiz": "n2"}]


def test_abgleich_zeigt_nur_abweichungen():
    a = [{"id": "x", "en": {"titel": "Measurable consequence of the table", "aussage": "Same text $x^2$ here for sure"}}]
    b = [{"id": "x", "en": {"titel": "A testable sequence of the table", "aussage": "Same text $y$ here for sure"}}]
    abw = U.abgleich(a, b, schwelle=0.85)
    assert [x[1] for x in abw] == ["en.titel"]                 # Formeln zählen nicht, Wortwahl schon
    assert "consequence" in abw[0][3] and "sequence" in abw[0][4]
    assert U.abgleich(a, a) == []
    assert U.abgleich(a, [])[0][1] == "(fehlt in B)"


def test_geprueft_je_gruppe():
    e = {"en": {"titel": "T", "beweis": {"text": "B"}}, "en_theorie": {"worum": "W"}, "geprueft": "zweitpruefung"}
    assert U.geprueft_je_gruppe(e) == {"beweis": "zweitpruefung", "kern": "zweitpruefung", "theorie": "zweitpruefung"}
    e["geprueft"] = {"kern": "zweitpruefung", "theorie": "mensch"}
    assert U.geprueft_je_gruppe(e) == {"beweis": "nein", "kern": "zweitpruefung", "theorie": "mensch"}
    assert U.schlechtester_stand(e) == "nein"
    assert U.schlechtester_stand({"en": {"titel": "T"}, "geprueft": "mensch"}) == "mensch"
    a = knoten("a")
    assert not [x for x in U.pruefe_eintrag(a, {"id": "g.a", "en": {"titel": "T"}, "geprueft": {"kern": "mensch"}})
                if x[0] == "FEHLER"]
    assert [x for x in U.pruefe_eintrag(a, {"id": "g.a", "en": {"titel": "T"}, "geprueft": {"quatsch": "mensch"}}) if x[0] == "FEHLER"]
    assert [x for x in U.pruefe_eintrag(a, {"id": "g.a", "en": {"titel": "T"}, "geprueft": {"kern": "vielleicht"}}) if x[0] == "FEHLER"]
