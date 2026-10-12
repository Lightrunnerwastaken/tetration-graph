"""Prüft die Übersetzung der Strukturtexte (wissen/uebersetzung_struktur.py)."""
import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from test_graph import knoten, regeln  # noqa: E402
from wissen import uebersetzung_struktur as U  # noqa: E402
from wissen.graph import Graph  # noqa: E402
from wissen.uebersetzung import flach, hole, quell_hash  # noqa: E402

KAPITEL = {"id": "k9", "titel": "Titel", "untertitel": "Unter", "teile": [
    {"abschnitt": "Abschnitt", "anker": "a"}, {"text": "Text $x^2$ siehe [[g.a]]"}, {"strang": "g.s1"},
    {"bild": "b", "unterschrift": "Bild"}]}


def graph(uebersetzung_struktur, **roh_extra):
    roh = {"gebiete": [{"id": "g", "name": "Grundlagen", "kurz": "Einstieg"}], "objekte": {"o": {"name": "Objekt", "tex": "x"}},
           "quellen": {"q": {"geprueft": "nein"}}, "experimente": [], "aussagen": [knoten("a")],
           "straenge": [{"id": "g.s1", "gebiet": "g", "name": "Werte", "verlauf": "Was a↑↑n ist", "mitglieder": ["g.a"]}],
           "ideen": [], "lehrbuch": {"kapitel": [copy.deepcopy(KAPITEL)]}, "uebersetzung_struktur": uebersetzung_struktur}
    roh.update(roh_extra)
    return Graph(roh).alles()


def eintrag(art, i, item, en, **kw):
    e = {"art": art, "id": i, "en": en, "geprueft": "nein", **kw}
    e["quelle_hash"] = quell_hash({p: hole(item, p) for p in flach(en)})
    return e


def test_pfade_und_vorlage_kapitel():
    assert U.pfade("kapitel", KAPITEL) == ["titel", "untertitel", "teile[0].abschnitt", "teile[1].text", "teile[3].unterschrift"]
    v = U.vorlage("kapitel", KAPITEL)
    assert v["titel"] == "Titel" and len(v["teile"]) == 4 and v["teile"][2] == {} and v["teile"][1] == {"text": KAPITEL["teile"][1]["text"]}
    assert U.vorlage("strang", {"name": "N", "verlauf": "V", "mitglieder": []}) == {"name": "N", "verlauf": "V"}
    assert U.pfade("gebiet", {"id": "x", "name": "N", "kurz": " "}) == ["name"]          # leere Texte zählen nicht


def test_aktuell_und_veraltet():
    s = {"id": "g.s1", "name": "Werte", "verlauf": "Was a↑↑n ist"}
    e = eintrag("strang", "g.s1", s, {"name": "Values", "verlauf": "What a↑↑n is"})
    g = graph([e])
    st = g.straenge[0]
    assert st["en"]["name"] == "Values" and st["_en_stand"] == "aktuell"
    g2 = graph([e], straenge=[{**s, "gebiet": "g", "name": "Werte (neu)", "mitglieder": ["g.a"]}])
    assert g2.straenge[0]["_en_stand"] == "veraltet"
    assert ("uebersetzung-veraltet", "strang:g.s1") in regeln(g2)


def test_fehler_in_eintraegen():
    s = {"id": "g.s1", "name": "Werte", "verlauf": "Was a↑↑n ist"}
    assert ("uebersetzung-struktur", "strang:g.nix") in regeln(graph([{"art": "strang", "id": "g.nix", "en": {"name": "x"}}]))
    assert ("uebersetzung-struktur", "quatsch:x") in regeln(graph([{"art": "quatsch", "id": "x", "en": {"name": "x"}}]))
    assert ("uebersetzung-struktur", "strang:g.s1") in regeln(graph([eintrag("strang", "g.s1", s, {"name": "V", "nonsense": "x"})]))
    assert ("uebersetzung-struktur", "strang:g.s1") in regeln(graph([eintrag("strang", "g.s1", s, {"name": "V"}), eintrag("strang", "g.s1", s, {"name": "V"})]))
    assert ("uebersetzung-struktur", "strang:g.s1") in regeln(graph([eintrag("strang", "g.s1", s, {"name": "V"}, geprueft="vielleicht")]))


def test_kapitel_teile_gleiche_laenge_und_invarianten():
    ok = {"titel": "Title", "teile": [{"abschnitt": "Section"}, {"text": "Text $x^2$ see [[g.a]]"}, {}, {"unterschrift": "Figure"}]}
    g = graph([eintrag("kapitel", "k9", KAPITEL, ok)])
    assert g.roh["lehrbuch"]["kapitel"][0]["en"]["titel"] == "Title"
    kurz = {"titel": "Title", "teile": [{}, {}]}
    assert ("uebersetzung-struktur", "kapitel:k9") in regeln(graph([eintrag("kapitel", "k9", KAPITEL, kurz)]))
    kaputt = {"teile": [{}, {"text": "Text $x^3$ see [[g.a]]"}, {}, {}]}
    g = graph([eintrag("kapitel", "k9", KAPITEL, kaputt)])
    assert ("uebersetzung-invariante", "kapitel:k9") in regeln(g)
    assert "en" not in g.roh["lehrbuch"]["kapitel"][0]                      # fehlerhafte Übersetzung wird nicht angehängt


def test_fortschritt_nur_wenn_uebersetzt():
    assert not [b for b in graph([]).befunde if b.regel == "uebersetzung-fortschritt" and b.id in U.ARTEN]
    s = {"id": "g.s1", "name": "Werte", "verlauf": "Was a↑↑n ist"}
    g = graph([eintrag("strang", "g.s1", s, {"name": "Values"})])
    assert any(b.regel == "uebersetzung-fortschritt" and b.id == "gebiet" for b in g.befunde)
