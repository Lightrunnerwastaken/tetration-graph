"""Prüft die Regeln des Prüfprogramms an kleinen künstlichen Graphen."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from wissen.graph import Graph  # noqa: E402

VOLL = {
    "pruefung": {k: "x" for k in ("f1_folgt_unmittelbar", "f2_benutzte_voraussetzungen",
                                  "f3_spezialfall", "f4_uebertragbar", "f5_neue_frage")},
    "bewertung": {"richtigkeit": "x", "neuheit": "klassisch", "nutzen": "x"},
    "beweis": {"art": "vollstaendig", "text": "x"},
}


def knoten(i, art="satz", status="bewiesen", **kw):
    d = {"id": f"g.{i}", "titel": i, "art": art, "ebene": "A", "status": status,
         "aussage": "x", "gebiet": "g", **VOLL}
    if status == "skizze":
        d["beweis"] = {"art": "skizze", "text": "x", "luecken": ["l"]}
    if status == "bedingt":
        d.setdefault("bedingungen", [])
    if art == "vermutung":
        d["gegenbeispielsuche"] = "x"
    d.update(kw)
    return d


def graph(*aussagen, straenge=None):
    roh = {"gebiete": [{"id": "g", "name": "G"}], "objekte": {"o": {"name": "o"}},
           "quellen": {"q": {"geprueft": "nein"}}, "experimente": [], "aussagen": list(aussagen),
           "straenge": straenge or []}
    return Graph(roh).alles()


def regeln(g, i=None):
    return {(b.regel, b.id) for b in g.befunde if i is None or b.id == i}


def test_bewiesen_auf_vermutung_wird_bedingt():
    g = graph(knoten("v", art="vermutung", status="vermutung"),
              knoten("s", stuetzt_sich_auf=["g.v"]))
    assert g.knoten["g.s"]["_status_eff"] == "bedingt"
    assert ("status-zu-stark", "g.s") in regeln(g)
    assert g.knoten["g.s"]["_bedingungen_auto"] == ["g.v"]


def test_bedingt_pflanzt_sich_fort():
    g = graph(knoten("a", art="annahme", status="annahme"),
              knoten("s", stuetzt_sich_auf=["g.a"]),
              knoten("t", stuetzt_sich_auf=["g.s"]))
    assert g.knoten["g.t"]["_status_eff"] == "bedingt"


def test_entlastete_annahme_traegt_beweis():
    g = graph(knoten("a", art="annahme", status="annahme"),
              knoten("p", kanten=[{"typ": "entlastet", "ziel": "g.a"}]),
              knoten("s", stuetzt_sich_auf=["g.a"]))
    assert g.knoten["g.a"]["_status_eff"] == "entlastet"
    assert g.knoten["g.s"]["_status_eff"] == "bewiesen"


def test_widerlegt_markiert_abhaengige_transitiv():
    g = graph(knoten("b", art="beobachtung", status="widerlegt", experimente=[], quellen=["q"]),
              knoten("v", art="vermutung", status="vermutung", stuetzt_sich_auf=["g.b"]),
              knoten("s", status="bedingt", bedingungen=["g.v"], stuetzt_sich_auf=["g.v"]))
    assert g.knoten["g.v"]["_neu_pruefen"] == [["g.b"]]
    assert g.knoten["g.s"]["_neu_pruefen"] == [["g.v", "g.b"]]


def test_nachgeprueft_hebt_markierung_auf():
    g = graph(knoten("b", art="beobachtung", status="widerlegt", quellen=["q"]),
              knoten("v", art="vermutung", status="vermutung", stuetzt_sich_auf=["g.b"],
                     nachgeprueft=[{"gegen": "g.b", "datum": "2026-10-08", "ergebnis": "ok"}]))
    assert g.knoten["g.v"]["_neu_pruefen"] == []


def test_kreis_wird_gemeldet():
    g = graph(knoten("a", stuetzt_sich_auf=["g.b"]), knoten("b", stuetzt_sich_auf=["g.a"]))
    assert any(b.regel == "kreis" for b in g.befunde)


def test_unbenutzte_voraussetzung_ist_hinweis():
    g = graph(knoten("s", voraussetzungen=[{"text": "f stetig", "benutzt": "nein"}]))
    assert ("abschwaechbar", "g.s") in regeln(g)


def test_verweise_und_vokabular():
    g = graph(knoten("s", status="fast", objekte=["fehlt"], stuetzt_sich_auf=["g.nix"]))
    r = regeln(g, "g.s")
    assert ("vokabular", "g.s") in r and ("verweis", "g.s") in r


def test_vermutung_ohne_gegenbeispielsuche():
    k = knoten("v", art="vermutung", status="vermutung")
    del k["gegenbeispielsuche"]
    assert ("gegenbeispiel", "g.v") in regeln(graph(k))


def test_wurzelbedingungen_durch_bedingte_kette():
    g = graph(knoten("a", art="annahme", status="annahme"),
              knoten("s", status="bedingt", bedingungen=["g.a"], stuetzt_sich_auf=["g.a"]),
              knoten("t", status="bedingt", bedingungen=["g.s"], stuetzt_sich_auf=["g.s"]))
    assert g.knoten["g.t"]["_wurzelbedingungen"] == ["g.a"]
    assert not any(b.regel == "bedingung-implizit" for b in g.befunde)


def test_ungenannte_wurzel_wird_gemeldet():
    g = graph(knoten("a", art="annahme", status="annahme"),
              knoten("v", art="vermutung", status="vermutung"),
              knoten("s", status="bedingt", bedingungen=["g.a"], stuetzt_sich_auf=["g.a", "g.v"]))
    assert ("bedingung-implizit", "g.s") in regeln(g)


def test_tragweite_zaehlt_transitiv():
    g = graph(knoten("a"), knoten("b", stuetzt_sich_auf=["g.a"]), knoten("c", stuetzt_sich_auf=["g.b"]))
    assert g.knoten["g.a"]["_tragweite"] == 2


def test_verwandte_nach_kanten():
    g = graph(knoten("a"), knoten("b", stuetzt_sich_auf=["g.a"]), knoten("c"))
    assert g.knoten["g.a"]["_verwandt"][0]["id"] == "g.b"


def test_straenge_zuordnung_und_doppelte():
    s = [{"id": "g.s1", "gebiet": "g", "name": "S1", "verlauf": "x", "mitglieder": ["g.a", "g.b"]},
         {"id": "g.s2", "gebiet": "g", "name": "S2", "verlauf": "x", "mitglieder": ["g.b", "g.nix"]}]
    g = graph(knoten("a"), knoten("b"), knoten("c"), straenge=s)
    r = regeln(g)
    assert g.knoten["g.a"]["_strang"] == "g.s1"
    assert ("strang", "g.b") in r and ("strang", "g.s2") in r and ("ohne-strang", "g.c") in r


def graph_mit_buch(aussagen, theorie=None, lehrbuch=None, straenge=None):
    roh = {"gebiete": [{"id": "g", "name": "G"}], "objekte": {"o": {"name": "o"}},
           "quellen": {"q": {"geprueft": "nein"}}, "experimente": [], "aussagen": list(aussagen),
           "straenge": straenge or [], "theorie": theorie or {}, "theorie_herkunft": {}, "lehrbuch": lehrbuch or {}, "bilder": {}}
    nach_id = {a["id"]: a for a in roh["aussagen"]}
    for tid, block in roh["theorie"].items():
        if tid in nach_id:
            nach_id[tid]["theorie"] = block
    return Graph(roh).alles()


def stufen(g, regel):
    return {(b.stufe, b.id) for b in g.befunde if b.regel == regel}


def test_theorie_toter_verweis_und_unbekannter_teil():
    block = {"worum": "siehe [[g.gibtsnicht]]", "genau": "x", "warum": "x", "einordnung": "x", "fazit": "x"}
    g = graph_mit_buch([knoten("a")], theorie={"g.a": block})
    texte = [b.text for b in g.befunde if b.regel == "theorie" and b.stufe == "FEHLER"]
    assert any("g.gibtsnicht" in t for t in texte)
    assert any("fazit" in t for t in texte)


def test_theorie_unvollstaendig_und_fremde_id():
    g = graph_mit_buch([knoten("a")], theorie={"g.a": {"worum": "x"}, "g.fremd": {"worum": "x"}})
    assert ("WARNUNG", "g.a") in stufen(g, "theorie")
    assert ("FEHLER", "g.fremd") in stufen(g, "theorie")


def test_theorie_verweis_mit_grossbuchstaben():
    block = {"worum": "[[g.Gibt_L]] und [[g.Fehlt_L]]", "genau": "x", "warum": "x", "einordnung": "x"}
    g = graph_mit_buch([knoten("a"), knoten("Gibt_L")], theorie={"g.a": block})
    texte = [b.text for b in g.befunde if b.regel == "theorie" and b.stufe == "FEHLER"]
    assert any("g.Fehlt_L" in t for t in texte)
    assert not any("g.Gibt_L" in t for t in texte)


def test_lehrbuch_strang_und_doppelte_aussage():
    s = [{"id": "g.s", "gebiet": "g", "name": "S", "mitglieder": ["g.a"]}]
    buch = {"kapitel": [{"id": "k", "titel": "K", "teile": [{"strang": "g.s"}, {"knoten": "g.a"}, {"strang": "g.fehlt"}]}]}
    g = graph_mit_buch([knoten("a")], lehrbuch=buch, straenge=s)
    assert ("FEHLER", "k") in stufen(g, "lehrbuch")
    assert ("WARNUNG", "g.a") in stufen(g, "lehrbuch")


def test_doppelter_yaml_schluessel_ist_fehler(tmp_path):
    from wissen.laden import _yaml
    p = tmp_path / "q.yaml"
    p.write_text("a:\n  x: 1\nb:\n  x: 2\na:\n  x: 3\n", encoding="utf-8")
    try:
        _yaml(p)
    except ValueError as e:
        assert "doppelter Schlüssel 'a'" in str(e)
    else:
        raise AssertionError("doppelter Schlüssel wurde nicht erkannt")
