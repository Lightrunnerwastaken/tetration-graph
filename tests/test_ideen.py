"""Prüft die Zuordnung von Aussagen zu Ideen der Landkarte."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from test_graph import knoten, regeln  # noqa: E402
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
    assert ("idee", "x.u") in r
    assert any(rr == "idee-doppelt" for rr, _ in r)


def test_ohne_ideendatei_keine_hinweise():
    g = graph(knoten("a"), straenge=[{"id": "g.s", "gebiet": "g", "name": "S", "verlauf": "x", "mitglieder": ["g.a"]}])
    assert not {r for r, _ in regeln(g) if r in ("ohne-idee", "idee", "idee-doppelt")}
