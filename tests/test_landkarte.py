"""Prüft die Rechnungen der Ideen-Landkarte an einem Mini-Graphen."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from test_ideen import graph, knoten  # noqa: E402
from wissen import landkarte as L  # noqa: E402

IDEEN = [{"id": "h1", "name": "H1", "kurz": "x"}, {"id": "h2", "name": "H2", "kurz": "x"},
         {"id": "h1.a", "eltern": "h1", "name": "A", "kurz": "x", "aussagen": ["g.a1", "g.a2"]},
         {"id": "h1.b", "eltern": "h1", "name": "B", "kurz": "x", "aussagen": ["g.b1"]},
         {"id": "h2.c", "eltern": "h2", "name": "C", "kurz": "x", "aussagen": ["g.c1"]},
         {"id": "h2.leer", "eltern": "h2", "name": "Leer", "kurz": "x"}]


def mini():
    a1 = knoten("a1", objekte=["o"])
    a2 = knoten("a2", stuetzt_sich_auf=["g.a1"])
    b1 = knoten("b1", stuetzt_sich_auf=["g.a1"], objekte=["o"])
    c1 = knoten("c1", art="vermutung", status="vermutung", kanten=[{"typ": "verwandt", "ziel": "g.a2"}])
    return graph(a1, a2, b1, c1, ideen=IDEEN)


def test_kopplung_von_hand():
    g = mini()
    w = L.kopplung(g, {k: a["_idee"] for k, a in g.knoten.items()})
    assert abs(w[("h1.a", "h1.b")] - 1.5) < 1e-9          # Kante 1.0 + Objekte 0.5
    assert abs(w[("h1.a", "h2.c")] - 0.25) < 1e-9         # verwandt
    assert ("h1.b", "h2.c") not in w


def test_lage_deterministisch_ohne_ueberlappung():
    g = mini()
    roh = {"ideen": IDEEN, "quellen": {"q": {"art": "literatur"}}}
    d1, d2 = L.berechne(g, roh), L.berechne(g, roh)
    assert d1 == d2
    H = d1["haupt"]
    for i in range(len(H)):
        for j in range(i + 1, len(H)):
            assert ((H[i]["x"] - H[j]["x"]) ** 2 + (H[i]["y"] - H[j]["y"]) ** 2) ** 0.5 >= H[i]["r"] + H[j]["r"] - 1e-6
    hp = {h["id"]: h for h in H}
    for u in d1["unter"]:
        p = hp[u["eltern"]]
        assert ((u["x"] - p["x"]) ** 2 + (u["y"] - p["y"]) ** 2) ** 0.5 + u["r"] <= p["r"] + 0.5
        assert u["r"] > 0
    assert next(u for u in d1["unter"] if u["id"] == "h2.leer")["n"] == 0


def test_hauptaussagen_hoechstens_ein_drittel_mindestens_eine():
    d = L.berechne(mini(), {"ideen": IDEEN, "quellen": {}})
    assert d["aussagen"]["g.a1"]["haupt"] and not d["aussagen"]["g.a2"]["haupt"]


def test_fuellstand_und_offen():
    d = L.berechne(mini(), {"ideen": IDEEN, "quellen": {}})
    c = next(u for u in d["unter"] if u["id"] == "h2.c")
    assert c["fuell"] == {"vermutung": 1} and c["offen"] == ["g.c1"]
    t = L.text_md(d)
    assert "Leer" in t and "g.c1" in t


def test_ohne_idee_landet_in_sammelidee():
    g = mini()
    g.knoten["g.c1"]["_idee"] = None
    d = L.berechne(g, {"ideen": IDEEN, "quellen": {}})
    assert "g.c1" in d["aussagen"]
    assert any(u["id"] == "_ohne" for u in d["unter"])
