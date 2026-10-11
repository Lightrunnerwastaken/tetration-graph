"""Baut den Wissensgraphen: prüfen, ableiten, Bericht + graph.json + eine HTML-Datei.

    python baue.py                 prüfen und bauen
    python baue.py --experimente   vorher alle Experimente neu rechnen
    python baue.py --nur-pruefen   nur Prüfbericht
Rückgabewert 1, wenn es FEHLER-Befunde gibt.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

WURZEL = Path(__file__).resolve().parent
sys.path.insert(0, str(WURZEL))

from wissen.ausgabe import bericht_md, html_bauen  # noqa: E402
from wissen.graph import Graph  # noqa: E402
from wissen.laden import lade  # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass


def experimente_rechnen(roh):
    for e in roh["experimente"]:
        skript = WURZEL / e["skript"]
        print(f"→ {e['id']} ({skript.name})")
        r = subprocess.run([sys.executable, skript.name], cwd=skript.parent, capture_output=True, text=True,
                           encoding="utf-8", env={**__import__("os").environ, "PYTHONIOENCODING": "utf-8"})
        letzte = (r.stdout.strip().splitlines() or ["(keine Ausgabe)"])[-1]
        print("   " + (letzte if r.returncode == 0 else "FEHLER: " + r.stderr.strip().splitlines()[-1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--experimente", action="store_true")
    ap.add_argument("--nur-pruefen", action="store_true")
    args = ap.parse_args()

    roh = lade(WURZEL)
    if args.experimente:
        experimente_rechnen(roh)
    g = Graph(roh, WURZEL).alles()
    stat = g.statistik()
    befunde = g.befunde_dict()
    aus = WURZEL / "ausgabe"
    aus.mkdir(exist_ok=True)
    (aus / "pruefbericht.md").write_text(bericht_md(g, stat, befunde), encoding="utf-8")
    if not args.nur_pruefen:
        daten = g_export(g, stat, befunde, roh)
        (aus / "graph.json").write_text(json.dumps(daten, ensure_ascii=False, indent=1), encoding="utf-8")
        html = html_bauen(WURZEL, daten)
        (aus / "tetration-graph.html").write_text(html, encoding="utf-8")
    zahl = {s: sum(1 for b in befunde if b["stufe"] == s) for s in ("FEHLER", "WARNUNG", "HINWEIS")}
    print(f"{len(g.knoten)} Aussagen, {len(g.kanten)} Kanten, {len(g.experimente)} Experimente — "
          f"{zahl['FEHLER']} Fehler, {zahl['WARNUNG']} Warnungen, {zahl['HINWEIS']} Hinweise")
    for b in befunde:
        if b["stufe"] != "HINWEIS":
            print(f"  {b['stufe']:8} {b['regel']:18} {b['id']}: {b['text']}")
    print(f"Bericht: ausgabe/pruefbericht.md" + ("" if args.nur_pruefen else "  ·  Ansicht: ausgabe/tetration-graph.html"))
    return 1 if zahl["FEHLER"] else 0


def _svg_schlank(svg):
    """XML-Kopf und Kommentare weg, damit das SVG direkt ins HTML passt."""
    import re
    svg = re.sub(r"<\?xml[^>]*>\s*", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg)
    return re.sub(r"<!--.*?-->", "", svg, flags=re.S)


def g_export(g, stat, befunde, roh):
    knoten = []
    for i, a in g.knoten.items():
        k = {key: v for key, v in a.items() if not key.startswith("_") or key in (
            "_status_eff", "_bedingungen_auto", "_neu_pruefen", "_tragweite", "_abhaengig", "_traegt", "_verwandt",
            "_teilweise_entlastet_durch", "_datei", "_wurzelbedingungen", "_strang", "_niveau")}
        k["befunde"] = [b for b in befunde if b["id"] == i]
        knoten.append(k)
    exp = []
    for e in g.experimente.values():
        d = {k: v for k, v in e.items() if not k.startswith("_")}
        d["zustand"] = e.get("_zustand")
        erg = e.get("_ergebnis") or {}
        d["ergebnis"] = {k: erg.get(k) for k in ("zusammenfassung", "pruefungen", "fehlgeschlagen", "gelaufen", "laufzeit_s")}
        d["svg"] = e.get("_svg")
        exp.append(d)
    from wissen import schema as S
    return {
        "gebaut": time.strftime("%Y-%m-%d %H:%M"),
        "eingabe_hash": roh["eingabe_hash"],
        "gebiete": roh["gebiete"],
        "objekte": roh["objekte"],
        "quellen": roh["quellen"],
        "knoten": knoten,
        "kanten": [k for k in g.kanten if k["nach"] in g.knoten],
        "experimente": exp,
        "straenge": roh.get("straenge") or [],
        "lehrbuch": roh.get("lehrbuch") or {},
        "bilder": {n: _svg_schlank(p.read_text(encoding="utf-8")) for n, p in (roh.get("bilder") or {}).items()},
        "befunde": befunde,
        "statistik": stat,
        "schema": {"arten": S.ARTEN, "ebenen": S.EBENEN, "status": S.STATUS, "kanten": S.KANTEN,
                   "neuheit": S.NEUHEIT, "prueffragen": S.PRUEFFRAGEN, "quellen_pruefung": S.QUELLEN_PRUEFUNG,
                   "niveau": S.NIVEAU, "theorie_teile": S.THEORIE_TEILE},
    }


if __name__ == "__main__":
    sys.exit(main())
