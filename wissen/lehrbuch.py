"""Prüfungen für das Lehrbuch: Theorie-Blöcke (daten/theorie/*.yaml) und Kapitelfolge (daten/lehrbuch.yaml)."""
import re
from collections import Counter

from . import schema as S

VERWEIS = re.compile(r"\[\[([a-z]{1,3}\.[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*)(?:\|[^\]]*)?\]\]")
BILD = re.compile(r"!\[\[bild:([A-Za-z0-9_.-]+)\]\]")


def _texte(x):
    if isinstance(x, str):
        yield x
    elif isinstance(x, dict):
        for v in x.values():
            yield from _texte(v)
    elif isinstance(x, list):
        for v in x:
            yield from _texte(v)


def _bild_ok(g, name):
    return name in g.roh.get("bilder", {}) or name in g.experimente


def pruefen(g):
    """Hängt Befunde an g an. g ist ein Graph nach aufbauen()."""
    roh = g.roh
    for tid in roh.get("theorie_doppelt", []):
        g.befund("FEHLER", "theorie", tid, "Theorie-Block doppelt")
    for tid, block in (roh.get("theorie") or {}).items():
        if tid not in g.knoten:
            g.befund("FEHLER", "theorie", tid, f"Theorie-Block für unbekannte Aussage (Datei {roh['theorie_herkunft'].get(tid)})")
            continue
        if not isinstance(block, dict):
            g.befund("FEHLER", "theorie", tid, "Theorie-Block muss ein Wörterbuch sein")
            continue
        for teil in block:
            if teil not in S.THEORIE_TEILE:
                g.befund("FEHLER", "theorie", tid, f"unbekannter Teil '{teil}' (erlaubt: {', '.join(S.THEORIE_TEILE)})")
        fehlt = [t for t in S.THEORIE_PFLICHT if not block.get(t)]
        if fehlt:
            g.befund("WARNUNG", "theorie", tid, "Theorie-Block unvollständig: " + ", ".join(fehlt))
        if g.knoten[tid].get("_niveau", g.knoten[tid].get("niveau")) == "einstieg" and not block.get("beispiel"):
            g.befund("HINWEIS", "theorie", tid, "Einstiegsknoten ohne Beispiel")
        for name in block.get("bilder") or []:
            if not _bild_ok(g, name):
                g.befund("FEHLER", "theorie", tid, f"Bild '{name}' gibt es nicht (lehrbuch/bilder oder Experiment-id)")
        for t in _texte({k: v for k, v in block.items() if k != "bilder"}):
            for ziel in VERWEIS.findall(t):
                if ziel not in g.knoten and ziel not in g.experimente:
                    g.befund("FEHLER", "theorie", tid, f"Verweis [[{ziel}]] zeigt ins Leere")
            for name in BILD.findall(t):
                if not _bild_ok(g, name):
                    g.befund("FEHLER", "theorie", tid, f"Bild '{name}' gibt es nicht")

    buch = roh.get("lehrbuch") or {}
    straenge = {s["id"]: s for s in roh.get("straenge") or []}
    gesehen_str, gesehen_kn = Counter(), Counter()
    for kap in buch.get("kapitel") or []:
        kid = kap.get("id", "?")
        for teil in kap.get("teile") or []:
            if "strang" in teil:
                if teil["strang"] not in straenge:
                    g.befund("FEHLER", "lehrbuch", kid, f"unbekannter Strang '{teil['strang']}'")
                gesehen_str[teil["strang"]] += 1
                for m in straenge.get(teil["strang"], {}).get("mitglieder", []):
                    gesehen_kn[m] += 1
            elif "knoten" in teil:
                if teil["knoten"] not in g.knoten:
                    g.befund("FEHLER", "lehrbuch", kid, f"unbekannte Aussage '{teil['knoten']}'")
                gesehen_kn[teil["knoten"]] += 1
            elif "bild" in teil:
                if not _bild_ok(g, teil["bild"]):
                    g.befund("FEHLER", "lehrbuch", kid, f"Bild '{teil['bild']}' gibt es nicht")
            elif "text" in teil:
                for ziel in VERWEIS.findall(teil["text"]):
                    if ziel not in g.knoten and ziel not in g.experimente:
                        g.befund("FEHLER", "lehrbuch", kid, f"Verweis [[{ziel}]] zeigt ins Leere")
                for name in BILD.findall(teil["text"]):
                    if not _bild_ok(g, name):
                        g.befund("FEHLER", "lehrbuch", kid, f"Bild '{name}' gibt es nicht")
            elif "abschnitt" not in teil:
                g.befund("FEHLER", "lehrbuch", kid, f"Teil ohne text/strang/knoten/bild/abschnitt: {list(teil)}")
    if buch:
        for sid in straenge:
            if not gesehen_str[sid]:
                g.befund("HINWEIS", "lehrbuch", sid, "Strang kommt im Lehrbuch nicht vor")
        for i, n in gesehen_kn.items():
            if n > 1:
                g.befund("WARNUNG", "lehrbuch", i, f"Aussage steht {n}-mal im Lehrbuch")
        for i in g.knoten:
            if not gesehen_kn[i]:
                g.befund("HINWEIS", "lehrbuch", i, "Aussage steht in keinem Kapitel")
    # Fortschritt je Gebiet als ein Hinweis
    for geb in g.roh.get("gebiete") or []:
        ids = [i for i, a in g.knoten.items() if a.get("gebiet") == geb["id"]]
        mit = sum(1 for i in ids if g.knoten[i].get("theorie"))
        if ids and mit < len(ids):
            g.befund("HINWEIS", "theorie-fortschritt", geb["id"], f"{mit} von {len(ids)} Aussagen mit Theorie-Block")
    return g


def fortschritt(g):
    return {"mit_theorie": sum(1 for a in g.knoten.values() if a.get("theorie")), "gesamt": len(g.knoten)}
