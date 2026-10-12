"""Erzeugt GESAMTVERZEICHNIS.md aus den Quellen des Graphen: jede Quelle mit Art, Prüfstand und den Gebieten/Aussagen, die sie zitieren.

Aufruf (aus dem Repo-Wurzelordner): C:/Users/janis/miniconda3/python.exe forschung/literaturkarte/verzeichnis.py
Kandidaten, die noch nicht im Graphen sind, stehen von Hand in KANDIDATEN.md und werden unten angehängt.
"""
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

WURZEL = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WURZEL))
from wissen.laden import lade  # noqa: E402

roh = lade(WURZEL)
Q = roh["quellen"]
zitiert = defaultdict(list)
for a in roh["aussagen"]:
    for q in a.get("quellen") or []:
        if isinstance(q, dict) and q.get("ref"):
            zitiert[q["ref"]].append(a["id"])

ART = {"literatur": "Literatur", "forum": "Forum", "eigene": "eigene Arbeit", "notiz": "Notiz/Wiki", "software": "Software"}


def jahr(q):
    m = re.search(r"(?<![\d.:])(1[6-9]\d\d|20[0-2]\d)(?![\d.])", f"{q.get('ort', '')} {q.get('titel', '')}")
    return m.group(1) if m else "—"


def zeile(k, q):
    ids = zitiert.get(k, [])
    geb = Counter(i.split(".")[0] for i in ids)
    gtxt = ", ".join(f"{g} ({n})" for g, n in sorted(geb.items())) or "— (nicht zitiert)"
    titel = str(q.get("titel", "")).replace("|", "/")
    autor = str(q.get("autor", "")).replace("|", "/")
    return f"| `{k}` | {autor} | {jahr(q)} | {titel} | {q.get('geprueft', '—')} | {gtxt} |"


teile = ["# Gesamtverzeichnis der Quellen (automatisch erzeugt)", "",
         "Erzeugt von `forschung/literaturkarte/verzeichnis.py` aus `daten/quellen.yaml` und `daten/quellen/*.yaml`. Spalte „zitiert in“: Gebiet und Zahl der Aussagen.", ""]
stat = Counter(q.get("art") for q in Q.values())
teile.append("**Stand:** " + ", ".join(f"{ART.get(a, a)} {n}" for a, n in stat.most_common()) + f" — zusammen {len(Q)} Quellen.")
teile.append("")
for art in ("literatur", "forum", "notiz", "software", "eigene"):
    eintraege = sorted(((k, q) for k, q in Q.items() if q.get("art") == art), key=lambda kq: (jahr(kq[1]), kq[0]))
    if not eintraege:
        continue
    teile += [f"## {ART[art]} ({len(eintraege)})", "", "| Schlüssel | Autor | Jahr | Titel | geprüft | zitiert in |", "|---|---|---|---|---|---|"]
    teile += [zeile(k, q) for k, q in eintraege]
    teile.append("")
kand = Path(__file__).with_name("KANDIDATEN.md")
if kand.is_file():
    teile += ["", kand.read_text(encoding="utf-8")]
Path(__file__).with_name("GESAMTVERZEICHNIS.md").write_text("\n".join(teile) + "\n", encoding="utf-8")
print(f"{len(Q)} Quellen geschrieben; nicht zitiert: {sum(1 for k in Q if k not in zitiert)}")
