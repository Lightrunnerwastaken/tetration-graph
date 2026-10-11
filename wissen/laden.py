"""Liest die YAML-Daten. Keine Logik ausser Zusammenführen und Herkunft merken."""
import hashlib
from pathlib import Path

import yaml


class _StrengerLader(yaml.SafeLoader):
    """Wie SafeLoader, aber ein doppelter Schlüssel ist ein Fehler (PyYAML überschreibt sonst stillschweigend)."""

    def construct_mapping(self, node, deep=False):
        gesehen = {}
        for k, _ in node.value:
            schluessel = self.construct_object(k, deep=deep)
            if schluessel in gesehen:
                raise ValueError(f"doppelter Schlüssel {schluessel!r} in {node.start_mark.name}, "
                                 f"Zeilen {gesehen[schluessel] + 1} und {k.start_mark.line + 1}")
            gesehen[schluessel] = k.start_mark.line
        return super().construct_mapping(node, deep=deep)


def _yaml(pfad: Path):
    with open(pfad, encoding="utf-8") as f:
        return yaml.load(f, Loader=_StrengerLader)


def lade(wurzel: Path) -> dict:
    daten = wurzel / "daten"
    roh = {
        "gebiete": _yaml(daten / "gebiete.yaml") or [],
        "objekte": _yaml(daten / "objekte.yaml") or {},
        "quellen": _yaml(daten / "quellen.yaml") or {},
        "experimente": _yaml(daten / "experimente.yaml") or [],
        "straenge": (_yaml(daten / "straenge.yaml") or []) if (daten / "straenge.yaml").is_file() else [],
        "aussagen": [],
        "dateien": [],
    }
    h = hashlib.sha256()
    for pfad in sorted((daten / "aussagen").glob("*.yaml")):
        inhalt = _yaml(pfad) or {}
        h.update(pfad.read_bytes())
        gebiet = inhalt.get("gebiet")
        for a in inhalt.get("aussagen") or []:
            a.setdefault("gebiet", gebiet)
            a["_datei"] = pfad.name
            roh["aussagen"].append(a)
        roh["dateien"].append(pfad.name)
    # Lehrbuch: Theorie-Blöcke je Gebiet (id -> Block) und Kapitelfolge
    roh["theorie"] = {}
    roh["theorie_herkunft"] = {}
    for pfad in sorted((daten / "theorie").glob("*.yaml")) if (daten / "theorie").is_dir() else []:
        h.update(pfad.read_bytes())
        for tid, block in (_yaml(pfad) or {}).items():
            if tid in roh["theorie"]:
                roh.setdefault("theorie_doppelt", []).append(tid)
            roh["theorie"][tid] = block
            roh["theorie_herkunft"][tid] = pfad.name
    nach_id = {a.get("id"): a for a in roh["aussagen"]}
    for tid, block in roh["theorie"].items():
        if tid in nach_id:
            nach_id[tid]["theorie"] = block
    roh["lehrbuch"] = (_yaml(daten / "lehrbuch.yaml") or {}) if (daten / "lehrbuch.yaml").is_file() else {}
    bilder = wurzel / "lehrbuch" / "bilder"
    roh["bilder"] = {p.stem: p for p in sorted(bilder.glob("*.svg"))} if bilder.is_dir() else {}
    for name in ("gebiete.yaml", "objekte.yaml", "quellen.yaml", "experimente.yaml", "straenge.yaml", "lehrbuch.yaml"):
        if (daten / name).is_file():
            h.update((daten / name).read_bytes())
    roh["eingabe_hash"] = h.hexdigest()[:16]
    return roh


def datei_hash(pfad: Path) -> str:
    return hashlib.sha256(pfad.read_bytes()).hexdigest()[:16]
