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


def _teile(daten: Path, name: str, h):
    """Hauptdatei daten/<name>.yaml (falls vorhanden) plus alle Ergänzungsdateien daten/<name>/*.yaml, in dieser Reihenfolge."""
    pfade = ([daten / f"{name}.yaml"] if (daten / f"{name}.yaml").is_file() else []) + \
            (sorted((daten / name).glob("*.yaml")) if (daten / name).is_dir() else [])
    for pfad in pfade:
        h.update(pfad.read_bytes())
        yield pfad, _yaml(pfad)


def _woerterbuch(daten: Path, name: str, h) -> dict:
    """Zusammenführen ohne stilles Überschreiben: ein Schlüssel in zwei Dateien ist ein Fehler."""
    alles, herkunft = {}, {}
    for pfad, inhalt in _teile(daten, name, h):
        for k, v in (inhalt or {}).items():
            if k in alles:
                raise ValueError(f"{name}: Schlüssel {k!r} doppelt, in {herkunft[k]} und {pfad.name}")
            alles[k], herkunft[k] = v, pfad.name
    return alles


def _liste(daten: Path, name: str, h) -> list:
    return [e for _, inhalt in _teile(daten, name, h) for e in (inhalt or [])]


def lade(wurzel: Path) -> dict:
    daten = wurzel / "daten"
    h = hashlib.sha256()
    roh = {
        "gebiete": _yaml(daten / "gebiete.yaml") or [],
        "objekte": _woerterbuch(daten, "objekte", h),
        "quellen": _woerterbuch(daten, "quellen", h),
        "experimente": _liste(daten, "experimente", h),
        "straenge": _liste(daten, "straenge", h),
        "ideen": _liste(daten, "ideen", h),
        "uebersetzung": _liste(daten, "uebersetzung", h),
        "aussagen": [],
        "dateien": [],
    }
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
    for name in ("gebiete.yaml", "lehrbuch.yaml"):
        if (daten / name).is_file():
            h.update((daten / name).read_bytes())
    roh["eingabe_hash"] = h.hexdigest()[:16]
    return roh


def datei_hash(pfad: Path) -> str:
    return hashlib.sha256(pfad.read_bytes()).hexdigest()[:16]
