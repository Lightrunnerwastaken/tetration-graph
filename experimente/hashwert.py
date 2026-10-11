"""Hash eines Experiments: Skript plus alle lokal importierten Hilfsmodule (rekursiv).

Ohne Abhängigkeiten, damit baue.py denselben Wert berechnen kann wie die Experimente.
"""
import hashlib
import re
from pathlib import Path

_IMPORT = re.compile(r"^\s*(?:from\s+([A-Za-z_]\w*)\s+import|import\s+([A-Za-z_]\w*))", re.M)


def lokale_module(skript: Path, gesehen=None) -> list[Path]:
    gesehen = gesehen if gesehen is not None else set()
    out = []
    for m in _IMPORT.finditer(skript.read_text(encoding="utf-8")):
        name = m.group(1) or m.group(2)
        pfad = skript.parent / f"{name}.py"
        if pfad.is_file() and pfad not in gesehen:
            gesehen.add(pfad)
            out.append(pfad)
            out += lokale_module(pfad, gesehen)
    return out


def experiment_hash(skript) -> str:
    skript = Path(skript)
    h = hashlib.sha256(skript.read_bytes())
    for p in sorted(lokale_module(skript), key=lambda p: p.name):
        h.update(p.name.encode())
        h.update(p.read_bytes())
    return h.hexdigest()[:16]
