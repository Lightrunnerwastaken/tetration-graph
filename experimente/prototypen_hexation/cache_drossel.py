"""Nur für die Prototypen: der gemeinsame Kneser-Zwischenspeicher wird nur gelesen; neue Werte gehen in kneser_cache_proto.json.

gemeinsam.cache_speichern schreibt sonst die ganze Datei (inzwischen ~190 MB) nach jedem Motoraufruf mit neuen Punkten; bei der
Blattverfolgung sind das hunderte Aufrufe. Hier: nur die neuen Einträge, höchstens alle 5 Minuten und beim Beenden.
gemeinsam.py selbst bleibt unverändert.
"""
import atexit
import json
import time
from pathlib import Path

import gemeinsam

_PROTO = Path(__file__).with_name("kneser_cache_proto.json")
_zuletzt = [time.time()]


def _speichern_jetzt():
    c = gemeinsam._cache
    neu = {k: v for k, v in c.items() if k not in _basis}
    _PROTO.write_text(json.dumps(neu), encoding="utf-8")


def _gedrosselt():
    if time.time() - _zuletzt[0] > 300:
        _speichern_jetzt()
        _zuletzt[0] = time.time()


if getattr(gemeinsam.cache_speichern, "__name__", "") != "_gedrosselt":
    _c = gemeinsam._cache_laden()
    _basis = set(_c.keys())
    if _PROTO.is_file():
        _c.update(json.loads(_PROTO.read_text(encoding="utf-8")))
    gemeinsam.cache_speichern = _gedrosselt
    atexit.register(_speichern_jetzt)
