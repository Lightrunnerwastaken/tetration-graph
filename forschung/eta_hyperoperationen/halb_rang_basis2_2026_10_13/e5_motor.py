"""Etappe 5: eigene Kopie des Hyperoperations-Motors der Vorrunde (nur gelesen, nicht veraendert):
engine_snapshot.gp + hyper.gp + even8.gp -> e5_engine.gp in diesem Ordner."""
import sys, time
from pathlib import Path
HIER = Path(__file__).resolve().parent
VOR = HIER.parent
sys.path.insert(0, str(VOR / '_backend'))
from fatou_backend import FatouGP, find_default_gp_exe

def motor(dps=30, timeout=600):
    f = HIER / 'e5_engine.gp'
    inhalt = (VOR / 'engine_snapshot.gp').read_bytes() + b'\n' + (VOR / 'hyper.gp').read_bytes() + b'\n' + (VOR / 'even8.gp').read_bytes()
    if not f.exists() or f.read_bytes() != inhalt:
        f.write_bytes(inhalt)
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=f, dps=dps, state_cache=False,
                   init_timeout=200, eval_timeout=timeout, n_workers=1)
