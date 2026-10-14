"""Etappe 3: Zugang zum Motor (über ../etappe2/kern.py) mit eigenem Zwischenspeicher etappe3/cache/."""
import sys
from pathlib import Path
HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER.parent / "etappe2"))
import kern  # noqa: E402
kern.CDIR = HIER / "cache"
from kern import Fix, gpzahl  # noqa: E402,F401
import mpmath as mp  # noqa: E402


def sexp(zs, dps, basis="exp(1)"):
    return kern.werte(basis, [f"sexp({gpzahl(z, dps + 10)})" for z in zs], dps)
