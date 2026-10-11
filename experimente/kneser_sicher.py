"""Kneser-Werte mit Zeitlimit: Aufruf in einem Unterprozess; hängt der Motor, wird der Stapel halbiert, bis der Punkt isoliert ist.

Erfolgreiche Werte landen über gemeinsam.kneser im Zwischenspeicher; Punkte ohne Ergebnis kommen als None zurück.
Aufruf als Skript (intern): python kneser_sicher.py <eingabe.json>
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import mpmath as mp

HIER = Path(__file__).parent


def _als_text(p):
    p = mp.mpc(p)
    return [mp.nstr(p.real, 45), mp.nstr(p.imag, 45)]


def _cache_neu():
    """Zwischenspeicher frisch von der Platte lesen (der Unterprozess hat ihn erweitert) und im Modul ersetzen."""
    import gemeinsam
    gemeinsam._cache = json.loads(gemeinsam._CACHE_DATEI.read_text(encoding="utf-8")) if gemeinsam._CACHE_DATEI.is_file() else {}
    return gemeinsam._cache


def kneser_sicher(art, basis, punkte, dps=40, zeit=45):
    c = _cache_neu()
    schl = [f"{art}|{basis}|{dps}|{mp.nstr(mp.re(p), dps)}|{mp.nstr(mp.im(mp.mpc(p)), dps)}" for p in punkte]
    erg = [mp.mpc(*c[s]) if s in c else None for s in schl]
    fehlend = [i for i, w in enumerate(erg) if w is None]

    def rechne(idx):
        if not idx:
            return
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as f:
            json.dump({"art": art, "basis": basis, "dps": dps, "punkte": [_als_text(punkte[i]) for i in idx]}, f)
            name = f.name
        proz = subprocess.Popen([sys.executable, str(HIER / "kneser_sicher.py"), name], cwd=HIER,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        try:
            proz.wait(timeout=zeit)
            ok = True
        except subprocess.TimeoutExpired:
            ok = False
            # ganzen Prozessbaum beenden — sonst bleiben die gp-Arbeiter des Motors als Waisen liegen
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proz.pid)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        Path(name).unlink(missing_ok=True)
        if ok:
            c2 = _cache_neu()
            for i in idx:
                if schl[i] in c2:
                    erg[i] = mp.mpc(*c2[schl[i]])
            return
        if len(idx) > 1:
            h = len(idx) // 2
            rechne(idx[:h])
            rechne(idx[h:])

    rechne(fehlend)
    return erg


if __name__ == "__main__":
    sys.path.insert(0, str(HIER))
    from gemeinsam import kneser
    daten = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    mp.mp.dps = daten["dps"]
    pkt = [mp.mpc(mp.mpf(r), mp.mpf(i)) if mp.mpf(i) != 0 else mp.mpf(r) for r, i in daten["punkte"]]
    kneser(daten["art"], daten["basis"], pkt, dps=daten["dps"])
