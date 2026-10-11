"""Gemeinsame Hilfen für alle Experimente: Ergebnisdatei, Prüfungen, Grafikstil."""
import hashlib
import json
import sys
import time
from pathlib import Path

import matplotlib

from hashwert import experiment_hash

matplotlib.use("Agg")
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass
import matplotlib.pyplot as plt  # noqa: E402

HIER = Path(__file__).resolve().parent
ERGEBNISSE = HIER / "ergebnisse"
BILDER = HIER / "bilder"
DATEN = HIER / "daten"
GP_REFERENZ = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\research\reference")

FARBEN = ["#3d6fb6", "#b5523b", "#2f8f83", "#b0832b", "#8c4fa3", "#6a8f2f", "#a3456e", "#5b6575"]

plt.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "axes.prop_cycle": matplotlib.cycler(color=FARBEN),
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})


class Lauf:
    """Sammelt Prüfungen und Daten eines Experiments und schreibt sie am Ende."""

    def __init__(self, exp_id: str, skript: str):
        self.id = exp_id
        self.skript = Path(skript).resolve()
        self.t0 = time.time()
        self.pruefungen = []
        self.daten = {}
        self.zusammenfassung = ""

    def pruefe(self, name: str, ok: bool, wert=None):
        self.pruefungen.append({"name": name, "ok": bool(ok), "wert": _plain(wert)})
        print(("  ok  " if ok else "  FEHL") + f"  {name}" + (f"  [{_plain(wert)}]" if wert is not None else ""))
        return ok

    def bild(self, fig, name=None):
        BILDER.mkdir(exist_ok=True)
        pfad = BILDER / f"{name or self.id}.svg"
        fig.savefig(pfad, format="svg", bbox_inches="tight")
        plt.close(fig)
        return pfad

    def schreibe(self):
        ERGEBNISSE.mkdir(exist_ok=True)
        out = {
            "id": self.id,
            "skript_hash": experiment_hash(self.skript),
            "gelaufen": time.strftime("%Y-%m-%d %H:%M"),
            "laufzeit_s": round(time.time() - self.t0, 1),
            "python": sys.version.split()[0],
            "zusammenfassung": self.zusammenfassung,
            "pruefungen": self.pruefungen,
            "fehlgeschlagen": [p["name"] for p in self.pruefungen if not p["ok"]],
            "daten": _plain(self.daten),
        }
        (ERGEBNISSE / f"{self.id}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), encoding="utf-8")
        print(f"{self.id}: {len(self.pruefungen) - len(out['fehlgeschlagen'])}/{len(self.pruefungen)} Prüfungen, "
              f"{out['laufzeit_s']} s")


def _plain(x):
    """mpmath/numpy-Werte in JSON-taugliche Werte."""
    try:
        import mpmath
        if isinstance(x, mpmath.mpf):
            return float(x) if abs(x) > mpmath.mpf(10) ** -300 or x == 0 else mpmath.nstr(x, 12)
        if isinstance(x, mpmath.mpc):
            return [_plain(x.real), _plain(x.imag)]
    except ImportError:
        pass
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    if hasattr(x, "item") and callable(x.item):
        try:
            return x.item()
        except Exception:
            return str(x)
    if isinstance(x, complex):
        return [x.real, x.imag]
    return x


def referenzwerte():
    """Kneser-Referenzwerte aus gp-tetration; lokale Kopie mit Herkunft als Rückfall."""
    kopie = DATEN / "kneser_referenz.json"
    quelle = GP_REFERENZ / "values.json"
    if quelle.is_file():
        roh = json.loads(quelle.read_text(encoding="utf-8"))["values"]
        auswahl = {k: (v["real"] if isinstance(v, dict) else v) for k, v in roh.items()
                   if k.split("|")[1] in ("e", "2") and not (isinstance(v, dict) and float(v.get("imag", 0) or 0) != 0)}
        DATEN.mkdir(exist_ok=True)
        kopie.write_text(json.dumps({
            "herkunft": str(quelle),
            "sha256": hashlib.sha256(quelle.read_bytes()).hexdigest(),
            "werte": auswahl}, indent=1), encoding="utf-8")
    return json.loads(kopie.read_text(encoding="utf-8"))


def gp_json(name):
    """phi_modes.json / mu_hub.json aus gp-tetration, lokal gespiegelt."""
    kopie = DATEN / name
    quelle = GP_REFERENZ / name
    if quelle.is_file():
        DATEN.mkdir(exist_ok=True)
        kopie.write_bytes(quelle.read_bytes())
    return json.loads(kopie.read_text(encoding="utf-8"))


# ------------------------------------------------------------------ Kneser-Werte über gp-tetration (mit Zwischenspeicher)
GP_SRC = Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
_CACHE_DATEI = DATEN / "kneser_cache.json"
_cache = None
_motor = None


def _cache_laden():
    global _cache
    if _cache is None:
        _cache = json.loads(_CACHE_DATEI.read_text(encoding="utf-8")) if _CACHE_DATEI.is_file() else {}
    return _cache


def cache_speichern():
    if _cache is not None:
        DATEN.mkdir(exist_ok=True)
        _CACHE_DATEI.write_text(json.dumps(_cache), encoding="utf-8")


def kneser(art, basis, punkte, dps=40):
    """art = 'sexp' oder 'slog'; basis als GP-Ausdruck ('exp(1)', '2'); punkte reell oder komplex (mpmath).

    Werte kommen aus dem Zwischenspeicher oder vom fatou.gp-Fork (gp-tetration). Ohne Motor und ohne
    Speichereintrag: Fehler mit Hinweis.
    """
    import mpmath as mp
    global _motor
    c = _cache_laden()
    schluessel = [f"{art}|{basis}|{dps}|{mp.nstr(mp.re(p), dps)}|{mp.nstr(mp.im(mp.mpc(p)), dps)}" for p in punkte]
    fehlend = [i for i, s in enumerate(schluessel) if s not in c]
    if fehlend:
        if _motor is None:
            if not GP_SRC.is_dir():
                raise RuntimeError("gp-tetration fehlt und Werte sind nicht im Zwischenspeicher")
            sys.path.insert(0, str(GP_SRC))
            from fatou_backend import FatouGP
            _motor = FatouGP(dps=dps, fatou_gp="fork")
        sitzung = _motor.session(basis)
        werte = (sitzung.sexp_batch if art == "sexp" else sitzung.slog_batch)([punkte[i] for i in fehlend])
        for i, w in zip(fehlend, werte):
            c[schluessel[i]] = [mp.nstr(mp.re(w), dps + 5), mp.nstr(mp.im(w), dps + 5)]
        cache_speichern()
    return [mp.mpc(*c[s]) for s in schluessel]


def cauchy_ableitungen(f_batch, x0, r, m_max=4, N=48):
    """Ableitungen f^(m)(x0), m = 0..m_max, aus Werten auf dem Kreis |z − x0| = r (exponentiell genau)."""
    import mpmath as mp
    zs = [x0 + r * mp.expjpi(2 * mp.mpf(j) / N) for j in range(N)]
    w = f_batch(zs)
    out = []
    for m in range(m_max + 1):
        s = mp.fsum(w[j] * mp.expjpi(-2 * m * mp.mpf(j) / N) for j in range(N)) / N
        out.append(mp.re(s * mp.factorial(m) / r ** m))
    return out


def gauss_legendre(f_batch, a, b, n=40):
    """∫_a^b f mit n-Punkt-Gauss–Legendre (f_batch bekommt alle Knoten auf einmal)."""
    import mpmath as mp
    x, w = _gl_knoten(n)
    h, m = (b - a) / 2, (a + b) / 2
    werte = f_batch([m + h * xi for xi in x])
    return h * mp.fsum(wi * mp.re(v) for wi, v in zip(w, werte))


_gl = {}


def _gl_knoten(n):
    import mpmath as mp
    key = (n, mp.mp.dps)
    if key not in _gl:
        x, w = [], []
        for i in range(1, n + 1):
            z = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
            for _ in range(100):
                p0, p1 = mp.mpf(1), z
                for k in range(2, n + 1):
                    p0, p1 = p1, ((2 * k - 1) * z * p1 - (k - 1) * p0) / k
                dp = n * (z * p1 - p0) / (z * z - 1)
                dz = p1 / dp
                z -= dz
                if abs(dz) < mp.mpf(10) ** (-mp.mp.dps + 5):
                    break
            x.append(z)
            w.append(2 / ((1 - z * z) * dp * dp))
        _gl[key] = (x, w)
    return _gl[key]
