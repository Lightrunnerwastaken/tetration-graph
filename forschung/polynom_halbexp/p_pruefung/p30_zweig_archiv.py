"""Strang P, Prüfung 30 (Meilenstein 6): Zweig-Zertifikat p24 NOCHMALS laufen lassen und dabei die Boxen archivieren.

Identische Mathematik wie p24_zweig_kohaerent.py (Funktion `uebergang` wird importiert, nichts geändert); nur
  * eigener Zustand (p30_zustand.json), Start wie p24 bei der p18-Box von p_45,
  * jede Box für N ≥ ARCHIV_AB wird als p30_boxen/box_N.json abgelegt (für (A): p_N gegen Π_N f_K, N = 140 … 193),
  * Fortsetzung über 193 hinaus möglich (nur auf aicgolling).
Aufruf: python p30_zweig_archiv.py NMAX [Zeitlimit_s] [Threads]
"""
import json
import sys
import time
from pathlib import Path

HIER = Path(__file__).resolve().parent
sys.path.insert(0, str(HIER))
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 193
ZEITLIMIT = float(sys.argv[2]) if len(sys.argv) > 2 else 1e9
THREADS = int(sys.argv[3]) if len(sys.argv) > 3 else 1
ARCHIV_AB = 100

import p24_zweig_kohaerent as z  # noqa: E402
from flint import ctx  # noqa: E402

z.ZUSTAND = HIER / "p30_zustand.json"
ARCHIV = HIER / "p30_boxen"
ARCHIV.mkdir(exist_ok=True)
sys.stdout.reconfigure(encoding="utf-8")


def schreibe(N, box, prot):
    z.ZUSTAND.write_text(json.dumps({"N": N, "box": [z.kodiere(t) for t in box], "protokoll": prot}, indent=1),
                         encoding="utf-8")
    if N >= ARCHIV_AB:
        (ARCHIV / f"box_{N}.json").write_text(json.dumps({"N": N, "box": [z.kodiere(t) for t in box],
                                                          "protokoll": prot[-1] if prot else None}),
                                              encoding="utf-8")


def main():
    try:
        ctx.threads = THREADS
    except Exception as e:  # ältere python-flint
        print("Threads nicht setzbar:", e)
    N, box, prot = z.lade_start()
    t_start = time.time()
    print(f"Start bei N = {N} (Box aus {'p30' if prot else 'p18'}), Threads {THREADS}", flush=True)
    while N < NMAX and time.time() - t_start < ZEITLIMIT:
        ok, neu, p = z.uebergang(box, N)
        print(json.dumps(p, ensure_ascii=False), flush=True)
        if not ok:
            print(f"ABBRUCH bei {N} → {N + 1}", flush=True)
            break
        prot.append(p)
        N, box = N + 1, neu
        schreibe(N, box, prot)
    print(f"Ende: N = {N}", flush=True)


if __name__ == "__main__":
    main()
