"""M26 — Gesamtlauf für KNESER: Satz M6-K und (a″)_∞ aus m16 (Eckmodell √L), m22 (Endstück), m23 (Kommutator), m24 (Mitte), m25 (C_E').
Mit quelle = "vorzertifikat" ein Vorzertifikat (Daten Strang A, Schwanz angenommen); mit quelle = <Pfad zur zertifizierten Datei von Strang P>
streng. Läuft die Teilskripte nacheinander (je < 20 min) und setzt zusammen:
  c ≥ min(c_mod − ε_L, c_m) − ‖R‖,   σ_min(J_∞|P_N) ≥ c / C_E'.
Aufruf: python m26_kneser_gesamt.py [quelle] [T1 T2]
"""
import re
import subprocess
import sys

sys.stdout.reconfigure(encoding="utf-8")
quelle = sys.argv[1] if len(sys.argv) > 1 else "vorzertifikat"
T1 = sys.argv[2] if len(sys.argv) > 2 else "6"
T2 = sys.argv[3] if len(sys.argv) > 3 else "14"
PY = sys.executable


def lauf(args, muster):
    out = subprocess.run([PY] + args, capture_output=True, text=True, encoding="utf-8").stdout
    print(out.strip().splitlines()[-1] if out.strip() else "(keine Ausgabe)")
    m = re.findall(muster, out)
    if not m:
        raise SystemExit(f"Muster {muster!r} nicht gefunden in {args[0]}")
    return float(m[-1])


c_mod = lauf(["m16_mellin_zert.py", "0.9898", "60", "sqrtL"], r"kleinste Unterschranke \[([0-9.]+)")
eps_L = lauf(["m22_kneser_ende.py", T1, quelle], r"ε_L ≤ ([0-9.]+)")
R = lauf(["m23_kneser_kommutator.py", T1, T2, "0.1", str(int(float(T2)) + 22), quelle], r"‖R‖ ≤ ([0-9.]+)")
c_m = lauf(["m24_kneser_mitte.py", T2, "0.5", "16", "2.8", "96", "8", "5", quelle], r"c_m ≥ ([0-9.\-]+)")
c = min(c_mod - eps_L, c_m) - R
c_s = f"{c - 1e-6:.6f}"
sig = lauf(["m25_kneser_carleson.py", "13", "0.02", c_s, quelle], r"σ_min\(J_∞\|P_N\) ≥ ([0-9.]+)\s*$")
print(f"\nKNESER ({quelle}): c_mod ≥ {c_mod}, ε_L ≤ {eps_L}, ‖R‖ ≤ {R}, c_m ≥ {c_m}  ⇒  Re⟨𝒯x,x⟩ ≥ {c:.4f}‖x‖²;  σ_min(J_∞|P_N) ≥ {sig}")
