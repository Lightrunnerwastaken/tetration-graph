"""Hals-Analyse: kleinster Abstand zwischen weit entfernten Parameterstücken der Kurve (aus daten/fein_*.json)."""
import json, sys
from pathlib import Path
H = Path(__file__).resolve().parent
for f in sys.argv[1:]:
    d = json.loads((H / "daten" / f).read_text())
    X = d["x"]; P = [complex(a, b) for a, b in zip(d["re"], d["im"])]
    best = []
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            dx = abs(X[i] - X[j]); dx = min(dx, 1 - dx)
            if dx < 0.1: continue
            best.append((abs(P[i] - P[j]), X[i], X[j]))
    best.sort()
    print(f, "kleinste Abstände (|x-y|>=0.1):", [(f"{a:.2e}", round(b, 4), round(c, 4)) for a, b, c in best[:4]])
    print("  Kurve:", " ".join(f"{X[i]:.2f}:{P[i].real:.3f}{P[i].imag:+.3f}i" for i in range(0, len(P), max(1, len(P)//24))))
