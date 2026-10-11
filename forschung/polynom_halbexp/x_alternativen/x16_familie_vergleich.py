"""X16 — Familien-Stabilität numerisch: Tangente J_∞(p_N) (echte Lösung, gemeinsam/pN_N100_bis_1000.json) gegen J_∞(f_K).
Vergleicht R, R̂, Rc auf Ps ∂K'-Gitter (x12/x13-Dateien gleicher N, d): max relative Abweichung, und Θ beider (x12b-Logik)."""
import json, sys
from pathlib import Path
import numpy as np
HIER = Path(__file__).resolve().parent
sys.stdout.reconfigure(encoding="utf-8")
for arg in sys.argv[1:]:
    a, b = arg.split("|")
    A = json.loads((HIER / a).read_text())["punkte"]; B = json.loads((HIER / b).read_text())["punkte"]
    out = []
    for key in ("R", "Rhat", "Rc"):
        rel = [abs(x[key] - y[key]) / max(abs(x[key]), 1e-300) for x, y in zip(A, B)]
        i = int(np.argmax(rel))
        out.append(f"{key}: max rel. Abw. {max(rel):.2e} (bei {A[i]['rand']} {A[i]['par']:.4g}), Median {np.median(rel):.1e}")
    print(f"{a} ↔ {b}: " + "; ".join(out))
