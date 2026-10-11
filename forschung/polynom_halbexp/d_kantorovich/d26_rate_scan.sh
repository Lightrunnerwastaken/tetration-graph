#!/bin/bash
# Asymptotische Eckmoden-Rate im kubischen Modell: je Konfiguration N = 400 ... 800
cd "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/d_kantorovich"
for cfg in "np.sqrt(L)|76.618" "1.25*np.exp(0.669j)|76.618" "1.35*np.exp(0.669j)|76.618" "np.sqrt(L)|105" "np.sqrt(L)|60"; do
  mu="${cfg%%|*}"; ph="${cfg##*|}"
  echo "== mu=$mu phi=$ph"
  C:/Users/janis/miniconda3/python.exe d25_kubik_modell.py 400,500,600,700,800 "$mu" "$ph" 2>&1 | grep "N="
done
echo FERTIG
