"""Fügt den Text aus einer Datei als neuen Eintrag oben in BRETT.md ein (nach dem Kopf)."""
import sys

p = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/BRETT.md"
t = open(p, encoding="utf-8").read()
e = open(sys.argv[1], encoding="utf-8").read().rstrip("\n") + "\n\n"
marker = "kennzeichnen.\n"
i = t.index(marker) + len(marker)
rest = t[i:].lstrip("\n")
t = t[:i] + "\n" + e + rest
open(p, "w", encoding="utf-8").write(t)
print("ok")
