"""Fügt einen Eintrag oben ins BRETT ein (direkt nach dem Kopf). Aufruf: python brett_eintrag.py datei.txt"""
import sys
pfad = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/BRETT.md"
neu = open(sys.argv[1], encoding="utf-8").read().strip() + "\n\n"
alt = open(pfad, encoding="utf-8").read()
kopfende = "`→ Strang Y:` kennzeichnen.\n\n"
i = alt.index(kopfende) + len(kopfende)
open(pfad, "w", encoding="utf-8").write(alt[:i] + neu + alt[i:])
print("eingefuegt")
