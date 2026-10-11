"""Fügt einen BRETT-Eintrag (Text aus Datei argv[1]) oben ein, direkt vor dem ersten '## ['-Eintrag."""
import os
import sys

HIER = os.path.dirname(os.path.abspath(__file__))
B = os.path.join(HIER, "..", "BRETT.md")
neu = open(sys.argv[1], encoding="utf-8").read().rstrip() + "\n\n"
alt = open(B, encoding="utf-8").read()
i = alt.index("\n## [") + 1
open(B, "w", encoding="utf-8").write(alt[:i] + neu + alt[i:])
print("eingefügt an Position", i)
