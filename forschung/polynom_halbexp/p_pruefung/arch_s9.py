"""Hilfsskript (Strang P): ersetzt ARCHITEKTUR.md §9.3–§9.4 und §9.6 durch den Text aus arch_s9_text.md (Abschnitte durch ===96=== getrennt)."""
import sys
from pathlib import Path

HIER = Path(__file__).resolve().parent
p = HIER.parent / "ARCHITEKTUR.md"
s = p.read_text(encoding="utf-8")
teil = (HIER / "arch_s9_text.md").read_text(encoding="utf-8")
neu34, neu6 = teil.split("===96===\n")
start = s.index("### 9.3 ")
end = s.index("### 9.5 Phase-2-Ergebnisse")
s = s[:start] + neu34 + s[end:]
st6 = s.index("### 9.6 ")
s = s[:st6] + neu6
p.write_text(s, encoding="utf-8")
print("ok", len(s))
