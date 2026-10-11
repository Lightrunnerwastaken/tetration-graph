pfad = "D:/Users/janis/Documents/Matura/tetration-graph/forschung/polynom_halbexp/BRETT.md"
t = open(pfad, encoding="utf-8").read()
neu = open("brett1.txt", encoding="utf-8").read().strip() + "\n\n"
kaputt = "Neue Einträge **oben** anfügen: `" + neu
assert kaputt in t
t = t.replace(kaputt, "Neue Einträge **oben** anfügen: `", 1)
kopfende = "`→ Strang Y:` kennzeichnen.\n\n"
i = t.index(kopfende) + len(kopfende)
t = t[:i] + neu + t[i:]
open(pfad, "w", encoding="utf-8").write(t)
print("ok")
