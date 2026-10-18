"""Fasst die E2-Protokolle zu einer Tabelle zusammen: pro Störung kleinstes x mit Verletzung, kleinste verletzende ungerade Ordnung, min. Radius."""
import re, glob, collections
d = collections.OrderedDict()
for f in sorted(glob.glob("daten/log_e2_*.txt")):
    for z in open(f, encoding="utf-8"):
        m = re.match(r"(\S+) (\S+) erste neg. ungerade: (\S+) \( (\d+) \)  rho~ (\S+)  x\+2 = (\S+)", z)
        if m:
            n, x, k, a, r, x2 = m.groups()
            d.setdefault(n, []).append((float(x), None if k == "-" else int(k), int(a), float(r)))
print("| Störung | Gitterpunkte | davon verletzt | kleinstes verletztes x | kleinste verletzende ungerade Ordnung | Wurzeltest-Radius ρ bei x = 0.3 (Kneser-Wert 2.3) |")
print("|---|---|---|---|---|---|")
for n, L in d.items():
    v = [t for t in L if t[1]]
    rr = min(L, key=lambda t: abs(t[0] - 0.3))
    print(f"| {n} | {len(L)} | {len(v)} | {min(t[0] for t in v) if v else '–'} | {min(t[1] for t in v) if v else '–'} | {rr[3]:.4g} (x={rr[0]:.1f}) |")
