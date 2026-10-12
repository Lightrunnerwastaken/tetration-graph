"""Mechanical invariants between a German source and its English translation (pilot version of the future check rules)."""
import json
import re
import sys
from collections import Counter

MATH = re.compile(r"\$\$.+?\$\$|\$.+?\$", re.S)
REF = re.compile(r"\[\[.+?\]\]")
IDENT = re.compile(r"\b(?:gr|fg|it|fp|kl|kn|bw|ar|br|exp)\.[a-z0-9_]+(?:\.[A-Za-z0-9_]+)*")
NUM = re.compile(r"(?<![\w.])\d+(?:\.\d+)?(?:\s?[×·]\s?10\^?[-−]?\d+)?")
# (german marker, english regex that must appear if the marker appears)
SIGNAL = [
    (r"\bnicht\b|\bkein", r"\bnot\b|\bno\b|\bnever\b|\bneither\b|n't|\bnon-|\bwithout\b|\bnor\b|cannot|can't|\bimpossible\b|\bnone\b|\bnowhere\b|\bun[a-z]+"),
    (r"\bnur\b", r"\bonly\b|\bmerely\b|\bsolely\b"),
    (r"\bgenau\b", r"\bexactly\b|\bprecisely\b|\bjust\b|\bif and only if\b|\bexact\b"),
    (r"\bjede[rsmn]?\b|\balle\b", r"\bevery\b|\bany\b|\ball\b|\beach\b"),
    (r"\bvermutlich\b|\bwohl\b|\bvielleicht\b", r"\bpresumably\b|\bprobably\b|\blikely\b|\bperhaps\b|\bsuppos|\bmay\b|\bmight\b|\bapparently\b"),
    (r"\bbedingt", r"\bconditional"),
    (r"\bVermutung", r"\bconjecture"),
]


def mathe(s):
    """Formulas as a multiset; the prose inside \\text{...} may be translated, so it is blanked."""
    return Counter(re.sub(r"\\text\{[^}]*\}", r"\\text{}", m) for m in MATH.findall(s))


def texte(o, pfad=""):
    """Yield (path, string) for every string leaf."""
    if isinstance(o, str):
        yield pfad, o
    elif isinstance(o, dict):
        for k, v in o.items():
            yield from texte(v, f"{pfad}.{k}" if pfad else k)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from texte(v, f"{pfad}[{i}]")


def pruefe_eintrag(quelle, ueb, name):
    fehler = []
    qt = dict(texte(quelle))
    ut = dict(texte(ueb))
    if set(qt) != set(ut):
        fehler.append(f"{name}: key structure differs: missing {sorted(set(qt)-set(ut))}, extra {sorted(set(ut)-set(qt))}")
    for pfad, q in qt.items():
        u = ut.get(pfad)
        if u is None:
            continue
        last = pfad.split(".")[-1]
        if last in ("id", "art", "status") or pfad.endswith("tex") or last == "art":
            if q != u:
                fehler.append(f"{name}:{pfad}: must be unchanged")
            continue
        if mathe(q) != mathe(u):
            a, b = mathe(q), mathe(u)
            fehler.append(f"{name}:{pfad}: formulas differ: lost {list((a-b).elements())[:3]} added {list((b-a).elements())[:3]}")
        if Counter(REF.findall(q)) != Counter(REF.findall(u)):
            fehler.append(f"{name}:{pfad}: [[refs]] differ")
        if Counter(IDENT.findall(q)) != Counter(IDENT.findall(u)):
            fehler.append(f"{name}:{pfad}: identifiers differ")
        qn = Counter(NUM.findall(MATH.sub(" ", q)))
        un = Counter(NUM.findall(MATH.sub(" ", u)))
        if qn != un:
            fehler.append(f"{name}:{pfad}: numbers outside formulas differ: lost {dict(qn-un)} added {dict(un-qn)}")
        if q.count("\n\n") != u.count("\n\n"):
            fehler.append(f"{name}:{pfad}: paragraph structure differs ({q.count(chr(10)*2)} vs {u.count(chr(10)*2)} blank lines)")
        for de, en in SIGNAL:
            if re.search(de, q, re.I) and not re.search(en, u, re.I):
                fehler.append(f"{name}:{pfad}: signal word '{de}' in German but nothing matching in English  [WARN]")
        if re.search(r"neighbour|behaviour|colour|analys(e|ing)\b|normalis|favour|centre\b|\bgrey\b|organis|recognis|honour", u, re.I):
            fehler.append(f"{name}:{pfad}: British spelling  [WARN]")
        if re.search(r"[äöüÄÖÜß]| (und|der|die|das|ist|wird|mit|für|nicht|auch) ", MATH.sub(" ", u).replace("van der", "van_der")):
            fehler.append(f"{name}:{pfad}: German left in translation  [WARN]")
    return fehler


def main(quell, ueb):
    Q = json.load(open(quell, encoding="utf-8"))
    U = json.load(open(ueb, encoding="utf-8"))
    if isinstance(Q, list):
        Qd = {e["id"]: e for e in Q}
        Ud = {e["id"]: e for e in U}
    else:
        Qd, Ud = Q, U
    alle = []
    for i, q in Qd.items():
        if i not in Ud:
            alle.append(f"{i}: missing")
            continue
        alle += pruefe_eintrag(q, Ud[i], i)
    hard = [f for f in alle if "[WARN]" not in f]
    warn = [f for f in alle if "[WARN]" in f]
    print(f"{len(Qd)} entries; {len(hard)} errors, {len(warn)} warnings")
    for f in hard + warn:
        print(" ", f[:300])


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
