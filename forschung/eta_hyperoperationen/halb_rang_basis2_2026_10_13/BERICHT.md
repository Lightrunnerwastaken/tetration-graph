# Halbrang H₃.₅(2,x) an Basis 2 — Etappen 1 und 2

13.10.2026 · Herleitungen in [MATHEMATIK.md](MATHEMATIK.md) · Reproduktion in [README.md](README.md)

## Ergebnis zuerst

**Etappe 1.** Der Halbrang-Jet H₃.₅(2,x)≈x+1+Σ_j u_j(x) für die echte Störung d=2ˣ−x−1 ist bis **Ordnung 7** berechnet, an reellen und komplexen Höhen, mit einer ausdrücklich gewählten Auswahl (perturbative Höhenklasse + Newtonfortsetzung im Rang). Die Glieder fallen von ~6·10⁻² auf 10⁻⁵–10⁻⁴ bei Ordnung 5–7. Damit:

| x | 2ˣ | **H₃.₅(2,x)** | Kneser 2↑↑x | H₄.₅(2,x) |
|---:|---:|---:|---:|---:|
| 0.25 | 1.189207 | **1.21246** | 1.224361 | 1.23115 |
| 0.5 | 1.414214 | **1.44302** | 1.458782 | 1.46827 |
| 0.75 | 1.681793 | **1.70164** | 1.713180 | 1.72118 |
| 1.25 | 2.378414 | **2.35338** | 2.336520 | 2.32438 |
| 1.5 | 2.828427 | **2.78231** | 2.748762 | 2.72266 |
| 1.75 | 3.363586 | **3.31673** | 3.278829 | 3.24697 |
| −0.5 | 0.707107 | **0.58923** | 0.544764 | 0.52576 |
| 0.5+0.5i | — | **1.3884+0.4720i** | 1.4208+0.4714i | 1.4395+0.4728i |

Unsicherheit der H₃.₅- und H₄.₅-Werte: heuristisch ~10⁻⁴ (siehe unten), bei 0.5+0.5i ~5·10⁻⁴. Die Leiter ist monoton im Rang: links von 1 steigend gegen x+1, rechts fallend.

**Aber:** Die ε-Reihe **divergiert**, und ab Ordnung 8 ist der Halbrang-Jet mit dieser Auswahl **gar nicht definiert**. Die ganzzahligen Rangjets c_n^{(j)} wachsen dort exponentiell im Rang n, mit höhenunabhängigen Raten (j=10: ρ≈117, j=12: ρ≈−2667). H₃.₅(2,x) ist so nur als asymptotische Entwicklung mit Genauigkeit ~10⁻⁴ bestimmt.

**Etappe 2.** Die Rangrekursion T(H₃.₅)=H₄.₅ gilt in jeder berechneten Ordnung auf der Genauigkeit der Rangfortsetzung (Ordnung 1: 10⁻³⁴; Ordnungen 2–7: ≤2·10⁻⁶ reell, ≤7·10⁻⁶ komplex). H₃.₅(2,1)=2 und H₃.₅(2,2)=4 gelten exakt (bewiesen). Zwei Halbschritte ergeben in der Konstruktion den perturbativen Rang-1-Jet T(2ˣ). Dessen optimal abgeschnittene Summe weicht von Kneser um **0.8·10⁻⁴ bis 3·10⁻⁴** auf den reellen Testpunkten ab, und Kneser liegt an allen reellen Punkten im Band der Partialsummen. **Das ist verträglich mit Kneser, aber kein Nachweis auf viele Stellen:** feiner als ~10⁻⁴ lässt die divergente Reihe nicht prüfen, Borel–Padé stabilisiert sich nicht.

## Die Auswahl (eine Wahl)

1. **Höhe:** In jeder ε-Ordnung muss der Koeffizient auf einer Halbebene holomorph sein, mit endlichem horizontalem Exponentialtyp und polynomialem vertikalem Wachstum (Klasse 𝒜 der Vorrunde). Das schliesst die periodische Freiheit θ(x) ordnungsweise aus. An ganzen Rängen sind alle Koeffizienten dann **exakte Exponentialpolynome** Σ_{m≤j}2^{mx}P(x) (Satz C).
2. **Rang:** Newtonreihe im Rang (polynomiales Wachstum in n, wie (A8) der Vorrunde), numerisch beschleunigt mit Levin-u.
3. **Summation in ε:** optimaler Abbruch (kleinstes Glied). Das ist eine Notlösung, keine Auswahl mit Eindeutigkeitssatz.

## Was bewiesen ist

- Lemma A/B, Satz C: Eindeutigkeit und Exponentialpolynom-Form aller ganzzahligen Rangjets in 𝒜, alle Ordnungen, mit expliziter Lösung der Kohomologiegleichung Δe_j=G_j (Verallgemeinerung von Δw=v(x+1), Δz=w·v′(x+1)).
- Satz D: H₃.₅(2,1)=2, H₃.₅(2,2)=4 für jede lineare Rangfortsetzung (also in die Konstruktion eingebaut, kein Test der Auswahl).
- Ordnung 1 in geschlossener Form (Vorrunde (A18)); hier numerisch auf 34 Stellen getroffen.
- Elementar: Newtonfortsetzung einer Folge ρⁿ divergiert bei |ρ−1|>1; bei |log|ρ||≥π fehlt die Eindeutigkeit über den Wachstumstyp; für ε<0 hat f_ε zwei reelle Fixpunkte.

## Numerisch (keine Intervallzertifikate)

Die Exponentialpolynome selbst sind mit Ballarithmetik (python-flint, 1000–3200 Bit) gerechnet; ihre Rundung ist kontrolliert (Kugelradien <10⁻³⁰⁰). Nicht zertifiziert sind die Rangfortsetzung (Levin) und alle ε=1-Summen.

### Halbrang-Glieder u_j(x) bei ε=1 (Fehler = Abstand Levin N=40 zu N=32)

| j | x=0.25 | x=0.5 | x=0.75 | x=1.5 |
|---:|---:|---:|---:|---:|
| 1 | −0.0391314 (exakt) | −0.0610330 (exakt) | −0.0526418 (exakt) | +0.302347 (exakt) |
| 2 | +1.0110e−3 | +3.5333e−3 | +4.1858e−3 | −2.30191e−2 |
| 3 | +5.8104e−4 | +7.3955e−4 | +3.3815e−4 | +2.93856e−3 |
| 4 | +1.144e−4 | −1.94e−5 | −1.433e−4 | +3.7283e−4 |
| 5 | −4.63e−5 | −1.294e−4 | −1.034e−4 | −1.8096e−4 |
| 6 | −6.59e−5 | −7.23e−5 | −1.34e−5 | −1.533e−4 |
| 7 | −3.16e−5 | +8.9e−6 | +4.18e−5 | −2.74e−5 |
| 8 | nicht definiert | nicht definiert | nicht definiert | nicht definiert |

Fehler der Rangfortsetzung: Ordnung 1 ≤10⁻²⁶, Ordnungen 2–3 ≤10⁻⁹, Ordnungen 4–7 ≤2·10⁻⁶. Komplexe Punkte (0.5+0.5i, 0.25+0.25i, 1.5+0.5i, 1.25+0.25i) verhalten sich gleich, mit etwas grösseren Gliedern; bei 0.5+1i fallen die Glieder nur auf ~4·10⁻³. Volle Daten: `data/auswertung_N40J8R7.json`, `data/summen.json`.

**Abfall der Ordnungen:** grob Faktor 10–20 pro Ordnung bis j≈4, danach Stagnation bei 10⁻⁵–10⁻⁴. Bei Rang 1 (52 Ordnungen gerechnet) sieht man, was danach kommt: Minimum ~3·10⁻⁵ bei j≈6–10, dann faktorielles Wachstum bis 10¹⁸ bei j=52 (log₁₀-Zuwachs 0.43 pro Ordnung bei j≈25, 0.58 bei j≈35; passend zu Gevrey-1 mit A grob 7–9, Fit unscharf). Abbildung: `figures/glieder_und_rangwachstum.png`.

### Etappe 2: Diskrepanzen

**Rang 1 (zwei Halbschritte) gegen Kneser**, ε=1, 52 Ordnungen:

| x | Summe bei optimalem Abbruch | Kneser | Abstand | kleinstes Glied | Band der Partialsummen 4…12 |
|---:|---:|---:|---:|---:|---|
| 0.25 | 1.224284 | 1.224361 | 7.7e−5 | 2.6e−5 | [1.224258, 1.224623] ∋ K |
| 0.5 | 1.459009 | 1.458782 | 2.3e−4 | 3.0e−5 | [1.458327, 1.459038] ∋ K |
| 0.75 | 1.713306 | 1.713180 | 1.3e−4 | 2.5e−5 | [1.712591, 1.713331] ∋ K |
| 1.25 | 2.336365 | 2.336520 | 1.5e−4 | 5.5e−5 | [2.336365, 2.336922] ∋ K |
| 1.5 | 2.748529 | 2.748762 | 2.3e−4 | 3.5e−6 | [2.748351, 2.749354] ∋ K |
| 1.75 | 3.279143 | 3.278829 | 3.1e−4 | 8.0e−5 | [3.277953, 3.279352] ∋ K |
| −0.5 | 0.544996 | 0.544764 | 2.3e−4 | 2.8e−5 | [0.544303, 0.545023] ∋ K |
| 0.5+0.5i | 1.420013+0.471276i | 1.420790+0.471367i | 7.8e−4 | 7.5e−4 | — |
| 0.5+1i | 1.312932+0.868181i | 1.323239+0.869340i | 1.0e−2 | 1.1e−2 | — |

Kneser-Werte über `experimente/gemeinsam.kneser('sexp','2',…)`, 40 Stellen, ein einziger Lauf. Der Abstand wächst mit |Im x| genauso wie das kleinste Glied; auch dort ist kein Widerspruch zu Kneser messbar. Zum Vergleich: die lineare Näherung (sexp=x+1 auf (−1,0]) liegt bei x=0.5 um 0.045 neben Kneser und ist klar ausgeschlossen.

**Rangrekursion H₄.₅(x+1)=H₃.₅(H₄.₅(x))** — Residuum je Ordnung (Jet-Ebene):

| x | j=1 | j=2 | j=3 | j=4 | j=5 | j=6 | j=7 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.25 | 2e−34 | 1.3e−7 | 8.9e−8 | 1.9e−7 | 6.4e−8 | 4.6e−8 | 6.9e−8 |
| 0.5 | 1e−34 | 7.9e−11 | 3.6e−8 | 5.4e−7 | 1.1e−7 | 2.9e−8 | 1.4e−7 |
| 0.75 | 4e−35 | 1.5e−9 | 1.2e−10 | 1.2e−8 | 2.0e−6 | 8.6e−9 | 7.6e−9 |
| 1.5 | 5e−36 | 2.9e−13 | 5.8e−10 | 1.1e−7 | 7.8e−7 | 9.3e−6 | 3.5e−5 |
| 0.5+0.5i | 4e−34 | 5.3e−12 | 3.4e−11 | 1.1e−6 | 5.1e−6 | 5.1e−6 | 6.9e−6 |

Die Residuen liegen auf der Höhe der Levin-Fehlerschätzungen; mit ganzzahligen Rängen (Kontrolle) ist die Gleichung auf 10⁻⁶⁰ erfüllt. Bei ε=1 mit abgeschnittenen Jets (Ordnung ≤7): Residuum 1.2·10⁻⁵ (x=0.25), 1.4·10⁻⁴ (0.5), 2.5·10⁻⁴ (0.75), 3.6·10⁻⁶ (−0.5), 8.2·10⁻⁴ (0.5+0.5i) — Grösse der weggelassenen Ordnungen.

An ganzen Höhen: |T(H₃.₅)(3) − H₄.₅(3)| je Ordnung ≤4·10⁻⁵ (j≤7); bei den Höhen 4–6 wird die Rangfortsetzung selbst ungenau (dort wachsen die Rangfolgen polynomial mit hohem Grad), das ist kein Gegenbefund.

### Welche Tetration reproduzieren zwei Halbschritte?

Die perturbative Tetration (der Rang-1-Jet) — und diese ist, soweit die divergente Reihe auflöst (~10⁻⁴), von Kneser nicht zu unterscheiden. Eine Tetration, die sich von Kneser durch ein 1-periodisches θ der Grösse unter ~10⁻⁴ unterscheidet, wäre ebenso verträglich. Eine Identifikation „Auswahl = Kneser auf viele Stellen“ ist mit diesem Weg nicht erreichbar.

## Was offen ist

1. **Rangfortsetzung ab Ordnung 8.** Die Rangfolgen wachsen wie ρ_jⁿ mit |ρ_j|: ≈7.6, 28, 117, 470, 2667 (j=8…12), unabhängig von x. Newton/Levin divergiert; ab |log ρ|≥π fehlt auch die Eindeutigkeit über den Wachstumstyp, und negative ρ (j=8, 12) liefern nichtreelle ρ^{1/2}. Ursache unbekannt.
2. **Summation in ε.** Gevrey-Ordnung und Borel-Singularitäten beweisen; eine Summationsvorschrift finden, die Kneser auszeichnet (oder nicht).
3. Damit ist eine Konstruktion von H₃.₅ mit vielen gesicherten Stellen auf diesem Weg **nicht** möglich. Der Engpass ist nicht mehr die Höhenklasse (die ist mit Satz C sauber), sondern die Entwicklung um den Successor selbst.

## Nächste Schritte (Vorschlag)

- Rangfortsetzung **im Borel-Raum**: B_n(t)=Σ_j c_n^{(j)}t^j/j! für kleine t in n fortsetzen und erst danach summieren; prüfen, ob die Raten ρ_j dort zu einer zahmen n-Abhängigkeit werden.
- Die Raten ρ_j erklären: Wachstum der Rangjets an negativen Höhen und in der Bernoulli-Komponente untersuchen (Komponenten heben sich um 10¹⁶⁸ auf).
- Alternativ eine Entwicklung um einen anderen Bezugspunkt, an dem die Fixpunktstruktur bei ε→0 nicht entartet (das ist beim Successor unvermeidlich: Fixpunkte wandern nach ∞).
- Die Werte H₃.₅(2,x) aus der Tabelle als Kontrollzahlen (±10⁻⁴) für jede künftige Konstruktion verwenden.

## Reproduzierbarkeit und Grenzen

- Längster Einzellauf: Rang 1 mit 52 Ordnungen, 730 s; `raenge.py --N 40 --J 12`: ~420 s. Keine Rechnung über 20 Minuten; ein Versuch mit 60 Ordnungen wurde nach 19 Minuten durch die Zeitgrenze beendet und nicht wiederholt.
- Ein einziger Kneser-Aufruf (12 Punkte, 40 Stellen, 5 s). **Hinweis:** `gemeinsam.kneser` hat dabei die fehlenden Punkte selbsttätig an `experimente/daten/kneser_cache.json` angehängt (Zeitstempel 13.10.2026 14:57:56, gleichzeitig mit `data/kneser_basis2.json`). Es wurden nur Einträge ergänzt, keine geändert; ich habe nichts zurückgesetzt. `gemeinsam.py`, `daten/`, `ansicht/` unverändert; keine Commits.
- Literatur: [LITERATUR.md](LITERATUR.md).
