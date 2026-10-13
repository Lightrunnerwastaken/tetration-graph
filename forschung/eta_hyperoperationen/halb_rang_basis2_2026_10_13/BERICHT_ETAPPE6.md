# Etappe 6: Einheitliche Kneser-artige Ränge — Ergebnis

13.10.2026 · gleicher Ordner · neue Dateien `e6_*.py`, Daten `data/e6_*`

## Ergebnis zuerst

1. **Der verallgemeinerte Cauchy-Löser trägt bis Rang 6.** H₄ (Kneser) stimmt mit den GP-Werten auf **1·10⁻¹²** überein, auch an komplexen Höhen. H₅ und H₆ sind einheitlich Kneser-artig berechnet. Die innere Kontrolle F(w+1)=f(F(w)), direkt im Rechteck ausgewertet, liegt bei ~4·10⁻¹⁰. Fixpunktpaare der unteren Stufen: L₄=1.469818+1.664342i (Multiplikator 0.154+0.691i), L₅=1.807677+1.596922i, L₆=2.008150+1.515636i.
2. **Die Ränge der Hauptrunde unterscheiden sich davon um eine θ-Welle.** Bei H₅ (Hauptrunde: reguläre Iteration am reellen Fixpunkt −1.7439) ist die Differenz eine **reine θ-Welle**: D/H′ = a·sin2πx + b(1−cos2πx) mit a=−4.40·10⁻⁶ und b=3.46·10⁻⁶. Mit zwei reellen Punkten angepasst, trifft das Modell auch 0.5+0.25i, 0.5+0.5i, 0.25+0.5i und 1.5+0.5i auf 3–4 Stellen. Bei H₆ (Hauptrunde: Zwei-Bogen-Auswahl) liegt die Differenz bei 1–3·10⁻⁵ und ist θ-artig, das Zwei-Parameter-Modell trifft aber nur auf 10–30 %. Das ist eine Folge höherer Harmonischer oder der Genauigkeit der alten Hexation. Die Differenz wächst also um etwa den Faktor 4–6 pro Rang. Hochgerechnet ergibt das ~10⁻³ bei N=8–9, passend zu den Knicken aus Etappe 5. Diese Erklärung ist plausibel, aber nicht nachgerechnet.
3. **Rang 7 ist mit diesem Löser nicht erreichbar.** Die Iteration konvergiert (Residuum 10⁻¹³), aber zu **unechten** diskreten Lösungen: Die Funktionalgleichung ist um 0.06–0.14 verletzt, und H₇(0.5) liegt je nach Startwert bei 1.41 oder 1.45. Ursache: Am linken Rand braucht man H₆⁻¹, und H₆ ist nicht injektiv. Die Newton-Inverse springt bei |y|≈3 auf einen fremden Urbildzweig (Steigung 49 statt ~1). Eine Zweigverfolgung pro Knoten ist in dieser Umsetzung zu langsam; zwei Läufe brachen an der 19-Minuten-Grenze ab.
4. **Die einheitliche Auswahl verbessert H₃.₅ nicht.** Mit H₃…H₆ einheitlich und H₃…H₆ der Hauptrunde unterscheiden sich die Newtonwerte bei 3.5 nur um 10⁻⁶–10⁻⁵, während der Abbruchfehler der Newtonreihe mit 4 Rängen ~10⁻³ beträgt. **Die Grenze ist die Zahl der Ränge, nicht die Auswahl**, zumindest bis Rang 6.
5. **Ein Abstieg im Rang (Fatou-artige Definition) ist numerisch instabil.** Getestet wurde H₃.₅=R^k(G), mit R(S)(y)=S(S⁻¹(y)+1) und G als Interpolation hoher Ränge. Bewiesen: R hängt nicht von der θ-Wahl der hohen Stufen ab, denn R(S∘ψ)=R(S) für ψ(x+1)=ψ(x)+1. Der Weg ist also selektionsfrei. Aber R verstärkt Fehler stark. Mit 3 Abstiegen streut H₃.₅(2,0.5) zwischen 1.44138 und 1.44699 je nach Interpolationsknoten, mit 4–5 Abstiegen kommen unsinnige Werte (0.99–1.56).

**Stand H₃.₅:** unverändert H₃.₅(2,0.5)≈1.4428±3·10⁻⁴ und H₃.₅(2,1.5)≈2.7824±2·10⁻⁴. Über 10⁻⁴ hinaus kommt in dieser Etappe keiner der Wege.

## Tabellen

**H₅ einheitlich gegen Hauptrunde** (D=H₅,K−H₅,alt; Modell aus x=0.25 und 0.75):

| z | D | D/H′ | Modell a·sin+b(1−cos) |
|---|---|---|---|
| 0.5 | 6.73e−6 | 6.914e−6 | 6.915e−6 |
| −0.5 | 6.90e−6 | 6.914e−6 | 6.915e−6 |
| 1.5 | 1.21e−5 | 6.914e−6 | 6.915e−6 |
| 1.25 | −1.30e−6 | −9.397e−7 | −9.397e−7 |
| 0.5+0.25i | 1.11e−5+1.02e−5i | 1.213e−5+1.012e−5i | 1.213e−5+1.012e−5i |
| 0.5+0.5i | 3.54e−5+4.91e−5i | 4.352e−5+5.078e−5i | 4.353e−5+5.078e−5i |
| 0.25+0.5i | −4.30e−5+3.48e−5i | −4.750e−5+3.993e−5i | −4.751e−5+3.993e−5i |

D/H′ ist 1-periodisch (Werte bei −0.5, 0.5 und 1.5 gleich) und wächst senkrecht wie e^{2π|y|}. Das ist die θ-Form aus Etappe 4, jetzt zwischen zwei *nichtperturbativen* Pentationen. Daten: `data/e6_vergleich.json`.

**H₆:** D=2.6·10⁻⁵ (0.5), 2.2·10⁻⁵ (0.75), 7.7·10⁻⁶ (0.25), 1.3·10⁻⁴+8.6·10⁻⁵i (0.5+0.5i); Modellabweichung 10–30 %.

**Rangabstieg** (`e6_abstieg.py`, Ränge 4–9 der Hauptrunde):

| Knoten, Start, Abstiege | H₃.₅(0.25) | H₃.₅(0.5) | H₃.₅(0.75) | H₃.₅(1.5) |
|---|---:|---:|---:|---:|
| 5..9, N=6.5, 3 | 1.22018 | 1.44699 | 1.70315 | 2.78168 |
| 6..9, N=6.5, 3 | 1.21057 | 1.44138 | 1.70082 | 2.78312 |
| 5..9, N=7.5, 4 | 1.58737 | 1.49894 | 1.72525 | 2.78737 |
| 4..9, N=7.5, 4 | 1.59398 | 1.56214 | 1.70875 | 2.75005 |
| 5..9, N=8.5, 5 | 0.98573 | 0.99298 | 1.19144 | 2.87903 |

## Gekoppeltes Problem und eine Gleichung für U allein: erste Einordnung

- Jede Wahl von u₃=H₃.₅ erzeugt mit dem Kneser-Operator eine eigene Kette u_{n+1}=T(u_n). Die Forderung U(u_n)=H_{n+1} ist auf Funktionsebene keine Gleichung in x, sondern eine Operatorgleichung. Eine **Funktionalgleichung in x für U allein habe ich nicht gefunden**, und ich sehe keinen Kandidaten. Die Rangrekursion verknüpft immer zwei verschiedene Funktionen (H₃.₅ und H₄.₅).
- Die einzige natürliche Zusatzbedingung ist **asymptotisch im Rang**: Die Kette u_n muss sich für n→∞ „halb zwischen“ H_n und H_{n+1} einordnen. Das ist die Fatou-Koordinaten-Idee am parabolischen Fixpunkt σ(x)=x+1 von T. Formal: H₃.₅ = lim_{n→∞} R^n(G_{n+3.5}), mit G aus dem nahezu linearen Bereich hoher Ränge (Gamma-Profil B_{m+1/2}). Diese Definition ist von θ unabhängig (Punkt 5 oben). Numerisch verlangt sie aber hohe Ränge mit einer Genauigkeit, die den Fehlerverstärkungsfaktor von R überwiegt. Den habe ich nicht bestimmt; er ist gross (Tabelle).
- **Bewertung:** Der begrifflich sauberste Weg ist der Fatou-artige Abstieg; er fixiert θ ohne Zusatzwahl. Mit den heute verfügbaren Rängen (bis 9, 11–22 Stellen) ist er aber nicht rechenbar.

## Was für Rang 7 und höher nötig wäre (Aufwandsschätzung)

Gebraucht wird eine schnelle, zweigtreue Inverse der unteren Stufe am linken Rand. Ein möglicher Weg: die Abel-Funktion der unteren Stufe direkt aus ihrer Streifendarstellung tabellieren (Chebyshev auf dem Bild des Streifens) statt Punkt-Newton.
- **Schätzung:** 1–2 Tage Entwicklung; danach pro Rang ~10–30 min Rechenzeit in doppelter Genauigkeit.
- Genauigkeit vermutlich ~10⁻⁹, mit einem Verlust von etwa einer Stelle pro Rang.
- Ränge bis ~12 wären damit denkbar. Ob die Newtonreihe im Rang dann die 10⁻⁴-Grenze unterschreitet, ist offen: Mit 7 Rängen lag die letzte Korrektur noch bei 2·10⁻⁴. Ich habe das nicht gestartet.

## Bewiesen / numerisch / offen

- **Bewiesen:** R(S∘ψ)=R(S) für ψ(x+1)=ψ(x)+1; der Rangabstieg ist von θ unabhängig.
- **Numerisch:** Löser gegen Kneser 10⁻¹²; θ-Welle zwischen den Pentationen (Amplitude ~5·10⁻⁶); H₆-Differenz ~3·10⁻⁵; Instabilität des Abstiegs.
- **Offen:** Rang ≥7 einheitlich; eine Funktionalgleichung für U allein; Stabilitätsanalyse von R.

## Reproduktion (aus diesem Ordner)

```powershell
python e6_raenge.py 5 30 500 4 2000      # H4, H5 (~20 s)
python e6_raenge.py 6 30 500 6 2000      # H6 auf H4, H5 aufbauend (~50 s); Rang 7 liefert unechte Lösungen
python e6_vergleich.py                    # theta-Welle gegen Hauptrunde
python e6_rangfolge.py 6                  # Rangfolgen und Newton bei 3.5
python e6_abstieg.py 5:3,6:3,5:4,4:4,5:5  # Rangabstieg (~7 min, Motor der Vorrunde)
```

`data/e6_H7_A30_N500.pkl` und der H₇-Eintrag in `data/e6_raenge_A30_N500.json` sind die **unechte** Rang-7-Lösung (Funktionalgleichung um 0.06 verletzt). Sie bleiben nur zur Dokumentation liegen; `e6_lade.lade_raenge(6)` lädt nur bis H₆. Keine `gemeinsam.kneser`-Aufrufe. Vorhandene Dateien unverändert; keine Commits.
