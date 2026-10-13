# Etappe 3: Rangfortsetzung im Borel-Raum — Ergebnis

13.10.2026 · gleicher Ordner · neue Dateien `e3_*.py`, Daten `data/e3_*.json`

## Ergebnis zuerst

**Das Ziel „H₃.₅(2,x) auf deutlich mehr als 4 Stellen“ ist auf diesem Weg nicht erreichbar, und der Grund ist jetzt genau benannt.** Die Divergenz der ε-Reihe und das exponentielle Wachstum im Rang sind **dasselbe Phänomen**, nämlich die periodische θ-Freiheit der Superfunktion:

- Die wachsenden Anteile der Jets (späte ε-Ordnungen bei festem Rang, hohe Ordnungen bei wachsendem Rang) haben in der Höhe die Form
  **A(x) ≈ a·sin(2πx) + b·(1−cos 2πx)**.
  Sie sind 1-periodisch, verschwinden an ganzen Höhen und wachsen wie cosh(2π Im x). Das sind genau die Funktionen, die zwei Tetrationen voneinander unterscheiden, also die θ-Moden.
- Die Borel-Singularitäten liegen bei Rang 1 bei |t|≈7–9 und bei Rang 2 bei |t|≈0.7–1.1. Sie rücken also von Rang 1 zu Rang 2 um etwa den Faktor 10 näher an den Ursprung (gemessen nur für diesen einen Schritt); bei Rang 0 gibt es gar keine. Daraus folgt das Ranggewicht ρ_jⁿ.
- Seitliche Borel-Summen unterscheiden sich genau um solche θ-Moden. **Eine Summationsvorschrift wählen heisst θ wählen.** Ohne eine nichtperturbative Zusatzbedingung (wie Knesers Bedingung für Im x→∞) kann die Entwicklung um den Successor nicht zwischen Kneser und einer θ-Variante entscheiden. Die Grenze liegt bei etwa 10⁻⁴ (Rang 1), und diese Zahl ist die Grösse der θ-Mehrdeutigkeit selbst.
- Einen nicht entartenden Bezugspunkt gibt es nicht. Jeder injektive Fixpunkt von T ist der Successor bei b=2 (Lemma 11 aus `rang_gamma`). Als unabhängige Gegenprobe wurde deshalb die volle (nicht entwickelte) Rangkette H₃…H₉ im Rang fortgesetzt. Sie stimmt mit dem perturbativen H₃.₅ auf 10⁻⁴ bis 5·10⁻⁴ überein, aber nicht genauer.

## 1. Woher kommen die Raten ρ_j? (numerisch, mit klarem Befund)

**Form in der Höhe.** Bei Rang n=40 wird die Ordnung j durch c_{40}^{(j)}(x)/c_{40}^{(j)}(0.5) normiert. Das Modell hat einen einzigen freien Parameter a, angepasst bei x=0.25 (b=½ folgt aus A(0.5)=1):

| x | j=12 gemessen | Modell a·sin2πx+½(1−cos2πx), a=0.7947 |
|---|---|---|
| 0.75 | −0.29450 | −0.29468 |
| 1.25 / −0.5 / 1.5 | 1.2941 / 1.0004 / 0.99957 | 1.2947 / 1 / 1 |
| 0.25+0.25i | 2.4940+1.1505i | 2.4940+1.1506i |
| 0.5+0.5i | 6.2940−9.1772i | 6.2960−9.1775i |
| 0.5+1i | 134.28−212.79i | 134.37−212.77i |
| ganze Höhen 0…6 | ≤4·10⁻²⁵ | 0 |

Für j=10 ist der Fit ebenso gut (Abweichung ≤2 %). Die **späten Glieder** der ε-Reihe zeigen dieselbe Form: bei Rang 2, Ordnung 36, mit a=1.67 auf ≤0.7 %, und bei Rang 1, Ordnung 36, auf den Punkten links von 1 auf ≤0.3 % reell und ≤8 % komplex, bei x=1.5 aber 23 % daneben. Dort ist die Periodizität in der Abel-Koordinate statt in x zu erwarten, weil die Funktion vom Successor abweicht. Daten: `data/e3_periodisch.json`.

**Lage der Borel-Singularitäten** (Padé der Borel-Transformierten, 36 Ordnungen, `data/e3_borel_raenge.json`):

| Rang | nächste Pole |t| | Zuwachs log₁₀|c_j| pro Ordnung bei j≈30 |
|---:|---|---:|
| 0 | keine (f_ε ist linear in ε) | — |
| 1 | 3.5–8.4 (Kranz um |t|≈8; die nächsten wechseln mit dem Padé-Grad) | ≈0.45 |
| 2 | 0.68–1.14 | ≈1.5 |

Das Verhältnis c₂^{(j)}/c₁^{(j)} wächst ab j≈14 um etwa den Faktor 10 pro Ordnung (bei j=36: 1.3·10¹⁹ bei x=0.5). Das passt zu A₁/A₂≈10.

**Raten bei hohen Rängen** (n=24…40, aus Etappe 1): |ρ_j|≈7.6, 28, 117, 470, 2667 für j=8…12; |ρ_{j+2}/ρ_j|≈15–23. Die Phasen folgen dem Muster arg ρ_j≈(j+2)π/2 (j=8 negativ, j=9 und 11 rein imaginäres Paar, j=10 positiv, j=12 negativ). Eine Beschreibung ρ_j≈−C·(i e^{π/2})^j mit C≈1.8·10⁻⁵ trifft die Grössenordnung auf ±25 %; das ist eine **Beobachtung ohne Herleitung**.

**Deutung (Heuristik, kein Beweis).** Bei einer divergenten Reihe tragen die späten Glieder die Höhenabhängigkeit der exponentiell kleinen Zusatzterme (Late-terms-Prinzip nach Dingle; hier nicht gelesen, nur als Stichwort). Diese Zusatzterme sind hier die periodischen Moden e^{±2πix}, also genau das, was die Klasse 𝒜 ordnungsweise ausschliesst. Weil die Moden im Kern von Δ liegen, kann sie keine perturbative Auswahl festlegen. Mit jedem Rangschritt rückt ihre Borel-Singularität näher, und so entstehen die Raten ρ_j.

## 2. Rangschritt im Borel-Raum: trägt nicht

- Rangfortsetzung und Borel–Laplace vertauschen nur, wenn die Borel-Transformierte B_n(t) als Funktion von n zur Auswahlklasse gehört. Bei Rang 0 hat sie keine Singularität, bei Rang 1 eine bei ≈8, bei Rang 2 eine bei ≈0.75. Die Singularität entsteht erst im ersten Rangschritt. Eine Rangfortsetzung A_s mit A_0=∞, A_1≈8, A_2≈0.75 ist durch die ganzzahligen Daten nicht eindeutig bestimmt.
- Noch grundsätzlicher: Selbst bei bekannter Singularitätslage hängt die Laplace-Summe von der Integrationsrichtung ab. Die seitlichen Summen unterschieden sich schon in Etappe 2 um 10⁻⁴–10⁻³, und der Unterschied besteht aus θ-Moden. Der Borel-Raum verschiebt die Wahl von θ, er beseitigt sie nicht.

## 3. Alternative: nicht entartender Bezugspunkt und volle Werte

**Bewiesen (Vorrunde, rang_gamma Lemma 11):** Ein injektiver Fixpunkt von T_b muss f(x)=x+1 mit b=2 sein. Es gibt also keinen anderen Koenigs-Punkt, und jede Entwicklung um einen Fixpunkt ist die Successor-Entwicklung, deren Grenze oben beschrieben ist.

**Nichtperturbative Gegenprobe:** Newtonreihe im Rang über die vollen Werte H_N(2,x), N=3…9, aus `../half_height_values.json`. Für N=4 ist das Kneser, für N≥5 die Fortsetzungskette der Vorrunde. Skript `e3_vollwerte.py`.

| Grösse | volle Rangkette (Newton M=7 bzw. 6 / Levin) | perturbativ (Etappe 1) | Abstand |
|---|---|---|---:|
| H₃.₅(2,0.5) | 1.44246 / 1.44266 (bei M=4…7 noch um 2·10⁻⁴ wandernd) | 1.44302 | 4–6·10⁻⁴ |
| H₃.₅(2,1.5) | 2.78254 / 2.78245 | 2.78231 | 1.4–2.3·10⁻⁴ |
| H₄.₅(2,0.5) | 1.46835 / 1.46839 | 1.46827 | ≈1·10⁻⁴ |
| H₄.₅(2,1.5) | 2.72319 / 2.72319 | 2.72266 | 5·10⁻⁴ |

Beide Wege stimmen im Rahmen ihrer Unschärfe überein. Der volle Weg ist durch die wenigen verfügbaren Ränge begrenzt: 7 Werte, Rang 9 hat nur 22 Stellen, und Rang ≥10 ist bei Basis 2 bisher nicht gerechnet. Mehr als ~3–4 Stellen liefert er auch nicht.

## Was bewiesen ist, was nicht

- **Bewiesen:** f_ε ist linear in ε, Rang 0 hat also keine Borel-Singularität. Es gibt keinen weiteren injektiven Fixpunkt von T (Lemma 11). Periodische Funktionen mit polynomialem vertikalem Wachstum sind konstant (Lemma A), die θ-Moden liegen also ausserhalb von 𝒜.
- **Numerisch, gut belegt:** die θ-Form der wachsenden Moden (Fit mit einem Parameter auf ~10⁻³ relativ, auch komplex); Borel-Pole bei Rang 1 ≈8 und bei Rang 2 ≈0.75; Übereinstimmung zwischen perturbativem und vollem H₃.₅ auf ≤6·10⁻⁴.
- **Offen:** eine Herleitung der Raten ρ_j und der Borel-Lagen A_n; die Stokes-Konstante bei Rang 1; welche seitliche Summe Kneser ist.

## Nächste Schritte (wenn weiter gewünscht)

1. **θ von Kneser bei Rang 1 messen.** Kneser minus perturbative Summe auf a·sin2πx+b(1−cos2πx) projizieren, mit hyperasymptotischer statt optimal abgeschnittener Summe. Wenn Kneser eine feste seitliche bzw. die Median-Summe ist, gibt das eine nichtperturbative Regel, die man auf den Halbrang übertragen kann. Das ist der einzige Weg, den ich sehe, damit „zwei Halbschritte = Kneser“ wirklich entscheidbar wird.
2. Alternativ müsste die volle Rangkette um Ränge 10–14 an mehreren Höhen verlängert werden. Das ist teuer, die Vorrunden scheiterten schon bei Rang 8–9.

## Reproduktion (aus diesem Ordner, Python `C:/Users/janis/miniconda3/python.exe`)

```powershell
python e3_raenge_hoch.py --J 36 --N 3 --prec 2000   # Ränge 1,2 mit 36 Ordnungen (~18 min; Rang 3 wurde wegen des Zeitbudgets nicht gerechnet)
python e3_borel_raenge.py                           # Borel-Pole, Quotienten, späte Glieder
python e3_periodisch.py                             # θ-Form der wachsenden Moden (nutzt data/raenge_N40J12.json)
python e3_vollwerte.py                              # Rangfortsetzung der vollen Werte N=3..9
```

Keine neuen Kneser-Aufrufe. Alle Läufe unter 20 Minuten; der längste dauerte 18 min (Rang 2, 36 Ordnungen). Vorhandene Dateien, BERICHT.md und MATHEMATIK.md, sind unverändert; keine Commits.
