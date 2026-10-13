# Etappe 5: Nichtperturbativer Halbrang aus den vollen Rechenarten — Ergebnis

13.10.2026 · gleicher Ordner · neue Dateien `e5_*.py`, Daten `data/e5_*.json`

## Ergebnis zuerst

1. **Nahe der reellen Achse trägt der Weg (|Im x| ≤ 0.5).** Mit dem Motor der Vorrunde habe ich die vollen Werte H_N(2,x+iy) für **N=3…9** an einem Raster berechnet. Der Motor wurde dafür nur in den eigenen Ordner kopiert; Rang 9 war in 74 s initialisiert. Die Newtonreihe im Rang gibt H₃.₅(2,x) mit einer letzten Korrektur von 0.4–5·10⁻⁴. Der Wert liegt 0.9–12·10⁻⁴ neben dem perturbativen Halbrang aus Etappe 1 (bei x=−0.5 5·10⁻³, dort konvergiert die Rangreihe schlecht); das Vorzeichen des Trends passt dazu, dass beide Wege gegen denselben Wert laufen (Tabelle unten).
2. **Die Rangrekursion H₄.₅(x+1)=H₃.₅(H₄.₅(x)) gilt auf 3·10⁻⁶ bis 1.2·10⁻⁴ (reell) und ≤4.4·10⁻⁴ (komplex).** Das ist dieselbe Grössenordnung wie beim perturbativen Weg bei ε=1 (Etappe 2: 1·10⁻⁵…3·10⁻⁴).
3. **Die Auswahlregel aus Etappe 4 (Beschränktheit auf senkrechten Geraden) lässt sich hier nicht anwenden**, aus zwei genau benennbaren Gründen:
   - **Sie wählt nichts aus (bewiesen, elementar).** Jede Rangfortsetzung dieser Art ist eine Summe Σ_k w_k(s)·H_{3+k}(z), deren Gewichte nicht von z abhängen. Jedes H_N ist auf senkrechten Geraden beschränkt: |2^{x+iy}|=2ˣ, Kneser strebt gegen L, die regulären Stufen sind in y periodisch. Also ist auch die Fortsetzung beschränkt. Die θ-Freiheit tritt erst über die Funktionalgleichung (Rangrekursion) auf, und die erfüllt die Fortsetzung nur bis ~10⁻⁴.
   - **Die Daten fehlen.** Für |Im x| ≥ 1 versagt die Hexation des Motors (Fehler „knhex: Newton bleibt nicht im Bogen“), damit auch alle Ränge darüber. Es bleiben nur N=3,4,5, und deren senkrechte Strukturen passen nicht zusammen: Bei x=0.5+2i ist H₃=0.26+1.39i, H₄=1.09+1.34i, H₅=1.34−1.28i (reguläre Iteration, in y periodisch mit eigener Periode). Die Newton-Korrekturen dort sind 0.3–0.4, also bedeutungslos.
4. **„Zwei Halbschritte = Kneser“ ist in dieser Konstruktion kein Test.** Eine Interpolation über ganzzahlige Ränge, an halbzahligen Stellen neu interpoliert, gibt H₄ per Konstruktion zurück. Der echte Test wäre, dass H₄.₅ die *Kneser-Superfunktion* von H₃.₅ ist. Dafür bräuchte man H₃.₅ auf einem ganzen Streifen bis |Im x|≈30 und seine Fixpunkte, und das liefert dieser Weg nicht.

**Stand der besten Werte (beide Wege zusammen):** H₃.₅(2,0.5) ≈ 1.4425–1.4430, H₃.₅(2,1.5) ≈ 2.7823–2.7825; Genauigkeit realistisch ±3·10⁻⁴. Mehr als 3–4 Stellen liefert keiner der beiden Wege.

## Tabellen

**H₃.₅ aus den vollen Werten (Newton durch N=3…9) gegen perturbativ (Etappe 1, Ordnungen 1–7):**

| z | voll (N=3…9) | letzte Korrektur | perturbativ | Abstand |
|---|---|---:|---|---:|
| 0.25 | 1.211889 | 1.9e−4 | 1.212431 | 5.4e−4 |
| 0.5 | 1.442465 | 2.1e−4 | 1.443028 | 5.6e−4 |
| 0.75 | 1.701332 | 1.2e−4 | 1.701664 | 3.3e−4 |
| 1.25 | 2.353489 | 9.1e−5 | 2.353236 | 2.5e−4 |
| 1.5 | 2.782541 | 1.0e−4 | 2.782278 | 2.6e−4 |
| 1.75 | 3.316862 | 4.1e−5 | 3.316772 | 9.0e−5 |
| −0.5 | 0.594111 | 1.6e−3 | 0.588696 | 5.4e−3 |
| 0.25+0.25i | 1.202587+0.218229i | 2.9e−4 | 1.203374+0.218388i | 8.0e−4 |
| 0.5+0.5i | 1.387425+0.472780i | 4.4e−4 | 1.388400+0.472330i | 1.1e−3 |
| 1.25+0.25i | 2.316433+0.383187i | 1.5e−4 | 2.316036+0.383017i | 4.3e−4 |
| 1.5+0.5i | 2.582991+0.906173i | 4.0e−4 | 2.581953+0.906814i | 1.2e−3 |

Zum Trend: Mit 5, 6 und 7 Rängen wandert der volle Wert bei x=0.5 über 1.44192 → 1.44226 → 1.44246 nach oben. Die Korrekturen fallen nur langsam (6·10⁻⁴ → 3.4·10⁻⁴ → 2.1·10⁻⁴). Ein Restbeitrag von einigen 10⁻⁴ ist also zu erwarten, und der perturbative Wert 1.44303 liegt in dieser Richtung. Ein sauberes Konvergenzgesetz gibt es mit 7 Werten nicht.

**Rangrekursion** (H₃.₅ aus N=3…9, H₄.₅ aus N=4…9; Kontrolle mit ganzzahligen Rängen ~10⁻²⁹):

| x | H₄.₅(x) voll | perturbativ (E1) | Residuum H₄.₅(x+1)−H₃.₅(H₄.₅(x)) |
|---|---|---|---:|
| 0.25 | 1.231209 | 1.231155 | −3.1e−6 |
| 0.5 | 1.468346 | 1.468269 | +6.8e−5 |
| 0.75 | 1.720546 | 1.721176 | +1.2e−4 |
| 0.25+0.25i | 1.228376+0.231486i | 1.228325+0.231293i | 6.6e−5 |
| 0.5+0.5i | 1.440525+0.472746i | 1.439519+0.472763i | 4.4e−4 |

## Was dabei auffällt

- **Die Ränge stammen aus verschiedenen Auswahlen:** Kneser bei N=4, reguläre Iteration an reellen Fixpunkten bei den ungeraden Stufen, Zwei-Bogen-Auswahl bei N=6 und N=8. An komplexen Höhen zeigt die Rangfolge deshalb Knicke. Beispiel Re H_N(0.25+0.5i) für N=7, 8, 9: 1.2451, 1.2461, 1.2439. Eine Rangfortsetzung mischt diese Auswahlen. Das begrenzt die Genauigkeit vermutlich stärker als die Zahl der Ränge.
- Die Grenze ~10⁻⁴ taucht in allen fünf Etappen auf: θ-Unschärfe der Reihe, Abstand zu Kneser bei Rang 1, Rangrekursion perturbativ und voll. Sie scheint die Grösse zu sein, bei der die Wahl der Auswahl sichtbar wird.

## Bewiesen / numerisch / offen

- **Bewiesen (elementar):** Eine Rangfortsetzung der Form Σ_k w_k(s)H_{3+k}(z) mit z-unabhängigen Gewichten ist auf senkrechten Geraden beschränkt, sobald alle H_N es sind. Die Regel aus Etappe 4 unterscheidet solche Fortsetzungen also nicht.
- **Numerisch:** alle Tabellenwerte; Motorgenauigkeit laut Motor 11 Stellen (Rang 9, 25 dps), Kontrolle der ganzzahligen Rekursion 10⁻²⁹.
- **Offen:** eine nichtperturbative Definition des Halbschritts, die eine Funktionalgleichung in x erfüllt; erst dann wird die θ-Freiheit greifbar und die Beschränktheitsregel wirksam.

## Nächster realistischer Weg

1. **Einheitliche Auswahl über alle Ränge:** jede Stufe Kneser-artig (Zwei-Bogen/Cauchy), auch an komplexen Höhen bis |Im|≈2. Der Cauchy-Löser aus Etappe 4 (`e4_kneser_eps.py`) lässt sich dafür auf beliebige Basisfunktionen f verallgemeinern: Er braucht nur f, f⁻¹ und ein Fixpunktpaar. Damit wären die Rangfolgen glatt, und eine Carlson-artige Fortsetzung im Rang hätte eine Grundlage.
2. **Gekoppeltes Problem statt Interpolation:** Unbekannt sind U=H₃.₅ auf einem Streifen und V=H₄.₅ als Kneser-Superfunktion von U (Cauchy-Löser). Gesucht ist die Lösung von V=T(U) und U·U-Konsistenz, mit der Rangfortsetzung von Etappe 5 als Startwert. Ob es eine geschlossene Gleichung für U allein gibt, ist die eigentliche mathematische Frage. Ohne sie bleibt jeder Halbrang eine Interpolation mit ~10⁻⁴ Spielraum.
3. Bis dahin sind die Kontrollzahlen H₃.₅(2,0.5)≈1.4428±3·10⁻⁴ und H₃.₅(2,1.5)≈2.7824±2·10⁻⁴.

## Reproduktion (aus diesem Ordner)

```powershell
python e5_werte.py 7 30                     # N=3..7, |Im| bis 4 (~15 s)
python e5_werte.py 9 25 0,0.25,0.5          # N=3..9, |Im|<=0.5 (~80 s)
python e5_fortsetzung.py e5_werte_N9.json   # Rangfortsetzung, Vergleich mit Etappe 1 (ergibt data/e5_fortsetzung.json; hier umbenannt in _N9)
python e5_rekursion.py 9 25                 # Rangrekursion (~80 s)
```

`e5_motor.py` setzt `engine_snapshot.gp`, `hyper.gp` und `even8.gp` der Vorrunde zu `e5_engine.gp` in diesem Ordner zusammen; die Quelldateien bleiben unverändert. Keine Kneser-Aufrufe über `gemeinsam.kneser`. Alle Läufe unter 2 Minuten. Vorhandene Dateien unverändert; keine Commits.
