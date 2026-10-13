# Etappe 4: Welches θ hat Kneser? — Ergebnis

13.10.2026 · gleicher Ordner · neue Dateien `e4_*.py`, Daten `data/e4_*.json`

## Ergebnis zuerst

1. **Kneser ist eine Resummation der perturbativen Reihe.** Für die ganze Familie f_ε(x)=x+1+ε(2ˣ−x−1) habe ich die Kneser-artige Superfunktion selbst berechnet: Cauchy-Integral auf einem Streifen, Grenzwerte L_ε und L̄_ε für Im x→±∞, normiert F(0)=1. Bei ε=1 stimmt sie mit der GP-Kneser-Tetration auf **2.9·10⁻¹¹** überein. Für ε=0.6…0.36 ist der Abstand zur Borel–Padé-Summe der Reihe **exponentiell klein**, D₀≈C(x)·e^{−B/ε} mit B≈7.6 (je nach Höhe 7.0–8.0). Er ist also kleiner als jede Potenz von ε (bei ε=0.36 etwa 5·10⁻¹⁰).
2. **Der Unterschied ist genau eine θ-Mode.** D₀ hat die Form a·sin2πx+b(1−cos2πx): 1-periodisch (D₀(−0.25)≈D₀(0.75)) und null an ganzen Höhen. Er wächst senkrecht wie e^{2π|y|}: bei ε=0.5 ist |D₀| = 1.9·10⁻⁶, 7.5·10⁻⁶, 3.2·10⁻⁵ für y=0.5, 0.75, 1, das sind die Faktoren 4.3 und 18 gegenüber e^{π/2}=4.8 und e^{π}=23. Kneser selbst bleibt beschränkt (Grenzwert L). **Unter allen Summen der Reihe ist Kneser also die eine, die in senkrechter Richtung beschränkt bleibt.** Das ist der Paulsen–Cowgill-Mechanismus (kn.satz.paulsen_cowgill), hier numerisch auf der Ebene der Resummation sichtbar gemacht.
3. **Kneser ist nicht die Borel–Padé-Summe längs der reellen Achse.** Der Abstand ist 6- bis 4000-mal grösser als das kleinste Glied der Reihe (je nach Höhe und ε). Ob Kneser die *Median*-Summe ist, lässt sich nicht entscheiden: Die massgebliche Borel-Singularität liegt bei t≈7.6±0.35i, fast auf der positiven Achse, und Padé löst sie nicht auf.
4. **Auswahlregel für den Halbrang:** Die übertragbare Form ist „beschränkt auf senkrechten Geraden“ (bzw. ein Grenzwert für Im x→±∞). Anwenden konnte ich sie **nicht**. Dafür bräuchte man einen nichtperturbativen Halbrang an komplexen Höhen, und der Halbrang-Jet existiert nur bis Ordnung 7 (Etappe 1). **Neue H₃.₅-Werte gibt es daher nicht**; die Werte aus Etappe 1 (±10⁻⁴) bleiben der Stand.

## 1. Die θ-Projektion bei Rang 1

**Verfahren.** Für jedes ε wird die Kneser-artige Lösung F_ε mit dem eigenen Cauchy-Verfahren berechnet (`e4_kneser_eps.py`: Rechteck [−1,1]×[−40,40], 700 Knoten, Residuum ≤10⁻¹³). Daneben steht S₀, die Borel–Padé–Laplace-Summe der 84 Ordnungen längs der reellen Achse (`e4_rang1_schnell.py`, `e4_scan.py`). Gemessen wird D₀=F_ε−S₀ zusammen mit dem kleinsten Glied τ. Als stabil gelten die Höhen, an denen Padé [24/24] und [27/27] auf ≤10 % übereinstimmen.

| ε | D₀(0.25) | D₀(0.75) | D₀(−0.25) | D₀(1.25) | τ(0.75) |
|---:|---:|---:|---:|---:|---:|
| 0.60 | +8.4e−7 | −3.9e−6 | −3.3e−6 | +1.1e−6 | 5.1e−8 |
| 0.50 | +7.2e−8 | −2.5e−7 | −2.2e−7 | +9.2e−8 | 4.0e−9 |
| 0.45 | +1.5e−8 | −4.2e−8 | −3.7e−8 | +1.9e−8 | 2.3e−10 |
| 0.40 | +2.2e−9 | −4.6e−9 | −4.1e−9 | +2.8e−9 | 4.2e−12 |
| 0.36 | +3.5e−10 | −5.1e−10 | −4.7e−10 | +4.3e−10 | 1.2e−13 |

Fit ln|D₀| gegen 1/ε (Restfehler ≤0.09): B=7.0 (x=0.25), 8.05 (0.75), 7.8 (0.875), 7.97 (−0.25), 7.1 (1.25). Die Abhängigkeit von x passt zu einem komplexen Exponenten. Mit a·e^{B/ε}≈0.57–0.73 und b·e^{B/ε} von −0.11 (ε=0.36) bis −0.47 (ε=0.6) dreht die Phase langsam; daraus ergibt sich grob A≈7.6±0.35i. Das ist eine Schätzung, kein Fit mit Fehlerbalken. Die Padé-Pole der Borel-Transformierten liegen in derselben Gegend (7.5±1.6i, 8.4±0.8i, je nach Ordnung), aber nicht stabil.

**Hochgerechnet auf ε=1:** D₀≈+8·10⁻⁵ (x=0.25), −8·10⁻⁴ (0.75), −3·10⁻⁴ (0.875). Das passt zur Unschärfe von 10⁻⁴–10⁻³ aus Etappe 2, wo Borel–Padé bei ε=1 nicht stabil war. Es ist eine Extrapolation über einen Faktor 1.7 in 1/ε und kein Messwert.

**Was „Kneser = Resummation“ hier heisst:** Die Differenz liegt jenseits aller Ordnungen (exponentiell klein in 1/ε) und hat die Form eines periodischen θ. Das ist numerisch über ε∈[0.36, 0.6] belegt, nicht bewiesen. Die Lösungen bei ε<1 sind Knesers Konstruktion übertragen auf f_ε (Grenzwerte gegen L_ε, L̄_ε); für ε<0.6 ist L_ε anziehend (|f′(L_ε)|<1). Die Eindeutigkeit dieser Wahl ist dort nicht bewiesen.

## 2. Kneser-Bedingungen und der Halbrang

| Kneser-Kennzeichen | bei Rang 1 | Sinn für den Halbrang? |
|---|---|---|
| Funktionalgleichung F(z+1)=f(F(z)), F(0)=1 | Grundlage | Nein, direkt. H₃.₅ hat keine eigene Gleichung in x; die Auswahl ist eine Operatorgleichung U∘U=T. Indirekt ja: H₄.₅ ist eine Superfunktion von H₃.₅, und deren θ wird von einer Kneser-Bedingung für H₃.₅ festgelegt. Das wählt T, nicht U. |
| Grenzwerte L, L̄ für Im z→±∞ (Paulsen–Cowgill) | zwingt θ zur Konstanz | Nur, wenn man weiss, wohin H₃.₅(2,x+iy) streben soll. Kandidat: ein Fixpunkt von H₃.₅ selbst (bei H₃ und H₄ gibt es komplexe Fixpunkte). Nicht hergeleitet. |
| Beschränktheit auf senkrechten Geraden | gilt für 2ˣ (|2^{x+iy}|=2ˣ) und für Kneser | **Ja**, das ist die natürliche Interpolation. Sie schliesst jede θ-Mode e^{±2πix} aus (Lemma A, Etappe 1) und wählt bei Rang 1 nachweislich Kneser unter den Summen (oben, Punkt 2). |
| Injektivität/Anfangsgebiet (Trappmann–Kouznetsov) | Eindeutigkeit der Abel-Funktion | Erst sinnvoll, wenn H₃.₅ als Funktion mit Fixpunktpaar vorliegt. |

**Vorschlag einer Regel:** H₃.₅(2,·) ist die auf senkrechten Geraden im Streifen 0≤Re x≤1 beschränkte Resummation des Halbrang-Jets.

Die Regel ist widerspruchsfrei mit Rang 0 (2ˣ) und Rang 1 (Kneser). Zwei Lücken bleiben offen:
- Der Halbrang-Jet existiert nur bis Ordnung 7, es gibt also nichts zu resummieren, was die θ-Grösse auflöst.
- Dass die beschränkte Wahl bei U auch U∘U=T_Kneser liefert, ist nicht gezeigt.

## 3. Anwendung auf H₃.₅

Nicht durchgeführt. Für eine Anwendung bräuchte es eine der beiden Grundlagen:
- den Halbrang-Jet in hoher Ordnung, was nach Etappe 3 an den Raten ρ_j scheitert;
- einen nichtperturbativen Halbrang an komplexen Höhen, etwa die Rangfortsetzung der vollen H_N(2,x+iy), N=3…9 oder mehr, die bisher nicht berechnet sind.

Die Werte aus Etappe 1 bleiben: H₃.₅(2,0.25/0.5/0.75/1.5) = 1.21246 / 1.44302 / 1.70164 / 2.78231, jeweils ±10⁻⁴.

## Bewiesen / numerisch / offen

- **Bewiesen** (übernommen): periodische Funktionen mit polynomialem vertikalem Wachstum sind konstant; daher ist eine vertikal beschränkte Resummation, wenn es sie gibt, bis auf eine Konstante eindeutig. Die Konstante legt F(0)=1 fest.
- **Numerisch:** Kneser-Löser gegen GP 2.9·10⁻¹¹; D₀∝e^{−B/ε} mit B≈7–8; θ-Form und senkrechtes Wachstum e^{2π|y|}; Kneser ≠ Borel–Padé-Summe der reellen Achse.
- **Offen:** Ist Kneser die Median-Summe? Lage und Typ der Borel-Singularität bei ≈7.6. Gilt die Regel „beschränkt“ beim Halbrang verträglich mit U∘U=T?

## Reproduktion (aus diesem Ordner)

```powershell
python e4_kneser_raster.py                              # 1 Kneser-Aufruf (47 Punkte), ~5 s
python e4_rang1_schnell.py --J 100 --prec 2200 --name J100 --budget 1000   # erreicht 84 Ordnungen (~18 min); Datei data/e4_rang1_J100.json danach in e4_rang1_J84.json umbenennen
python e4_kneser_eps.py --eps 1 --A 30 --N 400           # Kneser-Löser, Kontrolle gegen GP (~25 s)
python e4_scan.py 0.6,0.5,0.45,0.4,0.36 40 700          # D0 gegen eps (~13 min)
python e4_fit.py                                        # B, Form, Hochrechnung
python e4_senkrecht.py                                  # senkrechtes Wachstum (~4 min)
python e4_richtungen.py 0.4                             # seitliche Summen (Padé-Pole, Sprünge)
python e4_theta.py e4_rang1_J84.json 0.3                # seitliche Summen bei eps=1
python e4_stokes_folge.py                               # theta-Zerlegung der späten Glieder
```

Kneser-Aufrufe: genau einer (`e4_kneser_raster.py`). `gemeinsam.kneser` hängt dabei wie in Etappe 1 fehlende Punkte selbsttätig an `experimente/daten/kneser_cache.json` an. Alle Läufe unter 20 Minuten. Vorhandene Dateien unverändert; keine Commits.
