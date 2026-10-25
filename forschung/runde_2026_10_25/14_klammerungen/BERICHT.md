# Gebiet 14 — Andere Klammerungen stetig fortsetzen (Runde 25.10.2026)

Start 13:22, Ende etwa 13:50 (Uhr aus `date`).

**In einfacher Sprache.** Die ausgewogene Tetration (Iteration von x ↦ x^x) und die linke Pentation (Iteration von y ↦ y^(y^x)) sehen verschieden aus.
Nimmt man aber den Logarithmus, werden beide zur selben einfachen Abbildung F(m) = m·e^m. Deren Umkehrung ist die Lambert-W-Funktion.
F hat bei 0 einen «parabolischen» Fixpunkt: F(m) = m + m² + …. Dafür gibt es eine eindeutige, natürliche Fortsetzung auf gebrochene Schritte.
Man braucht also kein Kneser-Verfahren wie bei der gewöhnlichen Tetration. Alle linken Ränge ab 5 haben am Punkt 1 dieselbe Bauart. Sie
unterscheiden sich in einer Kennzahl, dem iterativen Residuum. Dafür habe ich eine geschlossene Formel gefunden. Nach dieser Kennzahl ist keine linke
Variante mit der Tetration zur Basis e^(1/e) verwandt.

## Befunde

Schreibweise wie in `br.satz.hyperop_varianten`: x[4']n = x^(x^n), x[5']0 = x, x[5'](n+1) = (x[5']n)[4']x, allgemein y[k+1']0 = y,
y[k+1'](n+1) = (y[k+1']n)[k']y. Ausgewogen: a[4]2^t = f^[t](a), f(x) = x^x. F(m) = m e^m, W = Lambert W (Hauptast).

| # | Aussage | Status | Beleg/Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| 1 | Log-Reduktion: ln(y^y) = F(ln y) und ln(z[4']x) = ζe^{xζ}, ζ = ln z; mit m = xζ wird der Schritt der linken Pentation zu F. Also a[4]2^t = exp F^[t](ln a), x[5']t = exp(F^[t](x ln x)/x) | bewiesen | eine Zeile (unten) | für x^x wohl bekannt (Folklore); für die linke Pentation vermutlich neu |
| 2 | Brücke: x[5']t = (x^x [4] 2^t)^{1/x} (linke Pentation = Wurzel der ausgewogenen Tetration von x^x) | bewiesen (ganze t; reelle t mit der kanonischen Fortsetzung aus 5) | aus 1 | vermutlich neu |
| 3 | Wikipedia-Konvention (a₄b = a^{a^{b-1}}): a₅b = exp(F^[b-1]((a−1) ln a)/(a−1)), dieselbe Normalform F | bewiesen | wie 1 mit c = a−1 | vermutlich neu |
| 4 | Fatou-Koordinate von F: Φ(m) = −1/m + ½ ln\|m\| − m/6 + m²/16 − 19m³/540 + m⁴/48 − 41m⁵/4200 + …, alle Koeffizienten rational | bewiesen (Rekursion exakt) | `fatou_reihe.py`, 22 Koeffizienten | neu als Liste (nicht in OEIS geprüft) |
| 5 | Kanonische Fortsetzung: Für x > 1 bzw. 0 < x < 1 liegt der Startwert x ln x im abstossenden bzw. anziehenden reellen Blatt. Die Fatou-Koordinate ist dort bis auf eine Konstante eindeutig, x[5']0 = x legt die Konstante fest | hergeleitet mit Lücke (Eindeutigkeitssatz für Fatou-Koordinaten zitiert, nicht nachgeprüft) | Standardtheorie | klassisch übertragen |
| 6 | Werte: 2[4]3 = 14.654511776029386087290142357920, 2[4]√2 = 2.5519515900664733159948710551, e[4]3 = 5164.0796335205401942184734139, ½[4]3 = 0.75719198254705640265505464700; 2[5']½ = 3.3049080007536475312904140452, 2[5']3/2 = 468247.12773186285106947604232, e[5']½ = 12.540371827877451627512874303, ½[5']½ = 0.56420394844110023442988843693; F^[½](1) = 1.5134281085001618745523784596 | numerisch | Einstellung A (dps 30, n = 300, K = 12) gegen B (dps 50, n = 1500, K = 18): 22 Stellen gleich; Gruppenprobe F^[½]∘F^[½](1) − e = 3·10⁻⁴⁶ | neu |
| 7 | Gleiche Bauart aller linken Ränge: T_k(ζ,s) := ln(e^ζ[k']s) = ζ + sζ² + (A_k s² + B_k s)ζ³ + O(ζ⁴) mit A₄ = ½, B₄ = 0 und A_k = 1, B_k = k − 9/2 für k ≥ 5 (Rekursion A_{k+1} = 1, B_{k+1} = A_k + B_k). Der Schritt jedes linken Rangs ≥ 5 ist bei y = 1 parabolisch mit zweitem Koeffizienten = Basis x | bewiesen für ganze s; für reelle s hergeleitet mit Lücke (Gleichmässigkeit der asymptotischen Entwicklung der gebrochenen Iterierten in s; Standardsatz, nicht nachgeprüft) | `rang_keime.py`: exakt für s = 1…4 bei k = 5, 6, 7 (Grad-2-Polynom, 4. Wert ist Gegenprobe) | vermutlich neu |
| 8 | Iteratives Residuum β = 1 − b/a² des Schritts von Rang k+1 (Basis x): β = ½ für die linke Pentation (unabhängig von x), β = (9 − 2k)/(2x) für k ≥ 5 (Hexation −1/(2x), Rang 7: −3/(2x)) | wie 7 | aus 7; für die Hexation zur Basis 2 direkt: G = ζ + 2ζ² + 5ζ³ + 65/6 ζ⁴ …, β = −1/4 | vermutlich neu |
| 9 | Keine linke Variante zu einer reellen Basis x > 0 ist formal konjugiert zur Tetration bei η = e^{1/e} (dort β = 1/3) oder zu e^z − 1 (β = 1/3) | bewiesen für Rang 5'; für Rang ≥ 6' wie 7 | β ist eine formale Invariante (Keime z + az² + …) | vermutlich neu |
| 10 | Die formale Halbiterierte F^[½](m) = m + m²/2 + m⁴/48 − m⁵/96 + m⁶/640 + 41m⁷/5760 − … divergiert faktoriell (Gevrey-1, Rate 1/(2π)) | numerisch | 35 Koeffizienten exakt, \|g_k\|(2π)^k/(k−2)! bleibt zwischen 0.2 und 47 (`halbiterierte_reihe.py`); dasselbe für die Fatou-Reihe | wohl erwartbar (Écalle–Voronin), für diesen Keim neu gerechnet |
| 11 | Folge aus 10: Es gibt keine bei 0 analytische Halbiterierte von F. Also ist a ↦ a[4]2^{½} bei a = 1 nur C^∞ und nicht analytisch, ebenso x ↦ x[5']½ bei x = 1. Die Fortsetzungen von beiden Seiten kleben glatt, aber nicht analytisch | Vermutung (hängt an der Divergenz in 10, die nicht bewiesen ist; C^∞-Kleben ist Standard) | wie 10 | vermutlich neu als Aussage über Hyperoperations-Varianten |
| 12 | Ausgewogene Pentation P(y) = y[4]y: ln P(e^ℓ) = ℓ + ℓ³/ln 2 − ℓ⁴/(2 ln 2) + O(ℓ⁵). Der Fixpunkt 1 ist dreifach (zwei Blätter, auf beiden reellen Seiten abstossend). Allgemein (formal): Der ausgewogene Schritt von Rang k ist ζ + ζ^{k−2}/(ln 2)^{k−4} + …, die Entartung wächst also mit dem Rang | numerisch für k = 5; hergeleitet mit Lücke für k ≥ 6 (nur formal) | (lnP − ℓ)/ℓ³ = 1.43932, 1.44637 bei ℓ = ±0.005; Mittel 1.44284, 1/ln 2 = 1.44270; Steigung −0.72 = −1/(2 ln 2) passt; A = B auf 22 Stellen | vermutlich neu |
| 13 | Komplexe Fixpunkte: F hat genau die Fixpunkte m = 2πik mit Multiplikator 1 + 2πik (abstossend, k ≠ 0). Für den Schritt der linken Pentation heisst das ln y = 2πik/x. Für x^x liegen alle über y = 1, auf verschiedenen Blättern des Logarithmus | bewiesen | m e^m = m ⇔ m = 0 oder e^m = 1 | elementar |
| 14 | Linke Hexation zur Basis 2: Der Schritt ist in Log-Koordinaten G(ζ) = F(F(F(ζ)))·e^{−ζ}, ganz ohne gebrochene Iterierte. Wert: 2[6']½ = 8.4762321835075403787413216014, ln = 2.1372660330827268961694455725 | numerisch | A (dps 30, n = 200, K = 12) gegen B (dps 45, n = 800, K = 16): 15 Stellen gleich (A druckt 15); Halbgruppenprobe 1.3·10⁻³⁷; Probe ln(2[6']1) = 1024 ln 2 | neu |
| 15 | Goodstein-artige Varianten (andere Startwerte, hereditäre Basisdarstellung) | negativ (nicht bearbeitet, Zeit) | — | — |

## Beweise und Details

**1 (Log-Reduktion).** ln(y^y) = y ln y = e^ℓ ℓ = F(ℓ) mit ℓ = ln y. Die ausgewogene Rekursion a[4]2b = (a[4]b)^(a[4]b) ist also ℓ ↦ F(ℓ) auf ℓ = ln(a[4]b),
daraus a[4]2^n = exp F^[n](ln a). Linke Pentation: ln(z^{z^x}) = z^x ln z = ζe^{xζ}. Mit m = xζ gilt xζe^{xζ} = F(m), also m_{n+1} = F(m_n),
m₀ = x ln x, und x[5']n = exp(m_n/x). **2:** Für a = x^x gilt ln a = x ln x = m₀, daraus x^x[4]2^n = exp F^[n](m₀) = (x[5']n)^x. **3:** Ebenso mit z^{z^{a−1}}.

**4 (Fatou-Reihe).** Weil ln F(m) = ln m + m exakt gilt und 1/F = e^{−m}/m, wird Φ(F) − Φ = 1 zu (1 − e^{−m})/m + βm + Σ c_k(F^k − m^k) = 1. Die Ordnung m¹
gibt β = ½. Die Ordnung m^{k+1} gibt c_k linear (Leitterm k c_k). Damit sind alle c_k rational.

**5 (Rechnung).** Für m > 0: Φ(m) = lim Φ_as(W^[n](m)) + n. Für m < 0: Φ(m) = lim Φ_as(F^[n](m)) − n. Dann F^[t](m) = F^[n](Φ_as⁻¹(Φ_as(W^[n]m) + t)),
Umkehrung von Φ_as mit dem Newton-Verfahren (`kanonisch.py`). Fehler etwa c_{K+1} n^{−K−1}.

**7 (Induktion).** Annahme: T_k(ζ,s) = ζ + sζ² + (A s² + B s)ζ³. Der Schritt von Rang k+1 zur Basis y = e^ζ ist u ↦ T_k(u, y) = u + yu² + (Ay² + By)u³.
Mit y = 1 + ζ + … und der Formel für die s-fache Iterierte eines Keims u + au² + bu³, nämlich u + sau² + (sb + s(s−1)a²)u³, folgt ab u = ζ:
ζ + sζ² + (s² + (A + B)s)ζ³. Startwert T₄ = ζe^{sζ} (A = ½, B = 0). Die Formel für die s-fache Iterierte ist für ganze s Induktion; für reelle s ist sie die
formale Flussformel. Dass die kanonische gebrochene Iterierte diese Entwicklung gleichmässig in s hat, ist die benannte Lücke. Gegenprobe exakt in Bruchzahlen
(`rang_keime_ausgabe.txt`): T5(z,n) z³ = 3/2, 5, 21/2, 18; T6 = 5/2, 7, 27/2, 22; T7 = 7/2, 9, 33/2, 26. Das sind genau n² + (k − 9/2)n.

**8–9.** Für f = z + az² + bz³ ist die Fatou-Koordinate −1/(az) + β ln z + …, mit β = 1 − b/a². β ist unter z ↦ cz invariant und eine formale Konjugationsinvariante.
exp_η um e: w + w²/(2e) + w³/(6e²), also b/a² = 2/3 und β = 1/3. Für die linke Pentation ist b/a² = ½, für Rang ≥ 6' ist β < 0.

**12.** ln P(e^ℓ) = F^[s](ℓ) mit s = ℓ/ln 2 und F^[s](ℓ) = ℓ + sℓ² + (s² − s/2)ℓ³ + O(sℓ⁴). Rang 6 formal: Der Fluss von ℓ + ℓ³/ln 2 + … mit Zeit ℓ/ln 2 gibt ℓ + ℓ⁴/ln²2.

**Zum Rechenaufwand.** Alle Läufe dauerten unter 30 s. Kein Prozess blieb stehen; den ersten sympy-Versuch (zu langsam) habe ich durch Bruchzahl-Reihen ersetzt
und danach kontrolliert, dass kein eigener Prozess mehr lief.

## Dateien
- `fatou_reihe.py` → `fatou_reihe_ausgabe.txt` (β, c₁…c₂₂, Wachstum)
- `rang_keime.py` → `rang_keime_ausgabe.txt` (Befund 7, exakt)
- `kanonisch.py 30 300 12` → `kanonisch_A.txt`, `kanonisch.py 50 1500 18` → `kanonisch_B.txt` (Befunde 6, 12)
- `halbiterierte_reihe.py 36` → `halbiterierte_ausgabe.txt` (Befund 10)
- `hexation_links.py 30 200 12` → `hexation_A.txt`, `hexation_links.py 45 800 16` → `hexation_B.txt` (Befund 14)

## Vorschläge für Grapheinträge
- «Linke Pentation und ausgewogene Tetration sind dieselbe Iteration von m e^m» — satz — bewiesen
- «x[5']t = (x^x [4] 2^t)^{1/x}» — satz — bewiesen
- «Fatou-Koordinate von m e^m: −1/m + ½ ln m − m/6 + m²/16 − …» — satz — bewiesen (Koeffizienten), numerisch (Divergenz)
- «Kanonische Fortsetzung linker und ausgewogener Hyperoperationen über den parabolischen Punkt 1» — satz — hergeleitet mit Lücke
- «Iteratives Residuum der linken Ränge: ½ bei Rang 5', (9 − 2k)/(2x) bei Rang k+1 ≥ 6'» — satz — hergeleitet mit Lücke
- «Keine linke Variante ist formal konjugiert zur η-Tetration» — satz — bewiesen (Rang 5'), mit Lücke (höhere)
- «Ausgewogene Ränge entarten: Schritt von Rang k ist ζ + ζ^{k−2}/(ln 2)^{k−4}» — vermutung/satz — numerisch (k = 5), formal (k ≥ 6)
- «a[4]√2 und x[5']½ sind bei 1 glatt, aber nicht analytisch» — vermutung — Belege: Divergenz der formalen Halbiterierten
- Experiment «Werte der kanonischen linken/ausgewogenen Hyperoperationen (2[4]3, 2[5']½, 2[6']½ …)» — experiment — numerisch, 15–22 Stellen

## Offene Fragen
1. Beweis der Divergenz der formalen Halbiterierten von m e^m. Ein Weg wäre zu zeigen, dass ein Écalle–Voronin-Invariant (Hornabbildung) nicht null ist; die Hornabbildung messen.
2. Gibt es eine geschlossene Form für die Fatou-Koeffizienten c_k (OEIS-Abgleich der Zähler/Nenner steht aus)?
3. Gleichmässigkeit in s (Lücke in 7/8) sauber zitieren oder beweisen; dann ist die Residuenformel für alle linken Ränge bewiesen.
4. Ausgewogene Ränge ≥ 6 numerisch (mehrblättrige Fatou-Koordinaten).
5. Goodstein-artige Varianten (Startwerte S(0) ≠ 1, hereditäre Basis) wurden nicht bearbeitet.
6. Hängt die kanonische Fortsetzung der linken Ränge stetig in einen gebrochenen *Rang* (zwischen 5' und 6') über? Wegen der einheitlichen Bauart (Befund 7)
   wäre ein Ansatz: B_k = k − 9/2 als Funktion von k ernst nehmen.
