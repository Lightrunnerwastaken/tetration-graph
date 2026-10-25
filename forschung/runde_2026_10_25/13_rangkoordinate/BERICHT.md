# 13 Rangkoordinate Φ — Konstruktion, Additivität, halber Rang

Forschungsrunde 25.10.2026, Gebiet 13. Start 13:22, Ende 13:47 (Uhr aus `date`). Baut auf Runde 18, Bericht 07 (B5–B7, `kn.beob.rangverschiebung`, `kn.verm.rangkoordinate`) auf.

## Kurz in einfacher Sprache

Auf Leitern, die beschränkt bleiben (Basis √2 und 1.5), kann man jeder Funktion f eine Zahl Φ(f) geben, ihre „Rangkoordinate“. Man wendet den Schritt „nächste Rechenart bauen“ immer wieder an und liest hoch oben ab, auf welcher Sprosse der normalen Leiter man gelandet ist. Diese Zahl existiert numerisch auf 5–6 Stellen, und sie hängt in erster Näherung linear von einer θ-Welle ab, aber nicht genau linear. Eine θ-Welle wirkt umso stärker, je höher der Rang ist, an dem man sie einbaut, und auch schnelle Wellen (Frequenz 2, 3) verschieben Φ genauso stark wie langsame. Für den halben Rang heisst das: Die Bedingung Φ = 3.5 legt die ganze Leiter darüber fest (zwei verschiedene Ansätze laufen nach wenigen Schritten zusammen). Den halben Rang selbst legt sie aber nicht fest, denn auf Höhe 3.5 unterscheiden sich die beiden Ansätze noch in der vierten Stelle. Bei Basis 2 geht die Konstruktion so nicht, weil die Leiter dort keinen Fixpunkt hat.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | Konstruktion: Φ(g) := lim_n [s_n − n], wobei ln(p(Tⁿg) − b) = L(s_n) und L die Interpolation der kanonischen Folge r ↦ ln(p_r − b) ist (p_r = Fixpunkt von H_r, H_3 = bˣ). Normiert ist Φ(H_r) = r. Falls der Grenzwert existiert, gilt Φ(Tg) = Φ(g)+1 exakt | **bewiesen** (Abel-Gleichung, wenn der Grenzwert existiert) / **numerisch** (Existenz) | √2, θ-Störungen an H_4: s_n − n stabil auf 5–6 Stellen ab n ≈ 6 (z. B. 3.969339, 3.969334, 3.969339, 3.969345) | neu (Präzisierung von kn.verm.rangkoordinate) |
| B2 | Gradient von Φ nach den θ-Moden an H_4(√2): ∂Φ/∂(sin 2πnx) = −6.125, +8.060, −2.84; ∂Φ/∂(1−cos 2πnx) = −8.042, −11.258, +9.91 (n = 1, 2, 3) | **numerisch** (n=1,2: 4 Stellen; n=3: ±0.01, läuft noch) | symmetrische Differenz ε = ±0.005, dps 20, Grad 30. Gegenprobe mit den Daten von Runde 18 (dps 30, Grad 40): Vorhersage c(sin, 0.01) = −0.06135, gemessen −0.06134; c(1−cos, 0.01) = −0.07753, gemessen −0.07739 | neu |
| B3 | Die Gradientenkomponenten fallen mit der Frequenz n **nicht** ab. dΦ ist also keine glatte Gewichtsfunktion auf dem Kreis, sondern verteilungsartig (Abhängigkeit von Punktwerten längs der Bahn von 1) | **numerisch** (n ≤ 3) / **Vermutung** (allgemein) | B2. Ein Punktmodell nur am Fixpunkt p_4 = 1.5492 (Wert, 1. und 2. Ableitung von η) passt **nicht** (Reste bis 7.7) | neu |
| B4 | Additivität: Φ ist differenzierbar (linear in erster Ordnung), aber nicht linear. Quadratische Terme: Q_ss ≈ −1.0 (sin1), Kreuzterm sin1·cos1 ≈ +17, sin1·sin2 ≈ −1.6 (in c = Σgᵢεᵢ + Σ Qᵢⱼεᵢεⱼ) | **numerisch** (2–3 Stellen) | c(s1, 0.02) = −0.12293, c(c1, 0.02) = −0.14921, zusammen −0.25834 statt −0.27214. c(s1, ±0.1) = −0.64164 / +0.62124 (ungerader Teil −0.63144 statt linear −0.6125) | neu |
| B5 | Rangabhängigkeit: dieselbe Welle sin 2πx wirkt stärker, je höher man sie einbaut: ∂Φ/∂ε = −6.125 (H_4), −6.526 (H_5), ≤ −6.82 (H_6), ≤ −7.06 (H_7), ≤ −7.24 (H_8). Die Zuwächse passen grob zu ≈ 1.6 ln r | **numerisch** (H_4, H_5 auf 3 Stellen; ab H_6 nur Schranken, weil die Leiter zu kurz ist) / **Vermutung** (logarithmisches Wachstum) | `exp.py sqrt2 20 30 15 rank` | neu |
| B6 | Basis 1.5 (H_4 = Kneser, gp-Werte aus dem Zwischenspeicher von Runde 18): ∂Φ/∂(sin 2πx) ≈ −4.89, ∂Φ/∂(1−cos 2πx) ≈ −10.74 an H_4 | **numerisch** (2–3 Stellen, letzte Glieder ändern sich noch um ~0.01) | dps 30, Grad 40, R = 13 | neu |
| B7 | Halber Rang H_{3.5}(√2): zwei Familien g_t = (1−t)H_3 + tH_4 und H_3^{1−t}H_4^t erfüllen Φ(g_t) = 3.5 bei t* = 0.62762 bzw. 0.62300. Die Leitern Tⁿg* laufen zusammen: relativer Unterschied von p − b nach n = 0…4 Schritten 1.2e−2, 9.9e−4, 1.7e−4, 4.5e−5, 2.0e−5, danach ein Boden von ~1e−5 (Rest-Unterschied in Φ ≈ 5e−5 aus der Nullstellensuche) | **numerisch** | `half2.py` | neu |
| B8 | Was Φ **nicht** festlegt: auf Höhe 3.5 selbst unterscheiden sich die beiden Lösungen deutlich: g*(0.5) = 1.223358 gegen 1.222821, g*(1.5) = 1.592052 gegen 1.591226, Fixpunkt 1.64203 gegen 1.63921. Die Bedingung Φ = 3.5 ist **eine** reelle Gleichung, und die Restfreiheit (die Niveaufläche Φ = 3.5) ist unendlichdimensional. Zum Beispiel liegt 8.060·sin 2πx + 6.125·sin 4πx in erster Ordnung im Kern von dΦ | **numerisch** (Unterschiede) / **bewiesen** (Zählung, gegeben B1 und B2) | B2, B7 | neu |
| B9 | Struktur: Φ ist die Fatou-artige Koordinate längs der langsamen Kurve. Die Querrichtungen werden überexponentiell zusammengezogen (B7: Faktoren ≈ 12, 6, 4, 2 pro Schritt bis zum Boden). Damit ist „die Leiter der halben Ränge“ hoch oben kanonisch, auf endlicher Höhe aber nur bis auf das Blatt Φ = const bestimmt. Rückwärts (T⁻¹) bläht sich das Blatt auf, deshalb gibt es keinen kanonischen Weg nach unten | **hergeleitet mit Lücke** (Lücke: Existenz der langsamen Kurve und der stabilen Blätterung im Funktionenraum ist nicht bewiesen) | B1, B7, Runde 18 B6 | neu |
| B10 | Genauigkeit von Φ bei **halbzahligen** Abständen nur 3–4 Stellen: 4- und 6-Punkt-Interpolation der Basisfolge geben 3.49888 bzw. 3.49913 (n = 7). Der Abstand fällt nur algebraisch (1.2e−3, 6.4e−4, 3.9e−4, 2.5e−4 für n = 4…7). Bei ganzzahlnahen Abständen (θ-Störungen) ist der Fehler viel kleiner | **numerisch** / negativ (Konvergenz langsam) | `half2.py` | — |
| B11 | Basis 2: die Konstruktion geht so nicht. Die Leiter hat keinen anziehenden Fixpunkt (H_r(2,2) = 4 für alle r, unbeschränkt), also gibt es keinen Koenigs-Schritt und keine Messgrösse p_r − b. Für eine Fassung über H_r(2,½) → 1.5 fehlte die Zeit (Kneser je Rang über `hyper(r,x)` ist teuer) | **negativ** (für diese Methode) / offen | — | — |

## Details

**Konstruktion (B1).** T ist der Koenigs-Schritt aus Runde 18: (Tf)(x) = σ⁻¹(λˣσ(1)) am kleinsten Fixpunkt p ≥ 1 von f, mit λ = f′(p) < 1. Für ein Startglied g mit nominellem Rang s₀ wird für jedes n der Wert s_n mit ln(p(Tⁿg) − b) = L(s_n) gesucht. L ist die Lagrange-Interpolation (4 oder 6 Punkte) der Basisfolge ln(p_r − b), r = 3…14 bei √2. Die Folge s_n − n konvergiert, denn die Querabweichungen sterben überexponentiell (Runde 18 B6). Was bleibt, ist der Fehler der Interpolation von L. Bei ganzzahlnahen Abständen ist er klein, bei halbzahligen grösser (B10). Er fällt mit r, weil ln δ_r ≈ −k ln r − ln ln r + … in r immer glatter wird. Dass die beiden Interpolationen zusammenlaufen, ist der numerische Hinweis darauf, dass es einen kanonischen Grenzwert gibt (Fatou-artig, durch die Asymptotik der Basisfolge festgelegt). Die Abel-Gleichung Φ(Tg) = Φ(g)+1 folgt direkt aus der Definition, denn die Folge für Tg ist die um eins verschobene Folge für g.

**Gradient (B2, B3).** Gestört wird H_4 ↦ H_4∘ψ mit ψ = x + εη, η 1-periodisch mit η(0) = 0. Die ganzen Werte bleiben also gleich. Die Komponenten für n = 1, 2, 3 haben die Beträge 6–11 und fallen nicht ab. Der Grund ist vermutlich, dass Tf von f längs der ganzen Bahn 1, b, f(b), … → p abhängt. Das sind Punkte, deren Phasen sich bei p_4 mod 1 = 0.549 häufen. Ein Modell, das nur η, η′, η″ am Fixpunkt benutzt, passt nicht (`punktmodell.py`).

**Halber Rang (B7, B8).** Bei √2 ist H_3 = bˣ (Fixpunkt 2) und H_4 die Kneser- bzw. Koenigs-Tetration (Fixpunkt 1.5492). Beide Familien verbinden H_3 und H_4 bei festem g(0) = 1 und g(1) = b. Φ(g_t) − 3.5 bei t = 0.3, 0.5, 0.7: lin −0.29635, −0.12755, +0.07922; geo −0.29320, −0.12296, +0.08406 (6-Punkt). t* stammt aus quadratischer Interpolation. Die Kandidaten für H_{3.5}(√2, 0.5) liegen bei 1.2228–1.2234 und hängen von der Familie ab. Kanonisch ist erst die Leiter darüber.

| n (Rang 3.5+n) | p(Tⁿg*_lin) | p_lin − p_geo | relativ zu p − b |
|---:|---:|---:|---:|
| 0 | 1.6420270 | 2.8e−3 | 1.2e−2 |
| 1 | 1.5080141 | 9.3e−5 | 9.9e−4 |
| 2 | 1.4694459 | 9.3e−6 | 1.7e−4 |
| 3 | 1.4517322 | 1.7e−6 | 4.5e−5 |
| 4 | 1.4418155 | 5.6e−7 | 2.0e−5 |
| 6 | 1.4313831 | 1.9e−7 | 1.1e−5 |
| 10 | 1.4230903 | 6.1e−8 | 6.9e−6 |

**Was Φ für den halben Rang festlegt.** Bezeichne mit 𝒜 die Menge der Funktionen „zwischen“ H_3 und H_4 mit g(0) = 1, g(1) = b, auf der T definiert ist. Φ: 𝒜 → ℝ ist (numerisch) glatt mit nichtverschwindendem Gradienten, und Φ⁻¹(3.5) ist eine Hyperfläche. Alle Punkte dieser Hyperfläche haben dieselbe Zukunft bis auf überexponentiell kleine Fehler. Festgelegt wird also die Halbrang-Leiter H_{r+½} für r → ∞, als Punkt auf der langsamen Kurve mit Φ = r + ½. Frei bleibt die ganze Hyperfläche, und für eine Auswahl auf endlicher Höhe braucht man eine zweite, unabhängige Bedingung. Kandidaten dafür wären Glattheit in r über alle Ränge oder eine Gamma-artige Interpolation der Querkontraktion (Runde 18, Frage 5).

## Vorschläge für Grapheinträge

- „Rangkoordinate Φ: Konstruktion über die Basisfolge ln(p_r − b), Abel-Gleichung exakt“ — satz/beobachtung — bewiesen (Abel) / numerisch (Existenz, 5–6 Stellen)
- „Gradient der Rangkoordinate nach θ-Moden fällt mit der Frequenz nicht ab“ — beobachtung — numerisch (n ≤ 3, b = √2)
- „Rangkoordinate ist differenzierbar, nicht linear (quadratische Kreuzterme)“ — beobachtung — numerisch
- „Dieselbe θ-Welle verschiebt Φ stärker, je höher der Rang (≈ log r)“ — vermutung — numerisch, Schranken
- „Φ = 3.5 legt die Halbrang-Leiter oben fest, nicht H_{3.5} selbst (√2: g*(0.5) = 1.2228 bzw. 1.2234)“ — beobachtung — numerisch
- „Rangkoordinate bei b = 2 mit Fixpunkt-Messgrösse nicht konstruierbar“ — gescheitert — negativ

## Offene Fragen

1. Beweis der Existenz von Φ: eine langsame Kurve und eine stabile Blätterung des Koenigs-Schritts im Funktionenraum.
2. Liegt dΦ als Verteilung auf dem Kreis an den Phasen der Bahn 1, b, H_4(b), …? Das liesse sich mit Moden n bis 10 und höherem Grad testen.
3. Wächst |∂Φ/∂ε| an H_r wirklich wie log r? Dafür braucht man längere Leitern (R ≥ 25).
4. Gibt es eine natürliche zweite Bedingung, die auf der Hyperfläche Φ = 3.5 einen Punkt auswählt?
5. Basis 2: eine Rangkoordinate über H_r(2,½) → 1.5 − 0.62/r², mit Kneser je Rang aus dem gp-Motor.

## Reproduktion (Ordner 13_rangkoordinate, Python = C:/Users/janis/miniconda3/python.exe)

```
python exp.py sqrt2 20 30 15 base        # Basisleiter r = 3..14, ~31 s
python exp.py sqrt2 20 30 15 grad        # B2 (12 Läufe, ~5 min, wiederaufnehmbar über cache/)
python exp.py sqrt2 20 30 15 add         # B4
python exp.py sqrt2 20 30 15 rank        # B5
python exp.py sqrt2 20 30 15 half        # Familien t = 0.3, 0.5, 0.7
python half2.py                          # B7, B8, B10
python exp.py 1.5 30 40 13 base grad1    # B6 (Kneser-Werte aus cache/sexp_1.5_30.json, Kopie aus Runde 18)
python punktmodell.py                    # B3 (negativer Test)
```
Die Rohdaten liegen in `cache/` (einzelne Läufe) und `daten/` (Auswertungen). Der gp-Motor wurde nicht gestartet, und es gab keine Commits.
