# Runde 18.10.2026 · 02 · Knesers θ-Welle bei höheren Rängen

## Kurz gesagt (einfache Sprache)

Bei der Tetration (Rang 4) kennt man die „Welle“ θ: Knesers Lösung ist die einfache Lösung am komplexen Fixpunkt, verschoben um eine kleine periodische Korrektur θ. θ hat scharfe Stellen genau bei den ganzen Zahlen.

Ich habe dasselbe für die höheren Rechenarten untersucht, mit dem gp-Motor auf 30–40 Stellen und mit dem Cauchy-Löser der Etappe 6 (doppelte Genauigkeit). Rechenbasis ist überall 2.

1. **Hexation (Rang 6) und Oktation (Rang 8).** Die Welle hat ihre scharfe Stelle **nicht bei den ganzen Zahlen**. Sie liegt an genau einer Stelle x* pro Periode, und zwar dort, wo die Hexation den reellen Fixpunkt p = −1.74391 von sexp erreicht. Die Hexation trifft diesen Wert bei x* = −2.7650220505…, die Oktation den Fixpunkt p₆ = −3.69941 von H₆ bei x₈* = −4.71388. Auch die **Art** der Stelle ist anders als bei der Tetration. Bei der Tetration wird θ dort unendlich, nur sehr langsam (superlogarithmisch). Bei der Hexation bleibt θ dort endlich und hat einen **Wurzelknick**: θ ≈ θ* + K·(z − x*)^α mit α = ln 2 / ln sexp′(p) = 0.4218. Deshalb fallen die Fourier-Koeffizienten wie eine **Potenz**: |a_n| ≈ 0.0634·n^(−1.4218). Bei der Tetration fallen sie dagegen wie 1/(n ln n …).
2. **Kneser-artige Pentation (Rang 5).** sexp hat nur einen *anziehenden* komplexen Fixpunkt L₄. Für die Welle bezogen auf diesen Punkt gilt: Ihre Fourier-Koeffizienten **wachsen**, etwa wie e^{2π·0.27·n}. Die Welle reicht also gar nicht bis zur reellen Achse. Sie endet auf einer fraktal aussehenden Grenze in der Höhe ≈ 0.27, dem Rand des Einzugsgebiets von L₄. Ob abstossender oder anziehender Fixpunkt, entscheidet also, ob es überhaupt eine reelle θ-Darstellung gibt.
3. **Welle zwischen regulärer und Kneser-artiger Pentation** (Amplitude ~3·10⁻⁶). Ihre Fourier-Koeffizienten fallen genau exponentiell, |c_n| ≈ 0.11·e^{−2π·1.6902·n}/n. Die Rate kommt von einem einzigen Verzweigungspunkt: der Stelle z = −0.61815 + 1.69019i, an der die Kneser-artige Pentation den Fixpunkt L₃ von 2^z annimmt. Damit ist auch erklärt, warum die Welle so klein ist.
4. **Basisunabhängigkeit.** |a_n(e)| ≈ |a_n(2)| ist kein Gesetz. Für n ≤ 4096 liegen die Werte zwar auffallend nahe beieinander: Über alle Basen von 1.45 bis 1000 ändert sich |a_1| nur um 5 %. Für astronomisch grosse n schwankt das Verhältnis aber zwischen 0.5 und 3.3, und das liegt an der Selbstähnlichkeit.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| B1 | Rang 6, Basis 2: Die Motor-Hexation ist Kneser an einem **abstossenden** komplexen Fixpunkt kp = −1.992948435932 + 1.445128496639i der regulären Pentation, mit Multiplikator μ = 1.368932747340 + 2.488144667623i. Die Pentation hat ausserdem einen **anziehenden** Fixpunkt 1.793730343239 + 1.573630243512i (μ = −0.02631 + 0.60814i). | numerisch | Newton auf der Taylorreihe, Fixpunkt auf 10⁻²⁸; μ stimmt mit `klam` des Motors überein | Der anziehende Fixpunkt und die Wahl sind in den Unterlagen nicht vermerkt |
| B2 | θ₆(z) = h(hex z) − z (h reguläre Abel-Funktion an kp) hat **keine negativen Moden** und ist genau bei **x* + ℤ** singulär, x* = ihex(p) = −2.76502205052560345111. Ganze Zahlen sind reguläre Punkte. | numerisch, mit Modell (B4) | a₋ₙ < 2.4·10⁻⁴⁹ (dps 40); arg(a_n/a_{n−1}) → −2πx* (mod 2π) = −1.47641, gemessen −1.476402 (n = 128); zwei Läufe (dps 40, y = 0.04; dps 30, y = 0.05, beide N = 512) stimmen auf 10⁻²⁹ (n ≤ 16) bzw. 2·10⁻¹³ (n = 128) überein | neu |
| B3 | **Lokale Form bei x\*:** θ₆(z) = θ* + K(z−x*)^α + O(\|z−x*\|^{2α}) + O(z−x*), mit α = ln b/ln λ₄ = 0.42182244029502478293 (λ₄ = sexp′(p) = 5.17180), θ* = 2 − x* + h(psing) = 5.17719763470 + 1.23143667244i und K = h′(psing)·ipent′(−2)·A₆^α·b^{−pcent}/sexp′(−1) = 0.2746753 + 0.4215492i. θ₆ bleibt also bei x* **endlich** (Hölder-stetig mit Exponent α), ist dort aber nicht analytisch. | hergeleitet mit Lücke (Zweigwahl, wie bei kn.satz.theta_singularitaet); numerisch bestätigt | θ(x*+iy) − θ* fällt wie y^α (Verhältnis pro Halbierung 1.237 → 1.299, Ziel 2^α = 1.340); Rest nach dem K-Glied fällt schneller als y^{2α} | neu |
| B4 | **Koeffizientenformel Rang 6:** a_n = e^{−2πinx*}·b_n·(1 + O(1/n)) mit b_n = [Qⁿ] h((c₆ + log(1−Q))/ln λ₄), c₆ = ln(A₆/2π) + iπ/2 − pcent·ln λ₄ (A₆ = hex′(x*) = 1.0325027039). Mit δ₁ = ½ + iκ₆/2π ist der Rest O(1/n²). | numerisch (Struktur hergeleitet) | n·Rest₁ = 0.33 → 0.309, n²·Rest₂ = 0.50 → 0.394 (n = 8…160); Rest₂ = 1.5·10⁻⁵ bei n = 160 | neu |
| B5 | **Asymptotik Rang 6:** a_n ≈ C₆·e^{−2πinx*}·n^{−1−α}, C₆ = K(i/2π)^α/Γ(−α) = 0.0053906 − 0.0631856i (\|C₆\| = 0.06342). Die Koeffizienten fallen also wie eine Potenz mit Exponent 1.4218, ohne iterierte Logarithmen. | hergeleitet aus B3 (Übertragungssatz für (1−Q)^α, Flajolet–Odlyzko) | a_n e^{2πinx*}/(C₆n^{−1−α}) = 1.050, 1.043, 1.035, 1.028, 1.021, 1.016 (n = 5…160); lokaler Exponent 1.4319 → 1.4289 | neu |
| B6 | **Rang 8:** θ₈ (Motor-Oktation, Kneser am abstossenden Fixpunkt kp₈ = −4.294772 + 1.283451i von H₇, μ = 3.41616 + 3.22021i) ist genau bei x₈* + ℤ singulär, x₈* = H₈⁻¹(p₆) = −4.7138756, p₆ = −3.6994078 = reeller Fixpunkt von H₆. Exponent vorhergesagt 1 + α₈, α₈ = ln λ₄/ln λ₆ = 0.730137 (λ₆ = 9.493087). | numerisch (Ort) / Vermutung mit Beleg (Exponent) | arg(a_n/a_{n−1}) = −1.8253 (n=2) … −1.798224 (n = 24), Ziel −2π·frac(x₈*) = −1.797773; lokaler Exponent 1.702 (n 6→12), 1.713 (12→24), steigt Richtung 1.730; negative Moden 2·10⁻¹⁰ (dps 20) | neu |
| B7 | **Gesetz für gerade Ränge r ≥ 4** (Motor-Familie: gerade Ränge Kneser am abstossenden komplexen Fixpunkt, ungerade Ränge regulär am reellen Fixpunkt p_{r−1} von H_{r−1}): θ_r ist genau bei x_r* + ℤ singulär, wobei H_r(x_r*) = p_{r−2}. Das ist der Wert, gegen den H_{r−1} für Re z → −∞ läuft. Lokal gilt θ_r ≈ θ* + K(z−x_r*)^{α_r} mit α_r = ln λ_{r−4}/ln λ_{r−2} und λ_k = H_k′(p_k). Dabei ist λ₂ = b (H₂(x) = b·x am Fixpunkt 0) und λ₀ = 1 (parabolisch). r = 4 ist der Grenzfall α = 0: x₄* = −1, der Exponent ist genau 1, dazu kommen iterierte Log-Korrekturen. | Vermutung, für r = 4, 6 hergeleitet und belegt, für r = 8 der Ort belegt | B2–B6; kn.satz.theta_koeffizienten | neu |
| B8 | **Kneser-artige Pentation (Rang 5, e6-Löser):** sexp₂ hat keinen für slog anziehenden komplexen Fixpunkt (9 Startwerte → p). Der Kneser-Fixpunkt L₄ = 1.469818 + 1.664342i ist **anziehend** (μ = 0.154239 + 0.690508i, \|μ\| = 0.7075). Für θ₅ = log_μ σ(H₅) − z (Koenigs-σ) gilt: Die Koeffizienten **wachsen**, \|a_1…a_10\| = 2.88, 9.42, 27.2, 154, 705, 3357, 1.42·10⁴, 5.25·10⁴, 1.76·10⁵, 6.53·10⁵. Die Wachstumsrate ist (ln\|a_n\| − ln\|a_{n−4}\|)/8π = 0.27–0.28 (n = 14–22). Die θ-Reihe konvergiert nur für Im z > y_c ≈ 0.27, θ₅ reicht nicht bis ℝ. | numerisch | Linien y = 0.30/0.35/0.5/0.7 stimmen überein (n ≤ 12 auf 2·10⁻⁴); keine negativen Moden (≤ 7·10⁻¹⁰) | neu |
| B9 | Die Grenze ist der Rand der Menge {z : H₅(z) im Einzugsgebiet von L₄ unter sexp}. Diese Menge ist 1-periodisch, weil sexp^n(H₅(z)) = H₅(z+n). Unterhalb der Grenze liegen fraktal durchmischt Punkte im Becken von L̄₄ und Punkte, deren Bahn entkommt. Der höchste Punkt ausserhalb liegt bei y ≈ 0.26 (Raster 0.02) und passt zur Rate in B8. **Vermutung:** Die Grenze ist eine natürliche Grenze von θ₅. | numerisch (Karte) / Vermutung | `s4_becken5.py` (100×25-Raster); Teilargument in Abschnitt 3 | neu |
| B10 | **Welle ψ = P⁻¹∘K₅ − z** zwischen regulärer (P) und Kneser-artiger Pentation (K₅): ψ ist auf ℝ analytisch, c₀ = 3.4572·10⁻⁶, c_{−n} = conj(c_n). Es gilt c_{−n} = (κ_n/n)·e^{2πinz_L3} mit langsam veränderlichem κ_n: \|κ_n\| = 0.1145, 0.1085, 0.1072 … 0.1211 (n = 1…14), arg κ_n = 1.647 … 1.881. Dabei ist z_L3 = −0.6181507204932638 + 1.6901945137048313i die Stelle mit **K₅(z_L3) = L₃** (Fixpunkt von 2^z, Verzweigungspunkt von slog). Der Abfall ist also genau e^{−2π·1.6902·n}/n. Typ: logarithmische Singularität mit superlog-artiger Korrektur. | numerisch (Ort auf ≥ 4 Stellen aus den Quotienten; Newton direkt auf 10⁻¹⁵) | Linien y = 0.8/1.2/1.5/1.6 stimmen auf 10⁻¹¹ (n=1) bis 4·10⁻⁴ (n=14) überein; Reihe bricht oberhalb y = 1.69 zusammen (y = 1.7: c₀ springt) | neu; erklärt die Amplitude 10⁻⁵ aus kn.beob.einheitliche_raenge |
| B11 | Rangskalierung der Welle aus B10: Die Amplitude ist ≈ C·e^{−2πd_r}, wobei d_r die Höhe des nächsten Verzweigungspunkts ist. Als grobe Schätzung dient die halbe Periode der regulären Lösung, π/ln λ_{r−1}: Rang 5 ergibt e^{−2π²/ln λ₄} = 6.1·10⁻⁶ (gemessen: \|c₁\| = 2.8·10⁻⁶), Rang 7 ergibt e^{−2π²/ln λ₆} = 1.6·10⁻⁴. | Vermutung (ein Datenpunkt plus Strukturargument) | passt zur Etappe-6-Schätzung „Faktor 4–6 pro Rang“ | – |
| B12 | \|a_n(e)\|/\|a_n(2)\| = 1.0065 (n=1) … 1.0050 (n=4096) aus den 44-stelligen Daten. Nach dem führenden Term liegt das Verhältnis bei ln ln n ≤ 3 nahe 1, schwankt dann aber: 0.83 (ln ln n = 10), 1.63 (15), 0.49 (10⁵), 2.1 (10⁶), 3.3 (10⁷), dazu Spitzen, wenn ln ln ln n einem Turmwert nahekommt. **Kein Gesetz.** | numerisch / aus kn.satz.theta_koeffizienten | `d1_universal.py` | beantwortet die Frage der Etappe 2 |
| B13 | Quasi-Universalität für kleine n: \|a_1(b)\| = 0.088953 (b = 1.45), 0.089238 (1.5), 0.090060 (2), 0.090647 (e), 0.091335 (4), 0.091955 (10), 0.092910 (100), 0.093248 (1000). Der Verlauf ist monoton und sättigt. Für n = 64 liegen die Werte zwischen 5.12·10⁻⁴ und 5.67·10⁻⁴, und arg a_n ≈ −π/2 + O(0.1) für alle b. | numerisch | dps 30, y = 0.1, N = 256; Basis 2 und e gegen die 44-stelligen Etappe-2-Werte auf ≥ 8 Stellen | Beobachtung neu; Erklärung offen |

## 1. Rang 6: Herleitung und Zahlen

**Aufbau.** f = P ist die reguläre Pentation des Motors (`pentx`), am reellen Fixpunkt p = −1.7439091761327370 von sexp mit λ₄ = sexp′(p) = 5.1717972635. Der Motor baut die Hexation `hexx` nach Kneser am Fixpunkt kp von P (`knhexinit`). kp ist für P abstossend (|μ| = 2.84) und für den Hauptzweig von ipent anziehend. Genau wie bei der Tetration (exp abstossend an L) ist deshalb h(w) = n + log σ(ipentⁿ w)/log μ eindeutig und rechenbar, mit `lib6.py`.

**Warum x\* und nicht die ganzen Zahlen.** Für die Tetration ist der Ausnahmewert 0, weil exp den Wert 0 als Grenzwert für Re → −∞ hat und sexp(−1) = 0. Eine Stufe höher übernimmt p diese Rolle: P(z) = Ψ_p(λ₄^z) → p für Re z → −∞, gleichmässig in Im z. Am Punkt x* = ihex(p) erreicht hex genau diesen Wert. Die Werte hex(k) = 0, 1, 2, … sind für h keine Ausnahmewerte. Daher ist θ₆ an den ganzen Zahlen analytisch.

**Lokale Form (B3).** Es gilt θ₆(z) = h(hex z) − z = 1 + h(V) − z mit V = ipent(hex z) = hex(z−1). Für z → x* gilt V = [log σ_p(hex z − p)]/ln λ₄ − pcent → −∞ (Im V → π/(2 ln λ₄)). Der wesentliche Unterschied zur Tetration: slog hat für V → −∞ einen **endlichen** Grenzwert, slog V = −2 + b^V/sexp′(−1) + …, weil sexp bei −2 logarithmisch singulär ist. Also gilt

ipent(V) = 1 + ipent(slog V) = psing + ipent′(−2)·b^V/sexp′(−1) + O(b^{2V}),  h(V) = 1 + h(psing) + h′(psing)·ipent′(−2)·b^V/sexp′(−1) + …

und b^V = (A₆(z−x*))^{ln b/ln λ₄}·b^{−pcent}·(1 + O(z−x*)). Daraus folgt B3. psing = ipent(−2) + 1 = −1.97859584473632 + 1.91185121532145i ist die Singularität der Pentation; der Motor gibt ihn als `psing` aus. Es gilt **K ≠ 0 streng**, denn h′ = 1/F′(h) ≠ 0, ipent′(−2) = 1/P′(·) ≠ 0 und A₆ ≠ 0. Alle drei Ableitungen sind nullstellenfrei nach dem S1-Argument aus kn.satz.theta_singularitaet (sexp′ ≠ 0 ⇒ Ψ_p′ ≠ 0 ⇒ P′ ≠ 0 ⇒ Ψ_kp′ ≠ 0). **Lücke:** h(w) = 1 + h(ipent w) und der Hauptzweig von ipent müssen längs hex(x* + iε) die analytische Fortsetzung sein. Das ist numerisch bestätigt (B4/B5 treffen die Koeffizienten), aber nicht allgemein bewiesen. Dieselbe Lücke gibt es bei Rang 4.

**Koeffizienten (B4/B5).** Q = q·e^{−2πix*}. Das Modell b_n ist der Cauchy-Koeffizient von h((c₆ + log(1−Q))/ln λ₄) auf |Q| = e^{−2π·0.05}, M = 512 (`s2_modell6.py`):

| n | \|a_n\| | Rest₁·n | Rest₂·n² | a_n e^{2πinx*}/(C₆n^{−1−α}) |
|---|---|---|---|---|
| 1 | 0.065905493 | 1.96 | 4.6 | – |
| 8 | 0.0034474438 | 0.329 | 0.499 | – |
| 16 | 0.0012768626 | 0.322 | 0.441 | – |
| 64 | 0.00017543917 | 0.313 | 0.405 | 1.028 (n=40), 1.021 (n=80) |
| 160 | 4.7359561·10⁻⁵ | 0.309 | 0.394 | 1.0164 |

Weitere Werte: a₀ = 5.16210149901427895 + 1.42703253050539489i, κ₆ = 0.335274554885, pcent = 2.07503075727351, ln λ₄ = 1.64322026129087539. Der reine führende Term −h′(v_n)/(n ln λ₄) trifft hier **nicht** (Verhältnis 1.63 bei n = 128, `s2b_fuehrend6.py`). Das ist zu erwarten, weil h′(v) ∝ b^v nicht langsam veränderlich ist. Richtig ist der Γ-Faktor aus B5, mit Grenzwert Γ(1−α) ≈ 1.55.

## 2. Rang 8 (B6) und das allgemeine Gesetz (B7)

Für die Oktation ist das Argument dasselbe, eine Stufe höher. H₇ ist die reguläre Iteration am reellen Fixpunkt p₆ von H₆ und läuft für Re → −∞ gegen p₆. Die Hexation H₆ hat bei x₆* − 1 eine logarithmische Singularität zur Basis λ₄, nämlich H₆(z) ≈ log(z − z_s)/ln λ₄. Daher hat H₆⁻¹(V) für V → −∞ den endlichen Grenzwert z_s mit Annäherung λ₄^V. Das ergibt α₈ = ln λ₄/ln λ₆. Bei der Tetration hat log_b für V → −∞ keinen endlichen Grenzwert. Die Verschachtelung bricht dann nie ab, und das ist die superlogarithmische Singularität (α₄ = 0). Den Ort habe ich für r = 8 auf 4·10⁻⁴ rad bestätigt, den Exponenten nur grob (1.713 gegen 1.730, die Konvergenz ist langsam). Rang 8 lief nur mit dps 20 und 64 Punkten. Mit ρ = 0.2 traf die Taylor-Kreisscheibe um kp₈ eine Singularität von H₇ und ergab einen falschen Multiplikator, mit ρ = 0.08 stimmt er mit `klam` überein.

## 3. Rang 5: anziehender Fixpunkt, θ reicht nicht bis ℝ (B8, B9)

Bei Rang 5 ist der Fixpunkt L₄ von sexp anziehend. Die Koenigs-Funktion σ ist auf dem ganzen Einzugsgebiet B(L₄) eindeutig und holomorph. Darum ist θ₅ auf der Komponente Ω von {z : H₅(z) ∈ B(L₄)}, die i∞ enthält, holomorph. Ausgenommen sind höchstens logarithmische Punkte mit H₅(z+n) = L₄; solche habe ich keine gesehen. Die Karte (`s4_becken5.py`, `.` = Becken L₄, `o` = Becken L̄₄, `#` = entkommt) zeigt unterhalb y ≈ 0.26 bei x ∈ [−0.1, 0.25] eine fraktal durchmischte Zone. Bei x ∈ [0.3, 0.9] reicht das Becken dagegen fast bis zur Achse. Die gemessene Wachstumsrate 2πy_c mit y_c ≈ 0.27 passt zum höchsten Punkt der Zone.

**Teilargument für die natürliche Grenze:** Angenommen, θ₅ liesse sich über einen Randpunkt z₀ ∈ ∂Ω fortsetzen und Ψ (die Umkehrung von σ) wäre längs s₀ + t, t ≥ 0, eindeutig fortsetzbar. Dann wäre H₅(z₀+n) = F(s₀+n) → L₄, also läge H₅(z₀) im Becken. Das widerspricht z₀ ∈ ∂Ω. Ψ ist hier aber nicht ganz, denn es entsteht aus Zweigen von slog. Die Annahme ist also nicht gesichert, und B9 bleibt Vermutung.

**Folge:** Eine Darstellung K₅(x) = F₅(x + θ₅(x)) mit konvergenter Fourier-Reihe gibt es auf der reellen Achse nicht. Die reelle Kneser-artige Pentation ist trotzdem wohldefiniert, denn der Löser baut sie direkt. Bei der Hexation des Motors wurde der abstossende Fixpunkt gewählt, und dort tritt das Problem nicht auf. Die Etappe-6-Hexation (am anziehenden Fixpunkt von K₅) dürfte dieselbe Grenze haben wie Rang 5; das ist ungeprüft.

## 4. Die Welle zwischen den Pentationen (B10, B11)

ψ(z) = ipent(K₅(z)) − z, mit ipent aus dem Motor (dps 30) und K₅ aus e6 (4·10⁻¹⁰). Auf ℝ ist ψ ≈ D/H′ aus Etappe 6, also a sin2πx + b(1−cos2πx) mit |c_{±1}| = √(a²+b²)/2 = 2.80·10⁻⁶. Auf Linien bis y = 1.6 lassen sich c_{−n} bis n = 14 bestimmen (`s5b_psi5_hoch.py`). Die Quotienten c_{−n−1}/c_{−n} streben gegen e^{2πiz_L3} = −1.800·10⁻⁵ + 1.651·10⁻⁵i (n = 13: −1.711 + 1.532, ·10⁻⁵, langsam wie bei Rang 4). Lokal gilt ψ ≈ 1 − z + ipent(slog(K₅(z))), und slog hat bei L₃ einen logarithmischen Verzweigungspunkt; daraus kommt c_n ~ 1/n mit langsam wanderndem Vorfaktor. Die Gegenrichtung χ = K₅⁻¹∘P − z hat ihren Punkt bei P⁻¹(L₃) = −0.95762 + 1.87156i (`s5c_hoehen.py`). Die Konsistenzbedingung ist χ(z + ψ(z)) = −ψ(z).

## 5. Basisabhängigkeit (B12, B13)

Führend gilt a_n ≈ −h_b′(v_n)/(n ln b), und h_b′(v) = h_b′(log_b v)/(v ln b) iteriert. Die ersten beiden Faktoren, 1/n und 1/ln n, hängen nicht von b ab. Danach geht b über ln ln b und über h_b′ an einer Stelle der Grössenordnung 1 ein. Deshalb ist das Verhältnis bei mässigen n nahe 1, für astronomische n aber eine beschränkte, schwankende Funktion von slog_b(ln ln n). Dass |a_1(b)| über drei Grössenordnungen von b nur um 5 % schwankt, ist eine Beobachtung, die ich nicht erklären kann. Grenzwerte für b → η⁺ und b → ∞ sind plausibel (Vermutung).

## Vorschläge für Grapheinträge

- `kn.satz.theta6_singularitaet`: θ der Kneser-Hexation ist genau bei x* = ihex(p) + ℤ singulär, mit algebraischem Knick vom Exponenten α = ln b/ln sexp′(p) | satz | bedingt (Zweigwahl)
- `kn.satz.theta6_koeffizienten`: a_n = C₆ e^{−2πinx*} n^{−1−α}(1+o(1)), Modell [Qⁿ]h((c₆+log(1−Q))/ln λ₄) mit Rest O(1/n²) | satz | bedingt / numerisch
- `kn.verm.theta_gerade_raenge`: Gesetz für gerade Ränge, Ort H_r⁻¹(p_{r−2}) und Exponent ln λ_{r−4}/ln λ_{r−2}; Tetration als Grenzfall α = 0 | vermutung | Ort für r = 4, 6, 8 belegt
- `kn.beob.hexation_zwei_fixpunkte`: Die reguläre Pentation hat einen abstossenden (kp, vom Motor benutzt) und einen anziehenden komplexen Fixpunkt; die Etappe-6-Hexation nutzt den anziehenden | beobachtung | numerisch
- `kn.beob.theta5_wachsende_koeffizienten`: θ der Kneser-artigen Pentation am anziehenden L₄ hat wachsende Koeffizienten (Rate 2π·0.27) und reicht nicht bis ℝ; Grenze = Rand des Beckens von L₄ | beobachtung | numerisch, natürliche Grenze als Vermutung
- `kn.beob.pentationswelle_singularitaet`: Die Welle zwischen regulärer und Kneser-artiger Pentation fällt wie e^{−2π·1.6902n}/n, bestimmt durch K₅(z) = L₃ | beobachtung | numerisch
- `kn.beob.theta_basis_quasiuniversal`: |a_n(b)| hängt für n ≤ 64 nur um ≤ 5 % von b ∈ [1.45, 1000] ab, ist asymptotisch aber nicht universell | beobachtung | numerisch
- Ergänzung zu `kn.satz.theta_singularitaet` (f4/f5): Die Übertragung auf Rang 6 gilt mit verschobenem Ort und anderem Typ; die Frage aus f5 ist damit beantwortet.

## Offene Fragen

1. Die Zweigwahl in B3 allgemein beweisen. Für Rang 6 reicht wie bei Rang 4 ein Punktvergleich plus Eindeutigkeit, nur nicht in Intervallarithmetik geprüft.
2. B9: Ist der Beckenrand wirklich eine natürliche Grenze von θ₅? Und hat die Etappe-6-Hexation (am anziehenden Fixpunkt) dieselbe Grenze?
3. Gibt es eine Kneser-artige Pentation an einem abstossenden komplexen Fixpunkt von sexp auf einem Nebenzweig? Auf dem Hauptzweig habe ich keinen gefunden.
4. Exponent α₈ genauer messen (höheres N, dps 30), Rang 10 prüfen.
5. B11 bei Rang 7 nachrechnen; dazu braucht es eine Kneser-artige Heptation. Achtung: ihyper(7, kp) ist entartet, weil kp Grenzwert von hex ist. Der Motor liefert dort Im = π/ln λ₆ exakt, aber das ist kein Treffer.
6. B13 erklären: Warum ist a_n für kleine n fast basisunabhängig, und gibt es Grenzfunktionen für b → η⁺ und b → ∞?

## Reproduktion (aus diesem Ordner, `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`)

```
python s1_theta6.py 40 0.04 512 ; python s1_theta6.py 30 0.05 512     # θ₆-Koeffizienten (je 1–1.5 min)
python s2_modell6.py daten/theta6_40_0.04_512.json daten/theta6_30_0.05_512.json   # B4 (~1 min)
python s7_lokal6.py ; python s2b_fuehrend6.py                          # B3, B5
python s6_theta8.py 20 0.03 64 0.08                                    # B6 (~1.5 min)
python s3_theta5.py ; python s3b_wachstum5.py ; python s4_becken5.py 100 25   # B8, B9 (Sekunden, e6-Daten nur gelesen)
python s5_psi5.py ; python s5b_psi5_hoch.py ; python s5c_hoehen.py    # B10, B11
python s8_fixpunkte_sexp.py                                            # B8 (keine komplexen slog-Fixpunkte)
python d1_universal.py ; python d2_basen.py 2 "exp(1)" 1.5 4 10 ; python d2_basen.py 1.45 100 1000   # B12, B13
```

Motor-Zwischenspeicher: nur `cache/` in diesem Ordner (`motor.py`). Die e6-Löserdaten und die Etappe-2-Koeffizienten werden nur gelesen. Es gibt keine Hintergrundprozesse und keine Commits.
