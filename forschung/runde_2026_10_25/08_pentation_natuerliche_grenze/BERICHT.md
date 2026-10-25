# Runde 25.10.2026 · 08 · Natürliche Grenze der θ-Welle der Pentation

## Kurz gesagt (einfache Sprache)

Bei der Kneser-artigen Pentation (Basis 2) ist der Fixpunkt L₄ von sexp anziehend. Deshalb reicht die Welle θ₅ nicht bis zur reellen Achse, und ihre Fourier-Koeffizienten wachsen. Die Vermutung aus Runde 18 hiess: Der Rand des Gebiets Ω, in dem H₅(z) ins Becken von L₄ fällt, ist eine **natürliche Grenze** von θ₅. Das habe ich **bewiesen**, und zwar für jede Superfunktion, die an einem anziehenden Fixpunkt gebaut ist. Der Trick: Man schiebt einen angeblichen Fortsetzungspunkt n Schritte nach rechts, bis alles in der kleinen Kreisscheibe um L₄ liegt, wo die Koenigs-Funktion umkehrbar ist. Dann sagt der Identitätssatz, dass der Punkt doch im Becken liegt. Ein Widerspruch. Zahlen dazu: Der höchste Randpunkt liegt bei y_max ≈ 0.2645 (x ≈ 0.03), und die Koeffizienten wachsen genau mit dieser Rate und mit der passenden Phase. Wichtig ist auch, was **nicht** gilt: Der Konvergenzkreis der Fourier-Reihe ist **keine** natürliche Grenze. θ₅ lässt sich unter ihm fortsetzen, überall dort, wo Ω tiefer reicht. Lückensätze (Hadamard, Fabry) sind hier fehl am Platz.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| N1 | **Satz (Durchschieben).** f sei holomorph auf D, und L sei ein anziehender Fixpunkt mit 0 < \|μ\| < 1 und Becken B. F sei holomorph auf der oberen Halbebene ℍ, mit F(ℍ) ⊂ D, F(z+1) = f(F(z)) und F(z) ∈ B für Im z gross. Ω sei die Komponente von F⁻¹(B) ∩ ℍ, die i∞ enthält, und θ = log σ(F)/log μ − z (σ = Koenigs-Funktion). Dann bleibt **jede** analytische Fortsetzung von θ längs Wegen in ℍ in Ω. Jeder Punkt von ∂Ω ∩ ℍ ist singulär, ∂Ω ∩ ℍ ist also eine natürliche Grenze. | **bewiesen** | Beweis in Abschnitt 1, benutzt nur den Satz von Koenigs, den Identitätssatz und F(z+n) = fⁿ(F z) | neu (die Vermutung B9 aus Runde 18 ist damit bewiesen, für K₅ unter den Annahmen N2) |
| N2 | Anwendung auf K₅ = H₅ (e6-Löser), f = sexp₂, L = L₄ = 1.469818 + 1.664342i, μ = 0.154239 + 0.690508i: Die Annahmen gelten auf dem Prüfraster. F(z+1) = f(F z) gilt auf 4·10⁻¹⁵ (innerhalb des Lösers), H₅ trifft den Schnitt (−∞, −2] von sexp auf y ∈ [0.02, 0.5] nicht. Damit ist ∂Ω eine natürliche Grenze von θ₅. | hergeleitet mit Lücke | Lücke: Die Existenz von K₅ mit diesen Eigenschaften ist nur numerisch gesichert (e6-Löser, ~4·10⁻¹⁰), nicht bewiesen | neu |
| N3 | **Konvergenzradius.** Für die Fourier-Reihe θ₅ = Σ aₙqⁿ gilt limsup \|aₙ\|^{1/n} = e^{2π y_c} mit y_c = max(y_max, Höhe der höchsten Nullstelle von σ∘H₅ in Ω). Dabei ist y_max = sup{Im z : z ∈ ℍ∖Ω}. Insbesondere gilt y_c ≥ y_max. | bewiesen (aus N1 und Cauchy-Hadamard) | {Im z > y_max} ⊂ Ω; dort hat θ höchstens logarithmische Punkte (Nullstellen von σ∘H₅) | neu |
| N4 | Höchster Randpunkt: y_max ≈ 0.2645 bei x ≈ 0.031. Eine zweite Spitze liegt bei 0.2635, x ≈ 0.0235. Im groben Raster liegt die Oberkante für x ∈ [0.01, 0.10] überall bei 0.26, ist also ein breites, flaches Dach. | numerisch | Raster 7.5·10⁻⁴ × 5·10⁻⁴, keine unentschiedenen Punkte; nur eine **untere** Schranke für y_max (dünne Spitzen könnten fehlen) | genauer als Runde 18 (0.26, Raster 0.02) |
| N5 | Die Koeffizienten bestätigen N3/N4 von zwei Seiten. (a) Die lokale Rate (ln\|aₙ\| − ln\|aₙ₋₈\|)/16π beträgt 0.264–0.277 für n = 18–32. (b) Die Phase −arg(aₙ/aₙ₋₁)/2π mod 1 beträgt 0.002–0.028 für n = 18–32, also liegt die dominante Singularität bei x ≈ 0.02–0.03. Das ist dieselbe Stelle wie die Spitze in N4. Also y_c = y_max: keine logarithmischen Punkte über dem Dach. | numerisch | Linien y = 0.30 und 0.33 (N = 512) stimmen bis n = 24 auf ≤ 2·10⁻³ überein; negative Moden ≤ 4·10⁻⁸; ab n ≈ 32 Rauschen | neu |
| N6 | **Der Konvergenzkreis ist keine natürliche Grenze.** θ₅ ist unterhalb von y_max holomorph, wo Ω tiefer reicht. Beispiel: z₁ = 0.6 + 0.12i (senkrechte Strecke 0.45i → 0.08i ganz im Becken). Dort gilt die Mittelwerteigenschaft auf 7·10⁻⁸ (ρ = 0.02) bzw. 1.7·10⁻⁷ (ρ = 0.03), während die Fourier-Glieder divergieren (\|a₆₀q⁶⁰\| = 1.8·10²²). Auf dem Kreis \|q\| = e^{−2πy_max} liegen nur die Dachspitzen als Singularitäten. | negativ (für die Kreis-Fassung) / numerisch | `s2_koeff.py`; Genauigkeit begrenzt durch e6 (~10⁻¹⁰), verstärkt durch log σ | neu |
| N7 | Lückensätze (Hadamard, Fabry, Pólya–Carlson) greifen nicht. Sie liefern den Konvergenzkreis als natürliche Grenze, aber hier sind alle aₙ ≠ 0 mit glattem Wachstum, und der Kreis ist nach N6 gar keine Grenze. Die natürliche Grenze kommt aus der Dynamik, nicht aus der Arithmetik der Koeffizienten. | negativ | N1, N6 | – |
| N8 | Unter dem Dach (x ∈ [0, 0.05], y ∈ [0.215, 0.265]) liegen 96 % der Rasterpunkte im Becken von **L̄₄** (2400 von 2500), nur 4 % in B(L₄). Die Grenze ist dort die Grenze zwischen den beiden Becken. Die Kastenzählung ergibt N = 34, 14, 7, 3, 2 für ε = 0.001 … 0.016, Steigung 1.04 (lokal 0.6–1.3). Auf Skalen ≥ 10⁻³ gibt es also **keinen** klaren Fraktalbeleg. | numerisch (grob) | `s3_kasten.py`, ein Fenster, 2500 Punkte | – |
| N9 | Die Bahn von L₃ (Fixpunkt von 2^z, Grenzwert von sexp für Im → ∞) läuft unter sexp nach L₄. Punkte mit grossem \|Im w\| liegen also im Becken von L₄. Das Abbruchkriterium \|w\| > 60 in `s4_becken5.py` (Runde 18) zählt solche Punkte fälschlich als „entkommt“. Hier ersetzt durch Re w > 60 oder nicht endlich. | numerisch | Bahn: 1.3245+1.2890i, 1.7005+1.4398i, … → 1.4804+1.6677i nach 12 Schritten (spiralförmig, wie μ) | Korrektur |
| N10 | Derselbe Satz gilt für die e6-Hexation (am anziehenden Fixpunkt von K₅) und für jede Kneser-artige Stufe an einem anziehenden Fixpunkt. Bei abstossenden Fixpunkten (Motor-Ränge 4, 6, 8) gibt es keine solche Grenze. Die Art des Fixpunkts entscheidet also: abstossend gibt eine Welle bis ℝ mit isolierten Singularitäten, anziehend gibt eine natürliche Grenze im Innern von ℍ. | bewiesen (Struktur, aus N1) / Vermutung (dass ∂Ω ∩ ℍ bei der e6-Hexation nicht leer ist) | nicht nachgerechnet | neu |

## 1. Beweis von N1

**Bezeichnungen.** σ: B → ℂ ist die Koenigs-Funktion, σ(w) = lim μ⁻ⁿ(fⁿ(w) − L). Sie ist holomorph auf dem ganzen Becken, mit σ∘f = μσ und σ′(L) = 1. Es gibt eine Umgebung N von L mit f(N) ⊂ N, sodass σ: N → D(0, r₀) biholomorph ist. Die Umkehrung heisse σ_N⁻¹. A = F⁻¹(B) ∩ ℍ ist offen und 1-periodisch, und Ω ist seine Komponente, die i∞ enthält. Auf Ω ist S(z) := σ(F(z)) holomorph mit S(z+1) = μS(z). Jeder Zweig von θ erfüllt exp(log μ·(z + θ(z))) = S(z), denn die Zweige unterscheiden sich nur um Vielfache von 2πi/log μ.

**Lemma.** U ⊂ ℍ sei eine Kreisscheibe und T holomorph auf U. Auf einer nichtleeren offenen Menge V ⊂ U ∩ Ω stimme T mit einem Zweig von θ überein. Dann gilt U ⊂ A.

*Beweis.* Es genügt, das für jede abgeschlossene konzentrische Teilscheibe U′ ⊂ U mit V ∩ int U′ ≠ ∅ zu zeigen. Auf U′ ist T beschränkt. Also strebt Eₙ(z) := μⁿ·exp(log μ·(z + T(z))) gleichmässig auf U′ gegen 0. Für n ≥ n₀ ist daher Wₙ := σ_N⁻¹∘Eₙ holomorph auf U′ mit Werten in N. Wähle v ∈ V ∩ int U′. Da F(v) ∈ B liegt und fⁿ → L lokal gleichmässig auf B gilt, gibt es eine Scheibe D_v ⊂ V ∩ U′ und ein n₁ ≥ n₀, sodass F(z+n) = fⁿ(F(z)) ∈ N für z ∈ D_v und n ≥ n₁. Dort ist σ(F(z+n)) = μⁿS(z) = Eₙ(z). Weil σ auf N injektiv ist, folgt F(z+n) = Wₙ(z) auf D_v. Beide Seiten sind holomorph auf int U′ (es gilt U′ + n ⊂ ℍ). Nach dem Identitätssatz gilt also F(z+n) = Wₙ(z) ∈ N ⊂ B auf ganz U′. Wegen F(z+n) = fⁿ(F(z)) folgt F(z) ∈ B, also U′ ⊂ A. ∎

**Satz.** γ: [0,1] → ℍ sei ein Weg, der in {Im z > Y₀} beginnt, und θ lasse sich längs γ fortsetzen. Dann gilt γ([0,1]) ⊂ Ω.

*Beweis.* Angenommen, γ verlässt Ω. Sei t* = inf{t : γ(t) ∉ Ω}, dann liegt z₀ = γ(t*) in ∂Ω. Das Funktionselement bei t* ist holomorph auf einer Scheibe U um z₀. Für t < t* nahe t* stimmt es auf einer kleinen Scheibe um γ(t) ⊂ Ω ∩ U mit einem Zweig von θ überein. Nach dem Lemma gilt U ⊂ A. Da U zusammenhängend ist und Ω trifft, folgt U ⊂ Ω, also z₀ ∈ Ω. Das widerspricht z₀ ∈ ∂Ω, denn Ω ist offen. ∎

Umgekehrt lässt sich θ längs jedes Weges in Ω fortsetzen, der die Nullstellen von S meidet (dort ist θ = log S/log μ − z lokal holomorph). Der maximale Bereich von θ in ℍ ist also genau Ω, ohne höchstens logarithmische Punkte. Gegenüber dem Teilargument von Runde 18 ist neu: Man braucht die globale Umkehrung Ψ = σ⁻¹(μ^s) gar nicht. Diese wäre wegen der slog-Zweige mehrdeutig. Es genügt die *lokale* Umkehrung σ_N⁻¹ nahe L, nachdem man n Schritte nach rechts geschoben hat.

**N3.** θ(z) = g(e^{2πiz}), und g ist auf |q| < e^{−2πy_max} holomorph bis auf logarithmische Punkte, denn {Im z > y_max} ⊂ Ω. Nach N1 ist g an jedem Punkt von ∂Ω̃ singulär, insbesondere an den Dachspitzen auf |q| = e^{−2πy_max}. Mit Cauchy-Hadamard folgt N3.

**Bemerkung zu Fraktalität.** N1 braucht keine Fraktalität. Ob ∂Ω fraktal ist (Hausdorff-Dimension > 1), ist eine eigene Frage. Lokal ist H₅ konform, deshalb hat ∂Ω dieselbe Dimension wie der Rand ∂B(L₄) ⊂ J(sexp). N8 gibt dazu nur einen groben Wert nahe 1.

## 2. Zahlen

- `s1_rand.py 100 0.01`: L₄, μ, Bahn von L₃, Prüfung der Funktionalgleichung. Grobes Raster (100 × 25, nmax 150): Klassen {L₄: 2288, L̄₄: 179, entkommt: 33, offen: 0}. Feinsuche um die zwei höchsten Stellen ergibt y_max = 0.26350 (x = 0.02350) und 0.26450 (x = 0.03100), ohne unentschiedene Punkte. Laufzeit 82 s.
- `s2_koeff.py 0.30 0.33 512 0.6 0.12`: Koeffizienten bis n = 120 (verlässlich bis n ≈ 30), Rate, Phase, Mittelwerttest bei z₁. Laufzeit ~1 s.

| n | \|aₙ\| (y = 0.30) | rel. Diff. y = 0.33 | lokale Rate | Phase −arg(aₙ/aₙ₋₁)/2π |
|---|---|---|---|---|
| 10 | 6.534·10⁵ | 1.2·10⁻⁵ | 0.2218 | 0.9405 |
| 14 | 6.915·10⁸ | 1.0·10⁻⁴ | 0.2434 | 0.9962 |
| 18 | 7.256·10¹¹ | 5.8·10⁻⁴ | 0.2769 | 0.0021 |
| 20 | 2.465·10¹³ | 1.1·10⁻³ | 0.2768 | 0.0230 |
| 24 | 1.350·10¹⁶ | 2.3·10⁻³ | 0.2670 | 0.0283 |
| 32 | 7.970·10²¹ | (Rauschen) | 0.2644 | 0.0230 |

- `s3_kasten.py 0.0 0.215 0.05 0.05 0.001`: Kastenzählung (N8). Die erste Ausgabeschleife im Skript (mit „N = 0“) ist fehlerhaft und wird nicht verwendet. Gültig ist nur die Zeile „Kasten (Seite, N)“.

## Vorschläge für Grapheinträge

- `kn.satz.theta_anziehend_natuerliche_grenze`: Bei einer Superfunktion an einem anziehenden Fixpunkt ist der Rand der Becken-Urbildkomponente Ω eine natürliche Grenze von θ (Durchschiebe-Lemma) | satz | bewiesen
- `kn.satz.theta5_natuerliche_grenze`: Anwendung auf die Kneser-artige Pentation, ∂Ω natürliche Grenze von θ₅ | satz | bedingt (K₅ numerisch)
- `kn.satz.theta5_konvergenzradius`: y_c = max(y_max, logarithmische Punkte) | satz | bewiesen
- `kn.beob.theta5_dach`: y_max ≈ 0.2645 bei x ≈ 0.03; Rate und Phase der aₙ treffen die Spitze | beobachtung | numerisch
- `kn.beob.theta5_unter_konvergenzkreis`: θ₅ setzt sich unter den Konvergenzkreis fort; der Kreis ist keine natürliche Grenze, Lückensätze greifen nicht | beobachtung | numerisch / negativ
- Änderung an `kn.beob.theta5_wachsende_koeffizienten`: Vermutung „natürliche Grenze“ → Satz (bedingt); Hinweis auf das falsche Entkommen-Kriterium in `s4_becken5.py`

## Offene Fragen

1. Gibt es Nullstellen von σ∘H₅ in Ω, also Punkte mit H₅(z+k) = L₄? N5 schliesst sie nur oberhalb des Dachs aus. Zur Prüfung: Argumentprinzip für S auf Teilgebieten von Ω.
2. Streng zeigen, dass ∂Ω ∩ ℍ ≠ ∅ ist, ohne Numerik. Dass ℝ ∩ A = ∅, ist klar (reelle Bahnen bleiben reell), aber daraus folgt noch nicht, dass Ω ≠ ℍ.
3. Hausdorff-Dimension von ∂Ω bzw. ∂B(L₄): grösseres Fenster und feinere Raster (≤ 10⁻⁴), dazu mehr Rechenzeit.
4. Genaues y_max als Supremum über dünne Spitzen. Gibt es eine geschlossene Kennzeichnung der Dachspitze, etwa H₅(z_tip) auf einem abstossenden Zyklus von sexp?
5. e6-Hexation: ∂Ω ∩ ℍ nachrechnen (N10).

## Reproduktion (aus diesem Ordner, `C:/Users/janis/miniconda3/python.exe -u`, `PYTHONIOENCODING=utf-8`)

```
python s1_rand.py 100 0.01                      # N2, N4, N9 (~80 s)
python s2_koeff.py 0.30 0.33 512 0.6 0.12       # N5, N6 (~1 s)
python s3_kasten.py 0.0 0.215 0.05 0.05 0.001   # N8 (~3 s)
```

Die e6-Daten (`eta_hyperoperationen/halb_rang_basis2_2026_10_13/data/e6_H*.pkl`) werden nur gelesen. Zwischenergebnisse liegen in `cache/`. Es laufen keine Hintergrundprozesse, und es gibt keine Commits.
