# Identitätssuche, Etappe 4 (14.10.2026): Wählt eine reelle Mittelwert-Gleichung allein Kneser aus?

## Kurz gesagt

- **Nein, keine einzelne der neuen Gleichungen reicht, und auch keine endliche Auswahl.** Ich habe für jede der vier Gleichungen
  (a) ∫θ² = (∫θ)², (b) σ-Produkt = 1, (c) ∫₁ᵉ slog′²/h′ = 1, (d) ∫ log(T′·h′∘T) = 0
  eine echte, glatte, wachsende Tetration T ≠ Kneser gebaut, die sie **exakt** erfüllt. Die Gleichung stimmt dann auf 10⁻³⁵ oder besser,
  die übrigen Kennzeichen sind um 10⁻⁴ bis 10⁻² verletzt. Dazu kommt ein Beispiel, das
  sechs Gleichungen zugleich erfüllt (a, c, d und die ersten drei „Fourier-Bedingungen“), und eines für a und b zugleich. Grund: Jede Gleichung ist nur eine
  (komplexe) Zahl, die Freiheit dagegen eine ganze periodische Funktion p. Man stellt zwei Zahlen in p nach, und die Gleichung stimmt wieder.
- **Unendlich viele Bedingungen reichen: „θ_T hat keine negativen Fourier-Moden“ ist gleichwertig zu Kneser** (Satz 3, bewiesen
  bis auf die bekannte Majoranten-Lücke). Reell ausgeschrieben: Kneser ist genau die Tetration mit
  **∫₀¹ h(T(x))·e^{2πimx} dx = 1/(2πim) für alle m = 1, 2, 3, …**
  Gleichwertig, als *eine* reelle Gleichung: Imaginärteil und Realteil von h(T(x)) − x sind ein Hilbert-Paar (Kramers–Kronig).
- **Eine einzige reelle Grösse mit Extremalprinzip gibt es:** der „Kneser-Defekt“
  E(T) = ¼∬|θ_T(x)−θ_T(y)|²/sin²π(x−y) dx dy − (1/π)·(von θ_T umlaufene Fläche) = 2Σ_{m≥1} m|θ̂_T(−m)|² ≥ 0.
  **Kneser ist die einzige Tetration mit E(T) = 0, also das eindeutige Minimum.** Das ist die klassische Ungleichung
  „Dirichlet-Energie ≥ Bildfläche, Gleichheit genau bei holomorph“, angewandt auf h∘T − id. Ehrlich gesagt ist das eine Umformulierung
  der Fourier-Bedingung, aber eine rein reelle und geometrische.
- **Kein Extremalprinzip mit den einzelnen Gleichungen:** Ihre erste Ableitung nach p ist bei Kneser nicht null. Kneser ist also
  weder Minimum noch Maximum von ∫ln T′ + ∫ln|h′∘T| oder ähnlichen Grössen. Jede kennzeichnende Gleichung muss entweder unendlich viele
  Bedingungen enthalten oder „von zweiter Ordnung“ sein, also eine Quadratsumme wie E.
- **Exakt lösbares Beispiel:** Verbiegt man Kneser mit einer Möbius-Abbildung des Kreises (Parameter a, |a| < 1), dann gilt exakt
  θ̂_T(−n) = i aⁿ/(2πn) und **E(T) = −ln(1−|a|²)/(2π²)**. Numerisch bestätigt auf 10⁻²⁹. Der Defekt misst hier genau den hyperbolischen Abstand von Kneser.
- **Momentenfamilie** M[θ_T^k] = M[θ_T]^k für alle k (gleichwertig: M[e^{cθ_T}] = e^{cM[θ_T]} für alle c):
  Für allgemeine Funktionen ist sie **schwächer** als „keine negativen Moden“. Gegenbeispiel: f = e^{−2πix}. Die Familie sieht nur die
  Werteverteilung von f, nicht die Reihenfolge. Unter den Tetrationen ist sie aber **infinitesimal starr** (Satz M, bewiesen für alle p, die in einem
  Streifen der Halbbreite 0.09 holomorph sind; numerisch für Halbbreite > 0.001). Ob sie Kneser auch nichtlinear festlegt, bleibt offen.
- **Lineare Tetration (Hooshmand)** verletzt alle Gleichungen (z. B. ∫₁ᵉ slog′²/h′ = 0.9675). Ihre negativen Moden fallen nur wie m⁻²,
  weil sie an den ganzen Zahlen nicht analytisch ist. Für jede reell-analytische Tetration fallen sie dagegen exponentiell (Lemma R).

**Antwort auf br.frage.bohr_mollerup (Teil):** Es gibt eine rein reelle Kennzeichnung, aber sie ist **nichtlokal** (ein Integral mit
cot- bzw. sin⁻²-Kern über eine Periode). Eine einzelne Mittelwert-Gleichung reicht nicht, eine endliche Liste auch nicht.
Eine *lokale* Bedingung wie Konvexität (Bohr–Mollerup) ist nicht gefunden. Ein Hinweis: Bohr–Mollerup ist selbst eine *unendliche Familie von
Ungleichungen*, und nach Abschnitt 4 muss jede kennzeichnende Bedingung von dieser Art sein.

## 1. Rahmen und Bezeichnungen

Basis e. Bezeichnungen wie in Etappe 3: L, σ, h, F, θ(z) = h(sexp z) − z, M[G] = ∫₀¹ G dx. Jede reell-analytische Tetration mit T′ > 0 und T(0) = 1 ist
T = sexp∘φ, φ(x) = x + p(x), p reell-analytisch, 1-periodisch, p(0) = 0, 1 + p′ > 0 (fg.satz.kohomologisch). Dann gilt

  θ_T(x) := h(T x) − x = θ(φ x) + p(x),  1 + θ_T′ = T′·h′(T) = (1 + θ′∘φ)·φ′.

Untersuchte Bedingungen (alle komplex, also je zwei reelle Gleichungen):

| | Bedingung | reelle Form |
|---|---|---|
| a | M[θ_T²] − M[θ_T]² = 0 | ∫θ_T² = (∫θ_T)² |
| b | M[e^{Lθ_T}]·M[e^{−Lθ_T}] − 1 = 0 | ∫σ(T)e^{−Lx}·∫e^{Lx}/σ(T) = 1 |
| c | M[1/(1+θ_T′)] − 1 = 0 | ∫₁ᵉ slog_T′(w)²/h′(w) dw = 1 |
| d | M[log(1+θ_T′)] = 0 | ∫ln T′ = −∫ln\|h′∘T\| und ∫arg h′∘T = 0 |
| e_m | θ̂_T(−m) := M[θ_T e^{2πimx}] = 0 (m ≥ 1) | ∫₀¹ h(T x) e^{2πimx} dx = 1/(2πim) |
| V_k | M[(θ_T − Mθ_T)^k] = 0 (k ≥ 2) | Momentenfamilie (V₂ = a) |

(e_m in reeller Form: ∫₀¹ x e^{2πimx} dx = 1/(2πim).)

## 2. Ergebnisse

| # | Aussage | Status | Beleg |
|---|---|---|---|
| S1 | **Linearisierung:** Für jede Bedingung C aus a–d, V_k, e_m gilt δC[p] = Σ_{k≥1} F̂_C(k)·conj(p̂(k)), mit F_C holomorph und periodisch in ℍ, F_C → 0 oben: F_a = 2gw, F_b = 2L sinh(Lg)·w, F_c = −2θ″/w², F_d = θ″/w, F_{V_k} = k g^{k−1} w, F_{e_m} = qᵐw (g = θ − a₀, w = 1 + θ′, q = e^{2πiz}). Konstanten in p wirken nicht. Es zählt nur der Teil von p mit positiver Frequenz, und zwar antilinear. | hergeleitet (B-S1) | numerisch: (∂C/∂ε_sin)/(∂C/∂ε_cos) = i auf 10⁻¹⁷ für alle C, m = 1, 2, 3; ∂a/∂ε_cos (m=1) = a₁ und (m=2) = a₂ + 2πi a₁², wie vorhergesagt (`log_s2_lin.txt`) |
| S2 | **Endlich viele Gleichungen legen Kneser nie fest:** Ist p ↦ (δC₁, …, δC_n) surjektiv auf ℂⁿ, so gibt es eine Kurve p_t ≠ konst. reell-analytischer Störungen (trigonometrische Polynome), längs der alle C_i exakt null sind. | bewiesen (Satz der impliziten Funktion in endlich vielen Parametern, B-S2); Surjektivität für jede einzelne Bedingung a–d gezeigt (F̂_a(1) = 2a₁ ≠ 0 explizit, sonst numerisch) | Gegenbeispiele Tabelle 3 |
| S2a | a und b sind **in erster Ordnung fast gleich**: b = L²V₂ + L⁴(V₄/12 + V₂²/4) + …; die Ableitungen stimmen auf den Moden 1 und 2 exakt überein (b = L²·a), verschieden erst ab Mode 3 (Faktor g³ ~ 10⁻³). | bewiesen (Reihenentwicklung) | Gegenbeispiel zu a erfüllt b bis auf 1.1·10⁻⁷ |
| S3 | **Kennzeichnung durch unendlich viele Gleichungen:** T reell-analytisch, T′ > 0, T(0) = 1. Dann: T = Kneser ⇔ θ̂_T(−m) = 0 für alle m ≥ 1 ⇔ N(T) := Σ\|θ̂_T(−m)\|² = 0 ⇔ E(T) = 0 ⇔ Im h(T x) = p.v.∫₀¹ (Re h(T y) − y) cot π(x−y) dy + konst. | **bewiesen** für „⇐“ mit Paulsen–Cowgill (B-S3); „⇒“ ist Z5 aus Etappe 3 (Lücke: Majorante). Beide Richtungen brauchen θ_T ∈ L¹ (Lokalform aus Etappe 2) | Kneser: \|θ̂(−m)\| ≤ 2·10⁻²⁹ (dps 30) für m ≤ 24 |
| S4 | **Extremalprinzip:** E(T) = ¼∬\|θ_T(x)−θ_T(y)\|²/sin²π(x−y) − (1/π)·Fläche(θ_T) = 2Σ_{m≥1} m\|θ̂_T(−m)\|² ≥ 0; Kneser ist das eindeutige Minimum. | bewiesen (Douglas-Formel + Flächenformel, B-S4) und S3 | Werte für Beispiele in Tabelle 3 |
| R | **Lemma R:** Für reell-analytisches T ist der Anteil negativer Frequenz von θ_T in einem Streifen um ℝ holomorph, also gilt \|θ̂_T(−m)\| ≤ C e^{−cm}. N(T) und E(T) sind endlich. | hergeleitet (B-R), Lücke wie Z0 (lokale Hardy-Klasse) | T = sexp(x+0.03 sin 4πx): \|θ̂_T(−m)\| fällt von 3·10⁻⁴ (m=3) auf 1.6·10⁻²⁰ (m=24). Lineare Tetration: nur ~ m⁻² (m²\|θ̂\| = 0.0017, 0.0012, 0.0009, 0.0008 bei m = 4, 8, 16, 24) |
| MÖ | **Möbius-Familie exakt:** φ mit e^{2πiφ} = e^{iα}(q−a)/(1−āq). Dann θ̂_T(−n) = i aⁿ/(2πn), N = Li₂(\|a\|²)/(4π²), E = −ln(1−\|a\|²)/(2π²). | bewiesen (B-MÖ) | a = 0.05, 0.3, 0.2+0.25i: θ̂_T(−n) für n ≤ 8 auf ≤ 5·10⁻²⁹ (dps 30) |
| M | **Momentenfamilie:** (i) Für allgemeine f gilt: Momentenfamilie ⇔ die Werteverteilung f_*(dx) stellt den Punkt Mf für Polynome dar; das impliziert **nicht** „keine negativen Moden“ (f = e^{−2πix}; f(x) = G(−x); f = G∘ψ mit maßtreuem ψ). (ii) Unter Tetrationen ist die linearisierte Familie starr: δV_k[p] = 0 für alle k ⇒ p konstant, für jedes p, das im Streifen \|Im z\| < δ holomorph ist mit e^{−2πδ} < ρ, wobei g = θ − a₀ auf \|q\| ≤ ρ univalent ist. (iii) Nichtlineare Eindeutigkeit offen. | (i) bewiesen; (ii) bewiesen (B-M) mit ρ ≥ 0.567 (Koeffizientenkriterium, δ > 0.090); numerisch univalent bis ρ = 0.995 (δ > 0.0008); (iii) offen | `log_s3.txt`; alle Gegenbeispiele zu a–d verletzen V₃ um ~2·10⁻⁴ |
| H | Lineare Tetration: a = −9.0·10⁻⁴ + 2.5·10⁻⁵ i, b = 1.5·10⁻³ − 8.1·10⁻⁴ i, ∫₁ᵉ slog′²/h′ = 0.96753, d = 0.0184, N ≈ 2.29·10⁻⁵, E ≈ 4.7·10⁻⁵ (m ≤ 24) | numerisch | `log_s4.txt` |
| X | Keine einzelne der Grössen a–d (Real- oder Imaginärteil) hat bei Kneser ein Extremum. | bewiesen (S1: Ableitung ≠ 0) | – |

### Tabelle 3: Gegenbeispiele (T = sexp(x + p(x)), alle mit 1 + p′ > 0.55, reell-analytisch, streng wachsend, T(0) = 1)

Newton bzw. Sehnenverfahren mit dps 40 und Quadratur k = 6; Nachprüfung unabhängig mit dps 50 und k = 7 (`s5_pruefen.py`, `log_s5_50_7.txt`).
Werte |·| bei der Nachprüfung.

| Fall | erfüllt (Restfehler, dps 50, k = 7) | p (gerundet) | verletzt (Auswahl) | N₈ = Σ_{m≤8}\|θ̂_T(−m)\|² | min(1+p′) |
|---|---|---|---|---|---|
| a | \|a\| = 3.5·10⁻³⁹ | 0.03 sin 4πx − 0.0352539 sin 2πx + 0.000461846(cos 2πx − 1) | c 0.031, d 0.039, \|θ̂(−1)\| 0.0090, V₃ 2.5·10⁻⁴, b 1.1·10⁻⁷ | 3.15·10⁻⁴ | 0.604 |
| b | \|b\| = 3.0·10⁻⁴⁰ | 0.03 sin 4πx − 0.0352530 sin 2πx + 0.000461637(cos 2πx − 1) | c 0.031, d 0.039, \|θ̂(−1)\| 0.0090, V₃ 2.5·10⁻⁴, a 5.8·10⁻⁸ | 3.15·10⁻⁴ | 0.604 |
| c | \|c\| = 3.9·10⁻³⁶ | 0.03 sin 4πx − 0.0460228 sin 2πx − 0.00167248(cos 2πx − 1) | a 7.3·10⁻⁴, b 1.4·10⁻³, d 0.025, \|θ̂(−1)\| 0.014 | 4.50·10⁻⁴ | 0.585 |
| d | \|d\| = 3.6·10⁻³⁸ | 0.03 sin 4πx − 0.0666380 sin 2πx − 0.00174313(cos 2πx − 1) | a 1.9·10⁻³, b 3.7·10⁻³, c 0.054, \|θ̂(−1)\| 0.025 | 8.80·10⁻⁴ | 0.554 |
| a+b | \|a\| = 1.1·10⁻³⁵, \|b\| = 2.1·10⁻³⁵ | 0.001 sin 8πx + 8.048·10⁻⁴ sin 2πx − 6.71·10⁻⁵(cos 2πx−1) − 1.8274·10⁻³ sin 6πx + 4.02·10⁻⁵(cos 6πx−1) | c 2.6·10⁻³, d 2.3·10⁻⁴, \|θ̂(−3)\| 6.3·10⁻⁴, V₃ 1.2·10⁻⁵ | 7.6·10⁻⁷ | 0.943 |
| alle | a 1.3·10⁻³⁸, c 3.0·10⁻³⁷, d 1.1·10⁻³⁷, \|θ̂(−1,−2,−3)\| ≤ 1.6·10⁻³⁹ | 10⁻⁷ sin 14πx + Moden 1–6 (größte: −2.87·10⁻⁷ sin 12πx), max\|p\| = 7.6·10⁻⁷ | \|θ̂(−4)\| 1.9·10⁻⁸, \|θ̂(−6)\| 1.2·10⁻⁷, V₃ 8·10⁻¹³ | 2.3·10⁻¹⁴ | 0.99997 |

Alle Restfehler liegen 30 Grössenordnungen unter den Verletzungen. Bei „alle“ ist die Störung absichtlich winzig: Das 6×6-System ist schlecht
konditioniert (kleinster Singulärwert 3.9·10⁻⁴, die Zeilen a, c, d, e_m sind auf hohen Moden fast abhängig), und schon bei t = 10⁻⁴ sind die nichtlinearen
Glieder von c zu gross. Nach S2 gibt es die Lösung für jedes kleine t; die Grösse der Störung spielt für die Aussage keine Rolle.

Bei den Einzelbeispielen ist p = 0.03 sin 4πx + ε₁ sin 2πx + ε₂(cos 2πx − 1), dabei sind ε₁ und ε₂ gelöst. Bei „a+b“ sind die Moden 1 und 3 gelöst (b ist von a erst ab Mode 3 unabhängig, S2a).
Bei „alle“ ist p = 10⁻⁷ sin 14πx + Σ_{k≤6}(α_k sin 2πkx + β_k(cos 2πkx − 1)) mit 12 gelösten Koeffizienten (Sehnenverfahren mit der theoretischen Ableitung aus S1, `s2b_alle.py`). Die vollständigen Koeffizienten (40 Stellen) stehen in `daten/gegenbeispiel_*_40.json`.

## 3. Beweise

**B-S1 (Linearisierung).** Exakt gilt θ_T(x) = θ(x+p) + p, also δθ_T = w·p mit w = 1+θ′. Für C = M[Φ(θ_T)] folgt δC = M[Φ′(θ)·w·p]. Für a ist das
2M[(θ−a₀)wp], für b L·M[(e^{Lg} − e^{−Lg})wp] (mit M[e^{±Lθ}] = e^{±La₀}), für V_k k·M[g^{k−1}wp] (es gilt M[g^j] = 0). Für die Grössen mit θ_T′ gilt
δθ_T′ = (wp)′. Partielle Integration ergibt δc = −M[(wp)′/w²] = −2M[p θ″/w²] und δd = M[(wp)′/w] = M[p θ″/w]. Die Randterme verschwinden, weil
p/w → 0 an den ganzen Zahlen (\|w\| ≍ 1/(x ln(1/x)), Etappe 2). Für e_m gilt δθ̂_T(−m) = M[w qᵐ p]. In allen Fällen ist F_C holomorph
und periodisch in ℍ mit F_C → 0 oben. Nach Z0 (Etappe 3) ist dann F̂_C(k) = 0 für k ≤ 0. Daraus folgt M[F_C p] = Σ_{k≥1}F̂_C(k)p̂(−k), und für reelles p gilt p̂(−k) = conj p̂(k). ∎
Folgerung: Für p = ε cos 2πkx ist δC = εF̂(k)/2, für p = ε sin 2πkx ist δC = iεF̂(k)/2. Die reelle 2×2-Jacobi-Matrix in diesen Richtungen ist die
Multiplikation mit F̂(k)/2, ihre Determinante ist \|F̂(k)\|²/4. Bei a ist F̂_a(1) = 2a₁ (denn (g²)′ hat keine Mode 1), also invertierbar.

**B-S2 (endlich viele Gleichungen).** Seien C₁, …, C_n Bedingungen, und die Vektoren (F̂_{C_i}(k))_k seien ℂ-linear unabhängig. Dann gibt es K, sodass die
Einschränkung auf die Moden 1…K surjektiv auf ℂⁿ ist. Wähle 2n reelle Richtungen e_j (sin/cos−1 dieser Moden), auf denen die Jacobi-Matrix invertierbar ist, und eine weitere
Richtung v ∉ span(e_j). Die Abbildung (t, ε) ↦ C(t v + Σε_j e_j) ist in endlich vielen Variablen C¹. Das braucht dominierte Konvergenz unter M, mit derselben Majorante wie Z0. Nach dem Satz
der impliziten Funktion gibt es ε(t) mit C(p_t) = 0 und p_t ≠ konst für kleine t ≠ 0. Für kleine t gilt 1 + p_t′ > 0. ∎
Achtung: Für {a, b} gilt Unabhängigkeit erst ab Mode 3 (S2a). Bei {a, b, c, d, e₁, e₂, e₃} auf 7 Moden ist der kleinste Singulärwert 1.6·10⁻⁷ (die b-Zeile ist klein). Deshalb habe
ich {a, b} und {a, c, d, e₁, e₂, e₃} getrennt gelöst. Auch {a, c, d, e₁, e₂, e₃} ist schlecht konditioniert (3.9·10⁻⁴ auf Moden 1–6). Lässt man Mode 1 weg, fällt der Wert auf 8·10⁻⁶. Die Gewichte F_C haben auf hohen Moden fast dieselbe Gestalt, weil alle von derselben Singularität bei q = 1 bestimmt sind.

**B-S3 (Kennzeichnung).** „⇒“: Z5 (Etappe 3). Alternativ: Ist Σ\|a_n\|² < ∞ (nach Etappe 2 gilt \|a_n\| ≍ 1/(n ln n)), so liegt θ in H²(Kreisscheibe). Seine Randwerte sind
fast überall die nichttangentialen Grenzwerte, also h(sexp x) − x, und sie haben keine negativen Moden.
„⇐“: θ_T ∈ L¹(0,1). Das folgt aus der Lokalform θ_T(x) = 2 − x + h(log(T′(0)x) + O(x)) wie in Etappe 2, S6, denn T ist bei 0 analytisch mit T′(0) > 0. Hat θ_T keine
negativen Moden, so ist G(z) := Σ_{k≥0} θ̂_T(k) e^{2πikz} das Poisson-Integral von θ_T, also in H¹, mit nichttangentialen Randwerten θ_T fast überall. Setze S(z) := F(z + G(z)) für z ∈ ℍ.
F ist ganz und F∘h = id, also hat S fast überall auf (0,1) die Randwerte F(x + θ_T(x)) = F(h(T x)) = T(x). T setzt sich holomorph in eine Umgebung von (0,1) fort, genannt T̃. S − T̃ ist auf
einem Halbkreis über einem Intervall holomorph und hat dort fast überall den nichttangentialen Grenzwert 0. Nach dem Eindeutigkeitssatz von Lusin–Privalov ist S ≡ T̃. Also setzt S die Funktion
T in ℍ fort. Ausserdem gilt S(z+1) = F(z+1+G(z)) = exp S(z), und G → θ̂_T(0) für Im z → ∞, also S → L (F(s) → L für Im s → ∞ bei beschränktem Re s). Spiegelung
liefert die untere Halbebene mit Grenzwert L̄. T ist auf (−2, ∞) reell-analytisch, also ist T holomorph auf Re z > −2, und die vertikalen Grenzwerte sind L und L̄.
Nach Paulsen–Cowgill (kn.satz.paulsen_cowgill) ist T = Kneser. ∎
Die Hilbert-Form: Eine L¹-Funktion f = U + iV hat genau dann keine negativen Moden, wenn V − M[V] die konjugierte Funktion von U ist (Kern cot π(x−y)).
Lücken: θ_T ∈ L¹ und die Majorante für „⇒“ (beide aus der Lokalform an den ganzen Zahlen); Paulsen–Cowgill ist Literatur, am Original nicht nachgeprüft.

**B-S4 (Defekt E).** Für f ∈ L² auf dem Kreis gilt ∬\|f(x)−f(y)\|²/sin²π(x−y) dx dy = Σ_k \|f̂(k)\|²∫₀¹4sin²πkt/sin²πt dt = 4Σ_k\|k\|\|f̂(k)\|² (Fejér-Kern; Douglas-Formel).
Für die umlaufene Fläche gilt A = ½∫Im(f̄ f′) dx = πΣ_k k\|f̂(k)\|². Also ist Douglas/4 − A/π = Σ(\|k\|−k)\|f̂(k)\|² = 2Σ_{k≥1}k\|f̂(−k)\|² ≥ 0. Gleichheit gilt genau dann, wenn keine negativen Moden
da sind, also bei Kneser (S3). Endlichkeit für Kneser: Σk\|a_k\|² < ∞, denn \|a_k\| ≍ 1/(k ln k), die Reihe konvergiert gerade noch. Für allgemeines reell-analytisches T folgt die Endlichkeit des negativen Teils aus Lemma R. ∎

**B-R (Lemma R).** θ_T ist auf ℝ∖ℤ reell-analytisch. Nahe jeder ganzen Zahl ist θ_T Randwert von oben der Funktion Θ(φ(z)) − z, Θ = h∘sexp, die auf einem oberen
Halbkreis holomorph ist. Denn φ ist dort reell-analytisch mit φ′ > 0, bildet also den oberen Halbkreis nach ℍ ab. Zerlege f = θ_T = F₊ + F₋ auf dem Kreis (Riesz-Projektion, F₋ holomorph für \|q\| > 1).
Lokal ist F₋ = (Θ∘φ − id) − F₊ auf dem Bogen. Links steht eine von innen holomorphe Funktion, rechts eine von aussen holomorphe, mit gleichen Randwerten fast überall. Beide liegen lokal in H¹. Nach dem Satz von Painlevé
bzw. Morera für H¹-Randwerte setzt sich F₋ holomorph über den ganzen Kreis fort. Also ist F₋ für \|q\| > ρ₀ < 1 holomorph, und die Koeffizienten fallen
exponentiell. ∎ Lücke: lokale H¹-Zugehörigkeit von Θ∘φ (dieselbe Majorante wie Z0).

**B-MÖ.** B(q) = e^{iα}(q−a)/(1−āq) ist ein Automorphismus der Kreisscheibe, und φ = (1/2πi) log B(q) ist auf ℝ ein reell-analytischer Kreisdiffeo mit φ′ = (1−\|a\|²)/\|1−āq\|². Es gilt
θ(φ(z)) = a₀ + Σ a_n B(q)ⁿ, also holomorph in der ganzen Scheibe, ohne negative Moden. Auf dem Kreis ist
p = φ − x = (1/2πi)[log(1 − a q̄) − log(1 − āq)] + konst. Der erste Summand hat nur negative Moden: −aⁿ/(2πi n) = i aⁿ/(2πn). Daraus folgen N und E durch Summation.
(Normierung p(0) = 0 über α.) ∎

**B-M (Momente).** (i) M[f^k] = M[f]^k für alle k ist gleichwertig mit ∫P dν = P(Mf) für alle Polynome P, wobei ν = f_*(dx). Das hängt nur von ν ab. e^{−2πix} hat alle Momente 0, aber eine
negative Mode. Ebenso erfüllt jede Umordnung G∘ψ (ψ maßtreu, z. B. x ↦ 2x mod 1 oder x ↦ −x) die Familie, wenn G sie erfüllt.
(ii) Nach S1 gilt δV_k = k·ℓ(g^{k−1}w) mit ℓ(F) = (1/2πi)∮_{\|q\|=ρ} F(q)Π(q) dq/q, Π(q) = Σ_{j≥1} conj(p̂(j)) q^{−j}. Das konvergiert für ρ > e^{−2πδ}, wenn p im Streifen
\|Im z\| < δ holomorph ist, und ℓ(qʲ) = conj p̂(j). Ist g auf einer Umgebung von \|q\| ≤ ρ univalent, so ist K = g(\{\|q\| ≤ ρ\}) eine abgeschlossene Jordan-Scheibe mit zusammenhängendem Komplement.
Die Funktion Φ = (qʲ/w)∘g⁻¹ ist auf einer Umgebung von K holomorph (w ≠ 0 nach Satz A, Etappe 3) mit Φ(0) = 0. Nach Runge gibt es Polynome P_n → Φ gleichmässig auf K, mit P_n(0) = 0. Dann
konvergiert P_n(g)·w gleichmässig auf \|q\| = ρ gegen qʲ. ℓ ist dafür stetig, also gilt ℓ(qʲ) = lim Σ_i c_{n,i}ℓ(gⁱw) = 0 und damit p̂(j) = 0 für alle j ≥ 1. ∎
Univalenzradius: Aus Σ_{n≥2} n\|a_n\|r^{n−1} < \|a₁\| folgt \|g′/a₁ − 1\| < 1, also Re(g′/a₁) > 0 auf der konvexen Scheibe, also Univalenz (Noshiro–Warschawski). Das ergibt
r = 0.567 (a_n auf 44 Stellen, n ≤ 4096; der Rest n > 4096 ist bei r = 0.567 kleiner als 10⁻⁹⁰⁰). Numerisch hat g′ keine Nullstelle in \|q\| < 0.995, und die Randkurve g(ρe^{it}) ist für ρ ≤ 0.995
doppelpunktfrei (`s3_univalenz.py`, Gleitkomma). Vermutlich ist g in der ganzen Scheibe univalent. Dann wäre die linearisierte Momentenfamilie für jedes reell-analytische p starr.
Warum (ii) nicht schon die nichtlineare Eindeutigkeit gibt: Die Familie hat unendlich viele Bedingungen, und das Bild der Ableitung ist dicht, aber nicht abgeschlossen. Der Satz der impliziten Funktion greift nicht.

**X (kein Extremum).** Ist Φ in Richtung p differenzierbar mit δΦ[p] ≠ 0, so wechselt Φ(p_t) − Φ(0) bei t = 0 das Vorzeichen. Nach S1 gilt das für Re und Im von a–d
(F̂(1) ≠ 0). Eine Grösse, die Kneser als Extremum auszeichnet, muss bei Kneser verschwindende erste Ableitung haben. E und N tun das (quadratisch in den negativen Moden). ∎

## 4. Einordnung

1. **Warum eine Gleichung nie reicht:** Die Menge der Tetrationen ist „eine Funktion gross“ (p), jede Mittelwert-Gleichung ist eine Zahl. Die Gleichungen schneiden also nur
   Flächen von endlicher Kodimension heraus. Kneser liegt auf allen, aber auch unendlich viele andere. S1 zeigt genauer: Jede Gleichung sieht in erster Ordnung nur
   den „antiholomorphen“ Teil von p, gewichtet mit einer holomorphen Funktion F_C. Erst unendlich viele Gewichte (e_m: F = qᵐw, ein vollständiges System) fangen alles ein.
2. **Bohr–Mollerup-Vergleich:** Für Γ ist log-Konvexität eine Familie von Ungleichungen an jeder Stelle. Hier ist die kennzeichnende Bedingung ebenfalls unendlich
   (alle m) oder als eine Quadratsumme geschrieben (E = 0). Sie ist aber nichtlokal (Hilbert-Kern). Ob es eine *lokale* Kennzeichnung gibt, zum Beispiel ein Vorzeichen von
   Ableitungen von slog, bleibt offen. Die Gegenbeispiele hier zeigen nur, dass analytisch + monoton + endlich viele Integralgleichungen nicht reicht. Das ist die
   Verschärfung von fg.kor.nichteindeutigkeit.
3. **Hooshmand:** Die lineare Tetration ist das Extrem auf der anderen Seite. Sie ist nur C¹ und nach Hooshmand die einzige konvexe auf (−1, 0). Gemessen am Defekt liegt sie
   nahe bei Kneser: N(T_lin) ≈ 2.3·10⁻⁵ entspricht einer Möbius-Verbiegung mit \|a\| ≈ 0.03. Die Konvexitätsbedingung zielt also in eine ganz andere Richtung als
   die Holomorphie-Bedingung. Kneser hat einen Wendepunkt in (−1, 0), Etappe 3.
4. **Neuheit (Vermutung):** Die Aussage „keine negativen Moden ⇔ Kneser“ ist im Kern Paulsen–Cowgill in Hardy-Raum-Sprache. Neu formuliert sind die reelle Form ∫h(T)e^{2πimx} = 1/(2πim),
   das Extremalprinzip E ≥ 0, die exakte Möbius-Formel, die Linearisierung S1 mit den expliziten Gewichten und die infinitesimale Starrheit der Momentenfamilie.
   Literatur nicht geprüft.

## 5. Offen

1. Nichtlineare Eindeutigkeit der Momentenfamilie (gibt es T ≠ Kneser mit M[e^{cθ_T}] = e^{cMθ_T} für alle c?). Mögliche Wege: Univalenz von g auf der ganzen Scheibe beweisen und
   ein Fixpunktargument im Raum der p mit festem Streifen. Oder gezielt nach einer Tetration suchen, deren θ_T eine maßtreue Umordnung einer holomorphen Randfunktion ist.
2. Reichen weniger als alle e_m, zum Beispiel nur gerade m, oder nur Re θ̂_T(−m) = 0? Nicht untersucht.
3. Lokale (Bohr–Mollerup-artige) Bedingung: An den Gegenbeispielen und an Möbius-Familien testen, welche Vorzeichenbedingungen an slog-Ableitungen Kneser von ihnen trennen.
4. Grapheinträge (Vorschlag, nicht eingetragen): S3/S4 als Satz zu br.frage.bohr_mollerup („nichtlokale reelle Kennzeichnung“), S2 als Korollar zu kn.satz.theta_mittelwerte
   (beantwortet dessen f5_neue_frage mit „nein“), MÖ als Beispiel, M als Satz mit offener Frage.

## 6. Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`. Benutzt `../etappe3/lib3.py` und `../etappe3/daten/taylor_130_256.json`, keinen Motor.

| Skript | Inhalt | Laufzeit | Ausgabe |
|---|---|---|---|
| `e4_lib.py` | T = sexp(x+p), θ_T, log(1+θ_T′), Kenngrössen a–d, e_m, V_k; Störungsarten trig, Möbius, linear | – | – |
| `s1_moebius.py 30 5` | Kneser-Kontrolle, Möbius-Formel | 10 s | `daten/log_s1_30.txt` |
| `s2_gegenbeispiele.py lin 30 5` | Linearisierung (S1) | 2 min | `daten/log_s2_lin.txt` |
| `s2_gegenbeispiele.py a\|b\|c\|d\|ab 40 6` | Newton für die Gegenbeispiele | je ~4 min | `daten/gegenbeispiel_*_40.json`, `daten/log_s2_*.txt` |
| `s2b_alle.py 40 6 1e-7` | Gegenbeispiel „alle“ (Sehnenverfahren mit Ableitung aus S1) | 2 min | `daten/gegenbeispiel_alle_40.json`, `daten/log_s2b_alle.txt` |
| `s3_univalenz.py` | Univalenz von θ − a₀ in q | 1 min | `daten/log_s3.txt` |
| `s4_hooshmand.py 30 6 24` | lineare Tetration, Abklingen der negativen Moden | 20 s | `daten/log_s4.txt` |
| `s5_pruefen.py 50 7` | unabhängige Nachprüfung | ~2 min | `daten/log_s5_50_7.txt` |

Hinweis zur Quadratur: Mit k = 5 (321 Knoten) sind die Ableitungsgrössen c und d bei gestörtem T nur auf ~10⁻¹⁸ genau. Ab k = 6 stimmen k = 6 und k = 7 überein. Deshalb wurde mit k = 6 gelöst und mit k = 7 geprüft.
