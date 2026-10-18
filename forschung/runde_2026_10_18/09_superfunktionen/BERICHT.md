# 09 Superfunktionen: Der Zoo der Kneser-Superfunktionen und welche Gesetze universell sind (18.10.2026)

## Kurz gesagt (einfache Sprache)

- **Worum es geht.** Knesers Tetration ist eine *Superfunktion* von eˣ: S(z+1) = e^{S(z)}. Dieselbe Bauweise geht für andere Abbildungen f
  ohne reellen Fixpunkt, aber mit einem Paar komplexer Fixpunkte L, L̄. Ich habe einen allgemeinen Löser gebaut (Cauchy-Integral über einen
  Streifen, Gauss–Newton) und damit die Kneser-Superfunktionen von **x²+1, x⁴+1 und cosh x** berechnet, auf rund 15 Stellen.
  Zur Kontrolle liefert derselbe Löser für eˣ die bekannten Werte (sexp(½), θ-Koeffizient a₁, Wendepunkt) auf 10⁻¹⁵ bis 10⁻¹⁴.
- **Das allgemeine Gesetz.** Die θ-Welle ist genau dort singulär, wo S einen **singulären Wert von f** trifft, also einen kritischen Wert
  (Stelle mit f′ = 0) oder einen asymptotischen Wert (wie 0 bei eˣ). Einen *ausgelassenen* Wert der regulären Superfunktion braucht es nicht.
  Bei Polynomen lässt die reguläre Superfunktion **keinen** Wert aus (bewiesen), und θ ist trotzdem singulär.
- **Die Art der Singularität bestimmt den Abfall der Koeffizienten.** Ein kritischer Wert vom Grad d erzeugt eine d-te-Wurzel-Singularität.
  Dann fallen die θ-Koeffizienten wie **aₙ ≈ C₁·n^{−1−1/d}**, mit einer geschlossenen Konstante C₁. Für x²+1 und cosh (d = 2) stimmt
  die Vorhersage mit zwei Gliedern auf 0.65/n² (relativ), für x⁴+1 (d = 4) mit drei Gliedern auf 0.24/n. Der Logarithmus bei eˣ ist der
  Grenzfall d = ∞, und dann gilt das bekannte Gesetz aₙ ≈ −h′(c − ln n)/n mit Exponent genau 1.
- **Universell (bei allen vier Abbildungen bestätigt):** keine negativen Fourier-Moden von θ, Multiplikativität der Periodenmittel
  (∫θ² = (∫θ)², σ-Produkt = 1), ∫ slog′²/h′ = 1, ∫ log(1+θ′) = 0, genau ein Wendepunkt. Die Gegenprobe mit einer um 10⁻⁴ verbogenen
  Superfunktion verletzt jedes dieser Gesetze um 10⁻⁶ bis 10⁻⁴.
- **Nicht universell:** „S′ hat in der oberen Halbebene keine Nullstelle“ gilt bei eˣ, aber **nicht** bei Polynomen. Die Superfunktion von
  x²+1 wird bei z = 2.17993 + 2.59729i null, und dann ist S′(z+1) = 0. Die Mittelwertgesetze überleben trotzdem, weil dafür nur 1+θ′ ≠ 0
  nötig ist, und das gilt weiter.
- **Neues Phänomen: Polynome haben eine zweite natürliche reelle Superfunktion**, die Böttcher-Superfunktion φ⁻¹(exp(G·dᶻ)). Ich beweise,
  dass sie nicht Knesers Lösung ist. Numerisch liegt sie aber unglaublich nah daran: Für x²+1 unterscheiden sich die beiden nur um eine
  Verschiebung von höchstens **4.5·10⁻⁷**, für x⁴+1 um 2.2·10⁻⁴. Der Kneser-Defekt ist E = 3.9·10⁻¹⁴ bzw. 1.5·10⁻⁸.

## Tabelle der Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| Z0 | Allgemeiner Kneser-Löser (Kouznetsov-Kontur + Gauss–Newton) für jede ganze reelle f mit Hauptast g: ℍ → ℍ | Methode | eˣ: sexp(½) auf 2·10⁻¹⁶ (Einst. B), a₁ von θ auf 2·10⁻¹⁵ gegen Etappe 1, Wendepunkt −0.52176479262332 gegen Zertifikat auf 1.7·10⁻¹⁴; zwei Einstellungen A/B ≤ 7·10⁻¹⁵ für alle Abbildungen | Verfahren bekannt (Kouznetsov 2009); Gauss–Newton-Form und Anwendung auf x²+c, x⁴+c, cosh vermutlich neu |
| P1 | Für f = x^d + c (c ≠ 0) lässt die Poincaré-Funktion Ψ (also die reguläre Superfunktion F) **keinen Wert aus** | **bewiesen** (B-P1) | – | elementar, wohl bekannt |
| U1 | **Singulärer Ort:** θ ist auf ℝ genau an den Stellen singulär, an denen S einen Wert der Bahn eines singulären Werts von f trifft; im Zoo ist das z₁ + ℤ mit S(z₁) = v | hergeleitet mit Lücke (Zweigwahl wie S4 der Vorarbeit) | negative Moden ≤ 1.1·10⁻¹⁵, Periodizitätsfehler ≤ 2·10⁻¹⁴; Puiseux-Rest glatt wie 1/n² (keine weiteren Singularitäten auf \|q\| = 1) | Verallgemeinerung von kn-Etappe 2 (dort nur eˣ); vermutlich neu |
| U2 | **Puiseux-Gesetz:** Bei einem kritischen Wert vom Grad d ist θ(z) = Φ((z−z₁)^{1/d}) mit Φ holomorph, also aₙ ~ Σ_{d∤k} C_k n^{−1−k/d} mit **C₁ = h′(w_c)·(iA/(2πa))^{1/d}/Γ(−1/d)** (A = S′(z₁), f(w_c + X) = v + aX^d + …) | **bewiesen** unter (V) (B-U2: Transfer nach Flajolet–Odlyzko, hier wörtlich anwendbar) | x²+1: a/(Vorhersage) = 0.99993 bei n = 1024, zweigliedrig Rest 0.65/n²; cosh gleich; x⁴+1: dreigliedrig Rest 0.24/n | vermutlich neu |
| U3 | Log-Fall (asymptotischer Wert) ist der Grenzfall d → ∞: aₙ ≈ −h′(c − ln n)/(n ln b) | bekannt (Etappe 2/Lücken), hier mit neuem Löser reproduziert | eˣ: \|a₁₀₂₄\| = 2.2013·10⁻⁵ wie dort | – |
| U4 | Multiplikativität der Periodenmittel, σ-Produkt, ∫slog′²/h′ = 1, M[log(1+θ′)] = 0, keine negativen Moden: **universell** | hergeleitet (Z0-Satz der Vorarbeit; Majorante bei d < ∞ trivial, weil θ beschränkt ist) | Tabelle 2: 10⁻¹⁵…10⁻¹⁰, Gegenprobe ε = 10⁻⁴ verletzt um 10⁻⁶…10⁻⁴ | Anwendung auf andere Abbildungen neu |
| A1 | **S′ ≠ 0 in ℍ ist nicht universell:** S_K(x²+1) hat die Nullstelle 2.1799322204542 + 2.5972861621226i, S_K(x⁴+1) die Nullstellen 1.0435923930395 + 1.3498560554128i, 2.15828 + 1.11638i, 2.40764 + 1.28052i; dort ist S(z+1) = v und S′(z+1) = 0 | numerisch (Newton, \|S\| < 10⁻¹⁵, A/B gleich auf 13 Stellen) | eˣ: keine Nullstelle von S oder S′ in [−6, 2] × (0, 8] (Argumentprinzip) | vermutlich neu |
| A2 | 1 + θ′ ≠ 0 in {Im z > 0.01} für alle vier Abbildungen | numerisch (Windungszahl 0 auf 10⁻¹⁵; \|1+θ′\| ≥ 0.69) | `w_windung.json` | neu formuliert |
| A3 | Re(1+θ′) > 0 in ℍ für x²+1, cosh, x⁴+1 (min Re ≈ 0.72…0.78 bei y ≥ 0.003; Lokalform beschränkt arg); **falsch für eˣ**: dort Re(1+θ′) < 0 für \|z\| ≲ 10⁻⁶ | Vermutung (d < ∞, numerisch + Lokalanalyse) / numerisch (eˣ, Lokalform) | arg(1+θ′) → [arg h₁ − (1−1/d)π, arg h₁]; eˣ: max\|arg\| = 2.06 bei \|z\| = 10⁻³⁰ | neu |
| W | Genau ein Wendepunkt auf (z_c, ∞) für alle vier Abbildungen | rechnergestützt (Vorzeichen von S″ auf (−1.99, 1], Fortpflanzung exakt) | x²+1: −1.26011657467806, cosh: −0.74758444374635, x⁴+1: −1.30858991524120, eˣ: −0.52176479262334 (A/B ≤ 4·10⁻¹⁵) | x²+1, cosh, x⁴+1 neu |
| B1 | **Böttcher ≠ Kneser:** Für f = x^d + c (d gerade, f(x) > x auf ℝ) ist die Böttcher-Superfunktion S_B = φ⁻¹(exp(G(c)dᶻ)) reell-analytisch, aber nicht die Kneser-Lösung | **bewiesen** (B-B1, Green-Funktion + Cantor-Menge) | numerisch: θ̂_B(−1) = 1.4·10⁻⁷ ≫ Rauschen 10⁻¹⁶ | vermutlich neu |
| B2 | S_B(x) = S_K(x − 1 + p(x)) mit winzigem p: x²+1: max\|p\| = 4.476914·10⁻⁷, p̂(1) = 8.379·10⁻⁸ − 1.1223·10⁻⁷i; x⁴+1: max\|p\| = 2.199189·10⁻⁴. Kneser-Defekt E(S_B) = 3.92·10⁻¹⁴ bzw. 1.46·10⁻⁸ | numerisch (A/B gleich auf 6·10⁻¹⁷, zwei Höhen gleich auf 10⁻¹⁶) | `boettcher.json` | neu |
| B3 | Grösse des Abstands ≈ e^{−2πδ} mit effektiver Streifenbreite δ = 2.39…2.63 (d = 2, c = 2…0.5) und 1.49 (d = 4); grob passend zu π/(2 ln d) = 2.27 bzw. 1.13 | Vermutung, schwache Belege (5 Fälle) | `log_boettcher_scan1/2.txt` | neu |
| K1 | f(x) = eˣ − 1 + u (u > 0) ist affin konjugiert zu bˣ mit ln b = e^{u−1}: kein neuer Fall | bewiesen (eine Zeile) | – | trivial |
| K2 | f(x) = x + eˣ hat überhaupt keinen Fixpunkt in ℂ, gehört also nicht zur Kneser-Klasse (parabolisch bei −∞) | bewiesen (trivial) | – | trivial |

## 1. Rahmen und Löser

**Voraussetzungen (V).** f ist ganz und reell. Es gibt einen abstossenden Fixpunkt L ∈ ℍ (λ = f′(L), |λ| > 1, arg λ > 0) und einen Ast g
von f⁻¹, holomorph auf ℍ mit g(ℍ) ⊂ ℍ. Nach Denjoy–Wolff ist L dann der einzige Fixpunkt von g in ℍ, und gⁿ → L. g ist auf einem Intervall
(v, ∞) reell und wachsend mit g((v, ∞)) = (w_c, ∞). Dabei ist v ein singulärer Wert: entweder ein kritischer Wert vom Grad d (f(w_c + X) = v + aX^d + …)
oder ein asymptotischer Wert (w_c = −∞). Die Kneser-Superfunktion ist S = F(z + θ(z)) mit F(s) = Ψ(λˢ), θ holomorph und 1-periodisch in ℍ.

| f | g (Hauptast) | v, d, a | L | λ |
|---|---|---|---|---|
| eˣ | log w | 0, ∞ (asymptotisch) | 0.3181315 + 1.3372357i | = L |
| x²+1 | √(w−1) | 1, 2, 1 | ½ + (√3/2)i | 1 + √3 i (\|λ\| = 2, arg = π/3) |
| x⁴+1 | (w−1)^{1/4} | 1, 4, 1 | 0.7271361 + 0.4300143i | −0.0756512 + 2.4102617i |
| cosh x | arcosh w | 1, 2, ½ (kritischer Punkt 0) | 0.9760545 + 0.8709515i | 0.7333614 + 1.1591777i |

Bei allen ist g univalent auf ℍ, also ist die Schröder-Funktion σ = lim λⁿ(gⁿ − L) auf ℍ schlicht (Hurwitz). Das ist dasselbe Argument wie Lemma U der Vorarbeit.

**Normierung:** S(−1) = v. Also gilt S(−2) = w_c (Verzweigungspunkt, bei eˣ der log-Pol), S(0) = f(v). Für eˣ ist das Knesers sexp.

**Löser (`zoo.py`):** S(z) = (1/2π)∫_ℝ [f(u(t))/(1+it−ζ) − g(u(t))/(−1+it−ζ)] dt mit u(t) = S(it), Spiegelung für t < 0 und einem Schwanzglied ab T.
Unbekannt sind u an 640 (A: T = 40, Paneele 1, 16 Knoten) bzw. 1072 (B: T = 50, Paneele 0.75) Gauss-Knoten. Gauss–Newton in reellen
Variablen mit der Normierungszeile. Residuum ≤ 5·10⁻¹⁵, Rechenzeit 10–50 s. Werte: x²+1: S(½) = 2.99455081706702, S(−½) = 1.41228567119652;
x⁴+1: S(½) = 4.11874746910240; cosh: S(½) = 1.90455288767422 (A gegen B ≤ 7·10⁻¹⁵).
h = log σ/log λ: g-Iteration je Punkt bis |w − L| < 0.02, dann Schröder-Reihe (24 Glieder) per Newton.
(Eine erste Fassung iterierte bis 10⁻⁹ an L heran und verlor dabei 7 Stellen. Behoben.)

## 2. Der singuläre Ort und das Puiseux-Gesetz

**B-P1 (Ψ lässt bei x^d + c nichts aus).** Ψ ist ganz, weil Ψ(u) = fⁿ(Ψ(u/λⁿ)). Angenommen, Ψ lässt a aus. Wegen Ψ(λu) = f(Ψ(u)) lässt Ψ dann
alle Urbilder f⁻¹(a) aus. Nach Picard ist das höchstens ein Wert, also a = c (einziges Urbild 0). Dann muss Ψ auch die d verschiedenen
Urbilder von 0 auslassen, die Wurzeln von −c. Das sind mindestens zwei Werte, ein Widerspruch. ∎
Bei eˣ lässt Ψ dagegen 0 aus. Für den θ-Abfall ist das aber nicht entscheidend, wie U1/U2 zeigen.

**U1.** θ(z) = h(S(z)) − z. h setzt sich längs jedes Weges fort, der die Bahn P der singulären Werte meidet (σ ist aus den Ästen von f⁻¹ gebaut,
und die sind nur über singulären Werten verzweigt). S ist auf (−2, ∞) reell und wachsend. Im Zoo gilt P ∩ ℝ = {v, f(v), f²(v), …}.
Bei cosh ist auch der Wert −1 singulär (kritischer Punkt iπ), aber S trifft −1 auf ℝ nie, und cosh(−1) = cosh(1) liegt schon auf der Bahn.
Also ist θ auf ℝ ∖ (z₁ + ℤ) analytisch.
Die Lücke ist dieselbe wie S4 der Vorarbeit: Dass der Hauptast die richtige Fortsetzung ist, ist nur numerisch belegt (stetige Abtastung, Periodizitätsfehler ≤ 2·10⁻¹⁴).

**B-U2 (Puiseux-Form).** In ℍ gilt S(z−1) = g(S(z)), also θ(z) = 1 − z + h(g(S(z))) mit dem von oben fortgesetzten g.
Nahe z₁ ist S(z) − v = A(z−z₁) + B(z−z₁)² + … mit A > 0. Dann ist g(S(z)) − w_c = ((S−v)/a)^{1/d}·(analytisch in (S−v)^{1/d}).
Also ist θ(z) = Φ((z−z₁)^{1/d}) mit Φ holomorph nahe 0, und h ist bei w_c analytisch, weil w_c ∉ P.
Mit z − z₁ = log q/(2πi) = (i/2π)(1−q)(1 + (1−q)/2 + …) ist θ eine konvergente Reihe in (1−q)^{k/d}, und ganzzahlige Potenzen tragen nicht bei.
Nach U1 ist θ Δ-analytisch bei q = 1. Der Transfersatz für (1−q)^α (Flajolet–Odlyzko) gilt hier wörtlich, anders als im log-Fall:

  **aₙ = Σ_{k ≥ 1, d ∤ k} E_k · [qⁿ](1−q)^{k/d},  [qⁿ](1−q)^α = Γ(n−α)/(Γ(−α)Γ(n+1)),  E₁ = h′(w_c)(iA/(2πa))^{1/d}.**

Für d = 2 kommt das zweite Glied von w^{3/2}: E₃ = h₁√(iA/2πa)·[¼ + (B/2A − δ)·i/2π] + h₃(iA/2πa)^{3/2}, mit δ = A/12 bei cosh
(aus arcosh(1+u) = √(2u)(1 − u/12 + …)) und δ = 0 bei x²+c. Für x^d + c ist g exakt die d-te Wurzel.
Für 1 ≤ k < d ist dann E_k = h_k (iA/2π)^{k/d} + O(nächste Ordnung).

Numerik (`s2_theta.py`, `s3_puiseux.py`; Koeffizienten aus FFT auf Höhe y = 0.001 und 0.002, 16384 Punkte, Einstellungen A/B):

| f | A = S′(z₁) | h′(w_c) | C₁ = E₁/Γ(−1/d) | a₁ |
|---|---|---|---|---|
| x²+1 | 0.735213815530658 | 0.456190785716 + 0.501269092189i | 0.0030758457867 − 0.0653307358217i | 0.0029141086356 − 0.0606079753466i |
| cosh | 0.529371357442456 | 0.389661625777 + 0.415092697509i | 0.0020823322071 − 0.0658944248945i | 0.0018743443167 − 0.0611957758682i |
| x⁴+1 | 0.522186950532817 | 0.108291291134 + 0.431203298705i | 0.0071163214164 − 0.0481774593987i | 0.0062299041222 − 0.0798442956031i |

| n | x²+1: aₙ/(C₁n^{−3/2}) | x²+1: n²·Rest(2 Glieder) | cosh: n²·Rest | x⁴+1: aₙ/(C₁n^{−5/4}) | x⁴+1: n·Rest(3 Glieder) | eˣ: aₙ/(−h′(c−ln n)/n) |
|---|---|---|---|---|---|---|
| 16 | 0.99522 | 0.745 | 0.740 | 1.328 | 0.299 | 0.750 |
| 64 | 0.99880 | 0.675 | 0.671 | 1.225 | 0.271 | 0.817 |
| 256 | 0.99970 | 0.659 | 0.655 | 1.154 | 0.254 | 0.858 |
| 512 | 0.99985 | 0.657 | 0.653 | 1.128 | 0.248 | 0.873 |

(Gemeint ist jeweils der relative Rest. Der Imaginärteil des Quotienten ist ≤ 10⁻⁴. Die Koeffizienten selbst sind bis n = 512 auf ≤ 10⁻⁹ relativ sicher: A/B und zwei Höhen.
Bei eˣ konvergiert das Leitglied nur wie 1/ln n, wie in Etappe 2.)
Der lokale Exponent ist 1.4999 bei d = 2. Bei d = 4 ist er 1.2734 bei n = 1024 und strebt langsam gegen 1.25 (Korrektur n^{−1/4}).
**Folgerung:** Die Drift des lokalen Exponenten (eˣ: 1.25 → 1.12, wie 1 + 1/ln n) ist eine Besonderheit des log-Falls.
Bei algebraischen kritischen Werten liegt der Exponent sofort bei 1 + 1/d.

**Einheitliche Form:** Für jeden singulären Wert v mit lokaler Umkehrung g_v gilt aₙ ≈ [qⁿ] h(g_v(v + A·i(1−q)/2π)).
Für g_v = d-te Wurzel ergibt das n^{−1−1/d}, für g_v = log das Superlog-Gesetz. Heuristisch ist (1 + x/d)^d → eˣ, und der kritische Wert 0
vom Grad d geht in den asymptotischen Wert über (n^{−1−1/d} → n^{−1}). Das habe ich nicht nachgerechnet.

## 3. Universelle Mittelwert-Gesetze (`s4_mittel.py`)

Periode [−1, 0]. Die Substitution x = −1 + sin²(πs/2) glättet die Wurzelenden. Gauss–Legendre mit 300 bzw. 200 Knoten, Lösung A bzw. B.
Gegenprobe: T(x) = S(x + 10⁻⁴ sin 2πx).

| Grösse (soll 0) | eˣ | x²+1 | cosh | x⁴+1 | Gegenprobe ε = 10⁻⁴ |
|---|---|---|---|---|---|
| M[θ²] − M[θ]² | 5·10⁻¹² | 7·10⁻¹⁴ | 8·10⁻¹³ | 1·10⁻¹³ | 6–9·10⁻⁶ |
| M[e^{iθ}]/e^{iM[θ]} − 1 | 3·10⁻¹² | 1·10⁻¹⁴ | 5·10⁻¹⁴ | 6·10⁻¹⁴ | 3–5·10⁻⁶ |
| σ-Produkt − 1 | 9·10⁻¹² | 5·10⁻¹⁵ | 6·10⁻¹⁵ | 3·10⁻¹³ | 7·10⁻⁶–3·10⁻⁵ |
| M[1/(1+θ′)] − 1 (= ∫slog′²/h′ − 1) | 3·10⁻¹⁵ | 6·10⁻¹⁵ | 1·10⁻¹⁴ | 1·10⁻¹⁴ | 2.4–3.6·10⁻⁴ |
| M[log(1+θ′)] | 6·10⁻¹⁰ | 3·10⁻¹⁰ | 3·10⁻¹⁰ | 5·10⁻¹⁰ | 1.2–1.8·10⁻⁴ |
| θ̂(−1) | 4·10⁻¹² | 2·10⁻¹⁵ | 8·10⁻¹⁵ | 1·10⁻¹³ | 5·10⁻⁵ |

Bei eˣ ist die Quadratur wegen des Superlog-Endes etwas schlechter, bei log(1+θ′) wegen des log-Endes ebenso. Das sind Grenzen des Verfahrens, keine Abweichungen.
Gegenprobe zur Vorarbeit: ∫₋₁⁰ ln S′ = −7.786·10⁻⁴ plus ∫₋₁⁰ sexp = 0.49921818283 ergibt 0.4984395595, wie ∫₀¹ ln sexp′ in Etappe 3.
**Begründung der Universalität:** Der Mittelwertsatz Z0 braucht nur θ ∈ Hol(ℍ), Periode 1, Grenzwert oben und eine integrierbare Majorante.
Bei d < ∞ ist θ beschränkt (Puiseux), die Majorante also trivial. Die Gesetze mit 1+θ′ brauchen zusätzlich 1+θ′ ≠ 0 in ℍ (A2).
Die Kennzeichnung durch keine negativen Moden und der Defekt E aus Etappe 4 übertragen sich wörtlich. Die Eindeutigkeit (Paulsen–Cowgill)
ist dabei für allgemeine f vorausgesetzt, nicht bewiesen.

## 4. Was nicht universell ist: Nullstellen von S′ (`s5*.py`)

S′(z+1) = f′(S(z))·S′(z). Eine Nullstelle von S in ℍ erzeugt also bei x²+c und cosh eine Nullstelle von S′ bei z+1.
Bei eˣ ist S ≠ 0, weil F den Wert 0 auslässt, und die Vorarbeit bewies S′ ≠ 0 in ℍ.
Für Polynome lässt F nichts aus (P1), und tatsächlich: S_K(x²+1) = 0 bei 2.1799322204542 + 2.5972861621226i.
Bei x⁴+1 gibt es drei solche Stellen in Re ∈ [0, 2.5]. Dort ist S(z+1) = 1 = v und |S′(z+1)| < 10⁻¹³.
Das Argumentprinzip für S und S′ in den Rechtecken [x₀, x₀+1] × [0.02, 8], x₀ = −5.9 … 1.1, ergibt bei eˣ, x²+1 und cosh nur Nullen.
Bei cosh fand ich in Re ≤ 3 keine Nullstelle; ob es weiter rechts welche gibt, ist offen.
**Warum die Mittelwertgesetze trotzdem gelten:** S′ = F′(z+θ)·(1+θ′). Die Nullstellen von S′ stammen von Nullstellen von F′, also von
kritischen Punkten von Ψ, die es bei eˣ nicht gibt. 1+θ′ hat in {Im z > 0.01} keine Nullstelle (Windungszahl 0) und erfüllt |1+θ′| ≥ 0.69.
Der richtige universelle Satz heisst also: **z ↦ z + θ(z) ist in ℍ lokal umkehrbar**, nicht „S′ ≠ 0“. Für eˣ sind beide gleichwertig.

**A3 (Re(1+θ′) > 0).** Nahe z₁ gilt 1+θ′ ≈ (1/d)·h₁(A/a)^{1/d}·ζ^{1/d−1}. Daher liegt arg(1+θ′) in [arg h₁ − (1−1/d)π, arg h₁].
Für x²+1 ist das [−0.74, 0.83], für cosh [−0.75, 0.82], für x⁴+1 [−1.03, 1.33], also immer innerhalb (−π/2, π/2).
Auf den Abtastlinien y ≥ 0.003 ist min Re(1+θ′) = 0.78 (x²+1), 0.77 (cosh), 0.72 (x⁴+1).
Bei eˣ wandert arg(1+θ′) ≈ arg h′(log(Az)) − arg z dagegen unbeschränkt langsam: max|arg| = 1.76 bei |z| = 10⁻⁸, 2.28 bei |z| = 10⁻²⁵⁰.
Re(1+θ′) wird also sehr nahe an den ganzen Zahlen negativ.
Vermutung: Bei kritischen Werten endlichen Grads ist z + θ(z) in ℍ „close-to-convex“ (Re(1+θ′) > 0), bei eˣ nicht.

## 5. Ein Wendepunkt (`s6_wendepunkt.py`)

Für alle Zoo-Abbildungen gilt S > 0 auf (−2, ∞), und f′, f″ > 0 auf (0, ∞). Aus S″(x+1) = f″(S)S′² + f′(S)S″ folgen zwei Regeln:
Ist S″(x) ≥ 0, so ist S″(x+1) > 0. Ist S″(x+1) ≤ 0, so ist S″(x) < 0, solange x > −2. Das ist dasselbe Lemma wie bei eˣ.
S″ wechselt auf (−1.99, 1] genau einmal das Vorzeichen, an der Stelle x_w, und ist auf (x_w + 0.05, 1] ≥ 0.019 > 0.
Zusammen mit der Fortpflanzung und S″ < 0 nahe −2 (Wurzel- bzw. log-Verhalten) gibt es genau einen Wendepunkt:

| f | x_w (Normierung S(−1) = v) | A gegen B |
|---|---|---|
| eˣ | −0.521764792623339 (Zertifikat: −0.521764792623322) | 2·10⁻¹⁴ |
| x²+1 | −1.260116574678058 | 2·10⁻¹⁵ |
| cosh | −0.747584443746352 | 3·10⁻¹⁵ |
| x⁴+1 | −1.308589915241198 | 4·10⁻¹⁵ |

Status: rechnergestützt in Gleitkomma, nicht mit Intervallen zertifiziert.
Die Fortpflanzung ist allgemein: Gilt f′, f″ > 0 auf S((z_c, ∞)), so reicht ein Vorzeichenwechsel auf einem Stück der Länge > 1.

## 6. Neues Phänomen: die Böttcher-Superfunktion (`s7_boettcher.py`, `s8_scan.py`)

Für f = x^d + c (d gerade, f(x) > x auf ℝ) gibt es die Böttcher-Koordinate φ mit φ(f(w)) = φ(w)^d und φ(w) ~ w bei ∞.
Damit ist S_B(z) := φ⁻¹(exp(G(c)·dᶻ)) eine reelle, reell-analytische Superfunktion mit S_B(0) = c und S_B(−1) = 0 (Wurzel-Verzweigung), genau wie S_K.
G ist die Green-Funktion; G(1) = 0.40735452273948 für d = 2 und 0.17707602108109 für d = 4.
Sie entspricht der regulären Iteration am superattraktiven Fixpunkt ∞. Bei eˣ gibt es nichts Vergleichbares (∞ ist wesentlich singulär).

**B-B1 (S_B ≠ S_K, Beweis).** Die kritische Bahn 0 → c → … läuft nach ∞, weil f(x) > x. Also ist die gefüllte Julia-Menge K eine Cantor-Menge mit K ∩ ℝ = ∅.
G ist stetig, ≥ 0, harmonisch auf ℂ∖K, gleich 0 auf K, und auf der Böttcher-Umgebung gilt G = log|φ|.
Angenommen S_B = S_K. Dann gilt G(S_K(z)) = G(c)·Re dᶻ auf einer Umgebung der reellen Halbgeraden.
Sei U₀ die Zusammenhangskomponente von {z ∈ ℍ : S_K(z) ∉ K}, die einen Streifen über (−2, ∞) enthält. Den gibt es, weil S_K(ℝ) ∩ K = ∅.
Auf U₀ sind beide Seiten harmonisch, also gleich. Wegen G > 0 auf U₀ folgt Re dᶻ > 0 auf U₀, also liegt U₀ im Band B₀ = {0 < Im z < π/(2 ln d)}.
An Randpunkten von U₀ in ℍ ist S_K ∈ K, also G∘S_K = 0, also Re dᶻ = 0. Im offenen Band B₀ gibt es solche Punkte nicht.
U₀ ist daher in B₀ offen und abgeschlossen, also U₀ = B₀. Auf der Oberkante y₀ = π/(2 ln d) gilt dann aus Stetigkeit G(S_K(x+iy₀)) = 0 für alle x.
Die zusammenhängende Kurve x ↦ S_K(x+iy₀) liegt also in der total unzusammenhängenden Menge K und ist konstant. Dann ist S_K konstant, ein Widerspruch. ∎
(Gebraucht wird nur: S_K holomorph auf ℍ, nicht konstant, und reell auf ℝ.)

**B2 (Numerik).** S_B(x) = S_K(x − 1 + p(x)) mit 1-periodischem, reell-analytischem p. Die Werte stimmen in A und B auf 6·10⁻¹⁷ überein.
Die negativen Moden von θ_B = h∘S_B − id wurden auf zwei Höhen Y = 0.03 und 0.06 bestimmt (Übereinstimmung 10⁻¹⁶).

| f | max\|p\| | \|p̂(1)\| | \|θ̂_B(−1)\|, \|θ̂_B(−2)\|, \|θ̂_B(−3)\| | N = Σ\|θ̂(−m)\|² | E = 2Σm\|θ̂(−m)\|² |
|---|---|---|---|---|---|
| x²+1 | 4.476914·10⁻⁷ | 1.4006·10⁻⁷ | 1.401·10⁻⁷, 8.2·10⁻¹⁴, (1.3·10⁻¹⁶ Rauschen) | 1.96·10⁻¹⁴ | 3.92·10⁻¹⁴ |
| x⁴+1 | 2.199189·10⁻⁴ | 8.5462·10⁻⁵ | 8.546·10⁻⁵, 4.57·10⁻⁸, 1.06·10⁻¹¹ | 7.30·10⁻⁹ | 1.46·10⁻⁸ |

Scan d = 2: |p̂(1)| = 6.65·10⁻⁸ (c = 0.5), 1.40·10⁻⁷ (c = 1), 3.02·10⁻⁷ (c = 2).
Die effektive Breite −ln|p̂(1)|/2π beträgt also 2.63, 2.51, 2.39 (d = 2) und 1.49 (d = 4).
Das passt grob zur Bandbreite π/(2 ln d) = 2.27 bzw. 1.13 aus dem Beweis, mit c-abhängigem Vorfaktor.
Ehrlich: Läufe mit c = 4 sowie d = 4, c = 2 und d = 6 scheiterten an Überlauf in der Böttcher-Reihe bzw. im Newton-Schritt.
c = 0.3 ist unzuverlässig, weil der Schwanz bei T = 40 nur auf ~5·10⁻⁸ genau ist (arg λ = 0.42). B3 bleibt daher eine schwache Vermutung.
**Deutung:** Zwei völlig verschieden gebaute natürliche Superfunktionen, eine von L/L̄ her, eine von ∞ her, stimmen bei x²+1 auf 6–7 Stellen überein.
Für Rechnungen auf der reellen Achse sind sie praktisch gleich, und doch ist nur eine davon „Kneser“.
Der Defekt E aus Etappe 4 misst den Unterschied direkt.

## 7. Vorschläge für Grapheinträge

- „Kneser-Superfunktionen anderer Abbildungen: θ singulär genau über den singulären Werten von f“ — Satz (hergeleitet mit Lücke Zweigwahl)
- „Puiseux-Gesetz aₙ ~ h′(w_c)(iA/2πa)^{1/d} n^{−1−1/d}/Γ(−1/d) für kritische Werte vom Grad d“ — Satz (bewiesen unter (V)), numerisch für x²+1, cosh, x⁴+1
- „Die reguläre Superfunktion eines Polynoms x^d + c lässt keinen Wert aus“ — Satz (bewiesen)
- „Mittelwert-Gesetze von Kneser gelten universell (x²+1, cosh, x⁴+1)“ — Beobachtung (numerisch 10⁻¹⁰–10⁻¹⁵)
- „S′ ≠ 0 in ℍ ist eine Besonderheit von eˣ; universell ist 1+θ′ ≠ 0“ — Beobachtung (numerisch; Nullstelle 2.17993+2.59729i bei x²+1)
- „Genau ein Wendepunkt bei den Superfunktionen von x²+1, cosh, x⁴+1“ — Beobachtung (rechnergestützt)
- „Böttcher-Superfunktion eines Polynoms ist nicht Knesers Lösung“ — Satz (bewiesen)
- „Böttcher und Kneser unterscheiden sich bei x²+1 nur um 4.5·10⁻⁷“ — Beobachtung (numerisch)
- „Re(1+θ′) > 0 bei algebraischen kritischen Werten, nicht bei eˣ“ — Vermutung
- „eˣ − 1 + u ist konjugiert zu bˣ mit ln b = e^{u−1}“ — Korollar (bewiesen, trivial)

## 8. Offene Fragen

1. Eindeutigkeitssatz (Paulsen–Cowgill) für allgemeine f mit (V): Ist jede in ℍ holomorphe Superfunktion mit S → L reell nach Normierung gleich S_K?
2. Konvergieren die Kneser-Superfunktionen von (1+x/d)^d für d → ∞ gegen sexp, und wie geht n^{−1−1/d} gleichmässig in das Superlog-Gesetz über?
3. Genaue Asymptotik des Böttcher-Kneser-Abstands: Wo liegen die Singularitäten von p, und gilt ln|p̂(1)| = −2πδ(c, d) + O(1) mit explizitem δ?
4. Hat S_K(cosh) Nullstellen in ℍ (weiter rechts)? Der Mechanismus (F′ = 0, weil Ψ′ = 0 über den Urbildern von iπk) lässt es erwarten.
5. Re(1+θ′) > 0 beweisen (d < ∞), etwa über Schlichtheit von σ und die Lokalform.
6. Abbildungen mit **mehreren** reellen singulären Werten auf S(ℝ) (etwa f mit zwei reellen kritischen Werten): Dann überlagern sich zwei Puiseux-Reihen mit verschiedenen d.
   Abbildungen ohne reellen singulären Wert (etwa x + 1 + ε sin x, ein Kreisabbildungs-Lift) ergäben ein θ ohne reelle Singularität mit exponentiellem Abfall.
   Das passt aber wegen unendlich vieler Fixpunkte L + 2πk nicht in (V). Nicht untersucht.

## Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, Ordner `09_superfunktionen/`, Daten in `daten/`.

| Schritt | Aufruf | Dauer |
|---|---|---|
| Löser + Kontrolle eˣ | `t_exp.py`; `s1_loesen.py q1 q4 cosh`; `s1_loesen.py exp` | je 1–2 min |
| θ-Koeffizienten | `s2_theta.py <abb> 16384 0.001 0.002` (abb = exp, q1, cosh, q4) | je ~1 min |
| Puiseux-Vorhersage | `s3_puiseux.py` | < 1 min |
| Mittelwerte | `s4_mittel.py` | ~3 min |
| Nullstellen | `s5_nullstellen.py`, `s5b_rechts.py` (bricht bei cosh am Überlauf ab, Ergebnisse für q4/q1 vorher gedruckt), `s5c_w.py`, `s5d_pruef.py` | je < 10 min |
| Wendepunkte | `s6_wendepunkt.py` | ~2 min |
| Böttcher | `s7_boettcher.py`; `s8_scan.py 1` (bricht bei c = 4 ab), `s8_scan.py 2` | 1–5 min |
