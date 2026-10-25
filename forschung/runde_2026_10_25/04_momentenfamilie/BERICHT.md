# Gebiet 04 (Runde 25.10.2026): Legt die Momentenfamilie Kneser auch nichtlinear fest?

## Kurz gesagt

Die Frage war: Wenn für eine Tetration T der Mittelwert jeder Potenz von θ_T gleich der Potenz des Mittelwerts ist, muss T dann Kneser sein?
Ich zeige: **Ja, für jede Tetration, deren Kurve θ_T sich nicht selbst schneidet.** Der Trick ist, die unendlich vielen Gleichungen in eine
einzige Funktion zu packen (die Cauchy-Transformierte). Diese muss ausserhalb der Kurve genau 1/(ζ − c) sein, und daraus folgt, dass θ_T der Rand
einer holomorphen Funktion ist, also Kneser. Ausserdem gibt es eine **ausdrückliche Umgebung**: Jede Tetration T = sexp(x + p(x)) mit
−0.23 < p′ < 0.43 erfüllt die Momentenfamilie nur, wenn sie Kneser ist. Die alten Gegenbeispiele aus Etappe 4 liegen ausserhalb dieser Umgebung,
zwei davon (a und c) haben tatsächlich einen Doppelpunkt. Offen bleibt nur der Fall, in dem sich die Kurve selbst schneidet.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| C | **Cauchy-Form:** Momentenfamilie ⇔ C_T(ζ) := ∫₀¹ dx/(ζ − θ_T(x)) = 1/(ζ − c) auf der äusseren Komponente Ω_T von ℂ∖Γ_T (c = M[θ_T]). | bewiesen, Lücke nur bei der Lage „Re θ_T nach unten beschränkt“ (aus der Lokalform, Etappe 2) | Kneser: \|C − 1/(ζ−c)\| ≤ 7·10⁻²² bei ζ = −5, 20+5i, 25+10i, 26−3i (dps 20); Gegenbeispiel a (nur V₂ = 0): 2.4·10⁻¹⁰ bei ζ = −5, passt zu V₃/(ζ−c)⁴ ≈ 2.7·10⁻¹⁰ | neu |
| B | **Bogen-Bild:** θ_T(x) = u(x) − x, wobei u = u_K∘φ den **festen** Bogen Γ₀ = h((1, e)) monoton durchläuft (u_K(s) = s + θ(s)). Alle Tetrationen haben denselben Bogen, nur die Geschwindigkeit ändert sich. | bewiesen (Definition, h∘T reell-parametrisiert) | – | Umformulierung |
| D | **Doppelpunkt-Kriterium:** θ_T(x) = θ_T(y) mit x > y ⇔ es gibt eine waagrechte Sehne (t < s, Im θ(t) = Im θ(s)) der Kneser-Kurve mit ψ(s) − ψ(t) = (s − t) − μ(s,t), μ := Re θ(t) − Re θ(s), ψ = φ⁻¹. Die Menge der ψ mit ψ(s) − ψ(t) > (s−t) − μ für alle Sehnen ist **konvex** und enthält id. | bewiesen | – | neu |
| K | Kneser-Kurve Γ_K: doppelpunktfrei, positiv umlaufen (Fläche +0.0762, Umlaufzahl +1 um c), jede waagrechte Gerade schneidet sie höchstens zweimal, μ > 0 auf allen Sehnen. | numerisch (N = 400 dps 15 und N = 800 dps 18 identisch); Spitzenbereich x < 10⁻¹² nur analytisch (Lokalform) | `log_s4_*.txt`, `log_s2_kneser_18_800.txt` | – |
| m* | Sehnenkonstante m* := inf μ(s,t)/min(s−t, 1−s+t) = **0.3010** (erreicht bei t = 0.2414, s = 0.7412, Höhe 5.4532). | numerisch, 4 Stellen in beiden Einstellungen gleich | `log_s4_15_400.txt`, `log_s4_18_800.txt` | neu |
| S | **Sehnensatz:** Lip(φ⁻¹ − id) < m* ⇒ θ_T doppelpunktfrei. Gleichwertig hinreichend: 0.7686 < φ′ < 1.4306, also −0.231 < p′ < 0.431. | bewiesen aus D, mit m* numerisch und Spitze analytisch (Lücke: Spitzenabschätzung nur qualitativ) | Gegenprobe: p = ε sin 2πx bleibt bis ε = 0.09 (p′ bis ±0.57) doppelpunktfrei – das Kriterium ist hinreichend, nicht scharf | neu |
| J | **Jordan-Satz:** Ist θ_T auf (0,1) injektiv (Γ_T Jordankurve durch ∞) und gilt die Momentenfamilie, dann ist θ_T(x) = G(e^{2πi(x+β)}) mit der konformen Abbildung G: 𝔻 → D_T, G(0) = c. Also hat θ_T keine negativen Moden, und T = Kneser. | hergeleitet mit Lücke (Hardy-/Smirnov-Klasse von qΦ(G)G′, Nullstellen von θ_T′, Orientierung an der Spitze; dazu Paulsen–Cowgill über Etappe 4 S3) | Beweis unten | neu |
| N | **Nichtlineare Eindeutigkeit in expliziter Umgebung:** Jede reell-analytische Tetration T = sexp∘φ mit 0.7686 < φ′ < 1.4306 und M[θ_Tᵏ] = M[θ_T]ᵏ für alle k ist Kneser. Ohne Streifenbedingung an p (anders als Etappe 4, M(ii)). | hergeleitet mit Lücke (Lücken von J und S) | – | neu |
| X | Die Gegenbeispiele a (erfüllt V₂ exakt) und c aus Etappe 4 haben je **einen Doppelpunkt**: a bei x ≈ 0.371 / 0.613 (Punkt ≈ 26.1923 + 5.4732i), c bei x ≈ 0.341 / 0.638 (≈ 26.1961 + 5.4640i), jeweils eine kleine Schleife oben; bei a ist max p′ ≈ 0.60 > 0.43. Sie liegen also ausserhalb des Jordan-Falls, und Satz J sagt über sie nichts. | numerisch (N = 400 dps 15; Gegenprobe mit anderen Knoten k = 4, dps 20: ebenfalls je 1 Schnitt) | `log_s2_ac_15_400.txt`, `log_s1_20_4.txt` | – |
| O | Im Fall mit Doppelpunkten ist das Argument **nicht** abgeschlossen: Auf einer sich schneidenden Kurve gibt es viele darstellende Masse für c (Mischungen harmonischer Masse verschiedener Teilgebiete). Ein Gegenbeispiel habe ich nicht gefunden und nicht gesucht. | negativ (Methode reicht dort nicht) | – | – |
| R | Die zwei naheliegenden Konstruktionen für Gegenbeispiele scheitern: (i) θ_T = G(β − x) (nur negative Moden) widerspricht der Orientierung an der Spitze; (ii) θ_T = G∘ψ mit Grad-n-Abbildung ψ (n ≥ 2) bräuchte n Spitzen pro Periode, θ_T hat nur eine. | hergeleitet (i: Lokalform an der Spitze; ii: bewiesen) | – | – |

## Beweise und Einzelheiten

Bezeichnungen wie Etappe 4: Basis e, θ = h∘sexp − id, T = sexp∘φ, φ = id + p, θ_T = h∘T − id, M = ∫₀¹, c = M[θ_T], Γ_T = θ_T((0,1)).
Für x → 0± läuft θ_T über die Lokalform θ_T(x) = 2 − x + h(log(A φ′(0) x) + O(x)) (Etappe 2, S6) sehr langsam nach ∞ („Spitze“), nach rechts:
Re θ_T → +∞, |θ_T| ≲ log log log(1/|x|). Damit ist M[e^{λ|θ_T|}] < ∞ für alle λ (Etappe 3) und Re θ_T ≥ −R_T.

**Satz C.** Für Re ζ < −R_T gilt 1/(ζ − w) = −∫₀^∞ e^{t(ζ−w)} dt für alle w auf Γ_T. Fubini (Integrand ≤ e^{t(Re ζ + R_T)}) gibt
C_T(ζ) = −∫₀^∞ e^{tζ} Λ(t) dt mit Λ(t) = M[e^{−tθ_T}] = Σ(−t)ᵏM[θ_Tᵏ]/k! (dominiert durch M[e^{t|θ_T|}]). Unter der Momentenfamilie ist Λ(t) = e^{−tc},
also C_T(ζ) = 1/(ζ − c), und durch Fortsetzung auf ganz Ω_T. Umgekehrt: Stimmt C_T auf einer linken Halbebene mit 1/(ζ − c) überein, so folgt aus der
Eindeutigkeit der Laplace-Transformation Λ(t) = e^{−tc} für t > 0, und da Λ ganz ist, für alle t; die Taylor-Koeffizienten sind die Momente. ∎
Folge: c ∉ Ω_T (sonst hätte die holomorphe Funktion C_T dort einen Pol).

**Satz D.** θ_T(x) = θ_T(y) ⇔ u(x) − u(y) = x − y. Mit s = φ(x), t = φ(y) heisst das: u_K(s) − u_K(t) ist reell und gleich ψ(s) − ψ(t).
Im u_K = Im θ, also liegen θ(s), θ(t) auf gleicher Höhe, und Re(u_K(s) − u_K(t)) = (s − t) − μ(s,t). Die Bedingung „> für alle Sehnen“ ist linear in ψ, also konvex. ∎

**Satz S.** r := ψ − id ist 1-periodisch, also |r(s) − r(t)| ≤ Lip(r)·min(s−t, 1−s+t) < m*·d ≤ μ(s,t). Eine Gleichheit r(s) − r(t) = −μ wie in D ist
damit unmöglich. An der Spitze (t → 0, s → 1) ist μ ≍ 1/(X ln X) mit X = ln(1/x), also viel grösser als d ≍ x; das ist der Grund, warum m* > 0 ist,
aber dieser Teil ist nur über die Lokalform begründet, nicht numerisch abgetastet (Lücke). Lip(ψ − id) = sup|1/φ′ − 1|. ∎

**Satz J.** D_T sei die Komponente von Ĉ∖Γ_T, die die Spitze von innen enthält (nicht Ω_T). *Orientierung:* Nahe x = 0 ist θ_T(x) = θ(φ′(0)x) + O(x),
und der Abstand der beiden Seiten der Spitze (≍ 1/(X ln X)) ist viel grösser als die Störung O(x); daher liegt D_T wie bei Kneser **links** der Laufrichtung.
*Plemelj:* C_T(ζ) = ∫_Γ f(z) dz/(ζ − z) mit f = 1/θ_T′. An glatten Punkten springt C_T von links nach rechts um −2πi f. Mit C_T = 1/(ζ − c) auf Ω_T hat
Φ(ζ) := 1/(ζ − c) − C_T(ζ) auf D_T einen einfachen Pol bei c (Residuum 1) und Randwerte Φ(θ_T(x))·θ_T′(x) = 2πi.
Mit der konformen Abbildung G: 𝔻 → D_T, G(0) = c, und F(q) := q·Φ(G(q))·G′(q) gilt F(0) = 1, und auf dem Kreis (q = e^{iτ}) ist
F dτ = 2π dx, also F reell und ≥ 0 (Orientierung). Ein holomorphes F mit reellen Randwerten ist konstant (Im F harmonisch mit Randwert 0), also F ≡ 1 und
dx/dτ = 1/2π: G⁻¹(θ_T(x)) = e^{2πi(x+β)}. Damit ist θ_T Randwert der holomorphen Funktion G(e^{2πiβ}q), hat keine negativen Moden, und nach
Etappe 4, S3 (Paulsen–Cowgill) ist T = Kneser. ∎
Lücken: (a) „Im F hat Randwert 0 ⇒ Im F ≡ 0“ braucht F in der Smirnov-Klasse; D_T ist unbeschränkt mit sehr dünner Spitze, G liegt in H^p für p < ½, F
habe ich nicht abgeschätzt. (b) Isolierte Nullstellen von θ_T′ (Spitzen der Kurve) sind ausgenommen, das ändert nur eine Nullmenge. (c) G ∈ H¹ für
„keine negativen Moden“; θ_T ∈ L¹ ist bekannt. (d) Orientierung an der Spitze nur über die Lokalform.

**Warum die naheliegenden Gegenbeispiele nicht gehen (R).** Für allgemeine f ist die Familie schwächer als „keine negativen Moden“ (f = G∘ψ, ψ maßtreu).
Für Tetrationen ist Γ_T aber eine Kurve mit genau einer Spitze pro Periode und fester Orientierung. Ist sie doppelpunktfrei, legt Satz J die Werteverteilung
als harmonisches Mass von D_T fest, und das harmonische Mass einer Jordankurve hat nur die Parametrisierungen G∘(Drehung) und G∘(Spiegelung). Die Spiegelung
kehrt die Orientierung um. Eine maßtreue Abbildung vom Grad n ≥ 2 hätte n Urbilder der Spitze.

**Was offen bleibt (O).** Tetrationen mit Doppelpunkt gibt es schon bei kleinen Störungen (Gegenbeispiel a: \|p\| ≤ 0.065, aber p′ bis 0.60). Dort zerfällt
ℂ∖Γ_T in mehrere Teilgebiete U_j. Nach Satz C ist C_T = 1/(ζ−c) auf Ω_T, auf jedem U_j ist C_T holomorph, und an jeder Kante gilt die Sprungformel
summiert über die Blätter. Für eine nichtlineare Eindeutigkeit dort müsste man zeigen, dass sich die Sprünge nicht zu einem zweiten darstellenden Mass
zusammensetzen lassen. Ein möglicher Weg: Umlaufzahlfunktion n(ζ) der Kurve und ein Argument wie in J auf der „Riemannschen Fläche über Γ_T“
(dem immersierten Gebiet, falls es eines gibt; Blank'sches Problem).

## Vorschläge für Grapheinträge

- `br.satz.momente_cauchy` — „Momentenfamilie ⇔ Cauchy-Transformierte der Werteverteilung ist 1/(ζ−c) ausserhalb der Kurve“ — Satz — bewiesen (Lage-Lücke klein).
- `br.satz.momente_jordan` — „Momentenfamilie + doppelpunktfreie θ_T-Kurve ⇒ Kneser“ — Satz — hergeleitet mit Lücke (Smirnov-Klasse, Spitze).
- `br.satz.sehnenkriterium` — „θ_T doppelpunktfrei, wenn sup\|1/φ′ − 1\| < m* = 0.3010; konvexe Menge der doppelpunktfreien ψ“ — Satz — bewiesen bis auf m* numerisch.
- `br.kor.momente_umgebung` — „Momentenfamilie legt Kneser in der Umgebung 0.7686 < φ′ < 1.4306 fest“ — Korollar — hergeleitet mit Lücke; beantwortet br.frage.momentenfamilie teilweise.
- `br.beob.gegenbeispiel_a_doppelpunkt` — „Gegenbeispiel a (V₂ = 0) hat einen Doppelpunkt“ — Beobachtung — numerisch.
- `br.frage.momente_doppelpunkt` — „Gibt es T ≠ Kneser mit Doppelpunkt-Kurve und Momentenfamilie?“ — Frage — offen.

## Offene Fragen

1. Lücke (a) in J schliessen: Smirnov-Klasse für F = q·Φ(G)·G′ bei der dünnen Spitze (vermutlich über die Lokalform und |θ_T| ≲ log log log).
2. Fall mit Doppelpunkt: zweites darstellendes Mass ausschliessen oder ein Gegenbeispiel bauen (z. B. Störungen mit sin 4πx und grossem p′, an denen
   man V₂, …, V_K gleichzeitig nullsetzt und schaut, ob die Lösung mit wachsendem K konvergiert).
3. Scharfe Konstante: Die konvexe Menge aus D ist viel grösser als die Lipschitz-Kugel (p = 0.09 sin 2πx ist noch doppelpunktfrei). Einseitige Version:
   für s − t ≤ ½ zählt nur inf r′, für s − t > ½ nur sup r′.
4. Gilt derselbe Satz für Basis 2 und andere Basen über η (gleiche Lokalform, eigenes m*)?

## Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`. Kein Motor; benutzt (nur lesend) `identitaeten_2026_10_14/etappe4/e4_lib.py`
mit den Taylor-Daten aus Etappe 3.

| Skript | Inhalt | Laufzeit | Ausgabe |
|---|---|---|---|
| `s1_kurve.py 20 4` | Kneser, Gegenbeispiele a, c, Möbius 0.3: Schnitte, Fläche, Cauchy-Transformierte gegen 1/(ζ−c) | 1 min | `daten/log_s1_20_4.txt`, `daten/kurven_20_4.json` |
| `s2_fein.py <dps> <N> <fall>` | feine Abtastung, Doppelpunkte, Umlaufzahlen (fall: kneser, a, c, d, moeb0.3, sin0.06, sin20.03) | 1–2 min | `daten/fein_*.json`, `daten/log_s2_*.txt` |
| `s3_neck.py <json>…` | engste Stellen zwischen fernen Kurvenstücken | Sekunden | Bildschirm |
| `s4_sehnen.py <json> <ny>` | waagrechte Sehnen, Schnittzahl je Höhe, m* | Sekunden | `daten/log_s4_*.txt` |
