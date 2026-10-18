# Agent 08 — Lokale, reelle Kennzeichnung von Kneser (Runde 18.10.2026)

## Kurz gesagt (einfache Sprache)

- Gesucht war eine Bedingung, die man **an einzelnen Stellen der Zahlengeraden** prüfen kann (Vorzeichen von Ableitungen, wie die
  Log-Konvexität bei der Gammafunktion), und die genau Kneser auswählt.
- **Hauptergebnis (bewiesen):** Die richtige lokale Grösse ist der **Konvergenzradius der Taylor-Reihe** ρ_T(x) an einer Stelle x, also
  wie schnell die Ableitungen dort wachsen. Für *jede* Tetration gilt ρ_T(x) ≤ x + 2. **Kneser erreicht an jeder Stelle den Höchstwert x + 2.**
  Gleichwertig: an jeder Stelle haben die hohen Ableitungen von Kneser streng abwechselnde Vorzeichen (wie bei log(x+2)).
  „ρ_T(x) = x + 2 für alle x“ ist **genau** dasselbe wie „T ist in der ganzen Halbebene Re z > −2 holomorph“.
- **mike3s alte Vermutung** (alle ungeraden Ableitungen positiv) ist schärfer, als man dachte: Wer sie erfüllt — sogar nur ab irgendeiner Ordnung —,
  ist automatisch in der ganzen Halbebene holomorph (bewiesen, ohne Analytizität vorauszusetzen).
- **Damit ist die Frage genau zurückgeführt:** Diese reellen Bedingungen kennzeichnen Kneser **genau dann**, wenn „holomorph auf Re z > −2“ allein schon
  Kneser erzwingt (Paulsen–Cowgill *ohne* die Grenzwertbedingung bei i∞). Das ist eine bekannte offene Frage; ich konnte sie nicht entscheiden.
- **Alle getesteten Verbiegungen** von Kneser (13 Fälle: Sinus-Störungen bis 10⁻⁶, alle Gegenbeispiele aus Etappe 4, vier Möbius-Verbiegungen)
  haben Singularitäten in der Halbebene; ihr Radius bleibt bei 0.2–1.7 statt x + 2, und schon die 3., 5. oder 15. Ableitung wird negativ.
  Schon eine Störung von 10⁻³ macht T⁽⁵⁾(−0.6) negativ (−0.408 statt +5.409).
- **Ehrliche Negative (bewiesen):** Totale Monotonie von slog′ auf irgendeinem Intervall ist für **jede** Tetration unmöglich; slog′ ist für jede Tetration
  weder log-konkav noch log-konvex (also keine Pólya-Frequenz-Funktion); Vorzeichen fester Ordnung an den Rändern sind für alle Tetrationen gleich;
  endlich viele Vorzeichenbedingungen trennen Kneser nie von nahen Verbiegungen. Kneser hat T′ log-konvex (neu gemessen, Abstand 0.92), aber das gilt
  auch für alle nahen Verbiegungen; die lineare Tetration liegt hier genau auf dem Rand.

## Tabelle der Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit (Vermutung) |
|---|---|---|---|---|
| A | **Pringsheim-Kennzeichnung der Halbebene.** Für eine reell-analytische Tetration T auf (−2,∞) sind gleichwertig: (a) T ist holomorph auf H = {Re z > −2}; (b) an jeder Stelle x gilt (−1)^{k−1}T^{(k)}(x) > 0 für alle grossen k, sogar (−1)^{k−1}k(x+2)^k T^{(k)}(x)/k! → 1; (c) (b) gilt auf einer Folge x_j → ∞; (d) ρ_T(x) = x + 2 für alle x; (e) ρ_T(x) − x → 2 für x → ∞. | **bewiesen** (B-A) | Kneser: Quotient (b) = 1.0000 bei k = 99 für x ≤ 0.15 (`log_e1_N100_dps60.txt`) | neu (Formulierung und Beweis); Levenstein hatte (b) für Kneser beobachtet |
| B | **Extremalprinzip lokal:** ρ_T(x) ≤ x + 2 für jede Tetration; D_T(x) := x + 2 − ρ_T(x) ist ≥ 0 und nicht fallend. Kneser hat D ≡ 0, ist also an jeder Stelle maximal. | **bewiesen** (B-B); Kneser holomorph auf ℂ∖(−∞,−2] ist Literatur | – | neu |
| C | **Ungerade Ableitungen ⇒ Halbebene:** Ist T eine C^∞-Tetration auf (−2,∞) mit T^{(2n+1)} ≥ 0 auf (−2,∞) für alle n ≥ n₀ (irgendein n₀), so ist T reell-analytisch und holomorph auf H. | **bewiesen** (B-C, Lemma von Bernstein über absolut konvexe Funktionen, hier selbst bewiesen) | – | neu |
| D | **Zurückführung:** Jede der Bedingungen A(b)–(e) und C(n₀) kennzeichnet Kneser unter den reell-analytischen Tetrationen **genau dann bzw. dann, wenn** Q* gilt: „jede reelle Tetration, die auf H holomorph ist, ist Kneser“. Für A ist es eine Äquivalenz. Q* ist Paulsen–Cowgill ohne Grenzwert bei i∞ (Halbebenen-Fassung; stärker als kn.verm.tk_c2). | **bewiesen** (Zurückführung); Q* **offen** | – | neu |
| E | Jede Tetration mit einer Singularität in H verletzt A und C(n₀) für jedes n₀; für x → ∞ gilt dann D_T(x) ≥ Re z₀ + 2 − o(1) > 0. | **bewiesen** (aus A, C) | – | neu |
| E1 | T = sexp_K(x + ε sin 2πx) hat logarithmische Verzweigungspunkte z₀ in H (Lösungen von φ(z₀) = −2, auf dem Hauptblatt erreicht); ε = 0.03: z₀ ≈ −0.318 + 0.766i, 0.702 + 0.834i, …; ε = 10⁻³: z₀ ≈ −0.358 + 1.329i, …; Höhe ≈ ln(4(x+2)/ε)/2π. | hergeleitet (Höhe), numerisch (Lage, Blatt: senkrechter Weg kreuzt den Schnitt nicht) | `log_e3_a.txt`, findroot 30 Stellen | – |
| E2 | 13 Verbiegungen (Tabelle 2) verletzen C(0) an 8–20 von 20 Gitterstellen; kleinste verletzende Ordnung 3, 5 oder 15; Radius bei x = 0.3: 0.20–1.74 statt 2.3. Einzige Ausnahme im Messfenster: ε = 10⁻⁶ (Singularität bei Höhe 2.42 liegt ausserhalb aller Kreise mit x ≤ 0.4) — erwartet. | numerisch (dps 50, k < 100); Gegenprobe zweier Wege auf 12 Stellen (`log_t2.txt`) | Tabelle 2 | – |
| E3 | Möbius-Verbiegung (Parameter a) hat Verzweigungspunkte im Abstand δ = ln(1/\|a\|)/2π von der Achse; ρ_T(x) ≤ δ für alle x, also beschränkt (Kneser: x + 2 → ∞). | hergeleitet (B-E3); numerisch: Wurzeltest 0.230 / 0.519 / 0.789 gegen δ = 0.192 / 0.477 / 0.733 (Überschätzung durch endliches k erwartet) | `log_e2_c.txt` | – |
| K | Kneser: alle ungeraden Taylor-Koeffizienten c_{2n+1}(x) > 0 für 2n+1 ≤ 99 auf x ∈ [−1.95, 0.5] (Schritt 0.05) und, soweit zwei Rechenwege übereinstimmen (Ordnung 75–99), auf (0.5, 1.5]. Kleinster normierter Wert c_k k(x+2)^k: 0.952 (k = 1, x ≈ −1.65). | numerisch (dps 60; Erweiterung von gr.verm.ungerade_ableitungen von Ordnung 25 auf 99) | `log_e1_N100_dps60.txt`, `log_e1b.txt` | Erweiterung |
| N1 | Für **keine** reell-analytische Tetration ist slog′ auf irgendeinem offenen Intervall vollständig monoton oder absolut monoton; dasselbe für T′ absolut monoton. T′ vollständig monoton auf (a,b) ⇒ schon auf (−2,b]. | **bewiesen** (B-N1: Bernstein + Pringsheim + „slog nicht ganz“) | – | vermutlich neu (Forum: „slog nicht total monoton“ nur beobachtet) |
| N2 | Für **jede** Tetration ist log slog′ bei +∞ konvex und kann nicht überall konvex sein; also ist slog′ nie log-konkav (keine PF₂-, keine Pólya-Frequenz-Funktion) und nie log-konvex. | bewiesen (Asymptotik B-N2, hergeleitet mit Standard-Induktion) | – | vermutlich neu |
| N3 | **Ränder sind universell:** Für jede feste Ordnung n und jede Tetration: T^{(n)} hat bei −2 das Vorzeichen (−1)^{n−1}, für grosse x ist es positiv; slog^{(n)} hat bei +∞ das Vorzeichen (−1)^{n−1}, bei −∞ ist es positiv. Positivität von T′,…,T^{(n)} auf einem Intervall der Länge 1 überträgt sich nach rechts. | bewiesen (B-N3) | – | teils bekannt (gr.kor.konvex_nach_oben) |
| N4 | **Offenheit:** Für jedes N haben alle hinreichend kleinen Verbiegungen (Möbius \|a\| < a_N, Sinus ε < ε_N) dieselben Vorzeichen und Nullstellenzahlen von T′,…,T^{(N)} auf (−2,∞) wie Kneser. Eine kennzeichnende Vorzeichenbedingung muss also unendlich viele Ordnungen enthalten, und Verletzungen dürfen erst bei Ordnungen → ∞ auftreten. | hergeleitet (Lücke: Einfachheit der Nullstellen von Knesers T^{(k)} nur numerisch) | E2: ε = 0.03 → Ordnung 3, 10⁻³ → 5, 10⁻⁴ → 15, 10⁻⁶ → nicht im Fenster | – |
| L | Kneser: (log T′)″ > 0 auf (−2, ∞), Minimum 0.9196 bei x ≈ −0.04 (gleichwertig: 1/slog′ konvex). Die lineare Tetration hat (log T′)″ = 0 auf (−1, 1). Nicht kennzeichnend (offene Bedingung). | numerisch auf (−2, 2] (Schritt 0.02), rechts davon bewiesen über (log T′)″(x+1) = (log T′)″(x) + T″(x) | `log_e4.txt` | vermutlich neu |
| X | Kleine Störung, grobe Verletzung: ε = 10⁻³ gibt T⁽⁵⁾(−0.6) = −0.408212081901 (Kneser 5.40941046332); ε = 0.03 gibt T‴(0) = −5.98442555504 (Kneser 1.2747). Grund: p^{(k)} ~ ε(2π)^k, bzw. komplex eine Log-Singularität der Stärke 1 im Abstand ~ ln(1/ε)/2π. | numerisch, zwei Wege (Reihenkomposition / mp.diff, 12 Stellen gleich) | `log_t2.txt` | – |

### Tabelle 2: Verbiegungen T = sexp_K(x + p(x)), Gitter x = −1.5, −1.4, …, 0.4, Ordnungen k < 100, dps 50

„sin:ε“ = ε sin 2πx, „cos:ε“ = ε(cos 2πx − 1), „gb:…“ = Gegenbeispiele aus Etappe 4 (`gegenbeispiel_*_40.json`), „moeb:Re:Im“ = Möbius-Parameter a.

| Störung | Gitterpunkte | davon C(0) verletzt | kleinstes verletztes x | kleinste verletzende ungerade Ordnung | Wurzeltest-Radius ρ bei x = 0.3 (Kneser 2.3) |
|---|---|---|---|---|---|
| Kneser | 20 | 0 | – | – | 2.375 |
| sin:0.03 | 20 | 18 | −1.3 | 3 | 0.797 |
| sin:0.001 | 20 | 12 | −0.7 | 5 | 1.353 |
| sin:0.0001 | 20 | 8 | −0.3 | 15 | 1.739 |
| sin:0.000001 | 20 | 0 | – | – | 2.375 |
| cos:0.01 | 20 | 16 | −1.1 | 3 | 1.019 |
| gb:a | 20 | 20 | −1.5 | 3 | 0.419 |
| gb:d | 20 | 20 | −1.5 | 3 | 0.414 |
| gb:alle (max\|p\| = 7.6·10⁻⁷!) | 10 (Zeitgrenze) | 10 | −1.5 | 5 | 0.402 (x = −0.6) |
| moeb:0.3 | 20 | 19 | −1.4 | 3 | 0.378 |
| moeb:0.05 | 20 | 18 | −1.3 | 3 | 0.600 |
| moeb:0.01 | 20 | 17 | −1.2 | 5 | 0.844 |
| moeb:0.2+0.25i | 19 | 19 | −1.5 | 3 | 0.255 |

Der Wurzeltest min_{80≤k<100} |c_k|^{−1/k} überschätzt den wahren Radius um einige Prozent (Faktor k^{1/k}); bei Kneser liefert er 1.048·(x+2).
„gb:alle“ zeigt den Mechanismus deutlich: Die Störung ist winzig, enthält aber die Mode 7 (10⁻⁷ sin 14πx); deren Singularitäten liegen in Höhe
≈ ln(1/ε)/(2π·7) ≈ 0.37.

## Beweise

Bezeichnungen: Basis e, Tetration = Lösung von T(x+1) = e^{T(x)}, T(0) = 1, T′ > 0 auf (−2,∞). H = {Re z > −2}. ρ_T(x) = Konvergenzradius der
Taylor-Reihe von T an x = sup{r : T setzt sich holomorph auf D(x,r) fort}. Für andere Basen b > η gelten A–E wörtlich mit log_b statt log.

**Lemma 0 (Struktur bei −2).** Für jede reell-analytische Tetration ist T(x) = log(x+2) + B(x+2) mit B analytisch in einer Umgebung von 0.
Ist T auf H holomorph, so gilt T(z) = log(z+2) + B(z+2) auch komplex nahe −2.
*Beweis.* T(x+1) = log T(x+2) mit T(x+2) = T(y), y = x+2, T(0) = 1, also T(x+1) = log T(y) = y·U(y), U analytisch, U(0) = T′(0) > 0.
Dann T(x) = log(yU(y)) = log y + log U(y). Komplex: T(z+1) = e^{T(z)} gilt auf H (Identitätssatz), also ist T(z) ein stetiger Logarithmus von
T(z+1) = (z+2)U(z+2); die Differenz zu log(z+2) + log U ist in 2πiℤ, stetig, auf ℝ null. ∎

**B-A.** (a)⇒(b): Fix x > −2. Die Kreisscheibe D(x, x+2) liegt in H, ihr Rand berührt ∂H nur in −2. Nach Lemma 0 ist B = T − log(z+2) holomorph auf
einer Umgebung der abgeschlossenen Scheibe (auf dem Rand ausser −2 liegt man in H, bei −2 ist B analytisch). Also sind die Taylor-Koeffizienten
von B an x O((x+2+η)^{−k}) für ein η > 0, während log(z+2) die Koeffizienten (−1)^{k−1}/(k(x+2)^k) hat. Daraus (−1)^{k−1}k(x+2)^k c_k → 1.
(b)⇒(c) trivial. (c)⇒(a): An x = x_j sei (−1)^{k−1}c_k > 0 für k ≥ k₀. Die Reihe von −T(x−u) in u hat dann ab k₀ positive Koeffizienten. Ist ihr
Radius R endlich, so ist nach Pringsheim u = R singulär, also die Stelle x − R. Die Summe der Reihe stimmt auf (x−R, x] ∩ (−2,∞) mit T überein
(T ist reell-analytisch, Identitätssatz auf der Komponente von D ∩ (Umgebung von ℝ), die x enthält). Wäre x − R > −2, so wäre T dort analytisch
und setzte die Reihe fort — Widerspruch. Also R ≥ x+2. Umgekehrt R ≤ x + 2, weil T(x) → −∞ für x → −2. Somit ist T holomorph auf D(x_j, x_j+2).
Die Scheiben D(x, x+2) sind in x wachsend geschachtelt, die Fortsetzungen stimmen auf den kleineren überein, und ⋃_j D(x_j, x_j+2) = H. ∎
(a)⇔(d): D(x,x+2) ⊂ H und ⋃ = H. (e)⇒(d): ρ ist 1-Lipschitz, also ρ(x) − x nicht wachsend; ≤ 2 nach B, Grenzwert 2 ⇒ ≡ 2.

**B-B.** ρ_T(x) ≤ x+2: sonst wäre T bei −2 analytisch, aber T(x) → −∞. Monotonie von D_T: wie oben, ρ(y) ≤ ρ(x) + (y − x) für y > x.
Kneser ist holomorph auf ℂ∖(−∞,−2] (Literatur, kn.konstr.kneser), also D ≡ 0. ∎

**B-C (ungerade Ableitungen).** *Lemma (Bernstein).* Ist g ∈ C^∞(a,∞) mit g^{(2k)} ≥ 0 für alle k ≥ 0, so ist g an jedem x analytisch mit Radius ≥ x − a.
*Beweis.* Sei δ = x − a, E(u) = (g(x+u)+g(x−u))/2 für |u| < δ. Taylor mit Integralrest der Ordnung 2N: der Rest ist ∫₀^u (u−t)^{2N−1}/(2N−1)! ·E^{(2N)}(t) dt ≥ 0
für u > 0, denn E^{(2N)}(t) = (g^{(2N)}(x+t)+g^{(2N)}(x−t))/2 ≥ 0. Also Σ_{k<N} g^{(2k)}(x)u^{2k}/(2k)! ≤ E(u), und g^{(2k)}(x) ≤ M(2k)!/r^{2k} für jedes r < δ,
gleichmässig für x in kompakten Teilmengen (M = sup g auf einem etwas grösseren Kompaktum). Ungerade Ordnungen: g^{(2k)} ist ≥ 0 und konvex, also
|g^{(2k+1)}(x)| ≤ max(g^{(2k)}(x±h))/h; mit h = r/(2k+1) folgt |g^{(m)}(x)| ≤ C·M·m!/r^m für alle m. Dann geht der Lagrange-Rest gegen 0 für |u| < r. ∎
*Satz C.* Mit g = T^{(2n₀+1)} und a = −2 ist g auf jeder Scheibe D(x, x+2) analytisch, also (geschachtelte Scheiben) holomorph auf H; T entsteht aus g
durch (2n₀+1)-fache Integration auf dem einfach zusammenhängenden H. ∎ Bemerkung: Gilt die Bedingung nur auf (X,∞), folgt Holomorphie auf Re z > X.

**D.** Aus C bzw. A: {T erfüllt C(n₀)} ⊂ Hol(H) und {T erfüllt A(b)} = Hol(H) (für reell-analytische T). Gilt Q*, so besteht Hol(H) nur aus Kneser.
Umgekehrt: Ist T* ≠ Kneser holomorph auf H, so erfüllt T* die Bedingung A(b), A kennzeichnet dann also nicht. Für C(n₀) ist offen, ob ein solches T* auch
die ungeraden Ableitungen positiv hätte (nur „⇐“). Dass Kneser selbst C(0) erfüllt, ist mike3s Vermutung (gr.verm.ungerade_ableitungen), hier bis Ordnung 99 gemessen;
Kneser erfüllt aber beweisbar A(b) an jeder Stelle (für jedes x ab einer Ordnung k₀(x); k₀(0.5) ≈ 406 für die geraden nach Levenstein).

**E.** Ein singulärer Punkt z₀ ∈ H der Fortsetzung verbietet eine holomorphe Fortsetzung auf H (Monodromie auf dem einfach zusammenhängenden H), also
verletzt T nach A die Bedingung A(b) an allen x mit |x − z₀| < x+2, und nach C jede C(n₀). D_T(x) ≥ x + 2 − |x − z₀| → Re z₀ + 2. ∎

**E1 (Sinus-Störung).** φ(z) = z + ε sin 2πz. Für z = x₀ + it, t ≫ 0, ist ε sin 2πz ≈ (iε/2)e^{−2πix₀}e^{2πt}: das Bild des senkrechten Wegs ist nahezu ein
Strahl in fester Richtung, der −2 von oben erreicht, ohne (−∞,−2] zu kreuzen. sexp_K ist auf ℂ∖(−∞,−2] holomorph und hat bei −2 die Log-Singularität
(Lemma 0), φ′(z₀) ≠ 0, also T ≈ log(φ′(z₀)(z−z₀)): Verzweigungspunkt. Höhe aus |ε sin| ≈ |z₀+2|: t ≈ ln(2|z₀+2|/ε)/2π. Lücke: „nahezu Strahl“ nur für kleine ε
quantifiziert; numerisch für ε = 0.03 und 10⁻³ geprüft (Weg in 400 Schritten, kein Kreuzen).

**E3 (Möbius).** p enthält −(1/2πi)log(1 − āq) (Etappe 4, B-MÖ); bei q = 1/ā (Im z = −δ) geht φ logarithmisch nach −i∞, und ein Umlauf ändert φ um ±1.
Also T ↦ exp^{±1}(T): echte Verzweigung, T → L̄ beim Anlauf. Die Punkte liegen periodisch in jeder Periode, darum ρ_T(x) ≤ δ für alle x.

**B-N1.** Ist f = slog_T′ auf (α,β) vollständig monoton, so ist u ↦ f(β′−u) (β′ < β) absolut monoton, nach Bernstein also eine Potenzreihe mit
Koeffizienten ≥ 0. Pringsheim: der Konvergenzradius R wird bei einer reellen Stelle angenommen; slog_T ist auf ganz ℝ reell-analytisch (slog_T(w) = slog_T(e^w) − 1),
also R = ∞ und slog_T ganz. Dann gilt slog_T(e^w) = slog_T(w) + 1 auf ganz ℂ, bei w = L (e^L = L) folgt 0 = 1. Absolut monoton analog (Pringsheim-Punkt rechts).
Für T′ absolut monoton liegt der Pringsheim-Punkt rechts, T wäre ganz, Widerspruch zu T → −∞ bei −2; für T′ vollständig monoton auf (a,b) kann der Punkt
links nur −2 sein, also gilt die Monotonie dann schon auf (−2,b]. ∎ (Für Kneser selbst: T″ < 0 nur auf (−2, −0.5218), also höchstens b ≤ −0.5218.)

**B-N2.** Aus slog(w) = 1 + slog(log w): (log slog′)″(w) = (1 − G′(ℓ) + G″(ℓ))/w² mit ℓ = log w, G = log slog′; G′(ℓ) → 0, G″ → 0 (Induktion über die
Iteration), also > 0 für grosse w. Bei −∞: slog′(w) = e^w slog′(e^w), Vorzeichen von (log slog′)″ ist das von slog″(0) = −T″(−1)/T′(−1)³. Eine auf ℝ konvexe
Funktion, die für w → +∞ gegen −∞ geht (log slog′ ~ −log w), wäre überall fallend, kann also für w → −∞ nicht gegen −∞ gehen (log slog′ ~ w). ∎

**B-N3.** Bei −2: T^{(n)}(x) = (−1)^{n−1}(n−1)!/(x+2)^n + B^{(n)} (Lemma 0). Rechts: (e^T)^{(n)} = e^T·Y_n(T′,…,T^{(n)}) mit dem vollständigen Bell-Polynom Y_n
(nur nichtnegative Koeffizienten), also überträgt sich Positivität von T′…T^{(n)} von [x₀,x₀+1] auf [x₀+1,x₀+2] usw.; positive Lage entsteht, weil T′ unbeschränkt wächst.
slog bei ±∞ über Faà di Bruno mit slog^{(j)}(w) ≈ (−1)^{j−1}(j−1)! slog′(log w)/w^j·(1+O(1/log w)) (hergeleitet).

**N4.** Für Möbius |a| → 0 gilt T_a → T_K mit allen Ableitungen gleichmässig auf Kompakta von (−2,∞); am linken Rand konvergiert B_a → B_K nahe 0, am rechten Rand
hält strikte Positivität auf einem Intervall der Länge 1 und überträgt sich (N3). Sind Knesers Nullstellen von T^{(k)}, k ≤ N, einfach (numerisch: jede gerade
Ableitung bis 26 wechselt genau einmal, gr.verm.ungerade_ableitungen), bleiben Zahl und Lage stabil. Lücke: die Einfachheit ist nur gemessen.

**L.** (log T′)(x+1) = (log T′)(x) + T(x) ⇒ (log T′)″(x+1) = (log T′)″(x) + T″(x). Für x ≥ 0 ist T″ > 0 (Wendepunkt bei −0.5218), also überträgt sich
Positivität von (0,1] nach rechts. Auf (−2, 2] gemessen (Schritt 0.02, dps 40). Lineare Tetration: T′ = 1 auf (−1,0), T′ = e^x auf (0,1), also (log T′)″ = 0 dort.

## Ehrliche Einordnung

- Die gesuchte lokale Bedingung **existiert genau dann, wenn Q* stimmt**; A macht das zu einer Äquivalenz. Die Frage „gibt es eine Bohr–Mollerup-artige Kennzeichnung“
  ist damit nicht beantwortet, aber auf eine präzise komplexe Frage zurückgeführt, die schon im Forum offen ist (Open-problems-Survey; tid 882).
- Der lokale Charakter ist echt, aber wie bei Bohr–Mollerup braucht man **x → ∞**: An einer festen Stelle oder auf einem beschränkten Intervall reicht keine
  Bedingung (N4; Möbius-Verbiegungen mit kleinem a sind auf jedem beschränkten Gebiet holomorph und erfüllen dort A(b)).
- Ich habe keine Verbiegung gefunden, die A oder C erfüllt. Alle Verbiegungen der Form sexp_K∘(id + p) mit ganzem, nicht konstantem p haben (soweit geprüft)
  Verzweigungspunkte in H, weil z + p(z) = −2 in jedem Periodenstreifen ~2N Lösungen hat (Argumentprinzip; N = höchste Mode). Ob alle davon auf Nebenblättern
  liegen können (dann wäre T* holomorph auf H), habe ich nicht ausgeschlossen.
- Rechenfehlerquelle: Für x > 0.5 sind hohe Taylor-Koeffizienten aus den gespeicherten Daten unzuverlässig (Auslöschung im exp-Weg, Faktor ((x+2)/(x+1))^k);
  `rk_lib.kneser` verweigert das ohne `erlaube_exp=True`. Alle Verbiegungs-Rechnungen bleiben bei φ(x) ≤ 0.5.

## Vorschläge für Grapheinträge

- „Kneser hat an jeder Stelle den grösstmöglichen Taylor-Radius x + 2; ρ_T(x) = x+2 ∀x ⇔ holomorph auf Re z > −2“ — satz — bewiesen
- „Positive ungerade Ableitungen (ab beliebiger Ordnung) erzwingen Holomorphie auf Re z > −2“ — satz — bewiesen (Kante zu gr.verm.ungerade_ableitungen)
- „Lokale reelle Kennzeichnung ⇔ Paulsen–Cowgill ohne Grenzwert (Halbebene)“ — satz (Zurückführung) — bewiesen; Q* als frage — offen (Kanten br.frage.bohr_mollerup, kn.verm.tk_c2, kn.satz.paulsen_cowgill)
- „slog′ ist für keine Tetration auf einem Intervall total monoton; nie log-konkav (keine Pólya-Frequenz)“ — satz — bewiesen
- „Endlich viele Vorzeichenbedingungen trennen Kneser nie von nahen Verbiegungen; Ränder universell“ — satz — hergeleitet
- „Kneser: T′ log-konvex, min (log T′)″ = 0.9196 bei x ≈ −0.04; lineare Tetration auf dem Rand“ — beobachtung — numerisch
- „Sinus- und Möbius-Verbiegungen: Verzweigungspunkte in Höhe ln(1/ε)/2π bzw. ln(1/|a|)/2π; ungerade Ableitungen ab Ordnung 3–15 negativ“ — beispiel — numerisch
- gr.verm.ungerade_ableitungen: Messung von Ordnung 25 auf 99 erweitern (x ∈ [−1.95, 1.5]).

## Offene Fragen

1. **Q*:** Ist jede reelle Tetration, die auf Re z > −2 holomorph ist, Kneser? (Ansatz: θ_T = h∘T − z ist dann auf der oberen Halbebene bis auf
   harmlose Log-Punkte (T = L) periodisch holomorph; zu zeigen wäre, dass q = 0 keine wesentliche Singularität sein kann.)
2. Erfüllt jedes T ∈ Hol(H) auch C(0)? Gilt C(0) für Kneser (mike3)? Für hohe Ordnungen folgt es aus A, offen sind mittlere Ordnungen bei x > 1.
3. Ist für jede Tetration ≠ Kneser sogar sup_x ρ_T(x) < ∞ (wie bei allen Beispielen: beschränkt oder O(log x))? Dann wäre schon „ρ_T unbeschränkt“ kennzeichnend.
4. Gibt es ein Gegenstück auf der slog-Seite (Radius von slog_T an w gegen |w − L|)? Für ganze Störungen p̃ der slog-Seite entstehen keine neuen Singularitäten, nur
   stärkere bei L — dort müsste man die Art, nicht die Lage der Singularität messen.

## Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`, im Ordner `08_reelle_kennzeichnung`. Nur gelesen: `../../identitaeten_2026_10_14/etappe3/daten/taylor_130_256.json`
und `../../identitaeten_2026_10_14/etappe4/daten/gegenbeispiel_*_40.json`. Kein Motor, kein Zwischenspeicher.

| Skript | Inhalt | Laufzeit | Ausgabe |
|---|---|---|---|
| `rk_lib.py` | Reihenarithmetik, Kneser-Taylor an reellen x (Verschiebung ≤ 0.1, log-Weg links), T = sexp_K∘φ | – | – |
| `stoerungen.py` | Katalog der Störungen, p komplex | – | – |
| `t0_kontrolle.py`, `t1_stoer_kontrolle.py`, `t2_gegenprobe.py` | Kontrollen (Zentren, Reihen von p gegen Cauchy-Taylor 10⁻⁴⁸, zwei Wege für Ableitungen) | je < 20 s | `daten/log_t2.txt` |
| `e1_kneser_ungerade.py 100 60 0.05` | Kneser, alle Vorzeichen k < 100 | ~3 min | `daten/log_e1_N100_dps60.txt`, `.json` |
| `e1b_rechts.py` | Kneser x ∈ (0.5, 1.5] über zwei Zentren | ~1 min | `daten/log_e1b.txt` |
| `e2_stoerungen.py 100 50 <namen>` | Verbiegungen | ~1 min je Störung | `daten/log_e2_{a,b,c}.txt`, `tabelle_e2.py` → `daten/tabelle_e2.md` |
| `e3_singularitaeten.py <namen>` | Lösungen von φ = −2, Blattprüfung | ~1 min | `daten/log_e3_a.txt` |
| `e4_logkonvex.py` | (log T′)″ für Kneser | 10 s | `daten/log_e4.txt` |

Hinweis: Der Lauf `log_e2_b.txt` wurde bei „gb:alle“ nach 10 Gitterpunkten durch die 10-Minuten-Grenze beendet (deshalb keine JSON-Datei dazu; die Protokollzeilen sind vollständig).
