# Runde 18.10.2026 — Agent 05: Kombinatorik und erzeugende Funktionen

## Kurz gesagt

- **Die ganze Störungsreihe hat eine geschlossene Bauform.** Die Reihe G(y) = e^{G(y−1)} − 1 + u ist (formal) der Fluss einer Differentialgleichung H′ = λ_u(H). Dabei ist λ_u der „Erzeuger“ der Abbildung eˣ − 1 + u. Seine zwei Nullstellen sind die Fixpunkte x± = u − 1 − W(−e^{u−1}) (Lambert-W). An beiden Fixpunkten ist das Residuum genau 1/x±. Daraus folgt: Jede Schicht ist ein Polynom in tan θ, ln cos θ und θ. Die höchsten Logarithmus-Anteile summieren sich zu 2·tan(θ + v·ln cos θ / 3). Ein nacktes θ kommt erst in Schicht 5 vor; es stammt aus der genauen Frequenz ω = v(1 − v⁴/270 + …).
- **Dritte Schicht geschlossen** (bewiesen, geprüft): H₂ = (u/2)^{3/2}·tan θ·[(2/9)sec²θ·ln²cos θ − 2 sec²θ·ln cos θ + (41 tan²θ + 39)/9].
- **Bernoulli-Zahlen:** Der in u lineare Teil des Erzeugers ist λ₀·Σ B_k⁺ (λ₀ d/dx)^k (1/λ₀)/k!. Dabei ist λ₀ die Julia-Funktion von eˣ − 1.
- **Taylor-Koeffizienten von Kneser-sexp bei 0:** Die Differenz sexp(z) − log(z+2) ist holomorph bis zum Radius 3. Also gilt c_n = (−1)^{n+1}/(n·2ⁿ) + O((3−ε)^{−n}), und der Rest hat eine genaue Asymptotik über log log bei −3. Die Vorzeichen sind für n ≤ 21 positiv und ab dann (geprüft bis n = 100) streng wechselnd. Überraschung: Bei 0 sieht man die Asymptotik des Rests erst für sehr grosse n. Der Grund ist das riesige Wachstum von sexp auf der positiven Achse. Vom Zentrum −0.4 aus stimmt die Vorhersage auf 0.7 % (n = 80).
- **Die Folgen von x·eˣ sind jetzt erklärt** (bewiesen und durch Abzählen geprüft für n ≤ 6): Das Dreieck aus Etappe 1 zählt Ketten im Verband der punktierten Partitionen von Chapoton–Vallette, die in einem Maximum enden. Die Möbius-Zahl ist (−1)^{n−1}n^{n−2}.
- **Pentation:** Der nächste Rang ist an dieser Stelle hyperbolisch und nicht mehr parabolisch. Es entsteht also keine neue Baum-Kombinatorik; die Bäume aus Bäumen liefern nur die Kettenformel. Das ist ein ehrliches Negativergebnis.
- **Neue ganzzahlige Folgen** (nicht in der OEIS, am 18.10.2026 geprüft): das Ketten-Dreieck 2; 3, 12; 4, 84, 144; 5, 510, 2760, 2880; 6, 3210, 41940, 117360, 86400, seine Zeilensummen 2, 15, 232, 6155, 248916 und die zweite Schicht 1, 3, 26, 468, 14256, 657568, …

Interpreter: `C:/Users/janis/miniconda3/python.exe`. Kein Motorlauf. Die Taylor-Daten werden nur gelesen, aus `identitaeten_2026_10_14/etappe3/daten/`.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit |
|---|---|---|---|---|
| S1 | H(y) = G(y−1) ist formal der Fluss H′ = λ_u(H). Dabei ist λ_u der graduierte Erzeuger mit exp(λ_u∂)x = eˣ − 1 + u (Grad x → 1, u → 2). λ_u = (u + x²/2) + (−ux/2 − x³/12) + (u²/12 + ux²/6 + x⁴/48) + … | **bewiesen** (B1) | Lie-Kontrolle exakt bis Grad 14 (`s2_generator.py`) | im Geist klassisch (Écalle); für diese Reihe neu notiert |
| S2 | Die Fixpunkte von f_u sind x± = u − 1 − W_{0,−1}(−e^{u−1}) = γ ± iβ, mit x₊/v = 2i + 2v/3 − 2iv²/9 − 8v³/135 + iv⁴/135 − … (u = 2v²). Es gilt λ_u′(x±) = x±, also hat 1/λ_u dort das Residuum 1/x±. | **bewiesen** (B2) | Residuum exakt bis v¹¹; Lambert-Form auf 10⁻⁴² (`s6_lambert.py`) | Residuum = 1/log(Multiplikator) ist klassisch; die Lambert-Form hier eingesetzt |
| S3 | **Bauform aller Schichten:** H = γ + β tan Θ mit Θ = ω(y + C − R_u(H)) + (γ/β)·ln cos Θ. Dabei ist ω = \|x₊\|²/(2 Im x₊), π/ω = Re(2πi/x₊), und R_u ist eine Potenzreihe mit polynomialen Koeffizienten. | **bewiesen** (B3) | Folgerungen S4–S6 exakt geprüft | vermutlich neu in dieser Form (Störungs-Fatou-Koordinate) |
| S4 | Jede Schicht h_k (H = Σ v^{k+1}h_k(θ), θ = vy) ist ein Polynom in t = tan θ, L = ln cos θ und θ, mit deg_L h_k ≤ k. | **bewiesen** (B3); Ansatz gelöst für k ≤ 4 | `s1_schichten.py`: alle Koeffizienten gegen die exakten g_n geprüft (True) | neu |
| S5 | **Höchste Logarithmen:** [L^k] h_k = 2 tan^{(k)}(θ)/(3^k k!), erzeugend: Σ_k v^k [L^k]h_k = 2 tan(θ + vL/3). Das sind Hoffmans Ableitungspolynome und damit Tangens- und Sekanszahlen. | **bewiesen** (B3, Ordnungszählung) | k ≤ 4 exakt | neu |
| S6 | **Säkularterm:** ω/v = 1 − v⁴/270 + v⁶/42525 + 137v⁸/1020600 − 2041v¹⁰/2652284250. Das erste nackte θ steht in h₄: −θ·sec²θ/135. Es passt genau zu 2·(−1/270). | **bewiesen** (S3) und exakt bestätigt | h₄ aus dem Ansatz, ω aus der Fixpunktreihe; ω numerisch gegen Lambert-W auf 10⁻¹⁹ (u = 10⁻³) | neu |
| S7 | **Dritte Schicht:** h₂ = t[(2/9)L²(1+t²) − 2L(1+t²) + 41t²/9 + 13/3]. Daraus folgt [y^{2n−3}] g_n(y−1) = 13/12, 7/8, 301/720, 9577/60480, … | **bewiesen** (Ansatz eindeutig, da t, L, θ algebraisch unabhängig) | exakt gegen g_n | neu |
| S8 | h₃ und h₄ explizit (siehe Details) | bewiesen | exakt | neu |
| S9 | **Bernoulli-Formel:** μ := ∂_uλ_u\|_{u=0} = λ₀·Σ_{k≥0} B_k⁺ (λ₀ d/dx)^k(1/λ₀)/k!, μ = 1 − x/2 + x²/6 − x³/18 + 5x⁴/288 − x⁵/240 + x⁶/2880 + 61x⁷/302400 + … | **bewiesen** (B4) | exakt bis x⁹ | vermutlich neu (Variationsformel klassisch) |
| S10 | λ_u(0) = G′(0) = u + u²/12 + u³/120 + u⁴/5040 − u⁵/10080 + u⁶/332640 − 323u⁷/86486400 | exakt | `s2_generator.py` | Nenner nicht in der OEIS |
| T1 | **sexp_e: D(z) = sexp(z) − log(z+2) ist holomorph auf ℂ∖(−∞,−3].** Also c_n = [zⁿ]sexp = (−1)^{n+1}/(n2ⁿ) + d_n mit limsup\|d_n\|^{1/n} = 1/3. | **bewiesen** (B5), gestützt auf die bekannte Holomorphie auf ℂ∖(−∞,−2] | d_n aus zwei Einstellungen (dps 130/M 256 gegen 100/128): Abweichung ≤ 6·10⁻⁵⁶ relativ (n ≤ 50) | Leitterm im Geist bekannt (kn.satz.singulaer_minus2); Holomorphie von D vermutlich neu notiert |
| T2 | **Asymptotik des Rests:** d_n ~ (−1/(c+3))ⁿ·(1/n)·Σ_{k≥1} (1/Γ)^{(k)}(0)/k!·G^{(k)}(ln n) mit G(s) = log(s − K′_c), K′_c = ln((c+3)A), A = sexp′(0). Führend gilt d_n ≈ (−1/3)ⁿ/(n(ln n − 1.1864 + γ)). | **hergeleitet** (Transfersatz) | Zentrum −0.4: n·d_n·(−2.6)ⁿ = 0.22564 gegen Vorhersage 0.22415 (n = 80), 0.2678 gegen 0.2376 (n = 60), Abstand fällt | neu |
| T3 | **Vorlauf bei 0:** Für n ≤ 100 ist d_n bei 0 von einem positiven Beitrag mit Quotient d_{n+1}/d_n ≈ +0.34 beherrscht (Sattelpunkt auf der positiven Achse, sexp(3) ≈ 3.8·10⁶). Die Asymptotik T2 greift erst für n ≳ 10³ (Abschätzung). | numerisch + Erklärung | Quotienten an 5 Zentren konsistent mit einem Beitrag aus z ≈ +2.9 | neu (Warnung für Extrapolationen) |
| T4 | **Vorzeichen:** c_n > 0 für 1 ≤ n ≤ 21, danach sign c_n = (−1)^{n+1} (geprüft bis n = 100). Für alle n ≥ 21 folgt das aus T1, wenn max_{\|z\|=2.5}\|D\| ≲ 200 gilt. | n ≤ 100 numerisch (100+ Stellen); für alle n **hergeleitet mit Lücke** (Lücke: die Schranke für \|D\| auf \|z\| = 2.5 ist nicht zertifiziert) | `s3_taylor.py` | neu |
| T5 | Ganzzahligkeit nach einer Normierung: Für c_n·n·2ⁿ und Ähnliches gibt es keinen Hinweis. Der Rest d_n ist transzendent aussehend. | negativ | – | – |
| P1 | **Iterierte von x·eˣ = Ketten punktierter Partitionen:** n![xⁿ](xeˣ)^{∘t} = Σ_j B′(n,j)·C(t,j). B′(n,j) ist die Zahl der Ketten 0̂ < σ₁ < … < σ_j mit σ_j maximal (ein Block), in der Ordnung von Chapoton–Vallette. | **bewiesen** (B6); durch Abzählen geprüft für n ≤ 6 | `s4_punktiert.py` | Identifikation vermutlich neu |
| P2 | μ(0̂, M) = (−1)^{n−1} n^{n−2} für jedes der n Maxima. Summe (−n)^{n−1} = n![xⁿ] W(x) (Lambert). | bewiesen (Hall + Lagrange); durch Abzählen geprüft für n ≤ 6 | `s4_punktiert.py` | passt zu dim PreLie(n) = n^{n−1}; wahrscheinlich bei Chapoton–Vallette enthalten (Volltext nicht geprüft) |
| P3 | **Neue Folgen:** Dreieck B′ (Zeile 6: 6, 3210, 41940, 117360, 86400), Zeilensummen 2, 15, 232, 6155, 248916, zweite Schicht (2n−2)!·2^{n−1}·[y^{2n−2}]g_n(y−1) = −1, −3, −26, −468, −14256, −657568, −42752256, …, dritte Schicht ×6: 13, 126, 2408, 76616, 3656688, … | ganzzahlig: Dreieck bewiesen; zweite Schicht ganzzahlig beobachtet für n ≤ 12 | OEIS-Suche 18.10.2026: keine Treffer | neu |
| Q1 | **Pentation an der parabolischen Stelle:** Die Abbildung y ↦ G(y;u) hat bei 0 den Multiplikator λ_u(0) = u + … Sie ist also hyperbolisch (anziehend), und die nächste Stufe ist eine Koenigs-Linearisierung ohne Baumstruktur. „Bäume aus Bäumen“ ergeben nur Ketten von Ketten: Z(st) = Z(s)∘Z(t). | bewiesen (elementar) | – | Negativergebnis |

## Beweise und Details

### B1 (Fluss)
Sei φ_y = exp(yλ_u∂) der formale Fluss. Die Grade von λ_u sind ≥ 2, und H(0) = log(1−u) = O(u). Deshalb tragen zu jeder u-Ordnung nur endlich viele Lie-Glieder bei. K(y) := φ_y(H(0)) hat also in jeder u-Ordnung ein Polynom in y als Koeffizienten. Weiter gilt K(y+1) = φ₁(K(y)) = f_u(K(y)) und K(0) = H(0). Damit ist K(k) = H(k) für alle k ∈ ℕ, und Polynome, die auf ℕ übereinstimmen, sind gleich. ∎

### B2 (Residuen)
x± sind die Fixpunkte von f_u: e^x − 1 − x = −u. Mit y = x + 1 − u wird daraus y e^{−y} = e^{u−1}, also x = u − 1 − W(−e^{u−1}). An einem Fixpunkt p des Zeit-1-Flusses gilt λ(p) = 0 und f′(p) = e^{λ′(p)}. Wegen f_u′ = eˣ folgt λ′(x±) = x±. Formal in v gilt das über die Hensel-Faktorisierung von Λ_v(ξ) = λ_u(vξ)/v² = (ξ²+4)/2 + O(v). ∎

### B3 (Bauform)
Λ_v ist in Q[ξ][[v]] nach Hensel gleich (ξ − ξ₊)(ξ − ξ₋)·Q_v(ξ), mit Q₀ = ½. Also ist 1/Q_v polynomial, und die Division durch das normierte Quadrat ergibt 1/λ_u = Σ±1/(x±(x − x±)) + ρ_u mit polynomialem ρ_u. Integration ergibt y + C = (2β/r²)·arctan((H−γ)/β) + (γ/r²)·ln((H−γ)²+β²) + R_u(H) mit r² = β²+γ². Setzt man H = γ + β tan Θ, folgt S3. Für die Schichten gilt γ/β = v/3 + v³/135 − 13v⁵/5670 + …, β/v = 2 − 2v²/9 + v⁴/135 + … und ω/v aus S6. Jedes Auftreten von L = ln cos Θ kostet mindestens einen Faktor v (über γ/β). Entwickelt man ln cos Θ um θ, kostet jeder Schritt ein v ohne neues L. Deshalb ist der L^k-Anteil von h_k gleich [v^k] 2 tan(θ + vL/3), mit eingefrorenem L. Das nackte θ entsteht nur aus ω/v − 1 = O(v⁴). ∎

Explizit (`daten/s1_schichten.json`):
- h₁ = (2/3)(L sec²θ − 4 sec²θ + 1). Das ist K2 aus Etappe 1.
- h₃ = (2/81)L³(1+t²)(1+3t²) − L²(29t⁴/27 + 38t²/27 + 1/3) + L(47t⁴/9 + 902t²/135 + 197/135) − 1147t⁴/135 − 1409t²/135 − 2
- h₄ = [60L⁴t⁵ + 100L⁴t³ + 40L⁴t − 1220L³t⁵ − 2000L³t³ − 780L³t + 9330L²t⁵ + 15024L²t³ + 5694L²t − 31758Lt⁵ − 49956Lt³ − 18198Lt + 40539t⁵ + 61761t³ + 21240t − 18θ(1+t²)]/2430

### B4 (Bernoulli)
Wir variieren f_ε = exp(λ_ε∂)x. Es gilt δf = λ(f(x))·∫₀¹ ψ(φ_s x) ds mit ψ = δλ/λ. Mit δf = 1 und der Fatou-Zeit τ (D = d/dτ = λ d/dx) wird daraus A(τ+1) = ((e^D − 1)/D)Ψ(τ) mit A = 1/λ∘φ_τ. Also ist Ψ = D e^D/(e^D − 1)·A = Σ B_k⁺ D^k A/k!. Die Reihe konvergiert x-adisch, weil jedes D den Grad um 1 hebt. ∎

### B5 (T1)
sexp ist holomorph auf ℂ∖(−∞,−2] und hat dort keine Nullstelle ausser −1: Aus sexp(z+1) = e^{sexp z} folgt, dass ausserhalb des Schnitts nur sexp(−1) = 0 möglich ist. Auf (−3,−2) springt sexp(z) = log sexp(z+1) um 2πi, weil sexp(z+1) ∈ (−∞,0) und sexp′ > 0 ist. log(z+2) springt genau gleich. D ist also stetig über (−3,−2) und damit holomorph. Bei −2 gilt D = log(sexp(z+1)/(z+2)) → log A, und das ist holomorph. Bei −3 gilt lokal D = log(log(1/(Aw)) − log φ(w)) − log(1−w) mit w = z+3 und φ(0) = 1 holomorph. Mit ζ = −z/3 ist das log(ℓ − K′) + O(w/ℓ), ℓ = log 1/(1−ζ). Die Δ-Analytizität ist gegeben (einziger Schnitt (−∞,−3]). Der Transfersatz (Flajolet–Sedgewick VI.2, β → 0) ergibt T2. ∎

Zu T3: Die Fehlerschranke im Transfer enthält max_{\|z\|=3+η}\|D\|, und das ist ≈ sexp(3+η) ≥ 10⁶. Deshalb wird die führende Ordnung erst für n ~ 10³ sichtbar. Vom Zentrum −0.4 aus reicht der Vergleichskreis nur bis 2.3, wo sexp(2.3) ≈ 50 ist; dort ist T2 schon bei n ≈ 60–80 deutlich zu sehen.

### B6 (P1)
xE (punktierte Menge) ∘ F heisst: Man wählt einen Block aus F-Strukturen und darin einen Punkt. (xE)^{∘t} sind deshalb die Multiketten 0̂ ≤ π₁ ≤ … ≤ π_t, bei denen π_t ein einziger Block ist. Die Punkte erben sich nach oben, und das ist genau die Ordnung von Chapoton–Vallette. Die verschiedenen Elemente oberhalb von 0̂ bilden eine strikte Kette, die im Maximum endet. Ihre Lage in t Stufen wählt man auf C(t,j) Arten. Bei t = −1 kommt Σ_j(−1)^j B′(n,j) = Σ_M μ(0̂,M) (Hall) heraus. Das ist gleich n![xⁿ](xeˣ)^{∘−1} = n![xⁿ]W(x) = (−n)^{n−1}, und aus Symmetrie folgt P2. Die Diagonale ist Π_{k=2}^{n} 2·C(k,2) = n!(n−1)!. ∎

## Reproduktion

| Skript | Inhalt | Dauer |
|---|---|---|
| `s1_schichten.py 4` | Schichten h₀…h₄ (Ansatz), Gegenprobe gegen exakte g_n | ~2 min |
| `s2_generator.py 14` | λ_u, Fixpunktreihe, Residuen, ω, γ/β, Bernoulli-Formel | ~1 min |
| `s3_taylor.py`, `s3b_zentren.py` | T1–T4 aus den etappe3-Taylor-Daten (nur lesend) | < 1 min |
| `s4_punktiert.py` | Verband der punktierten Partitionen, n ≤ 6 | ~1 min |
| `s5_folgen.py`, `s6_lambert.py` | Folgen für die OEIS, Lambert-Kontrolle | < 1 min |

Ausgaben liegen in `daten/*.json`.

## Vorschläge für Grapheinträge

- „Störungsreihe als Fluss des Erzeugers λ_u von eˣ−1+u“ — Satz — bewiesen
- „Residuen 1/x± und Lambert-Fixpunkte x± = u−1−W(−e^{u−1})“ — Satz — bewiesen
- „Bauform der Schichten: H = γ + β tan Θ, Θ = ω(y+C−R) + (γ/β) ln cos Θ“ — Satz — bewiesen
- „Höchste Logarithmen der Störungsreihe = 2 tan(θ + v ln cos θ/3)“ — Satz — bewiesen
- „Dritte Schicht der Störungsreihe geschlossen“ — Satz — bewiesen (exakt geprüft)
- „Säkularterm ab Schicht 5 und Frequenz ω = v(1 − v⁴/270 + …)“ — Satz — bewiesen
- „u-linearer Teil des Erzeugers über Bernoulli-Zahlen“ — Satz — bewiesen
- „sexp(z) − log(z+2) holomorph auf ℂ∖(−∞,−3]; Taylor-Asymptotik mit log log bei −3“ — Satz — bewiesen (Asymptotik hergeleitet)
- „Vorzeichen der Taylor-Koeffizienten von sexp_e bei 0“ — Beobachtung — numerisch bis n = 100, allgemein mit Lücke
- „Vorlauf der Taylor-Asymptotik durch Wachstum auf der positiven Achse“ — Beobachtung — numerisch
- „Iterierte von x eˣ zählen Ketten punktierter Partitionen; Möbius (−1)^{n−1}n^{n−2}“ — Satz — bewiesen
- „Neue OEIS-Kandidaten: Ketten-Dreieck, Zeilensummen, Schichtfolgen“ — Idee — offen

## Offene Fragen

1. Hat R_u (die regulären Teile) eine geschlossene Form? Bei u = 0 ist es der Fatou-Rest von eˣ−1 (K6: −x/36 + x²/540 + …).
2. Bewiesen ist deg_θ h_k ≤ ⌊k/4⌋? Das folgt nur, wenn ω/v − 1 = O(v⁴) gilt (bewiesen) und R_u kein θ erzeugt (ebenfalls); eine scharfe Gradschranke ist nur für k ≤ 4 geprüft.
3. Ist (2n−2)!·2^{n−1}·[y^{2n−2}]g_n(y−1) für alle n ganz? Dafür muss man durch 3 teilen; nachgewiesen ist es nur für n ≤ 12.
4. T4 zertifizieren: eine Intervallschranke für \|sexp\| auf \|z\| = 2.5, über die Funktionalgleichung aus den Taylor-Daten.
5. Die Kaskade weiter: Der Beitrag von −k ist ein (k−1)-fach iterierter Logarithmus, und man kann sein Transfergewicht bestimmen. Ausserdem wäre zu prüfen, ob die Sattelpunktnäherung den Vorlauf T3 quantitativ beschreibt.
6. Den Volltext von Chapoton–Vallette (arXiv math/0410051) auf die Möbius-Formel prüfen.
