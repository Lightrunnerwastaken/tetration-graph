# Strang E — Modellprobleme: das Abschneide-Verfahren an einfacheren Gleichungen f∘f = g

Stand: 10.10.2026, 22:53 (Runde 3 abgeschlossen). Ordner `e_modelle/`. Python: Miniconda 3.13, python-flint (arb), mpmath.
Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.

## 0. Ziel und Werkzeug

Dasselbe Verfahren wie bei exp — p_N ∈ P_N (Grad ≤ N), [x^k] p∘p = [x^k] g für k = 0..N, Zweig stetig von N = 1 aus
(Newton von p_{N−1} mit a_N = 0; bei Versagen Homotopie F(a) = (1−t)F(a_start)) — an Modellen untersuchen, in denen die
ausgezeichnete Lösung explizit ist. Kern: `kern.py` (`zweig(gk, nmax, bits, ordnung=2|3)`, auch f∘f∘f = g),
geprüft gegen `gemeinsam/polynome_bis_194.json` (Abweichung 5e-61 für N ≤ 40, `test_kern.py`).

Übersicht der Modelle und Ergebnisse (Details unten):

| Modell | ausgezeichnete Lösung | Ergebnis | Status |
|---|---|---|---|
| P: Fixpunkt im Entwicklungspunkt (g(0)=0) | reguläre (Koenigs-)Iterierte | p_N = T_N f_reg exakt; parabolisch (eˣ−1): Divergenz | bewiesen |
| T: lineares Translationsmodell u(x+s)+u(x)=h | Lösung kleinsten Exponentialtyps | konvergiert ⇔ Typ(h) < π/s, Rate Typ·s/π | bewiesen + numerisch |
| K: anziehender Fixpunkt, ‖f*‖_{W_r} < r | reguläre Iterierte am anziehenden Fixpunkt | p_N → f* geometrisch, J_N⁻¹ gleichmässig beschränkt | bewiesen (N ≥ N₀) + numerisch |
| K′: g = (3−x)/(8−3x), f = 1/(3−x) | Möbius-Wurzel = reguläre Iterierte | ganzer Zweig ab N = 1, ‖p_N − T_N f‖ ≤ 7.7·6^{−N−1} | **vollständig bewiesen** (arb) |
| b^x, b < e^{1/e} (1.2, √2) | reguläre Iterierte am unteren Fixpunkt | konvergiert, √2: Rate 0.40 | numerisch |
| Möbius-Drehung (Fixpunkte A±iB) | Möbius-Wurzel (einzige bei A±iB analytische) | konvergiert ⇔ Pol von f_M ausserhalb D(a₀, \|A+iB−a₀\|) (⇔ f³(0) ≠ ∞ zuerst; A=0: c < π/3); Rate \|A+iB−a₀\|/\|Pol−a₀\|; darüber **Versagen** (keine Lösung nahe f_M) | numerisch (präzise), Erklärung Strang C |
| hyperbolische Möbius tanh(artanh x + c) | Möbius-Wurzel | Lösung nahe f_M existiert, Rate T = tanh(c/2) (Miller); Zweig ab N = 1 springt bei c = 2 auf falsche Familie | numerisch |
| b^x, b = 1.5, 2, 10, e; f∘f∘f = exp | Kneser | konvergiert, **polynomial** ~ N^{1−2s_b}, s_b = Re(2πi/log(L_b ln b)) | numerisch |

## 1. Ersatzformen (bewiesen)

Sei p ∈ P_N mit p'(0) ≠ 0, a₀ = p(0). Dann sind äquivalent:
(i) p∘p ≡ g mod x^{N+1};
(A) p = T_N^{(a₀)}[g∘p⁻¹] (Taylor-Polynom vom Grad N an der Stelle a₀; p⁻¹ = lokale Umkehrung bei a₀, p⁻¹(a₀) = 0);
(B) p = T_N^{(0)}[p⁻¹∘g] (p⁻¹ = lokale Umkehrung bei g(0) = p(a₀), Werte nahe a₀).
*Beweis.* (i)⇒(B): Komposition mit der analytischen Funktion p⁻¹ erhält Kongruenzen mod x^{N+1}, also p = p⁻¹∘(p∘p) ≡ p⁻¹∘g;
da deg p ≤ N, ist p gleich diesem Taylor-Polynom. (B)⇒(i): p ≡ p⁻¹∘g ⇒ p∘p ≡ p∘p⁻¹∘g = g. (A) analog mit y = p(x) nahe a₀. ∎
Linearisierung an einer exakten Lösung f (L δ = δ∘f + (f'∘f)δ) in Form (B): K δ = δ + (δ∘f)/(f'∘f) =: (I + M⁻¹C)δ,
M = Multiplikation mit m = f'∘f, C = Komposition mit f. In Abel-Koordinaten (α(f) = α + ½): mit δ = V(α)/α'(f) gilt
L δ = (V(α + ½) + V(α))/α'(g) (bewiesen, direkte Rechnung) — der Kern sind genau die ½-antiperiodischen V.

## 2. Modell P: Fixpunkt im Entwicklungspunkt (bewiesen)

Ist g(0) = 0, g'(0) = λ > 0, und wählt man bei N = 1 den Zweig a₁ = √λ, a₀ = 0, so ist für alle N
**p_N = T_N h**, h = formale reguläre halbe Iterierte (h(0) = 0, h'(0) = √λ).
*Beweis.* Bei a₀ = 0 hängt [x^k] p∘p nur von a₁..a_k ab, und a_k tritt mit dem Faktor a₁^k + a₁ = λ^{k/2} + √λ ≠ 0 auf; das
dreieckige System hat die Lösung T_N h. Jacobi-Matrix dort untere Dreiecksmatrix mit Diagonale ≠ 0. Der Startpunkt
(p_{N−1}, 0) erfüllt alle Gleichungen ausser der letzten; der Newton-Schritt (J dreieckig) ändert nur a_N und landet exakt
auf T_N h. ∎
Folgerungen: (a) |λ| ≠ 1 (hyperbolisch): h konvergiert (Koenigs) ⇒ p_N → h auf der Konvergenzkreisscheibe.
(b) **Gegenbeispiel** g = eˣ − 1 (parabolisch): h divergiert (Baker 1958, Écalle: Gevrey-1), also divergiert p_N(x) für jedes
x ≠ 0. Numerisch (`m6_parabolisch.py`): max|a_k^(N) − h_k| = 2e-47, a₀ = 0 für alle N; h_40 = 12.4, h_50 = −3.3e6,
h_60 = 5.4e12 (≈ n!/(4π)^n); p_N(1) = 1.271, 1.271, 1.271, 8.5, −2.8e6, 5.2e12 (N = 10…60). Eine reelle analytische
Lösung nahe 0 gibt es hier nicht (f(0) = a ≠ 0 erzwingt einen weiteren Fixpunkt a von g, und der ist nicht reell).

## 3. Modell T: lineares Translationsmodell (bewiesen + numerisch)

Linearisierung der halben Iterierten von g = x + 2s + εh an f = x + s: gesucht u_N ∈ P_N mit T_N[u_N(x+s) + u_N(x)] = T_N h.
**Satz T.** (1) u_N = φ_s(D) T_N h exakt, φ_s(z) = 1/(1 + e^{sz}), D = d/dx (nilpotent auf P_N).
(2) Hat h Exponentialtyp τ < π/s, so gilt u_N → u := Σ_n c_n h^{(n)} (c_n = [z^n]φ_s) lokal gleichmässig, geometrisch mit
jeder Rate > τs/π; u ist die einzige Lösung von u(x+s) + u(x) = h vom Exponentialtyp < π/s.
(3) Für h = e^{λx}, |λ| > π/s, divergiert u_N(0).
*Beweis.* (1) Die Verschiebung erhält P_N, also ist T_N überflüssig; I + e^{sD} = 2I + nilpotent ist invertierbar mit Inverser
φ_s(D) (endliche Summe). (2) D^n T_N h = T_{N−n} h^{(n)}, also u_N − u = Σ_{n≤N} c_n (T_{N−n}h^{(n)} − h^{(n)}) − Σ_{n>N} c_n h^{(n)}.
φ_s hat Pole bei iπ(2m+1)/s, also |c_n| ≤ C_ε(s/π + ε)^n; mit |h_j| ≤ C τ'^j/j! ist |h^{(n)}(x)| ≤ Cτ'^n e^{τ'R} und
|T_{N−n}h^{(n)} − h^{(n)}| ≤ Cτ'^n Σ_{k>N−n}(τ'R)^k/k! auf |x| ≤ R. Aufteilen n ≤ N/2 (Fakultät) und n > N/2 (geometrisch)
gibt Konvergenz. Eindeutigkeit: die Differenz zweier Lösungen ist 2s-periodisch und s-antiperiodisch, also eine
Fourier-Reihe in e^{iπ(2m+1)x/s} — vom Typ ≥ π/s, falls ≠ 0. (3) u_N(0) = Σ_{n≤N} c_n λ^n, Glieder ~ (|λ|s/π)^n. ∎
Numerisch (`m7_translation_linear.py`, s = ½, π/s = 2π): |u_N(0) − u(0)| bei N = 160: 4e-52 (λ = 3), 2e-4 (λ = 6),
0.038 (λ = 6.2), **6.1 (λ = 6.4)**, 2e16 (λ = 8); λ = 5i: 2e-16, λ = 7i: 9e7. Schwelle exakt |λ| = 2π.
(4) *Gleicher Typ, verschiedene Stärke (bewiesen, aus der exakten Formel):* u_N(0) = Σ_{n≤N} c_n n! h_n exakt. Haben die
Kernfunktion (φ_s hat einfache Pole bei ±2πi für s = ½, also c_n ~ (2π)^{−n}·Schwingung) und h denselben Typ, aber h die schwächere Borel-Singularität
(n! h_n ~ (2πi)^n n^{−σ−1}, „Lücke“ σ + 1 in den Koeffizienten-Exponenten), so ist u_N(0) − u(0) = −Σ_{n>N} c_n n! h_n ~ N^{−σ}
= **N^{1 − Lücke}**. Das ist genau die Form des exp-Gesetzes c_N ~ N^{1−2s} (Lücke 2s zwischen Kneser k^{−s−2} und Konkurrenz
k^{s−2}; Strang A/B) samt dem „Zusatzfaktor N“: er entsteht hier sichtbar als Schwanzsumme über n > N.
**Mechanismus in Reinform:** Das Abschneiden wählt unter allen Lösungen diejenige, deren Wachstum (hier: Exponentialtyp, d. h.
Abfall der Taylor-Koeffizienten) kleiner ist als das der Kernfunktionen e^{iπ(2m+1)x/s}; die Rate ist das Verhältnis der
Wachstumsgrössen. Ist die Lösung „grösser“ als der Kern, versagt das Verfahren.

## 4. Satz K: der kontraktive Fall (bewiesen) und ein explizites Modell (numerisch)

W_r := {δ = Σ δ_k x^k : ‖δ‖_r = Σ|δ_k| r^k < ∞} (Banach-Algebra, T_N hat Norm 1).
**Satz K.** Sei f* eine Lösung von f*∘f* = g mit f* ∈ W_{r'} für ein r' > r, ‖f*‖_r =: r₁ < r, und m := f*'∘f* ohne Nullstelle
auf |x| ≤ r. Dann gibt es N₀ und ρ₀ > 0, so dass für N ≥ N₀ das abgeschnittene System genau eine Lösung p_N in der Kugel
‖p − T_N f*‖_r < ρ₀ hat, und ‖p_N − T_N f*‖_r ≤ C θ^N (θ < 1). Insbesondere p_N → f* gleichmässig auf |x| ≤ r. Liegt der
von N = 1 fortgesetzte Zweig für ein N₁ ≥ N₀ in dieser Kugel, so bleibt er für alle N ≥ N₁ darin.
*Beweis.* (a) *Fixpunkt.* f*(D̄_r) ⊂ D_{r₁}, also hat f* genau einen Fixpunkt x* ∈ D_r, |μ| := |f*'(x*)| < 1 (Schwarz–Pick);
μ ≠ 0, da m(x*) = μ.
(b) *Grenzoperator.* A δ = δ∘f* + m δ = M(I + K)δ, K = M⁻¹C. C ist kompakt auf W_r: ‖f*^k‖_r ≤ r₁^k, also
‖C − C_k₀‖ ≤ (r₁/r)^{k₀+1} mit endlichem Rang C_k₀. 1/m ∈ W_r (Wiener-Lemma für W_r: m ∈ W_r ohne Nullstelle auf D̄_r).
Kern von I + K: δ∘f* = −mδ. Entwicklung bei x*: δ = d_j(x − x*)^j + …, Vergleich der Ordnung j gibt d_j(μ^j + μ) = 0,
μ^{j−1} ≠ −1 wegen 0 < |μ| < 1 ⇒ δ ≡ 0. Fredholm ⇒ I + K invertierbar.
(c) *Abschneiden.* M ist untere Dreiecks-Toeplitz-Matrix, daher T_N M T_N = T_N M und (T_N M|_{P_N})⁻¹ = T_N M⁻¹|_{P_N}.
Somit A_N := T_N A|_{P_N} = M_N (I + T_N K)|_{P_N}. Weil K kompakt und T_N → I punktweise, gilt ‖(I − T_N)K‖ → 0, also
ist B_N := I + T_N K für N ≥ N₀ auf ganz W_r invertierbar mit ‖B_N⁻¹‖ ≤ 2‖(I+K)⁻¹‖; B_N bildet P_N in sich ab und ist dort
injektiv, also bijektiv mit derselben Schranke. Damit ‖A_N⁻¹‖ ≤ 2‖(I+K)⁻¹‖‖M⁻¹‖ =: β/2 gleichmässig.
(d) *Konsistenz.* q_N := T_N f*. F_N(q_N) = T_N[(q_N − f*)∘q_N] + T_N[f*∘q_N − f*∘f*]; der zweite Term ist 0 (Kongruenz
mod x^{N+1}), der erste hat Norm ≤ Σ_{k>N}|f*_k| r₁^k ≤ (r₁/r')^{N+1}‖f*‖_{r'}.
(e) *Lipschitz.* Auf der Kugel ‖p‖_r ≤ r₂ < r gilt ‖F_N'(p) − F_N'(q)‖ ≤ γ‖p − q‖ mit γ unabhängig von N (alle
Ableitungen werden durch die Komposition mit einem p der Norm ≤ r₂ < r gedämpft: sup_k k r₂^{k−1}/r^k < ∞ usw.).
‖F_N'(q_N) − A_N‖ ≤ sup_k(k r₁^{k−1}/r^k)‖q_N − f*‖_r + Σ_{k>N} k|f*_k| r₁^{k−1} + ‖f*''‖‖q_N − f*‖_r → 0
(gleiche Abschätzungen wie in Satz K′, Teil 2), also ‖F_N'(q_N)⁻¹‖ ≤ β für grosse N.
(f) Kantorovich: h_N = βγ‖F_N'(q_N)⁻¹F_N(q_N)‖ → 0 geometrisch ⇒ Behauptung; der aufgefüllte Vorgänger p_{N−1} liegt für
grosse N in der Eindeutigkeitskugel, also führt Newton von dort auf p_N. ∎
**Modell K** (`m5_kontraktiv.py`): f(x) = 1/(3 − x), g = f∘f = (3 − x)/(8 − 3x), Fixpunkte 0.382 (anziehend, μ = 0.146) und
2.618. ‖f‖_r = 1/(3 − r) < r für 0.382 < r < 2.618 ⇒ Satz K greift. Numerisch: |J_N(T_N f)⁻¹| in W_1 = 9.3728, 9.3729,
9.3729, 9.3729 (N = 10, 20, 40, 60), in W_2 = 11.76 … 11.91; |a₀^(N) − f(0)| = 2.3e-11, 1.6e-28, 3.1e-55 (N = 10, 30, 60),
Rate 0.131 ≈ x*/R_f = 0.382/3 = 0.127; p_N(1) − f(1) ~ 3^{−N} (Taylor-Schwanz von f bei 1).

**Satz K′ (vollständig bewiesen, rechnergestützt mit arb; `satz_k_beweis.py`).** Für g(x) = (3 − x)/(8 − 3x) existiert der von
p_1 = 1/3 + x/8 aus fortgesetzte Zweig für **alle** N, und ‖p_N − T_N f‖_1 ≤ r_−(N) ≤ 7.71·6^{−N−1} (N ≥ 16; für kleinere N
explizit, z. B. 0.0179, 0.0038, 7.7e-4 für N = 1, 2, 3); also p_N → f = 1/(3 − x) gleichmässig auf |x| ≤ 1.
Beweisgang (alle Zahlen sind obere Schranken aus Intervallrechnung):
1. Grenzoperator, N ≥ n0 = 16: ‖(3−f)²‖_1 ≤ 7.97222 (= ‖M⁻¹‖), Blockzerlegung niedrig (k ≤ 16, 17×17-Matrix I + PKP exakt
   invertiert) / hoch (‖K x^k‖ ≤ 2^{−k}·7.97, ε ≤ 6.08e-5): ‖(I + XP)⁻¹‖ ≤ 1.20755 ⇒ ‖A_N⁻¹‖ ≤ 9.62755 für alle N ≥ 16
   (gemessen 9.3729 — die Schranke ist scharf).
2. Auf U = {‖p − f‖_1 ≤ 0.1056}: γ ≤ 3.0367 (Lipschitz von F_N'), ‖F_N'‖ ≤ 1.3023, N-unabhängig (Komposition mit p, ‖p‖ ≤ 0.6056 < 1).
3. Kantorovich um q_N = T_N f: N = 1…17 direkt (h ≤ 0.267, 0.090, 0.021, …), N ≥ 16 über 1. (h ≤ 6.7e-12, monoton in N):
   genau eine Lösung p_N* in B(q_N, min(r_+, 0.05)).
4. Zweig: p_1 = 1/3 + x/8 liegt in der Eindeutigkeitskugel (‖p_1 − q_1‖ = 1/72); für N < 16 Intervall-Newton vom aufgefüllten
   Vorgänger (1 Schritt für N = 1, 2; danach direkt) bis Kantorovich greift, Grenzwert liegt in der Eindeutigkeitskugel um q_{N+1};
   für N ≥ 16 Störungsargument mit monotonen Schranken (h₀ ≤ 9.5e-7). ∎
Damit ist ein nichtlineares Modell **vollständig** bewiesen: das Abschneiden trifft die ausgezeichnete (reguläre) Wurzel.
**Grenzen:** Satz K braucht einen anziehenden Fixpunkt in einer Kreisscheibe um 0, die f* in sich abbildet. Bei exp gibt es
keinen (‖f_K‖_{|L|} = 2.24 > |L|; f_K bildet keine Scheibe um 0 in sich ab); die symmetrische hyperbolische Möbius-Abbildung
erfüllt die Bedingung ebenfalls nicht ((r − T)² < T² − 1 unmöglich) — und dort wird das Verfahren für grosse c unregelmässig.

## 5. Möbius-Drehung: ein exp-ähnliches, explizites Modell mit scharfer Schwelle (numerisch, präzise)

g(x) = A + B·tan(arctan((x−A)/B) + c): Drehung der projektiven Geraden um die Fixpunkte A ± iB (|g'| = 1 dort), keine reellen
Fixpunkte, wie exp. Möbius-Wurzel f_M = Drehung um c/2: einzige reelle Wurzel, die bei A ± iB analytisch ist (alle anderen
α⁻¹(α + ½ + θ(α)), α = arctan((x−A)/B)/c, sind dort singulär mit Stärke π/c). f_M hat einen reellen Pol.
Symmetrisch A = 0, B = 1: τ = tan(c/2) = f_M(0), Pol bei 1/τ.

**Befund 1 — Konvergenz und Rate.** p_N → f_M geometrisch; die Fehler (a₀, p_N(x)) sind ein komplexes Modenpaar
Re(C N^β z^N). Gleitende Prony-Anpassung, N bis 90 (`m1e_mode.py`): |z| = 0.208→0.219 (c = 0.4), 0.355→0.371 (0.6),
0.557→0.572 (0.8), arg z = 96.1°, 99.5°, 103.1°. Erster Newton-Schritt η_N = |J⁻¹F(T_N f_M)| (`m1f_eta.py`) hat dieselbe Rate.
**Gesetz (numerisch, 9 Werte c ∈ [0.4, 1.4] symmetrisch + 5 Werte mit Fixpunkten L, L̄):**
  ρ = χ(0, f(0)) / χ(g(0), ∞) = χ(0, f(0)) / χ(f(0), Pol von f),
χ = chordaler Abstand bezüglich der Drehung (in x' = (x−A)/B: χ(u,v) = |u−v|/√((1+u²)(1+v²))); symmetrisch ρ = sin(c/2)/cos c.
Vergleich lokale Rate (mit Vorfaktor N^β, β ≈ −1.5 für η) gegen Gesetz: c = 1.03: 0.947/0.957; 1.04: 0.970/0.982;
1.055: 1.006/1.020; 1.065: 1.036/1.049; 1.1: 1.10/1.15; 1.2: 1.48/1.56; 1.4: 3.58/3.79; mit Fixpunkten L, L̄:
c = 0.705: 0.393/0.388; 1.0: 0.666/0.666; 1.15: 0.879/0.894; 1.2029: 0.979/1.000; 1.25: 1.080/1.111.
*Nachtrag (verschobene Fixpunkte, `m2b_eta.py`):* A = −1, B = 1, c = 0.6: Pol von f_M (2.233) liegt in D(a₀, R₀)
(a₀ = 0.896, R₀ = 2.143) ⇒ vorhergesagt Versagen mit Rate 1.603; gemessen η = 66, 7.6e4, 3.3e6, 2.1e8 (N = 20, 40, 50, 60),
lokal 1.515 ✓ — **Versagen schon bei kleinem Drehwinkel**, entscheidend ist allein die Lage des Pols relativ zu W.
A = −1, c = 0.2: Pol weit weg, Rate 0.181 vorhergesagt, gemessen 0.136–0.146 (N^β-Vorfaktor, β ≈ −π/c + 4).
A = +1, B = 1, c = 0.6: hier gewinnt der Miller-Term R₁/|ζ| = √2/4.233 = 0.334 gegen den Pol-Term 0.301; gemessen lokal
0.333, 0.330, 0.338 (N = 60→70→80→90) ✓. Das χ-Gesetz oben ist also der Pol-Term; vollständig gilt das Max-Gesetz (§6 Ende).
**Befund 2 — Schwelle.** ρ = 1 ⇔ Pol von g = f(0) ⇔ **f(f(f(0))) = ∞** (symmetrisch: c = π/3). Darunter konvergiert p_N
(bei c = π/3 selbst nur noch polynomial, η ~ N^{−1.5}), darüber gibt es **keine** Lösung nahe T_N f_M: η_N wächst geometrisch,
Newton von T_N f_M aus scheitert für alle N, und der von N = 1 fortgesetzte Zweig läuft weg (c = 1.1, 1.15, 1.2: a₁ → 0.04,
a_N ~ 2.4^N; bei c = 1.05 noch langsame Konvergenz, Fehler 8e-4 bei N = 55). In Abel-Koordinaten: mit W = Abel-Abstand von
0 nach ∞ ist ρ = sin(c/2)/sin(c(W − 1)); Schwelle W = 3/2 unabhängig von c; für c → 0 gilt ρ → 1/(2(W − 1)).
*Knapp über der Schwelle* (c = 1.055, ρ = 1.0206; `m8_langlaeufe.py`, `m8b_falte.py`): der Zweig folgt f_M scheinbar konvergent
(a₀-Fehler 2.6e-2, 5.3e-3, 2.0e-3, 1.0e-3, 6.1e-4 bei N = 20, 40, 60, 80, 100), kehrt dann um (5.1e-4, 5.4e-4, 5.6e-4 bei
N = 101–103), und die Fortsetzung wird instabil: erster Newton-Schritt vom aufgefüllten Vorgänger 14.0, 9.2, 6.3 (Koeffizienten-
Maximum), kleinster gewichteter Singulärwert ≈ 3e-4. Vorhergesagte Umkehr (Minimum von N^{−1.5}ρ^N) bei N ≈ 74 — passt.
⇒ **Numerisch gesehen kann das Verfahren über ~100 Grade „konvergent aussehen“ und trotzdem nicht konvergieren.** Eine
Warnung auch für exp: Daten bis N = 232 allein beweisen nichts; die Raten-Gesetze (hier ρ > 1) entscheiden.
**Befund 3 — „Möbius-exp“.** Fixpunkte L, L̄ und g(0) = 1 (c = 0.70510, g'(0) = 1.1925, f_M(0) = 0.4780, Pol bei 3.953):
p_N → f_M mit ρ = 0.388 (Fehler a₀ 7.5e-26 bei N = 50, 2.4e-38 bei N = 80). Weit unter der Schwelle c* = 1.2029.
**Befund 4 — Stabilität ist nicht der entscheidende Begriff** (`m1g_stabilitaet.py`). In l¹ mit gleichen Gewichten 1^k
(Radius der Fixpunkte) ist |J_N(T_N f_M)⁻¹| ≈ 2 für ALLE c (1.83/2.13/2.31 bei c = 0.6; 1.85/2.15/2.31 bei 1.0;
2.05/2.03/2.19 bei 1.2 für N = 20/40/60) — auch oberhalb der Schwelle, wo das Verfahren scheitert. Der Rest |F(T_N f_M)| wächst
schon bei c = 1.0 (1e4 → 5e11), η = |J⁻¹F| fällt dort aber (0.39 → 0.0061). Entscheidend ist also allein die Auslöschung in
J⁻¹F (vgl. Strang B 20:03: derselbe Befund bei exp).
**Befund 5 — ein einziger Kernmodus, wie bei exp** (`m1h_abel_moden.py`, c = 0.8). θ_N(w) := α(p_N(α⁻¹(w))) − w − ½ hat auf
w ∈ [−½, ½] (Höhe v = 0) die Fourier-Beträge |c_{±1}| = 9.4e-10, 2.6e-15, 1.4e-20 (N = 20, 40, 60), |c_{±3}| = 3.6e-14, 1.4e-23,
5.1e-30, gerade Moden (Nicht-Kernrichtungen) 4e-14, 5e-26, 9e-38 (Rauschen/Schwanz). Bei Höhe v = 0.1 wächst |c_{−1}| um
e^{2π·0.1} = 1.87 und |c_1| fällt um denselben Faktor — also exakt p_N ≈ f_M + c_N δ_{−1} + c̄_N δ_1 mit
δ_{∓1} = e^{∓2πiα}/α'∘f. Genau das haben A/B bei exp gefunden; der Unterschied liegt NUR im Gesetz für c_N
(Möbius: geometrisch ρ^N; exp: polynomial N^{1−2γ̄}).
**Deutung.** Der naive „Miller“-Vergleich (Strang A/B: c_N ≈ K·N·u_N/(δ_{−1})_N mit den Taylor-Koeffizienten u_N ~ τ^N von f_M
und (δ_{−1})_N ~ N^{π/c−1} der Konkurrenz) sagt die Rate τ voraus — falsch für grössere c (c = 0.8: τ = 0.423, gemessen 0.557).
Die richtige Erklärung hat Strang C (BRETT 20:37, Satz 9.5 + Nachtrag): das Abschneiden „sieht“ f auf W = D̄ ∪ D̄(a₀, R₀),
R₀ = |Fixpunkt − a₀|; der Pol von f_M antwortet mit der Rate R₀/|Pol − a₀| = √(1+τ²)·τ/(1−τ²) = sin(c/2)/cos c (identisch mit
meinem χ-Gesetz, auch im unsymmetrischen Fall: 1.3468/3.4747 = 0.3876), und Pol ∈ W ⇔ c > π/3 (= meine Schwelle f³(0) = ∞).
Bei exp hat f_K in W̄ keine Pole, nur die schwachen Randsingularitäten L, L̄ — dann bleibt allein der Miller-Mechanismus (polynomial).

## 6. Hyperbolische Möbius-Abbildung (numerisch)

g = tanh(artanh x + c), Fixpunkte ±1 reell (f anziehend bei +1, abstossend bei −1), Wurzel f_M = tanh(artanh x + c/2) =
(x + T)/(1 + Tx), T = tanh(c/2), Pol bei −1/T. Kernfunktionen ((1+x)/(1−x))^{iπ(2m+1)/c}: rein imaginärer Exponent, auf der
Einheitsscheibe beschränkt (Stärke 0) — die schwächste Konkurrenz aller Modelle; Satz K greift nicht (beide Fixpunkte auf |x| = 1).
**Lösung nahe f_M existiert** (`m3b_eta_hyp.py`): η_N = |J⁻¹F(T_N f_M)|_{0.5} fällt für c = 0.5, 1, 2 mit Rate 0.252, 0.469, 0.761
(N = 30 → 60) ≈ **T** = 0.245, 0.462, 0.762 — das ist der Miller-Quotient |Fixpunkt|/|Pol| = 1/(1/T); Cs Pol-Antwort
R₀/|Pol − a₀| ist hier kleiner. Zusammen mit §5: **Rate = max(Miller-Quotient, R₀/|Pol − a₀|)** — genau Cs Wachstumsgesetz
max(ρ/|ζ|, R₀/|ζ − a₀|) am Pol ζ (elliptisch gewinnt der zweite, hyperbolisch der erste Term).
**Aber der Zweig von N = 1 aus findet sie nicht immer** (`m3_moebius_hyp.py`, nach Korrektur des Newton-Abbruchs wiederholt):
c = 0.5 konvergiert (a₀-Fehler 1e-52 bei N = 80), c = 1.0 konvergiert nach Sprüngen (2.7e-3 bei N = 10, 6e-9 bei 20, 3e-28 bei 80),
c = 2.0 bleibt bis N = 60 auf einer anderen Lösungsfamilie (a₀ − T = −0.0119, −0.0153, −0.0125, −0.0077 bei N = 40, 45, 50, 60;
Homotopie bei N = 11), obwohl bei N = 60 eine Lösung im Abstand η = 4.5e-7 von T_N f_M existiert. ⇒ **Gegenbeispiel-Kandidat
für die Zweigwahl** (nicht für die Existenz): bei schwacher Konkurrenz (Stärke 0) und langsamer Rate kann die Fortsetzung in
N auf eine falsche Lösung springen. Für exp ist die Konkurrenz stark (s = 4.45), und der Zweig liegt numerisch bis N = 232 bei Kneser (Krawczyk-Boxen bis N = 120).

**Einheitliches Ratengesetz (numerisch, drei Möbius-Familien; Form von Strang C):** Seien R₁ = Abstand von 0 und R₀ = Abstand
von a₀ = f(0) zur nächsten Singularität der Konkurrenzlösungen (Fixpunkte), W = D(0, R₁) ∪ D(a₀, R₀). Für jede Singularität ζ
der kanonischen Lösung trägt ζ die Rate max(R₁/|ζ|, R₀/|ζ − a₀|) bei; liegt ζ in W (Rate > 1), scheitert das Verfahren.
 elliptisch symmetrisch: R₁ = 1, R₀ = √(1+τ²), ζ = 1/τ ⇒ max(τ, sin(c/2)/cos c) = sin(c/2)/cos c ✓ (gemessen bis 4 Stellen);
 elliptisch mit Fixpunkten L, L̄ (c = 0.705): max(0.348, 0.388) = 0.388 ✓;
 hyperbolisch: R₁ = 1, R₀ = 1 − T (nächster Fixpunkt +1), ζ = −1/T ⇒ max(T, (1−T)/(T+1/T)) = T ✓ (0.245/0.462/0.762 gegen
 gemessen 0.252/0.469/0.761); mit R₀ = 1 + T (fernerer Fixpunkt) käme 0.849 bei c = 2 heraus — widerlegt.
 exp: f_K hat in W̄ = D̄(0,|L|) ∪ D̄(a₀,|L − a₀|) keine Singularität ausser L, L̄ auf dem Rand ⇒ kein geometrischer Beitrag, nur
 der polynomiale Miller-Anteil N^{1−2s}. (Cs gemessenes Konvergenzgebiet von p_N bei exp ist genau dieses W.)

## 7. g = b^x und f∘f∘f = exp (numerisch)

**b < e^{1/e}** (`m4_basen.py`): b = 1.2 und √2: p_N → reguläre halbe Iterierte am **anziehenden** Fixpunkt (2 bei √2),
unabhängig berechnet über F_a (Taylor-Reihe bei a) und f = g^{−n}∘F_a∘g^n: Übereinstimmung 4e-43 (b = 1.2, N ≥ 40) bzw.
2.4e-43 (√2, N = 80), Rate √2: 0.40/Grad.
**b > e^{1/e}** (`m4b_kneser_vergleich.py`, Kneser-Werte von Strang A): p_N → Kneser für b = 1.5 und 2:
b = 1.5: |p_N(0) − sexp_b(−½)| = 9.5e-12, 1.6e-23, 5.3e-31, 1.2e-35 (N = 10, 40, 80, 120); lokale Exponenten (N^{−β}) 22.6, 27.1,
26.1, 26.8, 27.7; s_b = Re(2πi/log μ_b) = 14.25, μ_b = L_b ln b.
b = 2: 1.2e-9, 2.9e-16, 3.3e-19, 5.9e-21 (N = 10, 40, 80, 120); Exponenten 12.3, 13.2, 7.4, 9.5, 11.0 (Vorzeichenwechsel bei
N ≈ 55); s_b = 5.584.
**Polynomiales Gesetz (numerisch):** |p_N − f_K| ~ N^{−2s_b}: 2s_b = 28.5 (b = 1.5, gemessen ≈ 27), 11.2 (b = 2, gemessen ≈ 11),
8.9 (b = e; Strang B misst 7.9 auf |x| = 1). Das bestätigt über drei Basen das Miller-Bild von Strang B: Kneser-Schwanz
~ k^{−s−2}|L|^{−k}, Konkurrenz ~ k^{s−2}|L|^{−k}, Auswahl ~ N^{−2s}.
b = 10 (s_b = 3.29, |L_b| = 0.760, f(0) = 0.3939 — f(0) liegt schon bei 0.52|L_b|): konvergiert gegen Kneser, langsam:
|p_N(0) − f(0)| = 1.6e-7, 5.4e-11, 2.6e-12, 5.2e-13, 6.7e-14 (N = 10, 40, 60, 100, 140), Vorzeichenwechsel bei N ≈ 55; lokale
Exponenten 4.7, 5.6, 6.8 (N ≈ 90, 110, 130) — verträglich mit 2s−1 = 5.58, aber noch nicht asymptotisch sauber
(`daten/m4b_b10_e3.log`).
**f∘f∘f = exp** (Ordnung 3, Start p_1 = 1/3 + x; Vergleich mit sexp(−2/3), sexp(1/3) von Strang A; `daten/m4b_e3.log`):
konvergiert gegen Knesers Drittel-Iterierte: |p_N(0) − sexp(−⅔)| = 1.1e-8, 2.9e-11, 9.9e-14, 2.5e-15, 7.7e-16, 1.8e-16
(N = 10, 20, 40, 60, 80, 100), Vorzeichenwechsel bei N ≈ 55 (halbe Iterierte: N ≈ 51), Grössenordnung wie bei der halben
(a₀-Fehler −6.2e-17 bei N = 100). Erwartet, da der Kern von V(w+⅔) + V(w+⅓) + V(w) = 0 wieder mit k = ±1 (Stärke s = 4.447)
beginnt — dasselbe Gesetz N^{1−2s}; lokale Exponenten (4.6, 6.8) noch zu verrauscht für eine eigene Bestimmung.
Anmerkung: der erste b = 10-Lauf blieb bei N = 99 „stecken“ — Ursache war mein Newton-Abbruch
(absolute Toleranz 2^{−prec+60} bei Kondition ~1e38 unerreichbar ⇒ Homotopie-Kriechen); behoben durch Stagnations-Abbruch in
`kern.newton` (Test gegen polynome_bis_194 unverändert 5e-61).

## 8. Was sich auf exp überträgt

**8.1 Ein gemeinsames Bild aller Modelle (numerisch gut gestützt, als Arbeitshypothese).** Das Abschneiden wählt unter allen
Lösungen die mit dem schwächsten Wachstum der Taylor-Koeffizienten (bzw. in Abel-Koordinaten: kleinster Typ / schwächste
Singularität), und der Fehler ist in Abel-Koordinaten **ein** Kernmodus c_N e^{∓2πiw} (Möbius-Drehung: Befund 5; exp: A/B).
Die Grösse von c_N bestimmen zwei Mechanismen:
 (i) *Miller/Wachstumsvergleich* (Satz T exakt; exp: N^{1−2s}, über b = 1.5, 2, e bestätigt): Verhältnis der Koeffizientenwachstümer
     von kanonischer Lösung und Konkurrenz. Gleicher Radius (exp: beide |L|) ⇒ nur polynomial.
 (ii) *Singularitäten ζ von f im „Sichtgebiet“ W = D(0, R₁) ∪ D(a₀, R₀)*, R₁, R₀ = Abstand von 0 bzw. a₀ zur nächsten
     Konkurrenz-Singularität (Strang C, 20:37; Max-Gesetz max(R₁/|ζ|, R₀/|ζ − a₀|) in fünf Möbius-Lagen bestätigt, §5–6).
     Liegt ζ in W, gibt es keine Lösung nahe f (Möbius symmetrisch c > π/3, verschoben A = −1 schon c = 0.6); sonst geometrisch.
 Für exp: f_K hat in W̄ nur die schwachen Singularitäten L, L̄ auf dem Rand (keine Pole, keine Bahn nach ∞) ⇒ (ii) liefert nichts
 Schlimmeres als (i). Das Möbius-Modell zeigt aber: **ohne** diese quantitative Information kann das Verfahren scheitern, obwohl
 die kanonische Lösung eindeutig und reell ist (exp-ähnliche Lage: keine reellen Fixpunkte, konjugiertes Paar).

**8.2 Übertragbare Beweisschritte.**
 - Ersatzformen (A)/(B) und die Abel-Darstellung L ≅ S + I (§1): allgemein, bewiesen.
 - Satz K (Kompaktheit der Komposition ⇒ gleichmässige Inverse ⇒ Kantorovich) überträgt sich **nicht**: er braucht eine Scheibe
   um 0, die f* streng in sich abbildet (anziehender Fixpunkt). Für exp ist C nicht einmal beschränkt auf W_{|L|}. Übertragbar ist
   aber das Rezept „Blockzerlegung niedrig/hoch + endliche Intervallrechnung für den niedrigen Block“, falls eine Norm gefunden
   wird, in der der hohe Teil kontrahiert — bei exp müsste diese Norm in Abel-Koordinaten leben (vgl. Satz T: Borel-Norm, in der
   (S + I)⁻¹ = φ(D) beschränkt ist).
 - Satz T (lineares Modell) ist das saubere Vorbild für die lineare Kernaussage (Strang B: „Vermutung ist im Kern linear“):
   Auswahl = Wachstumsschranke, die die Kernfunktionen ausschliesst. Für exp entspricht das Gewichten |L|^k(1+k)^s mit
   s im Fenster zwischen Konkurrenz (k^{s_b−2}) und Kneser (k^{−s_b−2}) — Cs Fenster 4 ≤ s < 5.45.
 - **Warnung aus Befund 4:** In Fixpunkt-Gewichten ist |J_N⁻¹| auch dort beschränkt, wo das Verfahren scheitert (Möbius c > π/3).
   Ein Beweis muss die Auslöschung in J⁻¹F(q_N) quantitativ fassen (Cs Cauchy-/Pol-Antwort-Darstellung, Bs Faktor H), nicht
   |J⁻¹|·|F| abschätzen. Die Möbius-Drehung ist ein expliziter Prüfstein: jedes Beweisgerüst für exp muss dort für c < π/3 gelten
   und für c > π/3 versagen.
 - Parabolischer Fall (Modell P) und b < e^{1/e}: das Verfahren ist mit der regulären Iteration verträglich; bei e^{1/e} selbst
   (parabolisch, Fixpunkt e ≠ 0) siehe §7 (offen).

**8.3 Konkrete Vorschläge an die anderen Stränge.**
 - B/C/D: Ein Beweisansatz sollte zuerst am Möbius-Modell (A ± iB = L, L̄, c variabel) getestet werden — dort ist alles explizit
   (f_M, Kernfunktionen ((1−ix')/(1+ix'))^{π/c}·(1 + f²)…), und die Schwelle c* = 1.2029 (für Fixpunkte L, L̄) ist bekannt.
 - A: Das Gesetz c_N ~ N^{1−2s_b} lässt sich mit b = 1.5 (s = 14.2, sehr schnelle Konvergenz) und b = 10 (s = 3.3, langsam)
   als Prüfstein für die Konstante K nutzen.
 - Daten: `gemeinsam/modell_moebius_N80.json` — p_N (N ≤ 80, 60 Stellen), Taylor-Koeffizienten von g und f_M für die Drehung
   c = 0.8, c = 1.0 und das „Möbius-exp“-Modell (Fixpunkte L, L̄, g(0) = 1).

---

# Runde 2 (Auftrag 21:23): Bs Gerüst (H*) ⇒ (LP), quadratische Kontrolle — im Möbius-Modell durchgezogen

**Kurzfassung Runde 2.** (1) Bewiesen: lineares Problem = eine Pol-Antwort (M1); volles Problem p = Π^{(p)}f, p − f = r·E^{(p)}(ζ_p,·)
(M4); explizite Green-Funktion (M3). (2) Bewiesen unter (H*)_M: (LP) für alle c < π/3 mit Rate e^{−g(ζ_p)} (Satz M-lin); unter (H*) für die
Sekanten-Familie auf einer festen Kugel: Existenz einer Lösung des NICHTLINEAREN Systems mit derselben Rate (Satz M-voll, Brouwer) — die
quadratische Kontrolle entfällt. Beweis bricht genau, wo ζ_p ∈ Ḡ' (c ≥ π/3). (3) Numerisch: (H*)_M mit C_* ≈ 1, κ = 0 bis N = 160, gleichmässig
über die Sekanten-Kugel; das wahre Gesetz ist e^{−Ng(ζ)} (Green-Funktion von ℂ∖Ḡ', konstanter Vorfaktor), nicht Cs Max-Formel; bei exp ist
T(p) = Π^{(p)}f_K auf einer festen Kugel eine Kontraktion (Lipschitz 3e-7, 2e-8, 1e-9…4e-9 für N = 40, 80, 160). (4) Nicht gelungen: (H*)_M selbst zu beweisen.

Modell: symmetrische Drehung, Fixpunkte ±i, f = f_M = (x+τ)/(1−τx), τ = tan(c/2), a₀ = τ, Pol ζ_p = 1/τ; R₁ = 1, R₀ = |i − a₀| = √(1+τ²);
Ḡ' = D̄(0,1) ∪ D̄(a₀,R₀) (die Randkreise schneiden sich in ±i); J_∞u = u∘f + m·u, m = f'∘f; Π_N = J_N⁻¹T_NJ_∞ (Strang B, §1.2);
E_N(ζ,z) = Π_N k_ζ(z) − k_ζ(z), k_ζ = 1/(ζ − z). Werkzeug: `moebius_lin.py` (J_N an der EXAKTEN Lösung, Pol-Antworten in arb).

## R2.1 Exakte Strukturaussagen (bewiesen)

**Lemma M1 (lineares Problem = EINE Pol-Antwort).** Π_N f_M − f_M = r·E_N(ζ_p, ·), r = (1+τ²)/τ² = −Res_{ζ_p} f_M.
*Beweis.* Π_N hängt nur von den N-Jets an 0 und a₀ ab, also Π_N f(z) = (1/2πi)∮_Γ f(ζ)Π_N k_ζ(z)dζ (Γ kleine Kreise um 0, a₀).
ζ ↦ Π_N k_ζ(z) ist rational mit Polen nur in 0, a₀ und = 1/ζ + O(ζ⁻²) (Bs Lemma 7.1). Γ auf |ζ| = R → ∞ aufblasen: dabei wird nur
der Pol ζ_p von f überstrichen; ∮_{|ζ|=R} → 2πi f(∞). Also Π_N f(z) = f(∞) − Res_{ζ_p}[f Π_N k_·(z)]; dieselbe Rechnung mit k_ζ
statt Π_N k_ζ gibt f(z). Differenz. ∎ Numerisch bestätigt (N = 10, 30; 12 Stellen).
**Lemma M2 (Cayley-Bild).** w = (1+ix)/(1−ix): f ↦ w ↦ λw (λ = e^{ic}), J_∞ ↦ 𝓛Δ = Δ(λw) + λΔ(w) (mal Einheit φ'(λ²w)),
P_N ↦ 𝒫_N = {(i/2)(1+λw)²b(w)/(1+w)^N : deg b ≤ N}, Abschneiden = N-Jet bei w = 1, f_M ↦ U = (λ²w² − 1)/2 (Polynom!), Ḡ' ↦ Sektor
{−π/2 ≤ arg w ≤ π/2 + c/2}, Komplement Σ = Sektor der Öffnung ω = π − c/2, Pol ζ_p ↦ w_p = −1/λ = e^{i(π−c)}.
Folgerung: **ζ_p ∈ Ḡ' ⇔ c > π/3** (genau die Schwelle aus Runde 1). Kerne: w^{1+π(2m+1)/c}. (Direkte Rechnung; numerisch gegen die
x-Form geprüft, `r2_01_cayley_symbolisch.py`.) Für kleine N ist die Lösung kein geschlossener Ausdruck (sympy, N ≤ 4: grosse
palindromische Polynome in λ), und die Folge E_N(ζ,z) erfüllt keine kurze lineare Rekursion in N (getestet bis Ordnung 6, Grad 3).
**Lemma M3 (Green-Funktion explizit).** g(ζ) := g_{ℂ∖Ḡ'}(ζ, ∞) = log|(u − ū_∞)/(u − u_∞)|, u = (e^{−i(π/2+c/2)}w)^{π/ω},
u_∞ = e^{iπ(π/2−c/2)/ω}; am Pol: **e^{−g(ζ_p)} = sin(πc/(2ω)) / sin(π(π−2c)/(2ω))**, ω = π − c/2. g(ζ_p) > 0 ⇔ c < π/3.
(Mit ω = π statt π − c/2 erhielte man sin(c/2)/cos c — Cs Max-Formel ist die Green-Funktion ohne Eckkorrektur.)
**Lemma M4 (volles nichtlineares Problem = Pol-Antwort eines Sekanten-Operators).** Für p ∈ P_N sei
𝒥_p u := u∘p + M_{p,f}·u, M_{p,f} = (f∘p − f∘f)/(p − f) = (1+τ²)/((1−τp)(1−τf)), Π^{(p)} := (T_N𝒥_p|_{P_N})⁻¹T_N𝒥_p,
E^{(p)}_N(ζ,z) := Π^{(p)}k_ζ(z) − k_ζ(z). Dann gilt: p löst T_N(p∘p) = T_N g ⇔ p = Π^{(p)}f_M, und in diesem Fall exakt
**p − f_M = r·E^{(p)}_N(ζ_p, ·)**. *Beweis.* p∘p − f∘f = (p−f)∘p + (f∘p − f∘f) = 𝒥_p(p − f); also T_N𝒥_p(p − f) = 0 ⇔ Π^{(p)}(p − f) = 0
⇔ p = Π^{(p)}f (Π^{(p)} reproduziert P_N). ζ ↦ Π^{(p)}k_ζ(z) ist rational mit Polen nur in 0 und p(0) (Beweis wie Bs Lemma 7.1),
dann Residuenrechnung wie in M1. ∎ Numerisch an p = p_N: max_{D̄}|p_N − f − rE^{(p_N)}(ζ_p,·)| = 1.6e-41, 1.4e-41 (N = 40, 80).
Bedeutung: **es gibt keine separate „quadratische Kontrolle“** — das nichtlineare Problem ist ein Fixpunkt p ↦ Π^{(p)}f, und dessen
Abstand zu f ist eine Pol-Antwort. Das gilt wörtlich auch für exp (mit Bs Konturintegral über die Schlitze statt eines Residuums).

## R2.2 Sätze (bewiesen unter der jeweils genannten linearen Voraussetzung)

Sei Λ_N := {ζ : g(ζ) = 1/N} die Niveaulinie; dort ist dist(ζ, Ḡ') ≥ c₀/N (Ecken: ≍ N^{−ω/π}, also grösser).
**(H*)_M:** |Π_N k_ζ(z)| ≤ C_*N^κ/dist(ζ,Ḡ') für z ∈ D̄, ζ ∈ Λ_N, N ≥ N₀.
**Satz M-lin.** Unter (H*)_M gilt für alle c < π/3: sup_{D̄}|Π_N f_M − f_M| ≤ C·N^{κ+1}·e^{−N g(ζ_p)} → 0.
*Beweis.* Lemma M1 und Bs Lemma 7.2 (Maximumprinzip in ζ auf {g > 1/N}, E_N(·,z) = O(ζ^{−N−2})): |E_N(ζ_p,z)| ≤
max_{Λ_N}|E_N(·,z)|·e^{−(N+2)(g(ζ_p) − 1/N)} ≤ (C_*N^κ + 1)(N/c₀)·e·e^{−N g(ζ_p)}. ∎
**Wo der Beweis über der Schwelle bricht:** für c ≥ π/3 ist ζ_p ∈ Ḡ' (Lemma M2), g(ζ_p) ≤ 0, das Maximumprinzip im Äusseren erreicht den
Pol nicht mehr — und numerisch wächst E_N(ζ_p,·) dann tatsächlich (Runde 1, η_N; Strang C, Wachstum innerhalb W).
**Satz M-voll (Existenz + Konvergenz für das NICHTLINEARE Problem).** Sei c < π/3, η₀ > 0, und für N ≥ N₀ und alle reellen p ∈ P_N mit
‖p − f_M‖_{D̄} ≤ η₀ gelte: T_N𝒥_p|_{P_N} invertierbar (dann stetig in p) und **(H*_sek):** |Π^{(p)}k_ζ(z)| ≤ C_*N^κ/dist(ζ,Ḡ') für
z ∈ D̄, ζ ∈ Λ_N. Dann hat das abgeschnittene System für alle grossen N eine reelle Lösung p_N mit
sup_{D̄}|p_N − f_M| ≤ C·N^{κ+1}·e^{−N g(ζ_p)}.
*Beweis.* B := {p ∈ P_N reell : ‖p − f_M‖_{D̄} ≤ η₀} ist kompakt, konvex, nicht leer (T_N f_M ∈ B für grosse N). T(p) := Π^{(p)}f_M ist
stetig und reell. Nach M4-Rechnung (gilt für JEDES p ∈ B, nicht nur Lösungen): T(p) − f_M = r·E^{(p)}_N(ζ_p,·); wie in Satz M-lin
(Lemma 7.2 gilt für E^{(p)}, da p(0) ∈ Ḡ'°) ist ‖T(p) − f_M‖_{D̄} ≤ ε_N := C N^{κ+1}e^{−Ng(ζ_p)} ≤ η₀. Also T(B) ⊂ B; Brouwer ⇒
Fixpunkt p_N = Π^{(p_N)}f_M ⇒ T_N(p_N∘p_N) = T_N g (M4) und ‖p_N − f_M‖ ≤ ε_N. ∎
Offen bleibt hier nur die Zweig-Identität (dass die Fortsetzung von N = 1 aus diesen Fixpunkt trifft): dafür genügt Eindeutigkeit in B,
z. B. aus einer Lipschitz-Schranke für p ↦ E^{(p)}_N(ζ_p,·) (numerisch winzig, s. R2.3).
**Übertrag auf exp (Schablone):** Bs Satz B-Max benutzt nur (i) Rationalität in ζ mit Polen in den Entwicklungspunkten, (ii) Lemma 7.2,
(iii) die Schlitzstruktur von f_K. Alles gilt für Π^{(p)} mit p(0) statt a₀. Daher: **(H*) gleichmässig für die Sekanten-Operatoren 𝒥_p,
‖p − f_K‖_{D̄} ≤ η₀ ⟹ (Brouwer) Existenz einer Lösung p_N mit sup_{D̄}|p_N − f_K| ≤ C N^{κ+1−4.69}.** Die quadratische Kontrolle (Q)
ist damit in die Gleichmässigkeit von (H*) über eine FESTE Kugel verlagert.

## R2.3 Numerik zu den Voraussetzungen (`r2_03_hstern.py`, `r2_04_green.py`, `r2_05`, `r2_06_sekante.py`, Logs in `daten/r2_*`)

**(H*)_M im Band** (dist·max_{D̄}|Π_N k_ζ|, ζ im Abstand t/N, an der Ecke t/N^{ω/π}; t = 0.5/2/8):

| Ort | c = 0.6, N = 40 | N = 160 | c = 1.0, N = 40 | N = 160 |
|---|---|---|---|---|
| Ecke i, Winkelhalbierende | 0.406/0.887/1.01 | 0.39/0.873/1.00 | 0.348/0.877/1.02 | 0.331/0.858/1.03 |
| links −1 | 0.38/0.847/0.999 | 0.375/0.846/0.999 | 0.363/0.83/0.999 | 0.358/0.828/0.999 |
| rechts a₀+R₀ | 0.034/0.123/0.36 | 0.0087/0.034/0.123 | 0.018/0.068/0.226 | 0.0045/0.018/0.068 |
| ∂D(a₀) oben 60° | 0.051/0.179/0.469 | 0.013/0.051/0.179 | 0.025/0.093/0.294 | 0.0063/0.025/0.093 |
| ∂D oben 120° | 0.367/0.833/0.998 | 0.343/0.794/0.986 | 0.342/0.805/0.997 | 0.321/0.771/0.985 |

⇒ (H*)_M mit C_* ≈ 1.03, κ = 0 bis N = 160 — dieselbe Tabelle wie Bs b25 für exp.
**Gesetz R ist die Green-Funktion, nicht Cs Max-Formel** (`r2_04_green.py`, c = 0.6): max_{D̄}|E_N(ζ,·)|·e^{N g(ζ)} für N = 40/80/120/160:
ζ_p = 3.233: 0.117/0.127/0.133/0.136; ζ = −2: 0.492/0.491/0.491/0.490; ζ = −1.3+0.8i: 1.23/1.23/1.23/1.23; ζ = 1.6i: 1.23/1.34/1.40/1.43;
ζ = 2: 0.30/0.33/0.34/0.35; ζ = 1.2+1.2i: 0.61/0.68/0.72/0.73 — **konstant**, während |E_N|/max(1/|ζ|, R₀/|ζ−a₀|)^N exponentiell wächst
(ζ_p: 0.46 → 32.7; ζ = −2: 2.65 → 411). Raten: ζ_p: e^{−g} = 0.3705 gegen Max-Formel 0.3581; ζ = −2: 0.5215 gegen 0.5. Auch die Prony-Moden
aus Runde 1 (c = 0.6: 0.3709, c = 0.8: 0.5715 extrapoliert) treffen e^{−g} (0.3705, 0.5726), nicht ρ (0.358, 0.559). c = 1.0: Pol-Quotient
0.55 → 0.71 (N = 40 → 160), langsamer Vorfaktor ~N^{0.18}. ⇒ **Bs Lemma 7.2 ist scharf**: das Abschneiden approximiert auf Ḡ' mit der
Bernstein–Walsh-Rate des Kompaktums Ḡ' (wie eine Interpolation mit gleichverteilten Knoten). Für exp passt dazu Cs Gitter-Befund
„systematischer Vorlauf +0.016“: Bs e^{−g}(2.2) = 0.805 gegen Max-Formel 0.793 ist genau +0.015 im Logarithmus.
**Quadratischer Anteil** (`daten/r2_05_quadratisch.log`): max_{D̄}|p_N − Π_N f_M| / (max_{D̄}|Π_N f_M − f_M|)² = 0.022, 0.045, 0.066 (c = 0.6,
N = 20, 40, 60; danach Rechengenauigkeit), 0.18, 0.25, 0.40, 0.48, 0.58 (c = 0.8, N = 20…100): zweite Ordnung mit Faktor ≈ N^{0.8…1}.
**(H*_sek) gleichmässig über eine Kugel** (`r2_06_sekante.py`, c = 0.6): p = p_N + η·q, q zufällig mit ‖q‖_{ℓ¹} = 1, η = 0, 1e-8, 1e-4,
1e-2, 0.1: Band-Grössen (Ecke t = 0.5/8, links, rechts) bei N = 80: 0.392/1.01, 0.377/0.999, 0.017/0.219 für η ≤ 1e-4; 0.376/1.01 bzw.
0.477/1.01 für η = 0.01/0.1 (Rest unverändert); max|E^{(p)}(ζ_p,·)| = 4.1e-36, 4.1e-36, 4.1e-36, 6.7e-36, 1.5e-37 — d. h. T bildet die
Kugel vom Radius 0.1 in eine Kugel vom Radius ~1e-35 ab. Die Voraussetzung von Satz M-voll gilt numerisch mit grossem Spielraum.
**Stabilität im Gewicht 1^k** (Fixpunkt-Radius): |J_N⁻¹|_{ℓ¹} = 2.13, 2.43, 2.61, 2.73, 2.82 (c = 0.6, N = 40…200), 2.15 … 2.80 (c = 1.0):
nur ≈ 0.44·ln N, **keine Eckspitze**. Erklärung (Vermutung, Bs/Ds Streckband [N/|μ|, N]): bei der Drehung ist |f'(±i)| = 1, das Band ist leer;
bei exp |μ| = |f_K'(L)| = 1.17 > 1. ⇒ Die Möbius-Drehung modelliert Auswahl, Schwelle und Green-Gesetz, aber NICHT Ds Eckspitze.

**Die Brouwer-Reduktion bei exp numerisch** (`r2_07_exp_sekante.py`, Kneser-Koeffizienten k ≤ 1000 und p_N aus `gemeinsam/`;
M_p = Σ_k f_k h_k ohne Division durch p − f): T(p) := Π^{(p)}f_K für p = p_N + η·q, q zufällig, Σ|q_k||L|^k = 1;
max_{|z|=0.999|L|}|T(p) − f_K| bei N = 40: 4.82e-8 (η = 0, = |p_N − f_K|), 4.82e-8 (1e-6), 4.85e-8 (1e-3), 5.1e-8 (1e-2), 6.4e-8 (0.05);
N = 80: 1.59e-9, 1.59e-9, 1.60e-9, 1.79e-9, 2.96e-9. Und |T(p) − p_N| = 1.8e-13, 2.8e-10, 3.1e-9, 1.6e-8 (N = 40) bzw. 1.6e-14, 5.5e-12,
2.0e-10, 1.4e-9 (N = 80) für η = 1e-6 … 0.05; N = 160 (`daten/r2_07_exp_sekante_gross.log`): |T(p) − f_K| = 5.06e-11 (η = 0),
4.97e-11 (1e-3), 2.25e-10 (0.05), |T(p) − p_N| = 1.0e-12 (1e-3), 1.8e-10 (0.05) — **T ist auf der Kugel vom Radius 0.05 eine Kontraktion
mit Lipschitz-Konstante ≈ 3e-7 (N = 40), 2e-8 (N = 80), 1e-9…4e-9 (N = 160)** (numerisch, Zufallsrichtungen). N = 300 ist mit den
Kneser-Daten nicht prüfbar (Koeffizientenfehler ~1e-47 gegen Kondition ~1e84: schon T(p_N) ≠ p_N, Abweichung 925). Damit wäre p_N der EINZIGE Fixpunkt in dieser festen Kugel (Banach statt Brouwer) — das erledigte
auch die Zweig-Identität. Für exp reduziert sich die volle Vermutung also auf zwei lineare Aussagen über die Sekanten-Familie:
(LP_sek) sup_{p∈B}‖Π^{(p)}f_K − f_K‖_{D̄} → 0 und (Lip) sup_{p∈B}‖D_p Π^{(p)}f_K‖ < 1 (D_p Π^{(p)}f = −Π^{(p)}-Antwort auf (δ∘p)'·(Π^{(p)}f − f)-artige
kohärente Eingaben, wieder Pol-Antworten).

## R2.4 Was NICHT gelungen ist (ehrlich)

(H*)_M selbst ist auch im Möbius-Modell **nicht bewiesen**. Versucht: (a) geschlossene Lösung (sympy bis N = 4: keine); (b) Rekursion in N
(keine kurze); (c) Kontraktion in Sup-Normen auf Kreisen um w = 1 (Bernstein–Walsh-Wachstum schlägt den Faktor h^{−N}, h = (1+w)/(1+λw)); (d) „Taylor-
Eingaben“ statt Pol-Eingaben für den quadratischen Term (Antworten auf k_ζ mit ζ IN Ḡ' wachsen exponentiell — unbrauchbar; darum der Weg über M4).
Struktur, die ein Beweis ausnutzen müsste: Fehler = ein Kernmodus w^{1±π/c} (Runde 1, Befund 5); Nullstellen von S (Cayley-Lösung) liegen auf
dem Einheitskreis bei gleichabständigen Winkeln (Abstand c, = Bahn von g) und paarweise r ↔ 1/r̄ auf zwei Strahlen aus 0 (Schnitte der
Kernfunktionen) — ein Gleichgewichtsproblem vom Stahl-Typ (`r2_02_nullstellen.py`). Ein Riemann–Hilbert-Ansatz scheint der realistische Weg.

---

# Runde 3 (Auftrag 22:07): (H*)_M beweisen — Ergebnis: nicht gelungen; was gezeigt ist und woran es hängt

**Kurzfassung.** (H*)_M ist weiterhin unbewiesen. Bewiesen bzw. sauber festgestellt: (1) es gibt keine geschlossene Form vom
Trappmann-Typ (det J_N irreduzibel, N ≤ 5), und der Grund ist strukturell; (2) J_N ist exakt die Resolvente eines „komprimierten
Drehoperators“ 𝒮_N an der Stelle −λ — die Stabilität ist eine Pseudospektral-Aussage, keine Eigenwert-Aussage; (3) die
Negativrichtung lässt sich auf eine Überkonvergenz-Aussage zurückführen (Residuen-Lemma), diese aber nicht beweisen.
Numerisch: (H*)_M mit C_* = 1.0, κ = 0 bis N = 320 (c = 0.6; Ecke 0.389/0.873/1.0, links 0.375/0.846/0.999, oben 120° 0.285/0.681/0.95 bei N = 320; `daten/r3_03_hstern_gross.log`), Pol-Antwort·e^{Ng(ζ_p)} = 0.117, 0.127, 0.133, 0.136, 0.142, 0.147 (N = 40 … 320), Max-Formel dagegen um den Faktor 8483 daneben (N = 320).

## R3.1 Keine geschlossene Form (bewiesen für N ≤ 5; Begründung allgemein)

Symbolisch (sympy, `r3`-Rechnung im Log): **det J_N = (1+τ²)^N·Q_N(τ²)/(1−τ²)^{2N+2}** mit Q_N ∈ ℤ[τ²] irreduzibel über ℚ vom Grad
(N+1)(N+2)/2 + 1 (4, 7, 11, 16, 22 für N = 1…5), Q_N(0) = 2^{N+1}. Im Cayley-Bild ist die Systemdeterminante ebenso ein irreduzibles
palindromisches Polynom in λ (Grad 8, 14, 22, 32). Die Koeffizienten der Pol-Antwort sind also rationale Funktionen mit irreduziblem
Nenner wachsenden Grades — kein Produkt-/hypergeometrischer Ausdruck. Die Folge E_N erfüllt auch keine kurze Rekursion in N (Runde 2).
**Warum Trappmanns Technik nicht greift (bewiesen):** Bei f(x) = bx liegt der Entwicklungspunkt im Fixpunkt — das System ist dreieckig.
Im Cayley-Bild ist die Möbius-Drehung zwar linear (w ↦ λw), aber der Pol des Polynomraums liegt bei w = −1 (x = ∞), und −1 ist kein
Fixpunkt. Legt man den Pol künstlich in einen Fixpunkt (Polynome in w, Pol bei w = ∞ ↔ x = −i), so wird J_N exakt diagonal:
𝓛w^k = (λ^k + λ)w^k, Π_N u = 𝓛⁻¹T_1[𝓛u] — explizit lösbar, aber mit kleinen Nennern |λ^{k−1} + 1| = 2|cos((k−1)c/2)| (Siegel-Typ).
Das echte Problem hat diese Resonanzen numerisch NICHT (det J_N(c) ohne Vorzeichenwechsel auf 199 Gitterpunkten c ∈ (0, π/3), N ≤ 80).

## R3.2 J_N als Resolvente des komprimierten Drehoperators (bewiesen) und sein Spektrum (numerisch)

Für Δ ∈ 𝒫_N bestimmt der N-Jet bei w = 1 das Element; 𝒮_N := (Jet bei 1 ↦ Jet bei λ, in den gedrehten Koordinaten w/λ − 1).
Aus E(λw) + λE(w) = O((w−1)^{N+1}) folgt exakt **J_N ≅ 𝒮_N + λI = λ(I + λ⁻¹𝒮_N)** (bis auf feste invertierbare Koordinatenwechsel).
Also: (H*) ⇔ Resolventenschranke von λ⁻¹𝒮_N bei −1 auf kohärenten Eingaben. Spektrum (`r3_02_spektrum_schnell.py`, acb-Eigenwerte):
c = 0.6: zwei Eigenwerte ≈ λ^{±1}, die übrigen N−1 **reell positiv** (N = 20, 40, 60), Spanne 6e-4 … 1.3e4 (N = 20) bis 1.6e-11 … 3.2e14
(N = 60) — also extrem nicht-normal; c = 0.9: ebenso ab N = 40; c = 1.0, 1.04 (noch unter π/3!): 14…40 negative/komplexe Eigenwerte,
das Paar weicht von λ^{±1} ab; kein Eigenwert nahe −1 (Abstand ≥ 0.6 in allen Fällen, auch c = 1.2). ⇒ Weder die Lage der Eigenwerte
noch Positivität erklären die Schwelle; entscheidend ist die Nicht-Normalität (Pseudospektrum). Das Paar ≈ λ^{±1} gehört zu
span{1, w²} ∋ U = (λ²w² − 1)/2 = Bild von f_M: es konvergiert gegen λ^{±1} genau in den Fällen, in denen Π_N f_M schnell konvergiert.
**Rang-Struktur** (`r2`-Test): der bivariate Zähler Ω_N(ζ,z) = −E_N(ζ,z)ζ^{N+1}(ζ−a₀)^{N+1}(ζ−z) ist numerisch fast Rang 1
(σ₂/σ₁ = 0.026, 0.020, 0.017 für N = 10, 20, 30): E_N ≈ −u_N(z)v_N(ζ)/(ζ^{N+1}(ζ−a₀)^{N+1}(ζ−z)) — eine „effektive Interpolation“ mit
Knotenpolynom u_N. Das ist die natürliche Ansatzform für einen Riemann–Hilbert-Beweis (Gleichgewichtsmass von ∂Ḡ' ⇒ Green-Gesetz).

## R3.3 Negativrichtung (c > π/3)

**Residuen-Lemma (bewiesen).** Konvergiert Π_N f_M → f_M gleichmässig auf einem Kreis |ζ − ζ_p| = ε, so ist das unmöglich:
(ζ − ζ_p)(f_M − Π_N f_M) ist im Kreisinneren holomorph, auf dem Rand ≤ ε·δ_N, also |Res_{ζ_p}f_M| = (1+τ²)/τ² ≤ εδ_N → 0. ∎
Also bricht (LP) für c > π/3, SOBALD man weiss, dass Konvergenz auf D̄ Überkonvergenz auf ganz W = D(0,1) ∪ D(a₀,R₀) erzwingt (dann läge
ζ_p ∈ W im Konvergenzgebiet). Diese Überkonvergenz ist numerisch klar (Strang C: Konvergenzgebiet = W; Runde 2: Wachstum der Pol-Antwort
im Inneren von W), aber nicht bewiesen: ein Zwei-Konstanten-Bootstrap über J_∞E = O(x^{N+1}) (Kontrolle auf V ⇒ auf f(V_klein)) erreicht
von D̄ aus nur eine Scheibe vom Radius ≈ (1+τ²)(1−τ)/(1+τ) um a₀ (τ = 0.65: 0.30), ζ_p liegt aber 1/τ − τ = 0.89 entfernt — die
Bernstein–Walsh-Abschätzung auf D̄ ist zu grob. Status: **Skizze, Lücke = Überkonvergenz-Lemma**.

## R3.4 Übertrag auf exp (Abel-Bild)

Beide Modelle sind in Abel-Koordinaten (f = Translation um ½, J_∞ ≅ S + I, Kerne e^{±2πi(2m+1)w}) gleich; verschieden sind
(a) der Polynomraum: Möbius α = arctan(x)/c, P_N∘α⁻¹ = {B(cw)/cos^N(cw)} ist **periodisch** in w (Periode π/c, elliptische Drehung);
bei exp ist P_N∘sexp nicht periodisch; (b) die Fixpunkt-Asymptotik: Möbius α ≈ ∓(i/2c)log(x ∓ i) — „gerade“ Halbstreifen (|μ| = 1);
exp α ≈ log(x − L)/L mit μ = L, |μ| = 1.37 > 1 und arg μ ≠ 0 — „schräger“ Halbstreifen, daraus Ds Eckband [N/|μ|, N] und Eckmode
(Rate ≈ 0.074·ln|μ|, Strang D). **Übertragbar (wörtlich):** M1/M4-Konturdarstellung, Bs Lemma 7.1/7.2, Green-Gesetz, Brouwer/Kontraktion
über die Sekanten-Familie (bei exp numerisch Lipschitz 3e-7 → 1e-9). **Nicht übertragbar:** dass (H*) mit polynomialem κ für ALLE N gilt —
bei exp kippt es ab N ≈ 300–420 (Strang C) und braucht Cs Poisson–Jensen-Variante plus Ds Eckanalyse. **Arbeitsteilung, die daraus folgt:**
Die Möbius-Drehung ist genau „regulärer Teil ohne Ecke“: jedes Argument, das Strang C für den regulären Teil von (H*) bei exp entwickelt,
muss im Möbius-Modell bereits ALLES liefern (dort gibt es keinen Eckteil) — es ist der Prüfstein; Ds kubisches Modell ist der Prüfstein für
den Eckteil. Für exp gibt es kein Gegenstück zur Schwelle: f_K hat in W̄ keine Singularität ausser L, L̄ (Ecken).

## R3.5 Prüfstein „Möbius-exp“ für Strang C (numerisch, `r3_04_moebius_exp_hstern.py`, `daten/r3_04_moebius_exp.log`)

Drehung um L, L̄ (Fixpunkte von exp), c = 0.705104 (g(0) = 1): a₀ = 0.47801 (exp: 0.49856), R₀ = 1.34676 (exp: 1.34935), Aussenwinkel
der Ecke ω = 159.8° (exp: 158.9°), Pol der Wurzel bei 3.95268 — also praktisch Knesers W-Geometrie, aber ohne Eckmode und mit einem
glatten Ziel. Allgemeiner Aufbau `moebius_lin.setup_allg` (geprüft gegen den symmetrischen Fall). dist·max_{|z|=|L|}|Π_N k_ζ| (t = 0.5/2/8):
Ecke L (Winkelhalbierende, Abstand t/N^{0.888}): 0.315/0.797/1.01 (N = 40), 0.312/0.800/1.01 (80), 0.311/0.798/1.01 (120), 0.308/0.792/1.01 (160);
links −|L|: 0.298/0.752/0.995 … 0.293/0.750/0.996; rechts a₀+R₀: 0.027/0.10/0.31 → 0.0069/0.027/0.10; ∂D bei 120°: 0.29/0.74/0.99 → 0.27/0.70/0.97.
Pol: max_{|z|=|L|}|E_N(ζ_p,·)|·e^{Ng(ζ_p)} = 0.120, 0.131, 0.137, 0.142 (N = 40…160), e^{−g(ζ_p)} = 0.4092 (Max-Formel 0.3876).
⇒ Auf ganz D̄ (|z| ≤ |L|) gilt (H*) mit C_* = 1.01, κ = 0, und das Green-Gesetz mit konstantem Vorfaktor — in exp-Geometrie.
(Die erste Rechnung dieses Abschnitts hatte einen Fehler: ζ_p ging einmal 25-stellig, einmal als Gleitkommazahl ein; dadurch sättigte
die Pol-Antwort bei 7.5e-26. Behoben, Werte oben neu.)

## 9. Offene Punkte (Stand Runde 3)

 - **(H*)_M beweisen** — der einzige fehlende lineare Kern; danach sind Satz M-lin und (mit Gleichmässigkeit über die Sekanten-Kugel)
   Satz M-voll unbedingt. Kein elementarer Weg gefunden (keine geschlossene Form, keine Rekursion, keine Normkontraktion, Spektrum nicht
   entscheidend). Realistisch: Riemann–Hilbert/Steepest Descent für die Cayley-Gleichung A(w)R(λw) + B(w)R(w) = (w−1)^{N+1}H(w)
   (R, H Polynome; A = (1+λ²w)²(1+w)^N(w−Z), B = λ(1+λw)^{N+2}(λw−Z)), Ansatz „effektive Interpolation“ (Rang-1-Struktur R3.2).
 - Überkonvergenz-Lemma (Konvergenz auf D̄ ⇒ auf W) — würde mit dem Residuen-Lemma die Negativrichtung c > π/3 beweisen.
 - Zweig-Identität in Satz M-voll: Eindeutigkeit in der Kugel (numerisch Kontraktion mit winziger Lipschitz-Konstante).
 - exp: Green-Gesetz für Cs R_N(ζ) direkt prüfen (Möbius und Möbius-exp: konstanter Vorfaktor); Lipschitz von T(p) = Π^{(p)}f_K bei
   N ≥ 300 (braucht Kneser-Koeffizienten mit > 90 Stellen).
 - Nachrangig: hyperbolische Familie (Zweig bei c = 2), b = e^{1/e}.
