# Strang B — Reformulierung & Beweisstrategie

Stand: 10.10.2026, 20:30. Kennzeichnung: **bewiesen** / **Skizze** (Lücken benannt) / **numerisch** / **Vermutung**.

Bezeichnungen: p = p_N = Σ a_k x^k, a_0 = p(0), h = p − a_0; f = f_K = T(A + ½) (Kneser), A = slog, T = sexp;
L = 0.3181 + 1.3372i, μ := f_K'(L) = e^{L/2} (|μ| = 1.1717, arg μ = 38.3°); γ := 2πi/L = 4.44695 + 1.05794i;
E_N := p∘p − exp (Nullstelle der Ordnung N+1 bei 0); J_∞^f u := u∘f + (f'∘f)·u; J_N := T_N J_∞|_{P_N}
(P_N = Polynome vom Grad ≤ N, T_N = Taylor-Abschnitt bei 0).

## 0. Kurzfassung

1. **Die Auswahl ist eine Ein-Parameter-Frage (numerisch, sehr sauber).** In Kneser-Koordinaten ist
   θ_N(w) := A(p_N(T(w))) − w − ½ im Korridor bis auf ~10⁻³⁴ (N = 194) gleich c_N e^{2πiw} + c̄_N e^{−2πiw}: Moden 0, ±2
   sind ~|c_N|², ±3 sind ≤ 10⁻¹⁵|c_N|. Im Inneren ist p_N also eine *exakte* halbe e-Funktion aus der Familie
   f_c = T(A + ½ + c e^{2πiA} + c̄ e^{−2πiA} + O(c²)) mit c = c_N. Gemessen (N = 18 … 194):
   **c_N = C·N^{1−2γ̄}·(1 ± 5 %)**, C ≈ 0.18 − 0.08i; Betrag N^{−7.894} und Phase 2·Im γ·log N passen beide.
2. **Nur polynomialer Abstand.** max_{|x|=1}|p_N − f_K| ~ N^{−7.9} (= |c_N|), |p_N(L) − L| ~ N^{−4.83}. Jede Beweisidee mit
   geometrischen Reserven ist damit ausgeschlossen; Abschätzungen dürfen nur polynomiale Faktoren verlieren.
3. **Geometrie:** p_N konvergiert auf G' := D(0,|L|) ∪ D(a_0,|L−a_0|) (Strang C; meine Punkte 1.45, 1.6 bestätigen).
   G' hat bei L, L̄ eine einspringende Ecke (Innenwinkel 201.07°, Aussenwinkel απ mit α = 0.8830). Polynome vom Grad N
   lösen an so einer Ecke die Skala N^{−α} auf (Green-Funktion ~ r^{1/α}); Bestapproximation von (z − L)^β gibt N^{−αβ}
   = N^{−0.883·5.447} = N^{−4.81} — passt zur Randrate N^{−4.83 ± 0.05} (b08, b13). (Früher hier fälschlich β/θ = 4.876.)
4. **Reformulierungen (bewiesen):** Jet-/Dreiecksform J_N = H(S + Γ) (§1.1); J-Projektion Π_N, die ker J_∞ exakt
   annulliert (§1.2); Kern von J_∞ auf H(D) explizit = {T'(A+½)·θ(A) : θ ganz, θ(w+½) = −θ(w)} und **trivial auf Funktionen
   mit Wachstum < |z−L|^{−3.447}** (§1.3, lineares Gegenstück des Eindeutigkeitssatzes); θ-Gleichung für p_N (§1.4);
   Zwei-Konstanten-Schranke für E_N (§1.5); quantitativer Bôcher (Fourier-Lemma, §1.6).
5. **Ideen 1–5:** (1) Lyapunov–Schmidt trägt in der Form „eine komplexe Amplitude c_N“; (2) Fixpunktform ist *keine*
   Kontraktion: sie ist exakt eine vorkonditionierte Richardson-Iteration a ↦ a − (HS)⁻¹F(a) mit Linearisierung −S⁻¹Γ,
   deren Spektrum bei −1 liegt (§3); (3) Kneser-Koordinaten tragen (Punkt 1); (4) Potentialtheorie: Dichotomie in Supremum-
   normen auf Kreisen (§2.6), Zwei-Konstanten-Schranke ist das richtige Werkzeug für E_N; (5) Analogie Miller-Algorithmus
   (minimale Lösung einer Rekursion), aber mit nur polynomialer Trennung.
6. **Rückführungskette.** Vermutung ⟸ (LP-B) + (Q), wobei
   (LP-B) [linear]: |Π_N f_K(z)| ≤ B(|z−L||z−L̄|)^{−a} auf D für ein a < α − 1 = 3.447 — nach **Satz B-lin (§1.7, bewiesen
   bis auf (G2))** folgt daraus Π_N f_K → f_K lokal gleichmässig, ganz ohne Korridor-Bedingung;
   (Q) [quadratisch]: p_N − Π_N f_K = −J_N⁻¹T_N Q(p_N − f_K) → 0 lokal gleichmässig (numerisch: 1.6e−22 bei Fehler 3.3e−11).
   (LP-B) ist numerisch mit a = 0 erfüllt. Ihr Kern ist ein **inneres Problem** bei L, L̄ (Skala N^{−0.883}, Ecke von G',
   = Strang Ds Eckschicht von J_N⁻¹ = Sprung des Toeplitz-Symbols), §4.4–4.6; aussen regelt Fourier-Lemma + Zwei-Konstanten.
7. **Auswahlregel (numerisch, Vermutung B1/B2'):** Die Konkurrenz-Amplitude, die eine Singularität (1 − z/L)^s erzeugt,
   ist ≍ N^{1 + s_c − s} (s_c = Exponent der Konkurrenzrichtung, = 1 − γ); bestätigt für 5 Exponenten s und für deformierte
   Probleme (λ = ½, 1, 2). Mit Knesers s = 1 + γ und D₁ (Strang A) ergibt das c_N ≍ N^{1−2γ̄}. Das „+1“ (statt des naiven
   Koeffizientenquotienten) ist strukturell (Abschneidung), ein affines Modell bei L reproduziert es nicht.

## 1. Reformulierungen mit Beweisen

### 1.1 Zwei-Stellen-Jet-Form und Dreiecksstruktur — **bewiesen**

Sei c_i := p^{(i)}(a_0)/i! (Taylor-Koeffizienten von p an a_0; p(a_0 + y) = Σ_{i≤N} c_i y^i exakt), S = S(a_0) mit
S_{ij} = C(j,i) a_0^{j−i} (j ≥ i; obere Dreiecksmatrix mit Einsen auf der Diagonale), also c = S a, und
H_{ki} := [x^k] h(x)^i (untere Dreiecksmatrix, H_kk = a_1^k). Wegen p∘p = Σ_i c_i h^i und h^i = O(x^i) gilt für k ≤ N

    F(a) := ([x^k](p∘p) − 1/k!)_{k≤N} = H(a) S(a_0) a − e,   e = (1/k!)_k.

Also F = 0 ⇔ c = H⁻¹e = N-Jet von exp∘h⁻¹ an 0 (Substitution x = h⁻¹(y)) — die Ersatzform aus dem README.
Ableitung: DF[δ] = T_N[δ∘p + (p'∘p)δ] = H S δ + M δ, M = untere Toeplitz-Matrix der Reihe T_N(p'∘p). Mit Γ := H⁻¹M:

    J_N = H (S + Γ),   Γ = N-Jet an a_0 des Operators C δ := p'·(δ∘p⁻¹)   (untere Dreiecksmatrix),

denn H⁻¹ bildet die Reihe g(x) auf die Koeffizienten von g(h⁻¹(y)) ab und (p'∘p·δ)∘h⁻¹ = p'(a_0+y)·δ(p⁻¹(a_0+y)).
Diagonale: Γ_ii = p'(a_0)·a_1^{−i} = a_1^{−i−1} (aus c_1 a_1 = 1). Folgerung: det J_N = a_1^{N(N+1)/2} det(S + Γ).
Die Fixpunktform „a = S(a_0)⁻¹H(a)⁻¹e =: Ψ_N(a)“ ist **exakt** Ψ_N(a) = a − (HS)⁻¹F(a) (für alle a, nicht nur an
Lösungen), ihre Linearisierung an einer Lösung ist DΨ_N = −S⁻¹Γ. (Bewertung in §3.)

### 1.2 J-Projektion — **bewiesen**

Für eine feste Funktion f (holomorph nahe 0 und nahe f(0)) mit J_N^f invertierbar sei

    Π_N u := (J_N^f)⁻¹ T_N J_∞^f u        (u holomorph nahe 0 und nahe f(0)).

(a) Π_N P = P für P ∈ P_N. (b) Π_N u = 0, sobald J_∞^f u = O(x^{N+1}); insbesondere Π_N(ker J_∞^f) = 0.
(c) Π_N u hängt nur von den N-Jets von u an 0 und an f(0) ab: Π_N u = T_N u + (J_N)⁻¹ S_{≤N,>N} u_{>N}
    (Schwanz von u, an a_0 neu entwickelt, siehe §1.1).
(d) Sekanten-Identität: Ist f = f_K exakt und p ∈ P_N, so ist p∘p − exp = 𝒥(p − f) mit 𝒥u := u∘p + M_{p,f}u,
    M_{p,f} := ∫_0^1 f'(f + s(p − f)) ds; die abgeschnittene Gleichung ist also genau p = Π_N^{𝒥} f (Projektion zum
    Sekantenoperator). Linearisiert: p_N ≈ Π_N^{f_K} f_K.
Beweis: (a) J_N P = T_N J_∞ P. (b) T_N J_∞ u = 0. (c) J_∞ u = J_∞ T_N u + J_∞(Schwanz); T_N[(p'∘p)·Schwanz] = 0.
(d) f∘p − f∘f = M_{p,f}(p − f) (Hauptsatz längs der Strecke), f∘f = exp. ∎
Bedeutung: Π_N ist das lineare Auswahlprinzip. Es reproduziert Polynome und annulliert *alle* Störrichtungen der anderen
halben e-Funktionen exakt; die Frage „warum Kneser“ wird zur Frage, wie gut Π_N Knesers Funktion reproduziert.

### 1.3 Kern von J_∞ an f_K — **bewiesen** (bis auf die numerisch geprüfte Form von Ω, wie im Eindeutigkeitssatz)

Für θ holomorph auf Ω = A(D) sei δ_θ := T'(A + ½)·(θ∘A). Dann (Ableitung der Familie T(A + ½ + εθ(A)) nach ε)

    J_∞ δ_θ = T'(A + 1)·(θ + θ(· + ½))∘A,   T'(A(z) + 1) = e^z / A'(z) ≠ 0.

Jedes δ ∈ H(D) ist ein δ_θ (θ := (δ/T'(A+½))∘T; T'(w+½) ≠ 0 auf Ω, denn f_K' = T'(A+½)·A' hat in D keine Nullstelle:
**numerisch** Windungszahl 0 auf |z| = 0.9, 0.99, 0.999·|L|, min|f_K'| ≈ 0.346). Also: δ ∈ ker J_∞ ⇔ θ(w) + θ(w + ½) = 0 auf
B = Ω ∩ (Ω − ½). Wie in Schritt 3 des Eindeutigkeitssatzes (B waagrecht konvex, Breite > 1) setzt sich θ zu einer ganzen,
½-antiperiodischen Funktion fort: θ = Σ_{n ungerade} t_n e^{2πinw}.

**Satz (lineare Eindeutigkeit).** Sei 0 ≤ a < α − 1 = 3.447 (α = Re γ). Ist δ ∈ ker J_∞ ∩ H(D) mit
|δ(z)| ≤ K·|z − L|^{−a}|z − L̄|^{−a} auf D, so ist δ = 0. Für α − 1 ≤ a < 3α − 1 (= 12.34) bleibt höchstens span{δ_{+1}, δ_{−1}}
(δ_n := δ_θ mit θ = e^{2πinw}; δ_{−1} wächst bei L wie |z − L|^{1−α}, δ_{+1} ebenso bei L̄).
*Beweis (Skizze mit allen Zutaten):* Im Korridor B gilt nach der Asymptotik A(z) = log(z − L)/L + c + o(1)
(kn.beob.slog_scheibe): |z − L| ≍ e^{−κ y} mit κ = Im L + Re L·(Re L/Im L) = |L|²/Im L = 1.413 (die Randgeraden von Ω
haben Steigung −Re L/Im L), und |T'(w + ½)| ≍ |f(z) − L| ≍ |z − L|. Also |θ(w)| ≤ K' e^{(1+a)κ|y|} auf B; wegen der
Periodizität auf allen Höhen. Fourier-Koeffizient: t_n e^{−2πny} = ∫_0^1 θ(x + iy)e^{−2πinx}dx, also
|t_n| ≤ K' e^{((1+a)κ − 2π|n|)|y|} mit y → ∓∞ passend zum Vorzeichen von n ⇒ t_n = 0, sobald 2π|n| > (1 + a)κ, d. h.
|n| > (1+a)κ/2π = (1 + a)/α. Für a < α − 1 sind alle ungeraden n verboten. ∎
(Zur Asymptotik: nach Knesers Konstruktion ist A = u + ϑ(u), u = log Ψ(z)/L, Ψ Schröder-Funktion an L, ϑ 1-periodisch,
holomorph und beschränkt für Im u ≥ c — Strang A hat die Entwicklung explizit, κ_0 = 0.4884 + 0.9223i, |κ_1| = 2.8e−4. Daraus
T(w) − L ≍ e^{L(w − κ_0)} und T'(w) ≍ L(T(w) − L) für Im w ≥ c, gleichmässig auf waagrechten Strecken beschränkter Länge;
das liefert |z − L| ≍ e^{−κy} auf B. Für die Fourier-Abschätzung braucht man ausserdem, dass B auf jeder Höhe eine waagrechte
Strecke der Länge ≥ 1 enthält: das ist (G2) aus kn.beob.slog_scheibe, numerisch geprüft — derselbe Status wie im
Eindeutigkeitssatz.) Das ist das lineare Gegenstück zum Eindeutigkeitssatz und der Grund, warum Strang C in
l¹ mit Gewichten |L|^k(1+k)^s (Kneser liegt darin für s < 5.45, δ_{±1} nicht) gleichmässige Stabilität sieht.

### 1.4 Kneser-Koordinaten für p_N — **bewiesen**

Auf V_N := {z ∈ D : p_N(z) ∈ D} sei Θ_N := A∘p_N − A − ½ und θ_N := Θ_N∘T. Dann gilt (wo alle Terme definiert sind)

    θ_N(w) + θ_N(w + ½ + θ_N(w)) = ε_N(w),   ε_N(w) := A(e^{T(w)} + E_N(T(w))) − w − 1,

denn A(p∘p(z)) = A(p(z)) + ½ + Θ_N(p(z)) und A(p(z)) = A(z) + ½ + Θ_N(z), A(e^z) = A(z) + 1. ε_N hat bei w = −1 eine
Nullstelle der Ordnung N+1. Mit φ = id + π und g = A∘p∘T = φ⁻¹∘(φ + ½) wird daraus π(w + 1 + ε_N(w)) − π(w) = −ε_N(w):
**p_N ist genau dann eine exakte Lösung, wenn π periodisch ist; die Abschneidung macht π „periodisch bis auf ε_N“.**

### 1.5 Zwei-Konstanten-Schranke für den Abschneidefehler — **bewiesen**

Sei V ∋ 0 ein Gebiet mit Green-Funktion g_V(·,0), E holomorph auf V, |E| ≤ M, E = O(x^{N+1}). Dann
|E(x)| ≤ M·exp(−(N+1) g_V(x,0)) auf V. (u = log|E| − log M + (N+1)g_V ist subharmonisch auf V∖{0}, nach oben
beschränkt bei 0, ≤ 0 am Rand.) Anwendung: Ist |p_N| ≤ B auf D und V ⊂ {z ∈ D : p_N(z) ∈ D'} (D' ⋐ D), so ist E_N
dort ≤ M e^{−(N+1)g_V}. Der Korridor V_0 ist bei L lokal ein Keil mit Öffnung π − arg μ = 141.7° (in w: Streifen der
asymptotischen Breite π/Im L − ½ = 1.85), also g_{V_0}(x,0) ≍ |x − L|^{ν}, ν = 180/141.7 = 1.27 (gleich π/(1.85·κ)).
**Die abgeschnittene Gleichung gilt super-genau bis zum Abstand ~N^{−1/ν} = N^{−0.79} von L, d. h. bis zur Höhe
y ≈ log N/(νκ) = 0.56·log N, und nicht weiter.** (κ = |L|²/Im L = 1.413, |z − L| ≍ e^{−κy}.)

### 1.6 Quantitativer Bôcher (Fourier-Lemma) — **bewiesen**

Sei θ holomorph auf einer Umgebung von R = [a, a+2] × [−H, H] mit |θ(w + 1) − θ(w)| ≤ η für Re w ∈ [a, a+1] und
|θ| ≤ M auf den Kanten Im w = ±H. Für t_n(y) := ∫_a^{a+1} θ(x + iy) e^{−2πin(x+iy)} dx gilt
t_n'(y) = i e^{−2πin(a+iy)} (θ(a + 1 + iy) − θ(a + iy)), also für n ≠ 0

    |t_n(0)| ≤ M e^{−2π|n|H} + η/(2π|n|).

(Für exakt periodisches θ ist t_n(y) konstant — das ist Bôcher/Liouville in endlicher Höhe.) Zusammen mit §1.5:
Hält man θ_N bis zur Höhe H_N ≈ (1 − ε)·0.56·log N beschränkt (|θ_N| ≤ M), so ist die Amplitude des Störmodus n
≤ M·N^{−(1−ε)·3.50|n|} + (super-kleiner Rest) (2π·0.56 = 3.50). Das erklärt, warum p_N **ein-modig** ist (Modus ±3 um
N^{−7} gegenüber ±1 unterdrückt) und gibt bedingt die Innen-Rate N^{−3.5}; die wahre Rate N^{−7.894} ist besser (§2.2).

### 1.7 Lineare Rückführung: eine gewichtete Schranke genügt — **bewiesen** (bis auf die Ω-Asymptotik aus §1.3)

**Satz B-lin.** Sei P_N := Π_N f_K (lineares Teilproblem (LP), §4.1). Gibt es B < ∞ und a < α − 1 = 3.447 mit

    |P_N(z)| ≤ B·(|z − L| |z − L̄|)^{−a}   für alle z ∈ D und alle N ≥ N_0,

so gilt P_N → f_K lokal gleichmässig auf D (und damit a_k(P_N) → a_k^K für jedes k).
*Beweis.* (1) Montel: jede Teilfolge hat eine auf D lokal gleichmässig konvergente Teilfolge P_N → P_∞, und P_∞ erfüllt
dieselbe gewichtete Schranke. (2) Sei V := D_{0.7}; f_K(V̄) ⊂ D kompakt (f_K(0.7) ≈ 1.24 < |L|). Für u_N := f_K − P_N ist
E_N := J_∞u_N = u_N∘f_K + (f_K'∘f_K)u_N auf V durch eine von N unabhängige Konstante M beschränkt (Schranke für P_N auf
dem Kompaktum f_K(V̄) ∪ V̄ ⊂ D) und hat nach Definition von Π_N eine Nullstelle der Ordnung N+1 bei 0. Zwei-Konstanten (§1.5):
|E_N| ≤ M e^{−(N+1)g_V(·,0)} → 0 lokal gleichmässig auf V. Also J_∞(f_K − P_∞) = 0 auf V. (3) δ := f_K − P_∞ ∈ H(D) liegt
im Kern von J_∞ und wächst höchstens wie (|z − L||z − L̄|)^{−a}, a < α − 1. Nach dem Satz in §1.3 ist δ = 0. (4) Alle
Teilfolgen-Grenzwerte sind f_K. ∎
Bemerkung: Anders als im nichtlinearen Fall (kn.satz.polynom_halbexp_reduktion) braucht man **keine Korridor-Bedingung**,
und es darf sogar ein Wachstum bis |z − L|^{−3.44} zugelassen werden. Numerisch gilt die Schranke mit a = 0 (B ≈ 2.24).
Für das volle Problem fehlt dann „nur“ die quadratische Kontrolle p_N − Π_N f_K = −J_N⁻¹T_N Q(p_N − f_K) (§1.2, §4.1).

## 2. Numerische Struktur

### 2.1 Polynomiale Raten (b08_exponent.py, Kneser bis k = 1000 aus gemeinsam/)

| N | 20 | 40 | 80 | 120 | 160 | 194 |
|---|---|---|---|---|---|---|
| max_{|x|=1}|p_N − f_K| | 1.4e−8 | 5.0e−11 | 2.2e−13 | 8.8e−15 | 9.0e−16 | 2.0e−16 |
| |p_N(L) − L| | 1.5e−6 | 4.5e−8 | 1.6e−9 | 2.3e−10 | 5.6e−11 | 2.2e−11 |

Fit N ≥ 80: Exponenten −7.905 (Kreis) und −4.848 (bei L); lokale Exponenten pendeln (Periode ≈ 40 in N).
Profil d_N = p_N − q_N: |d_k||L|^k wächst in k wie ~k^{2.7} (Form einer δ_{±1}-Richtung, b01).

### 2.2 Ein-Modus-Struktur und Amplitudengesetz (b09_theta_fourier.py, b10_amplitude.py)

θ_N auf waagrechten Strecken der Länge 1 im Korridor (y = 0, 0.4, 0.8), 64 Punkte; A aus 300 slog-Koeffizienten,
Kneser-Kontrolle max|θ_K| = 2e−41. Für N = 194, y = 0: |ĉ_{±1}| = 1.73e−19, |ĉ_0| = 1.6e−42, |ĉ_{±2}| = 9e−38,
|ĉ_{±3}| = 7e−34. ĉ_1 ist auf y = 0, 0.4, 0.8 gleich (Periodizität). Amplitude c_N := ĉ_1:

| N | 18 | 30 | 50 | 70 | 102 | 150 | 194 |
|---|---|---|---|---|---|---|---|
| |c_N| | 3.1e−11 | 4.5e−13 | 7.8e−15 | 5.8e−16 | 2.8e−17 | 1.3e−18 | 1.7e−19 |
| c_N/N^{1−2γ̄} | .220−.120i | .185−.091i | .171−.108i | .196−.087i | .182−.088i | .186−.078i | .184−.072i |

**Vermutung B1:** c_N = C N^{1−2γ̄}(1 + o(1)), C ≈ 0.18 − 0.08i. (Naiver „Schwanzquotient“ Kneser/Konkurrenz gäbe
N^{−2γ̄}; der zusätzliche Faktor N ist strukturell, §2.3/b21, aber nicht hergeleitet.)

### 2.3 Lineares Modell der Auswahl (b11_projektion_moden.py, b11b_vergleich.py)

Π_N (an f_K, Koeffizienten aus gemeinsam/) auf u_s := (1 − z/L)^s + (1 − z/L̄)^{s̄}; gemessen der Modus e^{−2πiw} von
θ_{Π_N u − u}. Fit-Exponenten (N = 30 … 190): s = 2.5: −5.32; s = 1 + γ: −8.10; s = 5.5: −8.17, also ≈ −s − 2.7
(oszillierend). Der Quotient (nichtlineares c̄_N)/(linear, s = 1 + γ) ist fast konstant: −0.018−0.010i (N=30),
−0.027−0.012i (N=70), −0.029−0.011i (N=130), −0.030−0.011i (N=190). **Die nichtlineare Auswahl ist also im
Wesentlichen die lineare Projektion von Knesers Leitsingularität bei L** (Faktor → D₁, Knesers Singularitätskoeffizient, s. u.).

**Generische Exponentenregel (b17_exponent_s.py, komplexer Fit N = 60 … 180):** für u_s = (1 − z/L)^s + konj. ist die
Modusamplitude c_{−1}(u_s) ~ N^{a(s)} mit a(1+γ) = −8.07−2.15i, a(5.447) = −8.06−1.13i, a(5.947+1.058i) = −8.49−2.13i,
a(4.947+1.058i) = −7.56−2.09i, a(5.447+2.116i) = −8.03−3.11i; Vorhersage 2 − γ − s: −7.89−2.12i, −7.89−1.06i, −8.39−2.12i,
−7.39−2.12i, −7.89−3.17i. Imaginärteile stimmen, Realteile liegen 0.1–0.18 darunter (Vorlauf). **Vermutung B2:**
c_{−1}(u_s) ≍ N^{2−γ−s}, d. h. N mal Quotient der N-ten Taylorkoeffizienten von u_s und δ_{−1}.
Normierung (b19_universell.py): Q(s,N) := c_{−1}(u_s)/(N^{2−γ}L^N[z^N](1−z/L)^s) wächst ∝ N; Q/N ist für den Kern-Exponenten
s = 1+γ bemerkenswert konstant: 0.00249+0.00452i (N=60), 0.00274+0.00424i (90), 0.00284+0.00460i (120), 0.00258+0.00472i
(150), 0.00249+0.00456i (180); für andere s hängt Q/N von s ab und driftet (z. B. s = 2.5: 0.0038+0.0017i → 0.0026+0.0012i).
Also **keine** universelle Regel „nur der N-te Koeffizient zählt“; der Kern-Exponent (Knesers eigener) ist der saubere Fall.
**Der Zusatzfaktor N ist strukturell (b21_lambda_deformation.py).** Deformiertes Problem J'_λδ = δ∘f + λ(f'∘f)δ: Kern
θ(w+½) = −λθ(w), Konkurrenz-Exponent bei L s_c(λ) = 1 − γ + 2 log λ/L. Für u = u_{5.447}, komplexer Fit N = 60 … 180:
λ = 1: −8.058−1.111i (Regel 1 + s_c − s: −7.894−1.058i); λ = 2: −7.836−1.999i (Regel −7.661−2.039i, ohne „+1“: −8.661);
λ = ½: −8.285−0.234i (Regel −8.127−0.077i, ohne „+1“: −9.127). Überall Regel minus ≈ 0.16 (gemeinsamer Vorlauf oder
log-Korrektur). **Vermutung B2':** Amplitude ≍ N^{1 + s_c − s} — das „+1“ gehört zur Abschneidung, nicht zum Kern-Exponenten.
Ein affines Modell f(z) = L + μ(z − L) (b20_affinmodell.py) reproduziert das NICHT (keine saubere Potenz) — der
Mechanismus braucht die globale Zwei-Kreis-Struktur, nicht nur die Linearisierung bei L.
Mit Strang A (20:08): Knesers L-Singulärteil ist D₁(1 − z/L)^{1+γ}, D₁ = −0.03200 − 0.01124i; mein Quotient
nichtlinear/linear driftet von −0.018−0.010i (N=30) auf −0.030−0.011i (N=190) gegen D₁ — die nichtlineare Amplitude ist
D₁ × (lineare Antwort der reinen Singularität) bis auf Terme höherer Ordnung in der Singularitätsentwicklung.

### 2.4 Konvergenzgebiet (b02, b03, b06, b07; vgl. Strang C)

p_N → f_K bei x = 1.45, 1.6 (|p_194(1.6) − f_K(1.6)| = 3.9e−19; beste Taylorteilsumme dort 1.2e−6), bei −1.3, −1.35,
1.3e^{i·135°} usw. (in D); keine Konvergenz bei 1.5i, −1.5, 2.0. Passt zu G' = D(0,|L|) ∪ D(a_0,|L−a_0|).
|a_N^{(N)}| ≈ 1.59^{−N} (b05), Nullstellen auf einer langsam schrumpfenden Kurve (b05_nullstellen.png).
Ecke von G' bei L: Tangentenwinkel 76.62° (Kreis um 0) und 97.68° (Kreis um a_0) ⇒ Innenwinkel 180° + 21.07°.

### 2.5 Lebesgue-Konstanten der J-Projektion (b04_projektion_lebesgue.py)

Λ_N(r,R) = max_{|z|=r}(R/2π)∫|Π_N[1/(t−·)](z)| dφ, t = Re^{iφ} (Π_N an p_N): wächst exponentiell für R ≤ 1.45
(R = 1.35, r = 1: 5.9, 1.5e3, 1.1e6, 1.3e9 für N = 20, 40, 60, 80), ≈ 1.09 für R = 1.8. Der grosse Anteil hat immer
dieselbe Form (Quotient r = 0.5 zu r = 1 konstant ≈ 63): eine einzige fast-singuläre Richtung, wie in §2.2.

### 2.6 Dichotomie in Kreis-Normen (Skizze)

Sup-Normen auf einem Gebiet G für Urbild und Bild: Liegt G̅ ⊂ D (weg von L), so sind die Kernelemente δ_n auf G
beschränkt und polynomial gut approximierbar ⇒ exponentiell kleine Singulärwerte (echter Fast-Kern). Reicht G bis L,
so muss das Bild die Komposition δ∘p auf p(G) sehen, und p(D) reicht bis 2.23 > |L| ⇒ Vorwärtsoperator und F''
wachsen exponentiell (Strang C: ω ~ 1.25^N; meine Rechnung: Konsistenzfehler ~ (|h|_{|L|}/|L|)^N = 1.27^N). Keine
Wahl von Kreisradien vermeidet beides; der Ausweg muss die Struktur der konkreten Vektoren benutzen (Kegel, Strang C)
oder die Geometrie G' mit einspringender Ecke.

## 3. Bewertung der vorgeschlagenen Ideen

1. **Lyapunov–Schmidt** — trägt. Numerisch ist der relevante Fast-Kern genau ein komplexer Parameter (Modus ±1);
   die reduzierte Gleichung bestimmt c_N. Offen: gleichmässige Gutgestelltheit des Komplements (Strang C sieht sie in
   gewichteten l¹-Normen) und die Asymptotik der reduzierten Gleichung (Vermutung B1).
2. **Jet-Fixpunktform** — verworfen als Kontraktion. Ψ_N(a) = a − (HS)⁻¹F(a) exakt; DΨ_N = −S⁻¹Γ; da J_N = HS(I + S⁻¹Γ)
   fast singulär ist, hat S⁻¹Γ einen Eigenwert nahe −1; im Grenzfall ist DΨ_∞ in Kneser-Koordinaten θ ↦ −θ(· − ½), also
   Eigenwert 1 auf allen antiperiodischen Moden und Betrag 1 auf allen anderen Fourier-Moden. Keine Monotonie gefunden.
   Gemessen an p_N (b12_spektrum_psi.py): Eigenwert von DΨ_N nächst −1: −0.9926±0.028i (N=10), −1.0003 (20), −1.00005 (30),
   −1.00002 (40); Spektralradius 5.9e2, 1.7e6, 5.5e9, 1.9e13 (N = 10 … 40) — die Iteration divergiert sogar.
   Nützlich bleibt die Dreiecksstruktur (§1.1).
3. **Kneser-Koordinaten** — trägt am besten (§1.4, §2.2): fast-lineare Gleichung für die Fourier-Koeffizienten mit
   kleiner Kopplung (Moden ≠ ±1 sind ~c_N² bzw. N^{−2α}-unterdrückt).
4. **Potentialtheorie** — als a-priori-Argument per Widerspruch nicht ausreichend (es gibt andere Lösungszweige; jede
   Schranke muss den Zweig benutzen). Richtige Werkzeuge: Zwei-Konstanten-Satz für E_N (§1.5) und Ecken-Approximation
   auf G' (§2.4).
5. **Eigene:** Miller-Analogie (minimale Lösung, Konkurrenz nur polynomial getrennt); Kern-Satz §1.3.

## 4. Rückführung

### 4.1 Die Vermutung ist im Kern linear (numerisch, b14_linear_fk.py)

Lineares Teilproblem **(LP)**: P_N := Π_N f_K = J_N(f_K)⁻¹ T_N[eˣ + (f_K'∘f_K)·f_K] (eindeutiges P ∈ P_N mit
T_N[J_∞(f_K − P)] = 0). Exakt gilt p_N − Π_N f_K = −J_N⁻¹ T_N Q(p_N − f_K) mit dem quadratischen Rest
Q(u) = f∘(f+u) − f∘f − (f'∘f)u + u∘(f+u) − u∘f (§1.2). Gemessen auf |z| = |L| (inkl. z = L):

| N | 30 | 60 | 90 | 120 | 150 | 180 |
|---|---|---|---|---|---|---|
| max|Π_N f_K − f_K| | 2.08e−7 | 7.18e−9 | 9.40e−10 | 2.37e−10 | 8.11e−11 | 3.27e−11 |
| max|Π_N f_K − p_N| | 1.2e−15 | 3.4e−18 | 6.8e−20 | 6.2e−21 | 8.6e−22 | 1.6e−22 |

Also: **Vermutung ⇔ (LP) + quadratische Kontrolle**, und die quadratische Kontrolle ist numerisch harmlos
(Strang C, 19:42: ω_eff ~ N^{1.9} längs des echten Zuwachses). Ausserdem nahe L ganz zahm: |p_N'(L) − μ| = 1.2e−6,
7.5e−8, 2.2e−9 (N = 40, 80, 194), p_N''(L) beschränkt, N|p_N(L) − L| → 0. Die „Grenzschicht“ bei L ist also eine kleine
Störung von f_K, kein wildes Verhalten.

### 4.2 Warum einfache Normargumente für (LP) nicht reichen (numerisch + exakte Struktur, b15_konsistenz.py)

Exakt: Π_N f_K − q_N = J_N⁻¹ r_N, r_N = T_N[Schwanz_N(f_K)∘f_K] = H·j_{a_0}(Schwanz), und J_N = H(S + Γ), also
Π_N f_K − q_N = (S + Γ)⁻¹ j_{a_0}(Schwanz_N f_K). In l¹_{|L|,s} wächst r_N ≈ 1.47^N (343 → 4e17 für N = 60 → 150), J⁻¹
dämpft um bis 1e−28. In der natürlichen a_0-Jet-Norm (Radius ρ_1 um a_0) ist der Konsistenzfehler klein nur für
ρ_1 ≤ |L| − a_0 = 0.876 — dort hat (S + Γ) aber echte Fast-Kern-Richtungen (δ_n sind auf D(a_0, 0.876) beschränkt).
Ebenso zeigen die Auswahlgewichte (b16): c_N = Σ_{k>N} a_k Λ_N(k) exakt, aber |Λ_N(k)L^k| wächst bis 1e76 (N = 120):
c_N ~ 1e−18 entsteht aus massiver Auslöschung. **Darstellungen über Knesers Taylor-Schwanz bei 0 sind die falsche Basis**:
die natürliche Geometrie ist G' (zwei Kreise, beide mit Rand durch L, L̄), nicht D.

### 4.3 Abdeckung von D und Jet-Transfer (bewiesen bzw. numerisch)

(a) **numerisch** (Gitter 2518 Punkte, |x| ≤ 0.995|L|): D ⊂ V_0 ∪ D(a_0, |L − a_0|), V_0 = {x ∈ D : |f_K(x)| < |L|}.
Beide Stücke berühren L, L̄ — nur dort fehlt Kontrolle.
(b) **Lemma (Jet-Transfer, bewiesen).** Ist p ∈ P_N Lösung, U ∋ 0 ein Gebiet, auf dem p univalent ist, D(a_0, r) ⊂ p(U)
und |exp∘(p|_U)⁻¹| ≤ M auf D(a_0, r), so gilt |p(z)| ≤ M r/(r − |z − a_0|) und
|p(z) − exp(p⁻¹(z))| ≤ M (ρ/r)^{N+1}/(1 − ρ/r) für |z − a_0| ≤ ρ < r.
Beweis: Nach §1.1 ist p(a_0 + y) das Taylor-Polynom vom Grad N von g := exp∘p⁻¹ an a_0 (beides Polynome bzw. Jets vom
Grad N, die bis zur Ordnung N übereinstimmen); Cauchy-Abschätzung |g_i| ≤ M r^{−i}. ∎
Bedeutung: Kontrolle von p_N auf dem Korridor (dort lebt p_N⁻¹(D(a_0,r))) überträgt sich automatisch auf die Kreisscheibe
um a_0 — mit geometrischer Genauigkeit p_N ≈ exp∘p_N⁻¹. Die Abschneidung „wirkt“ nur noch bei L, L̄.

### 4.4 Beweisprogramm (Skizze; was fehlt, ist benannt)

Schritt 1 (aussen, Kneser-Koordinaten): Hat man für z ∈ V_0 mit |z − L|, |z − L̄| ≥ C/N die Schranke p_N(z) ∈ D,
|Θ_N(z)| ≤ M („Korridor-Kontrolle bis Skala 1/N“, Hypothese (K)), dann gilt nach §1.5 (E_N super-klein bis zum Abstand
N^{−0.79}, Höhe H_N ≈ 0.56 log N) und §1.6: θ_N = c_N e^{2πiw} + c̄_N e^{−2πiw} + O(|c_N|² + N^{−∞}) im Inneren und
|c_N| ≤ M' N^{−3.5}. — **bewiesen bis auf Buchhaltung**, *bedingt* auf (K). (Wahr ist N^{−7.9}.)
Schritt 2 (Scheibe um a_0): Jet-Transfer §4.3(b) gibt Schranken auf D(a_0, r) aus Schritt 1. — **bewiesen**, bedingt.
Schritt 3 (Grenzschicht): Unter (K) ist p_N auf D ∖ (C/N-Umgebungen von L, L̄) beschränkt (V_0-Teil direkt, Rest über
Jet-Transfer — nahe L mit Vorsicht, dort laufen ∂D und ∂D(a_0,·) unter 21° zusammen); die C/N-Bisse kosten nach
Bernstein–Walsh nur einen Faktor e^{O(C)} (Green-Funktion des Komplements eines Bisses der Grösse ρ ist O(ρ), N·ρ = O(1);
Bisse der Grösse N^{−0.79} wären zu teuer: e^{N^{0.21}}). ⇒ sup_D|p_N| ≤ B und die Korridor-Bedingung ⇒ **Konvergenz nach
kn.satz.polynom_halbexp_reduktion**. — Skizze.
**Offen bleibt genau (K)** — eine Aussage über p_N auf der Skala 1/N bei L, L̄. Gleichwertige lineare Form: die
J-Projektion e_N := Π_N f_K − f_K ist auf V_0 ∩ {|z − L| ≥ C/N} durch o(|z − L|) beschränkt (Grenzschichtproblem).

### 4.5 Vermutete Form der Lösung (numerisch gestützt)

e_N = p_N − f_K ≈ c̄_N δ_{−1} + c_N δ_{+1} (aussen) + Grenzschicht bei L, L̄; c_N = C N^{1−2γ̄}; |e_N(L)| ≈ 3·N^{−4.83}
(b13: N^{4.876}|e_N(L)| = 2.94, 3.07, 3.11, 3.12, 3.16 für N = 40 … 194, also Exponent ≈ −4.83; Ecken-Vorhersage −αβ = −4.81).
Übergang zwischen beiden Formen bei |z − L| ≈ N^{−0.88} (Gleichsetzen N^{−7.894}r^{−3.447} = N^{−4.83}) — das ist genau die
Auflösungsskala N^{−α} = N^{−0.883} der Ecke. Lokale Abklingexponenten bei festem Abstand r von L
(radial): r = 0.2: −7.9, 0.05: −7.1, 0.01: −5.4, 0.002: −5.0 (b13).

### 4.6 Drei Bilder desselben Kerns (Abgleich mit Strang D)

Strang D (20:10/20:15): G = lim J_N⁻¹ ist im Inneren Toeplitz mit Sprungsymbol σ = (1/(f'∘f))·1_{Γin} auf |x| = |L|,
Sprung genau bei arg L; die Abweichung J_N⁻¹ − P_N G P_N sitzt in einer Eckschicht (Spalten j ≈ N). Mein Befund b18: der
Konsistenzrest r_N ist als Polynom auf dem GANZEN Kreis riesig (N = 120: 7e12 bei φ = 0°, 6e11 bei 180°), also kann
T(σ) allein ihn nicht klein machen — die Auslöschung J_N⁻¹ r_N ≈ 1e−10 passiert in der Eckschicht. Damit sind
(i) Ds Eckschicht von J_N⁻¹ (Indizes k, j ≈ N), (ii) meine Grenzschicht bei L auf der Skala ~1/N (§4.4 (K)) und
(iii) der Sprung des Symbols bei arg L **drei Bilder desselben inneren Problems**. Lokal bei L ist J_∞ auf
(z − L)^s die Multiplikation mit dem Mellin-Symbol μ^s + μ = 2μ^{(s+1)/2} cosh((s − 1)L/4), Nullstellen genau an den
Kern-Exponenten s = 1 + (2n+1)γ; cosh hat die Γ-Faktorisierung cosh(πζ) = π/(Γ(½+iζ)Γ(½−iζ)), ζ = (s−1)L/(4π) —
der natürliche Kandidat für eine Wiener–Hopf-Lösung des inneren Problems (Trennlinie Re s ∈ (−3.447, 5.447)).

## 5. Nächste Schritte

1. Inneres Problem als Grenzproblem N → ∞ formulieren: Mellin-Symbol μ^s + μ (Γ-Faktorisierung oben), Ecke von G' mit
   Aussenwinkel 0.883π, Abschneidung als Paley–Wiener-Bedingung; Ziel: Exponent 1 − 2γ̄ (insbesondere der Faktor N) und die
   Konstante Q/N ≈ 0.0025 + 0.0045i aus b19 vorhersagen.
2. Für (LP): Satz B-lin (§1.7) reduziert alles auf eine gewichtete Schranke für Π_N f_K. Diese aus Ds Toeplitz-Bild + Ecken-
   analyse gewinnen; dabei muss die Auslöschung in der Eckschicht explizit werden (Strang E: sup|J⁻¹| allein reicht nie).
3. Nichtlinear: p_N − Π_N f_K = −J_N⁻¹T_N Q(p_N − f_K) — numerisch 2. Ordnung klein (§4.1); Kegel-Kontrolle (Strang C).
4. Vorhersagen für Strang E prüfen lassen (Basen 2, 10, 1.5; dritte Iterierte), um den Mechanismus zu testen.
5. Werkzeug für die Teilchen-Lösung von θ(w) + θ(w + ½) = h mit Schranken (für (Q) und Schritt 1): Loray–Touzet–Voronin,
   arXiv:2001.03033, Thm 9.5 (Hinweis Strang C).
6. Bemerkung zu Modellen (**bewiesen**, elementar): Eine Möbius-Abbildung hat an ihren zwei Fixpunkten die Multiplikatoren
   λ und 1/λ; bei reellen Koeffizienten und konjugierten Fixpunkten ist der zweite zugleich λ̄, also |λ| = 1 — **jedes reelle
   Möbius-Modell mit Fixpunkten L, L̄ ist elliptisch**. Es kann die abstossenden Fixpunkte von exp (|L| = 1.37) und damit die
   komplexen Kern-Exponenten 1 ± γ samt polynomialer Auswahl nicht nachbilden (dort ist γ_M = 2π/φ reell). Ebenso versagt
   das affine Modell f = L + μ(z − L) (b20). Ein treues Modell braucht abstossende konjugierte Fixpunkte und eine kanonische
   Lösung mit schwacher Kern-Singularität dort.

## 6. Dateien (alle in b_reformulierung/, Python = C:/Users/janis/miniconda3/python.exe, nur Daten aus gemeinsam/)

| Datei | Inhalt |
|---|---|
| b01_fehlerprofil.py | Koeffizientenprofil von p_N − q_N, |p_N − f_K| auf Kreisen, E_N auf Kreisen |
| b02_ueberkonvergenz.py, b06_ueberkonv_check.py, b07_links.py | Konvergenz ausserhalb D (ℝ⁺) bzw. links in D, Vergleich Taylorteilsummen |
| b03_konvergenzgebiet.py (+png, npz) | Karte log10|p_N − p_194| |
| b04_projektion_lebesgue.py | Lebesgue-Konstanten der J-Projektion Π_N |
| b05_top.py, b05_plot.py (+png, npy) | Spitzenkoeffizient a_N^(N), Nullstellen von p_N |
| b08_exponent.py | Abklingexponenten im Inneren und bei L (Kneser bis k = 1000) |
| b09_theta_fourier.py | θ_N in Kneser-Koordinaten, Fourier-Moden (Ein-Modus-Befund) |
| b10_amplitude.py (+json) | c_N für N = 10 … 194, Test c_N/N^{1−2γ̄} |
| b11_projektion_moden.py, b11b_vergleich.py | lineares Modell: Π_N auf Singularitäten (1 − z/L)^s; Quotient zu c_N |
| b12_spektrum_psi.py | Spektrum der Jet-Fixpunkt-Linearisierung DΨ_N = −S⁻¹Γ |
| b13_innen.py | Grenzschicht bei L |
| b14_linear_fk.py | (LP): Π_N f_K gegen p_N und f_K |
| b15_konsistenz.py | Konsistenzfehler r_N in l¹_{|L|,s} (wächst ~1.47^N) |
| b16_auswahlgewichte.py (+json) | Auswahlgewichte Λ_N(k) (massive Auslöschung) |
| b17_exponent_s.py, b19_universell.py | komplexer Exponent a(s), Normierung Q(s,N) |
| b18_lokalisierung.py | r_N auf Γin/Γout (keine Lokalisierung) |
| b20_affinmodell.py | affines Modell bei L (reproduziert die Regel nicht) |
| b21_lambda_deformation.py | λ-Deformation: „+1“ ist strukturell |
| brett_eintrag.py, brett1–6.txt | BRETT-Einträge |

---

# Runde 2 (ab 20:31): Voraussetzung von Satz B-lin — Rückführung auf eine Pol-Antwort-Schranke im Randband

## 7.0 Kurzfassung Runde 2

1. **Satz B-Max (bewiesen, §7.3):** Gilt die *polynomiale* Pol-Antwort-Schranke (H*) in einem Band der Breite ~N^{−0.883}
   um Ḡ' = D̄(0,|L|) ∪ D̄(a₀,|L−a₀|), so ist sup_{D̄}|Π_N f_K − f_K| ≤ C·N^{κ+1−4.69} → 0, also gilt die **lineare Vermutung
   (LP) gleichmässig auf D̄ mit Rate**. Neu gegenüber Cs Satz 9.5: der exponentielle Zerfall der Pol-Antwort weg vom Rand
   ((H_R1)-Schwanz, (H_R2)) muss NICHT vorausgesetzt werden — er folgt aus dem Maximumprinzip in der Pol-Variablen ζ
   (Lemma 7.2) mit der **explizit berechneten Green-Funktion** von ℂ∖Ḡ' (Lemma 7.3; Kapazität 1.5169; bei L: g ≈ 0.655·u^{1.1326}).
2. **(H*) numerisch:** dist(ζ,Ḡ')·max_{D̄}|Π_N k_ζ| ≤ 1 an allen 7 getesteten Randstellen, t = 0.5…8, N = 60…300 (b25, b25b;
   fallend in N);
   der Anteil der von Strang D gefundenen exponentiell wachsenden Eckspitze von J_N⁻¹ an der Pol-Antwort wächst nur ≈ N^{1.0}
   (N = 160…360, b30): die Spitze ist für kohärente Eingaben fast unsichtbar (Überlapp mit der f_K-Eingabe 8.7e−9 bei N = 300).
3. **Satz B-Stab (Implikation bewiesen, Voraussetzung FALSCH):** gleichmässige Stabilität |J_N⁻¹|_{ℓ²_{|L|,s}} ≤ S würde
   (LP) in drei Zeilen liefern (die rechte Seite T_N[eˣ + (f_K'∘f_K) f_K] ist trivial beschränkt). Aber der Eck-Singulärwert
   wächst (b28: 1.54 → 6.13 für N = 160 → 360, s = 0), in Übereinstimmung mit Strang D. Man braucht Kohärenz, nicht Operatornormen.
4. **Inneres Problem / Faktor N (Leitungsfrage 1):** nicht hergeleitet. Klar ist jetzt: Pol-Antworten und Randfehler leben auf der
   Eckskala N^{−α}, α = 0.883 (Aussenwinkel/π); die Modenamplitude c_N dagegen hat exakt die Phase von N^{1−2γ̄}, d. h. Skala 1/N —
   eine reine Eckskalen-Anpassung gäbe N^{−2αγ̄} = N^{−7.85−1.87i} und ist durch die Phase widerlegt. Der Faktor N gehört zur
   Taylor-Koeffizientenebene (Darboux-Skala 1/N), nicht zur Ecke (§7.6).

## 7.1 Satz B-Stab — eine korrekte Implikation mit falscher Voraussetzung

**Satz B-Stab (bewiesen bis (G2)).** Sei ½ < s < 4.95 und A_s := {Σc_kx^k : |c|_s² = Σ|c_k|²|L|^{2k}(1+k)^{2s} < ∞} (Banach-
Algebra für s > ½). Ist S := sup_N |J_N(f_K)⁻¹|_{A_s→A_s} < ∞, so gilt Π_N f_K → f_K gleichmässig auf D̄.
*Beweis.* J_∞f_K = f_K∘f_K + (f_K'∘f_K)f_K = eˣ + g f_K, g = eˣ/f_K'. f_K ∈ A_s (Darboux: a_k = O(k^{−6.447}|L|^{−k}), s < 5.95),
f_K' ∈ A_s (s < 4.95), 1/f_K' ∈ A_s (Wiener-Lemma der gewichteten Algebra; f_K' ≠ 0 auf D̄, numerisch §1.3), also v := eˣ + g f_K ∈ A_s.
|Π_N f_K|_s = |J_N⁻¹T_N v|_s ≤ S|v|_s; sup_{D̄}|P| ≤ ζ(2s)^{1/2}|P|_s. Satz B-lin ⇒ Konvergenz; Beschränktheit in A_s und kompakte
Einbettung A_s ↪ A_{s'} (½ < s' < s) ⇒ gleichmässig auf D̄. ∎
**Aber:** b27/b28 (J_N an f_K, Gewichte |L|^k(1+k)^s, Spektralnorm): bis N = 160 konstant 1.4352 (s=1), 1.5922 (s=2), 2.2472 (s=3)
(= Strang D), dann taucht ein Eck-Singulärwert mit Schwerpunkt bei k ≈ 0.93N auf: s = 0: 1.537, 1.841, 2.062, 2.686, 3.649, 6.132
(N = 160, 200, 220, 260, 300, 360); s = 2 überholt 1.5922 bei N ≈ 240. prec 700 vs 1100 Bit identisch, gewichtete arb-Radien
≤ 1e−38 (b28c): **echt**. Strang D (20:53): ≈ 1.012^N bis N = 400. ⇒ sup_N|J_N⁻¹| = ∞ in diesen Normen; B-Stab ist leer.
Die rechte Seite T_N J_∞ f_K ist trotzdem harmlos: ihr Überlapp mit dem Spitzen-Eingangsvektor ist 4.2e−8 (N=200), 8.7e−9 (N=300)
(b29). Es kommt also auf **kohärente Eingaben** an — genau das leistet die Pol-Darstellung.

## 7.2 Bausteine (bewiesen)

**Lemma 7.1 (Struktur in der Pol-Variablen).** Für festes z ist ζ ↦ Π_N k_ζ(z) rational mit Polen höchstens in 0 und a₀, und
E_N(ζ,z) := Π_N k_ζ(z) − 1/(ζ − z) = O(ζ^{−N−2}) für ζ → ∞; E_N(·,z) ist holomorph auf ℂ ∖ {0, a₀, z}.
*Beweis.* T_N[1/(ζ − f(x))] = T_N Σ_m h(x)^m/(ζ − a₀)^{m+1} (endliche Summe m ≤ N, da h^m = O(x^m)); T_N[g(x)/(ζ − x)] =
T_N[g(x)Σ_m x^mζ^{−m−1}]. Beide Koeffizientenvektoren sind Polynome in 1/(ζ−a₀) bzw. 1/ζ; J_N⁻¹ ist konstant. Für grosse |ζ| ist
k_ζ = Σ_m z^m ζ^{−m−1} (konvergent mitsamt den Jets an 0 und a₀), Π_N stetig auf Jets und Π_N z^m = z^m (m ≤ N), also
E_N = Σ_{m>N}ζ^{−m−1}(Π_N z^m − z^m). ∎ (Dasselbe gilt komponentenweise für Cs Vektor A⁻¹B_N[(ζ−x)⁻¹]_{>N}.)

**Lemma 7.2 (Bernstein–Walsh in ζ).** Sei K ⊂ ℂ kompakt mit zusammenhängendem Komplement, 0, a₀, z ∈ K, g_K die Green-Funktion
von ℂ∖K mit Pol ∞, t > 0. Dann für g_K(ζ) > t: |E_N(ζ,z)| ≤ max_{g_K(ζ')=t}|E_N(ζ',z)| · exp(−(N+2)(g_K(ζ) − t)).
*Beweis.* log|E_N(·,z)| + (N+2)g_K ist subharmonisch auf {g_K > t}, beschränkt bei ∞ (g_K = log|ζ| − log cap K + o(1),
E_N = O(ζ^{−N−2})), Maximumprinzip. ∎ Mit K = Ḡ', t = 1/N: exponentieller Zerfall weg vom Rand **aus einer polynomialen
Schranke auf der Niveaulinie** — das ersetzt Cs (H_R2) und den Schwanz von (H_R1). (Für Cs Vektor: log|·| einer holomorphen
vektorwertigen Funktion ist subharmonisch, gleicher Beweis.)

**Lemma 7.3 (Green-Funktion von ℂ∖Ḡ' explizit, bewiesen; b31_green_explizit.py).** M(ζ) = (ζ−L)/(ζ−L̄) bildet beide Randkreise
auf Geraden durch 0 ab, das Äussere von Ḡ' auf den Keil {β < arg < β+ω}, β = arg M(a₀+R₀) = −82.316°, ω = 158.934° (= Aussen-
winkel bei L); W = (e^{−iβ}M)^{π/ω} bildet ihn auf die obere Halbebene ab; g(ζ) = log|(W − W̄_∞)/(W − W_∞)|, W_∞ = W(∞).
Kontrolle: g = 0 auf beiden Bögen (≤ 2e−31), g(ζ) − log|ζ| → −0.41689 ⇒ cap(Ḡ') = 1.51689; am Schlitz ζ = L + u e^{iψ}
(ψ = 87.15° Winkelhalbierende): g/u^{π/ω} = 0.6553 (u=1e−4), 0.6526 (1e−2), 0.6045 (0.2), π/ω = 1.1326 = 1/α. Also liegt die
Niveaulinie g = 1/N bei L im Abstand ≈ (0.655 N)^{−0.883}, entlang der glatten Bögen im Abstand ~1/N.
Vergleich mit Cs Zerfallsgesetz R_N^{1/N} → 1/r_W := max(|L|/|ζ|, R₀/|ζ−a₀|): stets e^{−g} ≥ 1/r_W (ζ = 2.2: 0.805 vs 0.793;
ζ = −2.2: 0.640 vs 0.625) — das Maximumprinzip ist etwas schwächer als das Beobachtete, reicht aber.

**Lemma 7.4 (f_K um Ḡ' bis auf zwei Schlitze; Skizze mit numerisch geprüften Inklusionen).** f_K ist holomorph auf einer Umgebung von
Ḡ' ∖ (S_L ∪ S_{L̄}), S_L ein kurzer Schlitz von L in den Aussenkeil, stetig bis an beide Ufer. Gründe: (i) auf D Kneser;
(ii) Ḡ'∖D̄ ⊂ {z : |log z| < |L|} (63 207 Gitterpunkte, max|log z|/|L| = 0.9947), dort f_K = exp∘f_K∘log; (iii) am linken
Bogen f_K = log∘f_K∘exp (|e^z| ≤ 0.784 auf |z| = 1.02|L|, 100°…260°); (iv) nahe L Knesers lokale Form f_K = φ(√LΨe^{Lε(Ψ^γ)})
(Strang A), konvergent für kleines |Ψ^γ|, beliebige Umläufe < 2π. Daraus auch **Lemma 7.5 (Sprung):** |Δf_K(L + ue^{iψ})| ≤ C_f u^{5.447}
(Strang A misst C_f ≤ 0.95 am radialen Schlitz, u ≤ 0.2).

## 7.3 Satz B-Max (Hauptergebnis Runde 2) — bewiesen unter Lemma 7.4/7.5

**Voraussetzung (H*)** [linear, polynomial]: Es gibt C_*, κ, N₀ mit: für N ≥ N₀, alle z ∈ D̄ und alle ζ ∉ Ḡ' mit g(ζ) ≤ 1/N gilt
|Π_N k_ζ(z)| ≤ C_* N^{κ}/dist(ζ, Ḡ').

**Satz B-Max.** Unter (H*) gilt für N ≥ N₀:
  sup_{z∈D̄}|Π_N f_K(z) − f_K(z)| ≤ C·(N^{κ+1−α(2+Re γ)} + N^{κ−α(1+Re γ)} + N^{κ+1}e^{−ηN}),
mit α(1+Re γ) = 0.883·5.447 = 4.81 und α(2+Re γ) = 5.69. Für κ < 4.69 folgt Π_N f_K → f_K gleichmässig auf D̄, also (LP) — sogar
ohne Satz B-lin.
*Beweis.* (1) Kontur: Niveaulinie Γ_η = {g = η} (η klein fest), an den Stellen, wo sie die Schlitze S_L, S_{L̄} schneidet,
aufgetrennt und längs beider Schlitzufer nach L bzw. L̄ geführt, kleine Kreise um L, L̄. Nach Lemma 7.4 Cauchy-Formel für f_K;
Π_N ist linear und stetig in den N-Jets an 0, a₀ (beide im Inneren), also Π_N f_K(z) − f_K(z) = (1/2πi)∮ f_K(ζ)E_N(ζ,z)dζ
(Kreise → 0, da E_N(·,z) bei L für festes N, z ≠ L beschränkt). (2) Auf Γ_η: (H*) auf der Niveaulinie g = 1/N gibt
|E_N| ≤ (C_* N^κ + 1)/dist ≤ C N^{κ+1} (dist ≥ c/N, Lemma 7.3); Lemma 7.2 ⇒ |E_N| ≤ C N^{κ+1}e^{−(N+2)(η − 1/N)}.
(3) Schlitz ζ_u = L + ue^{iψ}: im Band g(ζ_u) ≤ 1/N (⇒ u ≤ (0.65N)^{−α}·(1+o(1))) direkt (H*) mit dist(ζ_u, Ḡ') ≥ c u:
∫|Δf||E_N|du ≤ C∫_0^{cN^{−α}} u^{5.447}N^{κ}u^{−1}du = O(N^{κ−5.447α}); ausserhalb Lemma 7.2: |E_N| ≤ C N^{κ+1}e^{−N g(ζ_u)},
g(ζ_u) ≥ c u^{1/α}: ∫u^{5.447}N^{κ+1}e^{−cNu^{1/α}}du = O(N^{κ+1−6.447α}). L̄ symmetrisch. ∎
**Zusatz (für Strang C).** Dieselbe Rechnung zeigt: Cs (H_R1) mit Zerfall und (H_R2) folgen aus einer polynomialen Schranke für
R_N(ζ) auf {g ≤ 1/N} — Satz 9.5 braucht nur noch diese.

*Bemerkung:* Satz B-Max braucht weder (G2) noch den Kern-Satz §1.3 noch Satz B-lin — die Konvergenz auf D̄ folgt direkt aus der
Integraldarstellung. B-lin bleibt als zweiter Weg, falls man (H*) nur mit Gewichten hat.
**Kette (alle Pfeile bewiesen, modulo Lemma 7.4/7.5):**  (H*) ⟹ [Satz B-Max] sup_{D̄}|Π_N f_K − f_K| → 0 ⟹ (LP);
(LP) + (Q) [quadratisch, Strang C] ⟹ Vermutung. Ohne (H*) ist nichts davon bewiesen: (H*) ist der verbleibende lineare Kern
(„Pol-Antwort im Randband“, numerisch mit C_* = 1, κ = 0 bis N = 140 und mit Spitzenanteil ~N bis N = 360).

## 7.4 Numerische Prüfung (b23, b25, b29, b30)

**(H*) im Band (b25):** dist(ζ,Ḡ')·max_{D̄}|Π_N k_ζ| für ζ im Abstand t/N:

| Ort | N=60, t=0.5/2/8 | N=100 | N=140 |
|---|---|---|---|
| Ecke L, Winkelhalbierende | 0.208/0.616/0.983 | 0.196/0.593/0.979 | 0.188/0.578/0.976 |
| Ecke L, radial | 0.196/0.580/0.924 | 0.185/0.557/0.918 | 0.178/0.543/0.914 |
| Ecke L, 15° an ∂D(a₀) | 0.057/0.195/0.462 | 0.053/0.182/0.434 | 0.051/0.175/0.420 |
| links, −|L| | 0.298/0.753/0.996 | 0.296/0.753/0.996 | 0.295/0.752/0.996 |
| rechts, a₀+R₀ | 0.017/0.066/0.220 | 0.011/0.041/0.145 | 0.0075/0.029/0.108 |
| D(a₀) oben 60° | 0.029/0.106/0.325 | 0.018/0.067/0.223 | 0.013/0.048/0.170 |
| D oben 120° | 0.292/0.746/0.995 | 0.291/0.745/0.995 | 0.290/0.745/0.995 |

Fortsetzung (b25b, 720 Kreispunkte): N = 200 / 300, t = 0.5/2/8 — Winkelhalbierende 0.181/0.562/0.972 bzw. 0.174/0.545/0.967;
radial 0.172/0.528/0.910 bzw. 0.164/0.512/0.904; 15° 0.049/0.168/0.407 bzw. 0.047/0.161/0.394; links 0.295/0.752/0.996 (beide);
rechts 0.0053/0.021/0.078 bzw. 0.0035/0.014/0.053; D(a₀) 60° 0.0088/0.034/0.125 bzw. 0.0059/0.023/0.087; D 120° 0.290/0.744/0.995
bzw. 0.289/0.744/0.996.
Überall ≤ 1 und fallend oder konstant in N: **(H*) mit C_* = 1, κ = 0 bis N = 300** — obwohl J_N⁻¹ dort die Eckspitze σ₁ = 3.65 hat.

**Pol-Fehler am Schlitz (b23):** u·M_N(u), M_N = max_{D̄}|E_N|: für u ≪ N^{−0.9} → 1 (Pol nicht auflösbar), Halbwert bei
u = 0.026, 0.0165, 0.012, 0.0095 (N = 60, 100, 140, 180) ⇒ Skala N^{−0.90…0.95}; für u ≥ 0.1 Zerfall e^{−0.51 N u}. Äussere Kontur
(Abstand 0.2): 6.6e−3, 9.0e−5, 1.3e−6, 1.8e−8. **Schranke gegen Wahrheit:** I_N = (1/π)∫|Δf|M_N du = 1.66e−8, 1.33e−9, 2.49e−10,
7.15e−11 gegen sup_{D̄}|Π_N f_K − f_K| = 7.2e−9 (N=60), 3.3e−11 (N=180): Faktor ≈ 2.2 — die Dreiecksungleichung im Integral
verliert nichts Wesentliches (wie bei C: 2.06–2.28). Achtung Zweig: auf dem ψ-Schlitz ist ein Ufer des Singulärterms um
e^{1.058·2π} ≈ 770 grösser als das andere (komplexer Exponent); der Sprung ist ≈ |D₁|·22.9·(u/|L|)^{5.447}.

**Eckspitze und Kohärenz (b29, b30):** gewichtete SVD von J_N⁻¹ (w_k = |L|^k): Überlapp |⟨v₁, w̃⟩|/|w̃| des Spitzen-Eingangsvektors

| Eingabe | N = 200 | N = 300 |
|---|---|---|
| T_N J_∞ f_K | 4.2e−8 | 8.7e−9 |
| Pol L + 0.02e^{iψ} | 2.8e−2 | 1.2e−2 |
| Pol L + 0.2e^{iψ} | 1.1e−4 | 2.7e−5 |
| Pol −1.6 / Pol 2.0 | 1.4e−8 / 9e−18 | 3.7e−9 / 5e−18 |
| e_N (oberster Koeffizient) | 2.0e−7 | 1.2e−10 |
| Zufallsvektor | 7.0e−2 | 7.9e−2 |

Spitzenbeitrag σ₁·Überlapp für Pole bei u = v·N^{−0.883} (v = 1): 0.60, 0.80, 1.00, 1.17, 1.35 (N = 160, 200, 250, 300, 360) ≈ N^{1.0},
normale Antwort 4.8 → 6.8 ≈ N^{0.43}: polynomial, obwohl σ₁ exponentiell wächst (1.54 → 6.13). Der Überlapp fällt also exponentiell
genau so, dass das Produkt polynomial bleibt.

## 7.5 Toeplitz-Bild (Leitungsfrage 2)

Ds Grenzoperator G = T(σ), σ = (1/(f'∘f))·1_{Γin}: an L springt σ von 1/μ (|·| = 0.853) auf 0. Die klassische Theorie der Abschnitte
T_N(a) mit stückweise stetigem a (Gohberg–Krupnik, Böttcher–Silbermann: Stabilität ⇔ T(a), T(ã) invertierbar, Sprünge durch
Bögen gefüllt) setzt ein *invertierbares* Symbol voraus; σ verschwindet aber auf dem ganzen Bogen Γout (dort dominiert die
Komposition, J ist dort „unendlich“). Der Satz ist also nicht anwendbar, und Instabilität der Abschnitte ist im Einklang mit der
exponentiellen Eckspitze. Was übrig bleibt, ist die **schwache** Aussage (H*): nicht |J_N⁻¹| beschränkt, sondern J_N⁻¹ auf der
Familie der kohärenten Eingaben T_N J_∞ k_ζ (ζ im Randband) polynomial beschränkt. Diese Familie hat in |L|-Gewichten nur Norm
~1/dist für ζ ausserhalb K := D̄ ∪ f_K(D̄) (Pole von J_∞k_ζ bei ζ und f_K⁻¹(ζ), beide ausserhalb D̄); für ζ ∈ f_K(D̄)∖Ḡ' (rechter Rand
von D(a₀)) liegt f_K⁻¹(ζ) in D und die Eingabe ist exponentiell gross — dort braucht auch (H*) Auslöschung (gemessen ist sie da:
b25 „rechts“ 0.0075…0.22).

## 7.6 Inneres Problem und Faktor N (Leitungsfrage 1) — Stand

- **Zwei Skalen.** (a) Ecke N^{−α}, α = 0.883: Pol-Antworten (b23, Cs H_R1), Randfehler |p_N(L) − L| ~ N^{−αβ} = N^{−4.81}
  (gemessen −4.83), Green-Funktion g ≈ 0.655 u^{1/α}. (b) Taylor-Skala 1/N: Modenamplitude c_N = C N^{1−2γ̄} — Betrag UND Phase;
  Strang A: innere Schicht (p_N − f_K)(L(1 − v/N))N^γ → Φ(v) für v ≥ 4 ist die Fortsetzung des äusseren Modus.
- **Eckskalen-Anpassung reicht nicht:** c̄ ρ^{1−γ} ≈ ρ^{1+γ} mit ρ = N^{−α} gäbe c ~ N^{−2αγ̄} = N^{−7.85 − 1.87i}; gemessen −7.888 − 2.107i
  (A, freier Fit) bzw. C_N stabil mit N^{1−2γ̄} (b10). Die Phase schliesst die Eckskala aus. Dass 2(1−α)Re γ = 1.04 ≈ 1 ist, ist Zufall.
- **Pol-Moden-Antwort** m_N(ζ) (b22): K_N(v) := m_N δ_N/(N ζ^{−N−1}) bei ζ = L(1 + v/N) ist NICHT skaleninvariant (|K| ∝ N^{−0.75} bei
  v = 0.25, N^0 bei v = 8, N^{+0.84} bei v = 16), d. h. m_N fällt auf der Eckskala ab (≈ e^{−0.93N^{0.883}u}, Vorfaktor ≈ N^{−2.2}).
  Erst das Integral gegen Knesers Sprung ~u^{1+γ} ergibt die saubere 1/N-Phase — ein Zwei-Skalen-Effekt, den ich nicht herleiten kann.
- **Mellin:** lokal ist J_∞ auf (z−L)^s Multiplikation mit μ^s + μ = 2μ^{(s+1)/2}cosh((s−1)L/4) = 2πμ^{(s+1)/2}/(Γ(½+iζ_s)Γ(½−iζ_s)),
  ζ_s = (s−1)L/(4π). Ein lokales (affines) Modell bei L reproduziert die Regel aber nicht (b20, b20b: zusätzlicher Faktor e^{0.092iN}
  durch den komplexen Knoten f(0) = L(1−μ)) — das innere Problem ist nicht rein lokal, die zwei Knoten 0 und a₀ gehören hinein.
- **Konstanten:** A: K ≈ 0.002 − 0.051i (c_N = K·N·u_N/(δ_{−1})_N); B: C ≈ 0.18 − 0.08i (c_N/N^{1−2γ̄}), linear Q/N ≈ 0.0025 + 0.0045i.
  Nicht hergeleitet.
- **Eckschicht im Indexraum = Streckband [N/|μ|, N] (numerisch, b32; Deutung = Vermutung).** Bei L wirkt die Komposition mit f auf
  Koeffizientenprofile e_k = L^{−k}φ(k/N) als komplexe Streckung x ↦ x/μ (x = k/N, |μ| = 1.1717); die duale Gleichung koppelt x an
  μx und stösst für x > 1/|μ| = 0.853 an die Abschneidung x = 1. Gemessen: der Eingangsvektor v₁ der Eckspitze hat seine Masse
  in k/N ∈ [0.87, 0.96] (Quantile 5/50/95 % bei N = 300: 0.873/0.950/0.960; N = 200: 0.825/0.945/0.960), also im Band
  [1/|μ|, 1]. **Deutung des Faktors N:** das Auswahlfunktional summiert über ein Band von ~0.15·N obersten Koeffizienten, auf dem
  u_k ≈ u_N(k/N)^{−s−1} fast konstant ist ⇒ c ≈ N·u_N·K(s)/(δ_{−1})_N mit K(s) = ∫_{1/|μ|}^1 x^{−s−1}Λ̃(x)dx (schwach s-abhängig, wie
  in b19 gemessen). Das ist eine Plausibilitätsdeutung, kein Beweis; Λ̃ ist nicht bestimmt.

## 7.7 Dateien Runde 2

b22_polmoden.py (+json), b23_polfehler.py (+log, json), b24_toeplitz_bruecke.py (nicht ausgeführt), b25_band.py (+log),
b26_separabel.py, b27_stabilitaet_reicht.py (+log), b28_stab_wachstum.py, b28b_stab_gross.py, b28c_genauigkeit.py (+logs),
b29_spitze_kohaerenz.py (+log), b30_spitze_skalen.py (+log), b31_green_explizit.py, b20b_affin_weit.py,
b25b_band_gross.py (+log, N = 200/300), b32_eckschicht_mu.py (+log), brett7.txt.

---

# Runde 3 (ab 21:19): quadratische Kontrolle und Lemma 7.4/7.5

## 8.0 Kurzfassung Runde 3

1. **Newton um x₁ = Π_N f_K ist harmlos (numerisch, b33/b36).** F_N(x₁) ist exakt rein quadratisch, und
   sup_{D̄}|x₁ − p_N| = c_N·ε_N² mit ε_N := sup_{D̄}|x₁ − f_K| und c_N = 0.132, 0.171, ≈0.20, 0.256 (N = 150, 200, 232, 300; ≈ N^{0.95}).
   Der zweite Newton-Schritt ist 9e−43 … 8e−48 (quadratische Konvergenz); die Korrekturen sitzen im Inneren des Indexraums
   (Schwerpunkt k/N ≈ 0.55), nicht in der Eckschicht (≈ 0.93); Überlapp der Eingabe F(x₁) mit dem Spitzenvektor 1e−3 … 1e−2.
   Das Wachstum, das Strang D beim Newton-Schritt sieht, gehört zum Zentrum q_N, nicht zu x₁.
2. **Satz B-Sek (bewiesen, a posteriori):** p_N ist exakt die *Sekanten-Projektion* von f_K. Erfüllt deren Pol-Antwort (H*), so
   folgt p_N → f_K auf D̄ (Beweis wie Satz B-Max). Die Vermutung ist damit äquivalent zu einer linearen Aussage (über einen von p_N
   abhängigen Operator).
3. **Satz B-Q (bewiesen unter einer Hierarchie (H*_n)):** Jede Ordnung n der Lösung p = x₁ + d ist ein n-faches Cauchy-Integral
   eines *n-linearen Pol-Kerns* Ξ_N^{(n)}(ζ₁,…,ζ_n), der in jeder Variablen rational ist (Pole nur in 0, a₀, Abfall O(ζ^{−N−2}),
   Lemma 8.1) — das Maximumprinzip in ζ (Lemma 7.2) gilt variablenweise. Unter (H*) und (H*_n): |Ξ^{(n)}| ≤ Cⁿ N^{nκ}/Πdist(ζ_i)
   im Randband, κ < 4.69, gibt es für grosse N eine Lösung p*_N mit sup_{D̄}|p*_N − x₁| ≤ C N^{2κ−9.38}, also p*_N → f_K.
   (H*₂) numerisch mit kleinen Konstanten (b34). **Offen:** Beweis von (H*), (H*_n); Identifikation p*_N = Zweig p_N (numerisch ja).
4. **Lemma 7.4/7.5 bewiesen aus Knesers Struktursatz (K)** (§8.3): die lokale Formel f_K = φ(√L·Ψ·e^{Lε(Ψ^γ)}) mit ε holomorph,
   ε(0) = 0 folgt aus (K); daraus Fortsetzung auf die geschlitzte Kreisscheibe um L und die Sprungschranke |Δf| ≤ C|ζ−L|^{5.447};
   die Fortsetzung um ganz Ḡ' folgt aus den Funktionalgleichungen und zwei **streng (arb) geprüften** Inklusionen (b35).
   Numerisch bleiben nur: (K) selbst bzw. seine Konstanten (Literatur + Strang A) und eine Prüfung am linken Bogen (b35b').

## 8.1 Exakte Struktur um x₁ (bewiesen)

Φ(p) := p∘p, F_N(p) = T_N[Φ(p) − eˣ], J := J_N(f_K), x₁ := Π_N f_K, u₁ := x₁ − f_K, Q(u) := Φ(f+u) − Φ(f) − DΦ(f)u.
(i) F_N(x₁) = T_N Q(u₁), weil T_N[Φ(f) − eˣ] = 0 und T_N DΦ(f)u₁ = J x₁ − T_N J_∞ f = 0 (Definition von Π_N).
(ii) Für d ∈ P_N: F_N(x₁ + d) = J d + T_N Q(u₁ + d).
(iii) Q(u) = Σ_{j≥2} (f^{(j)}∘f) u^j/j! + Σ_{j≥1} (u^{(j)}∘f) u^j/j! (Taylor in u an der Stelle f(x)), also
Q(u) = ½D²Φ(f)[u,u] + O(u³), D²Φ(f)[a,b] = (a'∘f) b + (b'∘f) a + (f''∘f) a b, f''∘f = g'/f', g = f'∘f = eˣ/f'.
Für festes N hängt T_N Q(u) nur von endlich vielen Jet-Koeffizienten von u in 0 und a₀ ab.

**Satz B-Sek (bewiesen).** Sei p ∈ P_N eine Lösung, 𝒥v := v∘p + M·v mit M := ∫₀¹ f'(f + s(p − f)) ds, und 𝒥_N := T_N𝒥|_{P_N}
invertierbar. Dann p = Π^{sek} f := 𝒥_N⁻¹T_N𝒥 f, denn 𝒥(p − f) = p∘p − f∘p + f∘p − f∘f = p∘p − f∘f, also T_N𝒥(p − f) = 0.
Π^{sek}k_ζ ist rational in ζ mit Polen in 0 und p(0) (Lemma 7.1 wörtlich), O(ζ^{−N−2}). Gilt |Π^{sek}k_ζ(z)| ≤ C N^κ/dist(ζ,Ḡ')
im Band {g ≤ 1/N}, z ∈ D̄, so ist sup_{D̄}|p − f_K| ≤ C N^{κ−4.69} (Beweis von Satz B-Max §7.3 ohne Änderung; p(0) liegt im Inneren
von Ḡ', sobald |p(0) − a₀| klein ist). ∎ Bedeutung: a posteriori — die Hypothese spricht über p selbst.

## 8.2 Multilineare Pol-Kerne und Satz B-Q

**Konstruktion.** F_λ(d) := J d + T_N Q(λu₁ + d), λ ∈ ℂ, d ∈ P_N. F_λ ist analytisch nahe (0,0), F_0(0) = 0, D_dF_0(0) = J.
Satz über implizite Funktionen (endlichdimensional): es gibt eine analytische Lösung d(λ) = Σ_{n≥2} λⁿ d^{(n)} für kleines |λ|, und
F_1(d(1)) = F_N(x₁ + d(1)) — falls die Reihe bei λ = 1 konvergiert. Koeffizientenvergleich: d^{(2)} = −J⁻¹T_N ½D²Φ(f)[u₁,u₁], allgemein
d^{(n)} = 𝒟_n[u₁,…,u₁] mit einer symmetrischen n-linearen Abbildung 𝒟_n, die eine *endliche* Summe aus J⁻¹, T_N, Produkten,
Ableitungen und Komposition mit f ist (nur Terme von Q bis Grad n). Mit E_ζ := Π_N k_ζ − k_ζ setze
  Ξ_N^{(n)}(ζ₁,…,ζ_n) := 𝒟_n[E_{ζ₁},…,E_{ζ_n}] ∈ P_N.

**Lemma 8.1 (bewiesen).** Für feste übrige Variablen ist ζ_i ↦ Ξ_N^{(n)} rational mit Polen höchstens in 0 und a₀, und O(ζ_i^{−N−2}).
*Beweis.* 𝒟_n braucht nur endlich viele Jet-Koeffizienten in 0 und a₀. Die von Π_N k_ζ sind Polynome in 1/ζ, 1/(ζ−a₀) (Lemma 7.1),
die von k_ζ in 0 bzw. a₀ sind ζ^{−m−1} bzw. (ζ−a₀)^{−m−1}, die von k_ζ^{(j)}∘f sind Polynome in 1/(ζ−a₀). Multilinearität erhält die
Polynomstruktur. Für |ζ| gross ist E_ζ = Σ_{m>N}ζ^{−m−1}(Π_N x^m − x^m) als Jet. ∎
**Darstellung.** Nach Lemma 7.4 und Satz B-Max (1) ist u₁ = (1/2πi)∮_Γ μ(ζ)E_ζ dζ als Jet in 0 und a₀ (Γ: Niveaulinie mit
Schlitzufern; μ = f_K, auf den Schlitzen der Sprung). 𝒟_n ist stetig auf Jets, also
  d^{(n)} = (2πi)^{−n}∮…∮ μ(ζ₁)…μ(ζ_n) Ξ_N^{(n)}(ζ₁,…,ζ_n) dζ₁…dζ_n.

**Hypothese (H*_n).** |Ξ_N^{(n)}(ζ₁,…,ζ_n)(z)| ≤ Cⁿ N^{nκ}/Π_i dist(ζ_i,Ḡ') für alle ζ_i mit g(ζ_i) ≤ 1/N, z ∈ D̄, N ≥ N₀.
*[Korrektur Runde 4, §9.3(c): gleichmässig im ganzen Band numerisch falsch (Paare mit linkem Pol ~1.012^N); B-Q gilt mit der lokalen
Form (H*_n^{loc}).]*

**Satz B-Q (bewiesen unter (H*) und (H*_n) für alle n ≥ 2, mit κ < 4.69).** Für N ≥ N₁ konvergiert Σλⁿd^{(n)} für |λ| ≤ 1 + δ,
p*_N := x₁ + d(1) löst F_N(p*_N) = 0, und sup_{D̄}|p*_N − x₁| ≤ C N^{2κ−9.38}, sup_{D̄}|p*_N − f_K| ≤ C N^{κ−4.69}.
*Beweis.* (1) Variablenweises Maximumprinzip: wende Lemma 7.2 nacheinander in ζ₁, ζ₂, … an (Lemma 8.1 erlaubt das für jede Wahl
der übrigen Variablen). Ergebnis: |Ξ^{(n)}| ≤ Cⁿ N^{nκ}Π_i w(ζ_i) mit w = 1/dist im Band und w = C N e^{−(N+2)(g−1/N)} ausserhalb
(dort dist ≥ c/N auf der Niveaulinie, Lemma 7.3). (2) Das n-fache Integral zerfällt in ein Produkt von n Einfachintegralen
∫|μ|w|dζ|; jedes ist genau die Abschätzung im Beweis von Satz B-Max, also ≤ C N^{−4.69}. Damit sup_{D̄}|d^{(n)}| ≤ (C' N^{κ−4.69})ⁿ.
(3) Für N^{4.69−κ} > 2C' konvergiert die Reihe für |λ| ≤ 2; sup|λu₁ + d(λ)| bleibt klein, F_λ(d(λ)) ist dort analytisch in λ und
nahe 0 identisch null, also überall (Identitätssatz). ∎
**Lücken:** (α) (H*_n) unbewiesen (gleiche Natur wie (H*), Strang C); (β) p*_N ist *eine* Lösung nahe x₁; dass es der Zweig p_N
(Fortsetzung von p₁ = x + ½) ist, ist nur numerisch (b33: Newton ab x₁ landet für N = 150, 200, 232 auf p_N aus gemeinsam/);
(γ) Satz B-Q zusammen mit (β) ergibt die Vermutung mit Rate N^{κ−4.69}.
*Bemerkung:* Für n = 2 allein (ohne Hierarchie) gilt bewiesen: (H*) + (H*₂) ⇒ sup|d^{(2)}| ≤ C N^{κ₂−9.38}, d. h. die erste
nichtlineare Korrektur ist um N^{4.69−κ} kleiner als der lineare Fehler ε_N.

## 8.3 Lemma 7.4/7.5 aus Knesers Struktursatz

**(K) Knesers Struktur (Literatur: Kneser 1950, Paulsen–Cowgill 2017; hier übernommen, von Strang A numerisch bestätigt).**
Es gibt c > 0 und eine 1-periodische, auf {Im w > c} holomorphe, beschränkte Funktion θ mit θ(w) → θ_∞ (Im w → ∞), so dass
T(w) = φ(exp(L·(w + θ(w)))) dort; φ = Ψ⁻¹, Ψ die Schröder-Funktion von exp bei L (Ψ(eᶻ) = LΨ(z), Ψ(L) = 0, Ψ'(L) = 1);
für Im w < −c das Konjugierte. f_K = T(A + ½), A = T⁻¹ die Abelfunktion.

**Lemma 8.2 (lokale Formel; bewiesen aus (K)).** Es gibt r > 0 und eine bei 0 holomorphe Funktion ε mit ε(0) = 0, so dass auf
D ∩ D(L,r):  f_K(z) = φ(√L·Ψ(z)·e^{Lε(q)}),  q = Ψ(z)^γ := exp(γ log Ψ(z)),  γ = 2πi/L,  √L := e^{L/2}.
*Beweis.* χ := id + θ ist für Im w gross injektiv (θ' → 0), χ(w+1) = χ(w)+1, Inverse χ⁻¹ = id + ϑ mit ϑ 1-periodisch, ϑ → −θ_∞.
Für z ∈ D nahe L ist A(z) = χ⁻¹(u), u := log Ψ(z)/L, mit dem Zweig von log, für den das an einer Stelle gilt (e^{Lχ(A)} = Ψ legt u
modulo γ fest, Stetigkeit legt den Zweig fest); Im u → +∞ für z → L. Dann f_K = T(A + ½) = φ(e^{L·χ(χ⁻¹(u) + ½)}) und
ε(u) := χ(χ⁻¹(u) + ½) − u − ½ ist 1-periodisch, holomorph, → 0 (θ_∞ hebt sich weg), also ε = Σ_{m≥1}ε_m q^m, q = e^{2πiu} = Ψ^γ.
e^{L(u+½)} = Ψ·e^{L/2}. ∎ (Strang A: ε = 2κ₁q + 4πiκ₁²q² + …, κ₁ = −8.35e−6 − 2.757e−4i; Formel bis |ζ−L| ≈ 0.14 auf 5e−37 geprüft.)

**Lemma 8.3 (Fortsetzung um L, Sprung = Lemma 7.5; bewiesen aus Lemma 8.2).** Sei Σ_L ein Strahl von L in den Aussenkeil von Ḡ'
(z. B. Winkelhalbierende ψ = 87.15°). Für kleines r ist die rechte Seite von Lemma 8.2 auf D(L,r)∖Σ_L holomorph: arg Ψ läuft dort
(stetig von D aus) in einem Intervall der Länge 2π, also |q| = |Ψ|^{Re γ}e^{−Im γ·arg Ψ} ≤ (C r)^{4.447}e^{2π·1.058} < Konvergenzradius.
Sie stimmt auf D ∩ D(L,r) mit f_K überein und ist stetig bis an beide Ufer. Die Ufer gehören zu q₊ und q₋ = q₊e^{±2πiγ}, daher
|Δf_K(ζ)| ≤ sup|φ'|·|√L|·|Ψ(ζ)|·|e^{Lε(q₊)} − e^{Lε(q₋)}| ≤ C_J|ζ − L|^{1+Re γ} = C_J|ζ − L|^{5.447}. ∎
Konstante (numerisch, Strang A): |Δf|/u^{5.447} = 0.94 … 0.76 für u = 0.005 … 0.2 am radialen Schlitz.

**Lemma 8.4 (Fortsetzung um Ḡ' bis auf die Schlitze = Lemma 7.4; bewiesen bis auf (K) und b35b').** Es gibt eine Umgebung W von Ḡ',
so dass f_K auf W∖(Σ_L ∪ Σ_L̄) holomorph und bis an die Schlitzufer stetig ist.
(i) D: Kneser.
(ii) Halbmond: F₂(z) := exp(f_K(Log z)), Log := (exp|_D)⁻¹ (exp ist auf D injektiv, da |Im w| < |L| < π). Auf {w ∈ D: e^w ∈ D} =
D ∩ {Re w < Re L} (konvex, enthält ein reelles Intervall) gilt f_K(e^w) = e^{f_K(w)} (Funktionalgleichung, Identitätssatz); also
F₂ = f_K auf D ∩ exp(D). **Streng (b35a, arb, 200 Bit):** |Log z|² ist subharmonisch auf dem Halbmond D(a₀,R₀)∖D̄ (er liegt in Re z > 0);
auf ∂D ist |Log z|² = log²|L| + arg² z ≤ |Log L|² = |L|² mit Gleichheit nur in L, L̄ (denn Log L = L); auf dem äusseren Bogen und den
Bögen |z − L| = 0.02 gilt max|Log z|²/|L|² ≤ 0.98290 (20 000 bzw. 4 000 Intervalle, keine Verletzung). Also liegt der Halbmond ohne
D(L,0.02) ∪ D(L̄,0.02) kompakt in exp(D) (|Log z| ≤ 0.9914|L|), und F₂ ist auf einer Umgebung davon holomorph.
(iii) Linker Bogen (∂D ausserhalb D(a₀,R₀), ohne D(L,0.02)): F₃(z) := Log f_K(eᶻ). **Streng (b35b):** dort Re z − log|L| ≤ −0.01933,
also eᶻ ∈ D kompakt (|L| = e^{Re L}). **Numerisch (b35b', 2000 Kneser-Koeffizienten):** f_K(eᶻ) bleibt weg von (−∞, 0] (min|f_K(eᶻ)|
= 0.734, Abstand von arg zu ±π ≥ 1.82 rad), also ist Log dort holomorph. Übereinstimmung: auf dem (zusammenhängenden) Durchschnitt
mit D ist e^{F₃} = f_K(eᶻ) = e^{f_K(z)}, also F₃ − f_K ∈ 2πiℤ konstant, und = 0 im reellen Punkt nahe −|L|.
(iv) Bei L, L̄: Lemma 8.3. Konsistenz mit (ii): nahe L ist exp(D) der Halbkreis der Richtungen (243°, 423°); der Schlitz (87°)
liegt ausserhalb, der Halbmond-Keil innerhalb.
(v) Die Stücke überlappen paarweise in zusammenhängenden Mengen, die D treffen; dort stimmen sie mit f_K überein ⇒ eine holomorphe
Funktion auf W∖Schlitze (Überlappungsgeometrie elementar, nicht ausgeschrieben). ∎
**Was nur numerisch bleibt:** (K) bzw. seine Konstanten (θ, c, r, C_J; Strang A); b35b'. Alles andere in Lemma 7.4/7.5 ist bewiesen.

## 8.4 Numerik Runde 3

**b33 + b36 (Newton um x₁; ε_N := sup_{|z|=|L|}|x₁ − f_K|, Knesers Reihe mit 4000 Koeffizienten):**

| N | σ₁(J⁻¹) | Überlapp F(x₁)↔v₁ | sup|J⁻¹F(x₁)| = sup|x₁ − p_N| | Schwerpunkt k/N | 2. Schritt | ε_N | sup|x₁−p_N|/ε_N² |
|---|---|---|---|---|---|---|---|
| 150 | 1.49 | 1.1e−2 | 8.64e−22 | 0.58 | 9.4e−43 | 8.10e−11 | 0.132 |
| 200 | 1.84 | 3.3e−3 | 6.55e−23 | 0.56 | 6.9e−45 | 1.96e−11 | 0.171 |
| 232 | 2.22 | 6.7e−3 | 1.78e−23 | 0.55 | 5.8e−46 | ≈ 9.4e−12 (interp.) | ≈ 0.20 |
| 300 | 3.65 | 1.3e−3 | 1.83e−24 (p = x₃) | 0.53 | 7.6e−48 | 2.67e−12 | 0.256 |

p_N für N ≤ 232 aus gemeinsam/ (Abstand x₂ − p_N = zweiter Schritt, d. h. Newton ab x₁ trifft *diesen* Zweig). J(f)⁻¹ und J(x₁)⁻¹
geben dieselbe Korrektur auf 3 Stellen. ε_N fällt wie ≈ N^{−4.9} (verträglich mit N^{κ−4.69}, κ ≈ 0); c_N ≈ N^{0.95}.

**b34 (bilinearer Kern, (H*₂)): dist(ζ)·dist(ζ')·max_{D̄}|Ξ_N^{(2)}(ζ,ζ')|**, Pole im Abstand t/N:

| Polpaar (ψ = 87.15° Winkelhalbierende bei L) | N = 100 | N = 160 | lokaler Exponent |
|---|---|---|---|
| L, t=0.5 × L, t=0.5 | 9.2e−4 | 1.3e−3 | 0.7 |
| L, t=2 × L, t=2 | 2.6e−3 | 4.0e−3 | 0.9 |
| L, t=8 × L, t=8 | 3.9e−5 | 8.2e−5 | 1.6 |
| L, t=2 × L̄, t=2 | 4.6e−5 | 4.1e−5 | −0.2 |
| L, t=2 × links (−|L|), t=2 | 2.1e−4 | 4.2e−4 | 1.5 |
| links × links, t=2 | 1.6e−5 | 4.5e−5 | 2.2 |

(N = 220 nicht gerechnet: geschätzt ≈ 30 min, über der Laufzeitgrenze.) Also (H*₂) mit sehr kleiner Konstante und κ₂ ≲ 2.2 aus zwei Stützstellen
(grob) — weit unter der Schwelle 2·4.69 von Satz B-Q und verträglich mit Cs Randgrösse S_N ~ N^{1.12} je Variable. Das
integrierte Gegenstück ist c_N = sup|x₁ − p_N|/ε_N² ~ N^{0.95} (b33/b36).

## 8.5 Zu Strang E, Satz T (Faktor N)

Exakte Analogie (bewiesen, eine Zeile): Π_N f_K − f_K = Σ_{m>N} a_m(Π_N x^m − x^m) (Π_N reproduziert P_N und ist stetig auf Jets in
0, a₀), mit Π_N x^m = J_N⁻¹T_N[f_K^m] für m > N. Das ist die Schwanzsumme von Satz T(4) (u_N − u = −Σ_{n>N} c_n n! h_n); der
„Zusatzfaktor N“ entspricht der Zahl der wirksam beitragenden Schwanzterme m ∈ (N, N + O(N)) — **Vermutung**, nicht ausgerechnet.
In der Pol-Darstellung ist es derselbe Mechanismus: E_ζ = Σ_{m>N}ζ^{−m−1}(Π_N x^m − x^m). Die Übertragung auf exp über
Abel-Koordinaten habe ich nicht versucht (Hauptaufgabe war die quadratische Kontrolle).

## 8.6 Stand der Kette nach Runde 3

(K) [Literatur] ⟹ Lemma 7.4/7.5 [bewiesen bis b35b'] ; (H*) ⟹ [B-Max] Π_N f_K → f_K ;
(H*) + (H*_n) ⟹ [B-Q] Lösungen p*_N → f_K ; p*_N = p_N [numerisch] ⟹ Vermutung.
Alternativ a posteriori: (H*) für den Sekantenoperator ⟹ [B-Sek] p_N → f_K.
Offen bleiben: (H*), (H*_n) (Kern: Pol-Antwort im Randband, Strang C), die Zweig-Identifikation, b35b' streng.
Dateien: b33_newton_x1.py (b33.log), b34_bilinear.py (b34.log), b35_lemma74_streng.py (b35.log), b36_epsN.py (b36.log).

---

# Runde 4 (ab 21:45): Rückführung der multilinearen Kerne, Zweig-Identifikation

## 9.0 Kurzfassung Runde 4

1. **Ableitungsidentität (bewiesen, exakt; b40 symbolisch geprüft).** Für beliebige Keime a, b mit r_a := J_∞a, r_b := J_∞b gilt
   D²Φ(f)[a,b] = −J_∞[(ab)'/(2f')] + (1/(2f'))·[b r_a' − b' r_a + a r_b' − a' r_b − (g'/g)(a r_b + b r_a)] + (r_a r_b)'/(2gf').
   Sie gilt für jedes f (rein algebraisch). Der quadratische Term ist also bis auf Terme mit J_∞a, J_∞b ein **J_∞-Bild**.
2. **Folgerung für den bilinearen Kern (bewiesen):** Ξ^{(2)}(ζ,ζ') = ¼·Π_N[(E_ζE_ζ')'/f'] − κ_top·t_N mit t_N := J_N⁻¹x^N und
   κ_top = (N+1)(E_ζ'(0)ρ_ζ + E_ζ(0)ρ_ζ')/(4f'(0)). Der bilineare Kern ist also die **lineare** Projektion Π_N einer expliziten
   Funktion plus ein Kopfterm. Numerisch (b41, N = 60…400) stimmt das bis auf Rundung (relativ ≤ 1e−20). Der Kopfterm ist
   ≤ 1.4e−9, die Kopfantwort T_N := sup_D̄|t_N|/|L|^N fällt wie 0.937^N (4.2e−9 bei N = 300, 6.6e−12 bei N = 400): **keine Spitze**.
3. **(H*) ⇒ (H*₂) lässt sich so nicht zeigen.** Der Cauchy-Weg braucht die Pol-Fehler auf dem Halbmond Ḡ'∖D̄. Dort wachsen die Pol-Antworten
   exponentiell (linker Pol 0.03 → 311, N = 60 → 400; Eckpol ab N ≈ 300: 0.57 → 3.5), und Π_N f_K konvergiert auf dem Halbmond
   nicht (b42: Fehler bei a₀+R₀ ab N ≈ 300 wachsend, 1.4e−10 bei N = 500). Der Kern selbst bleibt auf D̄ für Eckpol-Paare polynomial
   (N^{0.92}: 1.6e−3 → 9.2e−3, N = 60 → 400). Paare mit einem linken Pol wachsen mit der Spitzenrate (1.012–1.021 pro Grad), was
   in Satz B-Q harmlos ist, da ferne Variablen mit e^{−Nη} gedämpft werden. Die nötige Zusatzaussage ist **linear im Operator**:
   „Π_N angewandt auf Produkte von Eck-Pol-Fehlern ist polynomial beschränkt“ ((H*_Π), §9.3). Exponenten-Bilanz: gemessen
   κ₂ ≈ 0.9 ≈ α = 0.883 (eine Ableitung auf der Eckskala), Schwelle κ₂ < κ + 4.69.
4. **Normalform für alle Ordnungen (Skizze, §9.4):** Ist J_∞a = 0, so ist W := a/(2f') exp-invariant (W∘exp = eˣ·W; bewiesen).
   Daraus ergibt sich eine explizite formale Reihe S(a) = a + aa'/(2f') + … mit Φ(f + S(a)) − eˣ = R(a; J_∞a), wobei jeder Term
   von R einen Faktor (J_∞a)^{(i)} hat. Für Pol-Fehler (J_∞a = O(x^{N+1})) ist damit jeder Kern n-ter Ordnung gleich Π_N einer
   expliziten n-linearen Ableitungsform plus Kopfband-Terme.
5. **Zweig-Identifikation (Satz B-Zweig, §9.5): bewiesen unter den Kern-Hypothesen, ohne Eindeutigkeitskugel.**
   - Die Induktion N → N+1 läuft über eine Zwei-Parameter-Familie (λ, σ) mit analytischer Fortsetzung.
   - Exakt gilt: F_{N+1}(pad p_N) = c_{N+1}x^{N+1}, x₁(N+1) − pad x₁(N) = −ρ^f_N·t_{N+1}, und der Zweigschritt ist ≈ −c_{N+1}t_{N+1}.
   - Numerisch (b43, N = 200…299) trifft die Vorhersage den echten Schritt (4e−14 … 6e−13) auf 1e−22, mit J(pad p_N) auf 1e−25 … 1e−27.
6. **Bilineare Messung bei N = 220 (Aufgabe 3):** erledigt, zusätzlich N = 300 und 400. Die neue flint-Maschine braucht 25 s, 53 s bzw. 122 s.
7. **Ordnungen an echten Daten (b44, N = 100…500):** |d^{(2)}|/ε² ~ N^{0.94}, |p_N − x₁ − d^{(2)}|/ε³ ~ N^{1.82}. Es kommt also genau eine
   Eckskalen-Ableitung N^{≈0.9} pro Ordnung hinzu, und die Reihe aus Satz B-Q konvergiert mit dem Verhältnis ≈ N^{−3.8}.

## 9.1 Die Ableitungsidentität

**Satz B-Id (bewiesen).** Seien f, a, b holomorph (Keime), f' ≠ 0, g := f'∘f ≠ 0, r_a := a∘f + g·a, r_b := b∘f + g·b. Dann:
  f'·D²Φ(f)[a,b] = r_a'·b + r_b'·a − (g·ab)',  und  (g·c)'/f' = J_∞[c'/(2f')] + [(g(a r_b + b r_a))' − (r_a r_b)']/(2g f') für c = ab,
also die Formel in 9.0.1.
*Beweis.* (J_∞a)' = a'∘f·f' + g'a + g a' ⇒ a'∘f·b = (r_a' b − g'ab − g a'b)/f'; ebenso für b; Summe plus (f''∘f)ab = g'ab/f' gibt die
erste Gleichung. Für die zweite: J_∞[c'/(2f')] = (c'∘f)/(2g) + g c'/(2f'), und c∘f = (−ga + r_a)(−gb + r_b) = g²c − g(a r_b + b r_a) + r_a r_b,
c'∘f = (c∘f)'/f'. Einsetzen. ∎ (b40: drei zufällige rationale Beispiele, Differenz exakt 0.)
*Herkunft:* χ = ½ löst χ∘f + χ = 1; das ist der Grund, warum der Ansatz v = c'/(2f') ohne Kompositionsterme aufgeht.

**Korollar B-Ξ₂ (bewiesen).** Für a = E_ζ, b = E_ζ' ist r = J_∞E = O(x^{N+1}) (Definition von Π_N), also T_N[r_a' b] = (N+1)ρ_a b(0)x^N,
die übrigen Restterme verschwinden unter T_N, und
  Ξ^{(2)}(ζ,ζ') = −½J_N⁻¹T_N D²Φ(f)[E_ζ,E_ζ'] = ¼Π_N[(E_ζE_ζ')'/f'] − (N+1)(E_ζ'(0)ρ_ζ + E_ζ(0)ρ_ζ')/(4f'(0))·J_N⁻¹x^N.
Ebenso für d^{(2)}: d^{(2)} = Π_N[u₁u₁'/(2f')] − (N+1)ρ^f u₁(0)/(2f'(0))·t_N (u₁ = x₁ − f_K).

## 9.2 Numerik zur Zerlegung (b41; flint, 800–1400 Bit)

| N | T_N = sup|J⁻¹x^N|/|L|^N | (H*) auf D̄, Eckpol t=2 | (H*_G) äusserer Bogen: Eckpol t=2 / linker Pol / rechter Pol |
|---|---|---|---|
| 60 | 2.2e−2 | 0.610 | 0.649 / 0.034 / 0.747 |
| 100 | 1.7e−3 | 0.588 | 0.634 / 0.060 / 0.746 |
| 160 | 3.5e−5 | 0.568 | 0.603 / 0.255 / 0.746 |
| 220 | 7.3e−7 | 0.555 | 0.552 / 1.36 / 0.745 |
| 300 | 4.2e−9 | 0.543 | 0.568 / 14.6 / 0.745 |
| 400 | 6.6e−12 | 0.531 | 3.47 / 311 / 0.745 |

Bilinearer Kern dist·dist'·sup_D̄|Ξ^{(2)}| (Π_N-Teil = Gesamt auf alle Stellen; Kopfteil ≤ 1.4e−9 überall):

| Paar | N=60 | 100 | 160 | 220 | 300 | 400 | Verhalten |
|---|---|---|---|---|---|---|---|
| L,t=2 × L,t=2 | 1.6e−3 | 2.6e−3 | 4.0e−3 | 5.4e−3 | 7.1e−3 | 9.2e−3 | ≈ N^{0.92} |
| L,t=0.5 × L,t=0.5 | 6.5e−4 | 9.2e−4 | 1.3e−3 | 1.6e−3 | 2.1e−3 | 2.6e−3 | ≈ N^{0.73} |
| L × L̄ (t=2) | 4.9e−5 | 4.6e−5 | 4.1e−5 | 4.1e−5 | 3.9e−5 | 3.8e−5 | fallend |
| L × links | 1.3e−4 | 2.0e−4 | 4.2e−4 | 8.6e−4 | 2.2e−3 | 7.3e−3 | 1.012^N |
| links × links | 1.0e−5 | 1.6e−5 | 4.5e−5 | 1.4e−4 | 7.0e−4 | 5.7e−3 | 1.021^N |
| rechts × L | 1.7e−5 | 7.2e−6 | 2.0e−6 | 5.4e−7 | 9.6e−8 | 1.1e−8 | fallend |

**b42 (x₁ = Π_N f_K auf dem Halbmond; f_K dort exakt als exp∘f_K∘Log):** sup|x₁ − f_K| auf |z| = |L| bzw. auf dem äusseren Bogen
(Maximum bei a₀+R₀): N = 100: 5.7e−10 / 8.3e−10; 200: 1.8e−11 / 5.2e−11; 300: 2.5e−12 / 1.5e−11; 400: 6.2e−13 / 4.9e−11;
500: 2.1e−13 / 1.4e−10 (Rundung der Daten ≤ 1e−15 geprüft). Auf D̄ Konvergenz ~N^{−4.9}, auf dem Halbmond ab N ≈ 300 Wachstum ≈ 1.01^N.

## 9.3 Was folgt — und was nicht (Aufgabe 1)

**(a) Der Cauchy-Weg und seine Grenze.** |Π_N v(z)| ≤ (1/2π)∮_Γ|v(η)||Π_N k_η(z)||dη| mit Γ um Ḡ' (Pol-Antworten sind nur ausserhalb Ḡ'
kontrolliert; innen wachsen sie, Cs Gesetz R). Für v = (E_ζE_ζ')'/f' braucht das E_ζ(η) für η auf Γ, also Π_N k_ζ **auf dem Halbmond**.
*Satz (Skizze, Zählung grob):* Gilt (H*_G): sup_{z∈Ḡ'}|Π_N k_ζ(z)| ≤ C N^κ/dist(ζ,Ḡ') für ζ im Band, dann (H*₂) für Eckpaare mit κ₂ ≤ 3κ + 2
(Γ auf Niveau g = 2/N, Bernstein–Walsh-Faktor e², je Ableitung ein Cauchy-Faktor ≲ N, Länge O(1)). Bilanz: d^{(2)} ≪ ε_N verlangt
κ₂ < κ + 4.69, also κ < 1.35. Satz B-Q verlangt κ̃ = (3κ+2)/2 < 4.69.
**Aber (numerisch):** (H*_G) ist falsch — siehe Tabelle 9.2 (linker Pol ab N ≈ 150, Eckpol ab N ≈ 300 exponentiell), passend zu b42.
Die Spitze von J_N⁻¹ ist auf D̄ unsichtbar, aber ausserhalb D̄ wird sie wie (|z|/|L|)^k verstärkt. **(H*) ⇒ (H*₂) über Cauchy um Ḡ'
ist damit nicht verfügbar**; (H*) selbst (auf D̄) bleibt gut (0.53 bei N = 400, fallend).
**(b) Was stattdessen nötig ist (präzise, linear in Π_N):**
  (H*_Π) Für Eckpole ζ, ζ' (feste Umgebung U von L, L̄, ausserhalb Ḡ', g ≤ η₀) gilt
         |Π_N[(E_ζE_ζ')'/f'](z)| ≤ C N^{κ₂} e^{−N(g(ζ)+g(ζ')−2/N)}/(dist(ζ)dist(ζ')), z ∈ D̄;
  für Pole ausserhalb U dasselbe mit Zusatzfaktor e^{cN}, c < η₀ (Wahl des fernen Konturniveaus η₀; Strang A: f_K beschränkt auf Γ_{0.1}).
  Plus (T'): der Kopfterm |κ_top|·sup_D̄|t_N| ≤ C N^τ/(dist·dist') (numerisch ≤ 1.4e−9 an allen Stellen; die Kopfantwort selbst
  fällt wie 0.937^N, der Koeffizient ρ_ζ|L|^{N+1} wächst, das Produkt mit E_ζ(0) bleibt winzig).
Nach Korollar B-Ξ₂ ist (H*_Π) + (T') ⇔ (H*₂^{loc}). Die Aussage spricht über **denselben** linearen Operator Π_N wie (H*), nur mit
Eingaben „Produkt zweier Pol-Fehler, abgeleitet“ statt „Cauchy-Kern“. Numerisch: κ₂ ≈ 0.92, c ≈ 0.012 (ein ferner Pol) bzw. 0.021
(zwei), also c < η₀ = 0.1 mit grossem Abstand.
**(c) Korrektur zu Satz B-Q (§8.2).** Die dort angesetzte (H*_n) „gleichmässig polynomial im ganzen Band“ ist für n = 2 numerisch
falsch (Paare mit linkem Pol wachsen ~1.012^N). Satz B-Q bleibt richtig mit (H*_n^{loc}): polynomial mit lokalem Bernstein–Walsh-
Abfall e^{−N(g−1/N)} für Variablen in U, Zusatzfaktor e^{cN} je Variable ausserhalb U, c < η₀. *Beweisänderung:* Die Kontur ist
Γ_{η₀} ausserhalb U und die Schlitzufer in U. Variablen auf Γ_{η₀} geben je e^{(c−η₀)N}|Γ|max|f_K| → 0, Variablen auf den Ufern je
C N^{κ−4.69} wie bisher. Also |d^{(n)}| ≤ (C'N^{κ−4.69} + C''e^{(c−η₀)N})ⁿ, Rest wie gehabt.
**(d) Exponenten-Bilanz (Vermutung, gestützt).** Gemessen κ₂ ≈ 0.92 und c_N = sup|x₁ − p_N|/ε_N² ~ N^{0.95} (§8.4). Beides ≈ α = 0.883:
**eine Ableitung auf der Eckskala N^{−α} je Ordnung** (s₂ = aa'/(2f') hat genau eine Ableitung). Erwartet daher κ_n ≈ nκ + (n−1)α.
Die Reihe in B-Q hat dann das Verhältnis N^{κ+α−4.69} ≈ N^{−3.8}.
**Test an echten Daten (b44; x₁ von Strang A, p_N aus pN_N100_bis_1000):** d^{(2)} nach Korollar B-Ξ₂ (direkt und zerlegt),
Rest := p_N − x₁ − d^{(2)} ≈ d^{(3)}:

| N | ε_N | sup|p_N − x₁| | sup|d^{(2)}| | sup|Rest| | c₂ = |d²|/ε² | c₃ = |Rest|/ε³ |
|---|---|---|---|---|---|---|
| 100 | 5.71e−10 | 2.961e−20 | 2.961e−20 | 9.8e−30 | 0.091 | 0.053 |
| 200 | 1.89e−11 | 6.546e−23 | 6.546e−23 | 1.3e−33 | 0.184 | 0.196 |
| 300 | 2.64e−12 | 1.826e−24 | 1.826e−24 | 7.3e−36 | 0.261 | 0.394 |
| 400 | 6.49e−13 | 1.407e−25 | 1.407e−25 | 1.8e−37 | 0.334 | 0.642 |
| 500 | 2.22e−13 | 2.027e−26 | 2.027e−26 | 1.1e−38 | 0.410 | 0.980 |

c₂ ~ N^{0.94}, c₃ ~ N^{1.82}: **genau eine Potenz N^{≈0.9} je Ordnung**, wie in 9.3(d) vorhergesagt (α = 0.883). Die Entwicklung
x₁ + d^{(2)} trifft p_N auf 12 Stellen. Der Kopfanteil von d^{(2)} ist ≤ 1e−25 (N ≤ 300); für N = 400 und 500 ist er aus den 60-stelligen
Daten nicht auflösbar (ρ^f verliert Stellen), c₂ und c₃ stammen dort aus der direkten Formel.

## 9.4 Höhere Ordnungen: formale Normalform (Skizze)

**Lemma 9.1 (bewiesen).** Ist J_∞a = 0, so ist W := a/(2f') exp-invariant: W∘exp = eˣ·W. Folglich vertauscht der Fluss τ_t von W mit exp,
und p_t := τ_t⁻¹∘f∘τ_t ist eine exakte Lösung von p∘p = exp mit ∂_t p_t|₀ = f'W − W∘f = a.
*Beweis.* a∘f = −ga ⇒ a∘f∘f = (g∘f)·g·a, f'∘f∘f = g∘f, also W∘exp = g a/2; eˣ = (f∘f)' = g f' ⇒ eˣW = g a/2. Für die Ableitung:
W∘f = (a∘f)/(2g) = −a/2. ∎ (Die Kernelemente sind also genau die Tangenten an die Familie der Nicht-Kneser-Halbiterierten, und
diese Familie ist **lokal explizit** in a.)
**Satz B-NF (Skizze).** Es gibt eine formale Reihe S(a) = a + aa'/(2f') + Σ_{k≥3}s_k(a), s_k homogene Differentialpolynome in a mit
festen, bei 0 holomorphen Koeffizienten aus f, so dass für **jeden** Keim a gilt: Φ(f + S(a)) − eˣ = R(a; r_a), r_a = J_∞a, wobei jeder
Term von R einen Faktor r_a^{(i)} (i ≥ 0, keine Komposition) enthält und der lineare Teil r_a ist.
*Konstruktion:* S(a) := Entwicklung von τ_1⁻¹∘f∘τ_1 − f mit der Ersetzungsregel a^{(i)}∘f ↦ C_i(a) (C_0(a) = −ga, C_i aus
(a^{(i−1)}∘f)' = f'·(a^{(i)}∘f)).
*Beweisidee:*
- Für echte Kernelemente ist alles exakt (Lemma 9.1).
- Die Jets von Kernelementen an einem Punkt sind beliebig: a = f'W − W∘f, W = θ(α_e)/α_e' mit 1-periodischem θ, und die Jets in
  α_e(x₀) und α_e(x₀)+½ sind unabhängig. Also gilt die Identität formal in jeder Ordnung.
- Für allgemeines a ist a^{(i)}∘f = C_i(a) + (linear in r_a,…,r_a^{(i)}), und Φ(f + S) enthält nur einfache Kompositionen mit f.
*Kontrolle:* Ordnung 2 ist Satz B-Id. Im Translationsmodell f = x + s (Strang E) ergibt sich s₂ = aa'/2 und s₃ = 0 (direkt nachgerechnet).
**Folgerung (Skizze).** Mit dem Fixpunktansatz a = u₁ + (Π_N − I)(S(a) − a) ist p := f + S(a) = x₁ + Π_N(S(a) − a) ∈ P_N. Wegen
J_∞a = O(x^{N+1}) bleibt in T_N[Φ(p) − eˣ] = T_N R nur ein **Kopfband** (Indizes ≥ N+1−i). Es wird durch einen endlichdimensionalen
Schritt in den Richtungen J_N⁻¹x^{N−i} (Kopfantworten, numerisch exponentiell klein) korrigiert. Jeder Kern n-ter Ordnung ist also Π_N
einer expliziten n-linearen Ableitungsform in Pol-Fehlern plus Kopfband-Terme: **(H*_n) ⇔ (H*_{Π,n}) + (T_Band)**.
**Lücken:** die Konvergenz der formalen Reihe S, die Schranken (H*_{Π,n}) selbst und die Fixpunktnorm. Letztere darf nicht auf dem
Halbmond leben (9.3a); sie muss eine Kernnorm auf D̄ sein.

## 9.5 Zweig-Identifikation (Aufgabe 2)

**Exakte Formeln (bewiesen).**
(i) F_{N+1}(pad p_N) = c_{N+1}x^{N+1}, c_{N+1} := [x^{N+1}](p_N∘p_N) − 1/(N+1)!, denn alle Gleichungen bis Grad N gelten schon.
(ii) x₁(N+1) − pad x₁(N) = −ρ^f_N·t_{N+1}, ρ^f_N := [x^{N+1}]J_∞(x₁(N) − f_K), t_M := J_M(f_K)⁻¹x^M. Beweis:
J_{N+1}(x₁(N+1) − pad x₁(N)) = T_{N+1}J_∞(f − x₁(N)) = −ρ^f_N x^{N+1}, da T_N J_∞(x₁(N) − f) = 0.
(iii) Die Homotopie T_{N+1}[p∘p − eˣ] = σx^{N+1}, σ: c_{N+1} → 0, beginnt in pad p_N. Solange J_{N+1}(p(σ)) invertierbar ist, gilt
p(0) − pad p_N = −c_{N+1}∫₀¹J_{N+1}(p(s))⁻¹x^{N+1}ds. Diese Homotopie ist die präzise Form von „stetig fortgesetzt“; Newton ab pad p_N
ist ihr numerisches Verfahren.

**Satz B-Zweig (bewiesen unter (A)–(C)).** Für M = N+1 sei 𝔊_M(λ,σ,e) := T_M[Φ((1−λ)f + λx₁(M) + e) − eˣ] − σx^M, e ∈ P_M.
(A) [Start] Für ein N₀ ist p_{N₀} = x₁(N₀) + d_{N₀}(1), wobei d_N(λ) die analytische Familie aus Satz B-Q ist.
(B) [Kerne] Für alle N ≥ N₀ ist d_N(λ) analytisch auf |λ| ≤ 1+δ. Die Lösung e_M(λ,σ) von 𝔊_M = 0 mit e(0,0) = 0 ist analytisch auf
    |λ| ≤ 1+δ, |σ| ≤ S_M, mit S_M ≥ sup_{|λ|≤1+δ}|σ_N(λ)|, σ_N(λ) := [x^M](Φ(P_N(λ)) − eˣ), P_N(λ) := (1−λ)f + λx₁(N) + d_N(λ).
(C) [Regularität] J_M(p) ist längs des σ-Pfads e_M(1,σ), σ ∈ [0, c_M], invertierbar (das heisst: der Zweig ist wohldefiniert).
Dann gilt p_N = x₁(N) + d_N(1) = p*_N für alle N ≥ N₀.
*Beweis (Induktion).*
1. P_N(λ) löst T_N[Φ(P_N(λ)) − eˣ] = 0, denn F_λ(d) = T_N[Φ(f + λu₁ + d) − eˣ] (T_N J_∞u₁ = 0).
2. Mit e_N(λ) := λx₁(N) + d_N(λ) − λx₁(M) ∈ P_M ist P_N(λ) = (1−λ)f + λx₁(M) + e_N(λ), also 𝔊_M(λ, σ_N(λ), e_N(λ)) = 0.
3. Für λ → 0 gehen e_N(λ) und σ_N(λ) gegen 0 (σ_N(0) = [x^M](f∘f − eˣ) = 0). Die lokale Eindeutigkeit des Satzes über implizite
   Funktionen bei (0,0,0) (endlichdimensional, D_e𝔊 = J_M(f) invertierbar) gibt e_N(λ) = e_M(λ, σ_N(λ)) für kleine λ.
4. Nach dem Identitätssatz gilt das bis λ = 1, also pad p_N − x₁(M) = e_M(1, c_M).
5. σ ↦ e_M(1,σ) ist ein stetiger Lösungspfad von pad p_N nach x₁(M) + e_M(1,0) = p*_M. Nach (C) ist er die eindeutige Fortsetzung,
   also p_M = p*_M. ∎
**Wozu keine Eindeutigkeitskugel nötig ist:** Die Identifikation läuft über analytische Fortsetzung in (λ,σ). Normen von J_N⁻¹
(mit Eckspitze) kommen nicht vor; gebraucht werden nur kohärente Kerne. (B) folgt aus (H*_n^{loc}) zusammen mit Kernen, in denen
einige Eingaben die Kopfantwort t_M sind. Nach B-Id reduzieren sie sich wieder auf Π_M von expliziten Funktionen plus Kopfterme;
die lineare σ-Antwort ist σ·t_M. **Die Lücke ist genau (H*)/(H*_n^{loc}) + (T)**, dazu (A) und (C).
**(A) praktisch:** Kantorovich mit Zentrum x₁ statt q_N hat η = |J⁻¹F(x₁)| ≈ c_Nε_N² (1.8e−24 bei N = 300, §8.4) statt ε_N. Mit der
gemessenen Spitzenrate β ~ 1.012^N trägt das (Schätzung, Vermutung) bis N von einigen Tausend. Als arb-Zertifikat bei einem festen
N₀ ≤ 300 sollte es machbar sein; Strang A hat bei N = 300 die Gleichheit „Newton ab Π_N f_K“ = „fortgesetzter Zweig“ auf 6e−61.

**Numerik b43 (Zweigschritte, Daten polynome_195_bis_300 und lin_PiN):**

| N → N+1 | |c_{N+1}||L|^{N+1} | T_{N+1} | sup|p_{N+1} − pad p_N| | Abw. zu −c·t_{N+1} (J(f)) | Abw. mit J(pad p_N) |
|---|---|---|---|---|---|
| 200 → 201 | 2.3e−7 | 2.5e−6 | 5.6e−13 | 1.1e−22 | 1.6e−25 |
| 250 → 251 | 5.8e−7 | 9.9e−8 | 5.8e−14 | 2.2e−23 | 2.0e−27 |
| 296 → 297 | 1.4e−5 | 5.1e−9 | 6.9e−14 | 1.1e−22 | 3.4e−27 |
| 299 → 300 | 9.9e−6 | 4.2e−9 | 4.2e−14 | 3.5e−22 | 1.3e−27 |

x₁-Identität (ii) bei 296 → 297 und 299 → 300 bestätigt (Abweichung 1e−22 = Abstand p_N ↔ x₁). Es gilt c_{N+1} = ρ^f_N auf 4 Stellen.
Der Zweigschritt ist rund 60-mal kleiner als ε_N und genau ein Vielfaches der (exponentiell kleinen) Kopfantwort.

## 9.6 Stand der Kette nach Runde 4

- (K) [Literatur] ⟹ Lemma 7.4/7.5.
- (H*) auf D̄ ⟹ [B-Max] Π_N f_K → f_K auf D̄. Auf dem Halbmond gilt das **nicht** (b42).
- Satz B-Id ⟹ Kern 2. Ordnung = Π_N(explizit) + Kopfterm.
- Satz B-NF (Skizze) ⟹ alle Ordnungen ebenso.
- (H*_{Π,n}) + (T) ⟹ (H*_n^{loc}) ⟹ [B-Q] p*_N → f_K.
- (A) + (B) + (C) ⟹ [B-Zweig] p_N = p*_N ⟹ Vermutung.

Offen: (H*), (H*_{Π,n}), die Konvergenz von S, (A) als Zertifikat und (C).

Dateien: b40_identitaet.py, b41_zerlegung.py (b41.log, b41b.log, b41_<N>.json), b42_halbmond.py (b42.log),
b43_zweigschritt.py (b43.log), b44_ordnungen.py (b44.log, b44b.log).

---

# Runde 5 (ab 22:23): konsolidierte Beweisarchitektur — siehe ../ARCHITEKTUR.md

Stand 23:04. Die Architektur aller Stränge steht in `polynom_halbexp/ARCHITEKTUR.md`. Neu von Strang B dort:
- **Satz H1 (§4.2, bew.|(F),(H)):** Unter der Familien-Annahme (H) — Cs (E) für die Linearisierungen an Fast-Lösungen in einer
  Kugel B_N vom Radius ~N^{κ−4.69}, plus Invertierbarkeit und Kopf-Stabilität — ist die Lösungsmenge über dem λ-Pfad von f_K nach
  x₁ + d eine einblättrige Überlagerung. Also hat das System genau eine Lösung in B_N, und sie konvergiert gegen f_K. Benutzt werden nur
  die a-priori-Schranke aus B-Sek, der Satz über implizite Funktionen und die Überlagerungstheorie; keine Operatornormen.
- **Satz H2 (§4.3, bew.|(F),(H),(A)):** Zweig-Identifikation per σ-Homotopie mit der exakten Identität
  c_M(T_M𝒥_{P₀})⁻¹x^M = P₀ − Π^{P₀}_M f.
- **b45:**
  - Bs „(H*) mit κ = 0“ (gemessen mit Kopf) und Cs „Spitze ab N ≈ 420“ (ohne Kopf) widersprechen sich nicht. Der kopffreie
    Anteil beim linken Pol wächst wie 1.0122^N: sup_D̄|y_ζ| = 1.5, 4.6, 15.6 bei N = 200, 300, 400.
  - **Korrektur zu §9.0 Punkt 3 / BRETT 22:22:** Π_N f_K konvergiert im Inneren des Halbmonds (x = 1.6, 1.8 bis N = 500), nur
    nicht gleichmässig bis zum äusseren Rand a₀ + R₀.
- **b46:** Es Sekanten-Abbildung T(P) = Π^P_N f_K hat bis N = 400 auch in Eckmoden-Richtung einen Lipschitz-Quotienten von etwa
  ε_N·(1…10). Das Risiko „feste Kugel“ ist damit bis N = 400 ausgeräumt (ARCHITEKTUR §6.5).
