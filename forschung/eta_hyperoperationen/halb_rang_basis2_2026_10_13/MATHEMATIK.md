# Halbrang H₃.₅(2,x): Herleitungen, Auswahl, Beweise und Grenzen

Stand 13.10.2026. Konvention wie in den Vorrunden: H₃(b,x)=bˣ, Rangschritt T(f)=S mit S(0)=1, S(x+1)=f(S(x)), H_{r+1}=T(H_r). „Bewiesen“ heisst unten: unter den genannten Voraussetzungen vollständig hergeleitet. Alles andere ist als **numerisch** oder **offen** markiert.

## 1. Die Familie und der Jet

Bei Basis 2 ist der Successor σ(x)=x+1 ein Fixpunkt von T. Wir schreiben

\[
f_\varepsilon(x)=x+1+\varepsilon\,d(x),\qquad d(x)=2^x-x-1,
\]

also f₀=σ und f₁=2ˣ. Für ganzzahligen Rang n und für den Halbrang wird der Jet

\[
T^s(f_\varepsilon)(x)=x+1+\sum_{j\ge1}\varepsilon^j c_s^{(j)}(x)
\]

betrachtet. Gesucht ist c_{1/2}^{(j)}, und dann die Auswertung bei ε=1.

Für einen Jet g=x+1+δ, δ=Σε^i g_i, lautet die Superfunktionsgleichung für S=x+1+e:

\[
e(x+1)-e(x)=\delta\bigl(x+1+e(x)\bigr)
=\sum_{r\ge0}\frac{\delta^{(r)}(x+1)}{r!}\,e(x)^r,\qquad e(0)=0. \tag{1}
\]

In Ordnung ε^j ergibt das

\[
\Delta e_j=G_j:=\sum_{i=1}^{j}\sum_{r=0}^{j-i}\frac{g_i^{(r)}(x+1)}{r!}\,[\varepsilon^{j-i}]\,e^r. \tag{2}
\]

G_j enthält nur e_1,…,e_{j−1}. Für j=1,2 und g=f_ε sind das die bekannten Gleichungen Δw=v(x+1), Δz=w·v′(x+1) der Vorrunde; (2) ist ihre Verallgemeinerung auf alle Ordnungen und auf Eingaben mit beliebig vielen Ordnungen (nötig ab Rang 2).

## 2. Die Auswahl (ausdrücklich eine Wahl)

**Höhenklasse 𝒜.** 𝒜 besteht aus den Funktionen, die auf einer Halbebene Re x>−ρ (ρ>0) holomorph sind und dort

\[
|h(x+iy)|\le C\,e^{a\max(x,0)}(1+|x|+|y|)^p
\]

erfüllen (endlicher horizontaler Exponentialtyp, polynomiales vertikales Wachstum). Das ist die Halbebenenklasse aus `halb_rang_analytisch_2026_10_12`, §2.

**Lemma A (bewiesen).** Ist G∈𝒜, so hat Δe=G höchstens eine Lösung e∈𝒜 mit e(0)=0.

*Beweis.* Die Differenz h zweier Lösungen ist 1-periodisch, setzt sich also zu einer ganzen periodischen Funktion fort, und erfüllt die Wachstumsschranke auf jeder vertikalen Geraden. Ihre Fourierkoeffizienten a_k=∫₀¹h(x+iy)e^{−2πik(x+iy)}dx hängen nicht von y ab und sind höchstens C′(1+|y|)^p e^{2πky}. Für k≥1 lässt man y→−∞, für k≤−1 y→+∞; also a_k=0 für k≠0, h ist konstant, und h(0)=0. ∎

**Lemma B (bewiesen).** Exponentialpolynome E=Σ_{m=0}^{M}2^{mx}P_m(x) (P_m Polynome) liegen in 𝒜 und sind abgeschlossen unter Verschiebung x↦x+1, Ableitung und Produkt. Zu jedem E gibt es genau eine Lösung von Δe=E in 𝒜 mit e(0)=0, und sie ist wieder ein Exponentialpolynom mit denselben Exponenten m:

- m=0: Σ_k p_k x^k ↦ Σ_k p_k (B_{k+1}(x)−B_{k+1})/(k+1) (Bernoulli-Polynome);
- m≥1: 2^{mx}q(x) mit λq(x+1)−q(x)=P(x), λ=2^m; q ist das eindeutige Polynom gleichen Grades, gefunden durch Rückwärtseinsetzen von der höchsten Potenz an (Diagonalelement λ−1≠0);
- danach eine Konstante, sodass e(0)=0.

*Beweis.* Einsetzen; Eindeutigkeit nach Lemma A. ∎

**Satz C (bewiesen).** Für jedes ganze n≥0 und jedes j≥1 gibt es genau einen formalen Jet von T^n(f_ε), dessen Koeffizienten c_n^{(j)} in 𝒜 liegen und e(0)=0 in jeder Stufe erfüllen. Jeder Koeffizient ist ein Exponentialpolynom

\[
c_n^{(j)}(x)=\sum_{m=0}^{j}2^{mx}P_{n,j,m}(x).
\]

*Beweis.* Induktion über n und j mit (2) und Lemma B. Exponentenschranke: g_i^{(r)}(x+1) hat Exponenten ≤i, [ε^{j−i}]e^r hat Exponenten ≤j−i; Produkte addieren Exponenten. ∎

Damit ist die **Auswahl an ganzen Rängen** festgelegt: polynomiales vertikales Wachstum in jeder ε-Ordnung („perturbative Auswahl“). Für Rang 1 entspricht sie formal der Forderung, die periodische Freiheit θ(x) in jeder Ordnung auszuschliessen. Ob sie mit Kneser übereinstimmt, ist eine Frage, keine Annahme (§5).

**Rangauswahl.** Für festes x (und festes j, r) wird n↦∂_x^r c_n^{(j)}(x) durch die Newtonreihe im Rang

\[
c_s=\sum_{p\ge0}\binom sp\Delta_n^pc_0 \tag{3}
\]

fortgesetzt; das ist die Auswahl „polynomiales Wachstum im Rang“ (Carlson-Typ <π) aus der Vorrunde, (A8). Numerisch wird (3) mit der Levin-u-Transformation der Partialsummen beschleunigt. Die Ableitung nach x vertauscht mit (3), weil (3) linear in der Folge ist.

**Satz D (bewiesen).** Für jede lineare Rangfortsetzung, die konstante Folgen konstant fortsetzt (also auch (3)), gilt H₃.₅(2,1)=2 und H₃.₅(2,2)=4 exakt, in jeder Ordnung.

*Beweis.* Für jede Funktion g mit g(1)=2 folgt T(g)(1)=g(T(g)(0))=g(1)=2; induktiv T^n f_ε(1)=2 für alle n, also c_n^{(j)}(1)=0. Ebenso T^n f_ε(2)=T^{n−1}f_ε(T^nf_ε(1))=T^{n−1}f_ε(2)=…=f_ε(2)=3+ε, also c_n^{(1)}(2)=1 und c_n^{(j)}(2)=0 für j≥2. Konstante Folgen bleiben konstant. Bei ε=1: 2 und 4. ∎

Die Bedingung prüft also die Auswahl nicht; sie ist in der Konstruktion eingebaut.

**Ordnung 1, geschlossen (übernommen, bewiesen in der Vorrunde (A18)).** c_s^{(1)}=L_sd=B_{s+1}(x)[½·₂F₁(1,x+s+1;s+2;½)−1]. Sie dient hier als Kontrolle der Numerik.

## 3. Prüfidentitäten

**(P1) Superfunktionsgleichung an ganzen Höhen.** Mit (1) und den Taylordaten von c_{n−1} bei k+1 folgt c_n(k+1) aus c_n(k). `pruefe_ganz.py` prüft das für k=1…6, n≤24, j≤10.

**(P2) Rangrekursion an beliebigen Höhen.** Mit U=T^{1/2}f, W=T^{3/2}f muss W(x+1)=U(W(x)) gelten. Mit w=W−(x+1):

\[
W(x+1)-(x+2)=w(x)+\sum_{r\ge0}\frac{u^{(r)}(x+1)}{r!}w(x)^r ,
\]

ausgewertet als ε-Reihe (`auswertung.py`, `rekursions_residuum`). Das ist die Jet-Form von T(H₃.₅)=H₄.₅. Für ganzzahlige Ränge (U=f, W=Tf) ist sie exakt erfüllt (Kontrolle im Programm, 10⁻⁶⁰).

**(P3) Ganze Höhen.** T(Uf)(k) wird aus den Taylordaten von U an den Höhen 1…k direkt rekursiv berechnet und mit der Fortsetzung von c_n(k) bei s=3/2 verglichen.

Die beiden Halbschritte T^{1/2}(T^{1/2}f) sind in dieser Konstruktion als Rangfortsetzung von n↦T^n(T^{1/2}f) definiert. Gilt T^n(T^{1/2}f)=T^{n+1/2}f (geprüft für n=1 durch P2, P3), so ist T^{1/2}(T^{1/2}f) die Fortsetzung von s↦T^{s+1/2}f bei s=1/2, also T(f). Der eigentliche Test von T^{1/2}∘T^{1/2}=T ist daher P2/P3 zusammen mit dem Vergleich T(f) gegen Kneser.

## 4. Zwei strukturelle Hindernisse

### 4.1 Die ε-Reihe divergiert (numerisch)

Bei Rang 1 wachsen die Glieder c_1^{(j)}(x) ab etwa j≈10 schneller als geometrisch: log₁₀|c_j| steigt pro Ordnung um 0.43 (j≈25) und 0.58 (j≈35); bei j=52 ist |c_j|≈10¹⁸. Das passt zu Gevrey-1-Wachstum |c_j|≈C·j!·A^{−j} mit A grob 7–9; ein sauberer Fit von A ist wegen der Vorzeichenwechsel nicht möglich. Die geometrische Alternative (endlicher Konvergenzradius) passt schlechter, weil der Quotient |c_j/c_{j−1}| im Mittel weiter steigt (≈2.6 bei j≈22, ≈3.7–4 bei j≈36).

*Heuristischer Grund (kein Beweis).* Die Fixpunkte von f_ε erfüllen 2^L−L−1=−1/ε. Für ε→0⁺ wandern sie nach Re L≈log₂(1/ε)→∞, also nicht analytisch in ε. Für jedes ε<0 dagegen hat f_ε **zwei reelle** Fixpunkte (bewiesen: 2ˣ−x−1 ist stetig, null bei 0 und strebt nach ±∞ gegen +∞, nimmt also jeden positiven Wert 1/|ε| zweimal an). Die Dynamik ändert bei ε=0 ihren Typ. Eine konvergente Reihe um ε=0 ist deshalb nicht zu erwarten.

Borel–Padé ([24/24] bis [26/26], Laplace längs arg t=0, ±0.15) stabilisiert sich nicht: Streuung 10⁻⁴–10⁻³, die seitlichen Summen haben Imaginärteile bis 10⁻³. Die Padé-Pole der Borel-Transformierten liegen etwa auf einem Kreis |t|≈8 in vielen Richtungen. Eine eindeutige Borel-Summe auf der positiven Achse ist damit **nicht** belegt.

Folge: Der Jet bestimmt die Funktion bei ε=1 nur bis auf eine Unschärfe von der Grössenordnung des kleinsten Glieds bzw. der Streuung der Partialsummen, auf dem reellen Abschnitt etwa 10⁻⁴.

### 4.2 Exponentielles Wachstum im Rang ab Ordnung 8 (numerisch)

Für festes x wachsen die exakten Koeffizienten c_n^{(j)}(x) ab j≈8 exponentiell in n (Ballarithmetik, Radius <10⁻³⁰⁰; bei 2400 und 4000 Bit identisch; die Superfunktionsgleichung P1 gilt). Fit log|c|=a+n·log|ρ|+k·log n über n=24…40:

| j | ρ_j bei x=0.5 | ρ_j bei x=1.5 | Vorzeichen |
|---:|---:|---:|---|
| 8 | ≈ −7.6 (Quotient, noch driftend) | ≈ −7.6 | alternierend |
| 9 | Mischung, geom. Mittel ≈ 28 | ≈ 28 | gemischt |
| 10 | 117.33 | 117.33 | positiv |
| 11 | Mischung, geom. Mittel ≈ 470 | ≈ 470 | gemischt |
| 12 | −2667.3 | −2667.3 | alternierend |

Die Raten hängen nicht von x ab (auch x=0.25 gleich), also nicht von der Höhe, sondern von der Ordnung. Ordnungen 6 und 7 beginnen bei n≈30 zu steigen (Quotienten 1.3–2.4, fallend); Ordnungen 1–5 bleiben bis n=40 zahm.

**Folgerung (bewiesen, elementar).** Für a_n=ρ^n ist Δ^pa_0=(ρ−1)^p; da |binom(s,p)|≍p^{−Re s−1} für nichtganzes s, divergiert (3), sobald |ρ−1|>1. Ausserdem setzen e^{s log ρ} und e^{s(log ρ+2πik)} dieselben ganzzahligen Werte fort; für |log|ρ||≥π gibt es also keine durch den Wachstumstyp (<π) eindeutige Fortsetzung, und bei ρ<0 ist schon die natürliche Fortsetzung ρ^{1/2} nicht reell. Die termweise Rangfortsetzung ist deshalb ab Ordnung 8 **nicht definiert**. Die Ursache des Wachstums ist offen; eine Vermutung ist, dass die effektive Borel-Weite der ε-Reihe mit jedem Rang schrumpft (dann wüchse c_n^{(j)} wie q^{jn}, Raten ≈ q^j; beobachtet wächst |ρ_j| um den Faktor 4–5 pro Ordnung).

## 5. Was die Konstruktion über Kneser sagt

Satz C legt T(f_ε) als formalen Jet eindeutig fest. Kneser-Tetration ist nicht in 𝒜 (Vorrunde: keine vollständige Tetration hat endlichen Exponentialtyp bei b>e^{1/e}), und der Jet ist divergent. Ein Satz „Kneser = perturbative Auswahl“ kann daher nur in der Form „Kneser ist eine (Borel-)Summe der Reihe“ stehen. Numerisch liegt Kneser im Band der Partialsummen 4…12 an allen reellen Testpunkten (Tabelle in BERICHT.md). Andere Tetrationen, die sich von Kneser um eine 1-periodische Störung θ von weniger als ~10⁻⁴ unterscheiden, lassen sich mit diesem Werkzeug nicht ausschliessen; die lineare Näherung (sexp(x)=x+1 auf (−1,0], bei 0.5 um 0.045 von Kneser verschieden) ist klar ausgeschlossen.

## 6. Offene Beweisfragen

1. Gevrey-Ordnung und Lage der Borel-Singularitäten des Rang-1-Jets beweisen (Vermutung: Gevrey-1, A in der Nähe von 8–9).
2. Ursache und genaue Form der Rangraten ρ_j; insbesondere ob sie auf eine Borel-Weite A_n mit A_{n+1}/A_n→const zurückgehen.
3. Eine Rangfortsetzung, die nicht termweise in ε arbeitet (etwa im Borel-Raum oder über eine Koordinate, in der T linear wird), und ein Beweis, dass sie mit (3) für die Ordnungen 1–7 übereinstimmt.
4. Ein Kriterium, das Kneser unter den Summen der divergenten Reihe auszeichnet.
