# Lücken S3, S4, S8 und die Z0-Majorante (14.10.2026)

## Kurz gesagt

- **Der Schlüssel:** Die Schröder-Funktion σ am Fixpunkt L ist auf der ganzen oberen Halbebene ℍ **schlicht** (umkehrbar eindeutig). Grund: Der Hauptzweig-Logarithmus bildet ℍ umkehrbar eindeutig in sich ab und zieht alles zu L. Für schlichte Funktionen gibt es klassische, scharfe Schranken (Koebe, Grunsky). Man wendet sie nicht direkt bei z an, sondern beim **dritten iterierten Logarithmus** u₃ = log log log sexp z. Dieser Punkt liegt immer mindestens π/2 über der reellen Achse, und dort sind die Schranken gut.
- **Majorante (Ziel 1, Lücke von kn.satz.theta_mittelwerte): geschlossen.** Auf dem Halbkreis D = {|z| < 0.1, Im z > 0} gilt (Basis e)
  - |θ(z)| ≤ 31.7 + 2.91·ln W(z), wobei W(z) ≈ ln ln(1/|z|). Also wächst |θ| höchstens wie **2.91·ln ln ln(1/|z|)**.
  - |1+θ′(z)| ≤ 2.24 / (|z|·ln(1/(1.13|z|))) (Gitterwert; rein analytisch ≤ 19 statt 2.24),
  - |1+θ′(z)| ≥ c/(|z|·ln(1/|z|)·(ln ln(1/|z|))⁴), und arg(1+θ′) wächst höchstens wie 2·ln ln ln(1/|z|).
  
  Wegen der Periode gilt dasselbe an jeder ganzen Zahl. Damit sind Z1–Z5 und V3 (Etappe 3) vollständig bewiesen, ebenso die Majoranten-Teile von Etappe 4. Eine kleine Präzisierung: „Φ(θ) für jede ganze Funktion Φ“ gilt nur für Φ mit log M_Φ(t) = O(e^{ct}). Dazu gehören alle tatsächlich benutzten Φ, auch e^{e^w}.
- **S3 (zitiert): nicht mehr nötig.** Ersetzt durch Lemma U (bewiesen): σ ist schlicht auf ℍ, h′ = σ′/(Lσ) ist auf ℍ eindeutig, und σ setzt sich über jeden reellen Punkt ausserhalb von {0, 1, e, eᵉ, …} analytisch fort.
- **S4 (Zweigwahl): auf einen einzigen Punkt zurückgeführt.** Stimmen die beiden Zweige an *einem* Punkt des Halbkreises überein, so stimmen sie überall überein (Eindeutigkeit von Hochhebungen). Das Zertifikat bei z₀ = 0.05i: Abstand 4·10⁻⁴⁸, nötig wäre nur < 0.12. Ergebnis: Der **Hauptzweig** ist richtig, ohne Zusatzkonstante (Zweigzahl k = 0 an vier Punkten). Daraus folgt: θ ist auf ℝ∖ℤ analytisch, und exakt gilt θ(z) = 2 − z + h(log log sexp z) auf D.
- **S8 (Transfer): bewiesen**, sogar etwas schärfer als behauptet:
  - aₙ = bₙ + O(1/(n² ln n)), mit bₙ = [qⁿ] h((c + log(1−q))/ln b),
  - mit der Korrektur δ₁ gilt aₙ = b2ₙ + O(1/(n³ ln n)),
  - führend gilt aₙ = −h′(vₙ)/(n ln b) + O(1/(n ln² n)) mit vₙ = (c − ln n)/ln b. Mit K Gliedern der Transfer-Reihe ist der Rest O(1/(n ln^{K+1} n)).
  - Daraus folgt c/(n ln n (ln ln n)⁴) ≤ |aₙ| ≤ C/(n ln n). Der **Exponent ist also genau 1** (bewiesen).
  
  Die Daten bestätigen das: n² ln n·|aₙ−bₙ| = 0.019 → 0.012, n³ ln n·|aₙ−b2ₙ| = 0.024 → 0.012, n ln²n·|bₙ − Leitglied| ≈ 0.14 (n = 16…4000, Basis e und 2 gleich).
- **Nebenbefund zu Etappe 4 (B-S1):** Die Gewichte F_C enthalten w = 1+θ′, und w ist **nicht** integrierbar (∫|θ′| = ∞, weil θ unbeschränkt ist). F_d = θ″/w verhält sich wie −1/z. Die Aussage von S1 bleibt richtig, wenn man die Koeffizienten regularisiert: F̃_C(j) = M[F_C·(e^{−2πijx} − 1)]. Das ist erlaubt, weil p(0) = 0. Siehe 6.
- **Was bleibt:** (1) Das Ein-Punkt-Zertifikat für S4 ist in Gleitkomma gerechnet (50 Stellen, 46 Grössenordnungen Spielraum), nicht in Intervallen. (2) Die Konstanten in den O-Termen des Transfers existieren, sind aber nicht ausgeschrieben. (3) Kneser-Darstellung und Paulsen–Cowgill bleiben Voraussetzung wie bisher. (4) Die genaue Grössenordnung |aₙ| ≍ 1/(n ln n ln ln n ⋯) bleibt heuristisch. Bewiesen ist nur der Bereich oben, unten mit (ln ln n)⁴ Verlust.

## 1. Voraussetzungen und Bezeichnungen

Basis e (allgemeine Basis b > η: Abschnitt 7). L = Fixpunkt von e^w mit Im L > 0, L = 0.3181 + 1.3372i, |L| = 1.3745. Log ist der Hauptzweig. log₊ ist der Zweig mit arg ∈ (0, 2π), holomorph auf ℂ∖[0,∞).

Vorausgesetzt (wie in Etappe 2–4):
- (V1) Kneser: sexp ist holomorph auf ℂ∖(−∞,−2], reell und mit sexp′ > 0 auf (−2,∞). Es gilt sexp(z) = F(z+θ(z)) für z ∈ ℍ, mit θ holomorph und 1-periodisch in ℍ und θ → a₀ für Im z → ∞.
- (V2) Ψ ist ganz, mit Ψ(Lu) = e^{Ψ(u)}, Ψ(0) = L, Ψ′(0) = 1, F(s) = Ψ(e^{Ls}). Nach S1 aus Etappe 2 gilt Ψ′ ≠ 0, Ψ ist also lokal umkehrbar.

Numerik: Taylor-Daten von sexp bei 0 (`etappe3/daten/taylor_130_256.json`, 130 Stellen) und θ-Koeffizienten n ≤ 4096 (Etappe 2, 44 Stellen). Kein Motorlauf.

## 2. Lemma U (ersetzt S3)

**Lemma U.**
- (a) Logⁿ(w) → L lokal gleichmässig auf ℍ.
- (b) σ(w) := lim Lⁿ(Logⁿ w − L) existiert lokal gleichmässig auf ℍ, ist holomorph und **schlicht** auf ℍ, mit σ(Log w) = σ(w)/L, σ(L) = 0 und σ′(L) = 1. Es ist dieselbe Funktion wie `lib3.Basis.sigma`.
- (c) Ψ∘σ = id auf ℍ. σ hat in ℍ nur die Nullstelle L.
- (d) h′ := σ′/(Lσ) ist eindeutig und meromorph auf ℍ, mit einem einzigen Pol bei L. Es gilt h′(w) = h′(Log w)/w. Jeder Zweig h von (1/L)·log σ erfüllt h(w) = 1 + h(Log w) mod 2πi/L.
- (e) σ setzt sich von oben über jeden Punkt x ∈ ℝ∖O, O = {0, 1, e, eᵉ, …}, holomorph fort, mit σ ≠ 0 und σ′ ≠ 0. Für Re v < 0 ist die Fortsetzung über die negative Achse σ(v) = L·σ(log₊ v).

*Beweis.*
- (a) Log bildet ℍ in den Streifen {0 < Im < π} ⊂ ℍ ab, ist also eine holomorphe Selbstabbildung von ℍ, aber kein Automorphismus. Ausserdem ist Log L = L, weil Im L ∈ (0,π). Überträgt man das auf die Kreisscheibe mit L ↦ 0, so gilt nach Schwarz |g(ζ)| ≤ c_R|ζ| auf |ζ| ≤ R mit c_R < 1. Daraus folgt (a).
- (b) Das ist Koenigs' Linearisierung im Becken (Multiplikator 1/L, |1/L| = 0.73). Jede Funktion Lⁿ(Logⁿ − L) ist schlicht, weil Log auf ℂ∖(−∞,0] injektiv ist. Der nichtkonstante Grenzwert ist nach Hurwitz schlicht.
- (c) Ψ(σ(w)) = Ψ(Lⁿσ(uₙ)) = expⁿ(Ψ(σ(uₙ))) = expⁿ(uₙ) = w, mit uₙ = Logⁿ w nahe L, wo Ψ die lokale Umkehrung von σ ist. Gilt σ(w) = 0, so ist σ(uₙ) = 0, also uₙ = L und damit w = L, denn Log ist injektiv.
- (d) Leite σ(w) = Lσ(Log w) ab.
- (e) Sei x ∈ ℝ∖O, x > 0, und x₀ = x, x_{k+1} = ln x_k. Wegen ln y < y−1 gibt es ein erstes m mit x_m < 0 (x_m = 0 hiesse x ∈ O). w ↦ Log^m w ist nahe x holomorph und bildet ℍ in ℍ ab. Setze σ(w) = L^{m+1}σ(log₊ Log^m w). log₊(x_m) = ln|x_m| + iπ ∈ ℍ ist ≠ L, weil Im L ≠ π. Die Zusammensetzung aus lokal schlichten Abbildungen gibt σ′ ≠ 0. Für x < 0 ist m = 0. ∎

Damit wird der zitierte Satz über singuläre Werte von Poincaré-Funktionen nicht mehr gebraucht. Ob die Fortsetzung an den Punkten von O wirklich verhindert ist, wird hier nicht benutzt. Dass θ an den ganzen Zahlen singulär ist, zeigt S5 ohnehin unabhängig.

## 3. Satz Z (Zweigwahl, S4)

**Halbkreis.** Sei r = 0.1 und D = {|z| < r, Im z > 0}. Rechnergestützt mit den Taylor-Daten gilt: Σ_{k≥2} k|c_k| r^{k−1} = 0.0610 < A = 1.0918 (der Rest ab k = 128 ist ≤ 4·10⁻¹²⁵). Nach Noshiro–Warschawski ist sexp also auf |z| < r schlicht. Ausserdem gilt auf |z| ≤ r:
- sexp z = 1 + Az(1+ε₀) mit |ε₀| ≤ 0.027,
- Log sexp z = Az(1+ε₁) mit |ε₁| ≤ 0.034,
- |sexp′/A − 1| ≤ 0.056.

Folgerungen für z ∈ D:
- sexp(D) ⊂ ℍ. Denn wäre sexp(z) reell, so wäre sexp(z̄) = sexp(z) mit z̄ ≠ z, im Widerspruch zur Schlichtheit. Nahe der Achse ist Im sexp > 0.
- sexp(z−1) = Log sexp z, und sexp(z−2) = Log Log sexp z. Beide Seiten sind auf D stetig, unterscheiden sich um 2πiℤ und stimmen am Rand (0, r) überein, wo alles reell ist.

**Lemma H (Hochhebung).** Auf einer zusammenhängenden Menge D stimmen zwei stetige Funktionen p₁, p₂ mit Ψ∘p₁ = Ψ∘p₂ überall überein, wenn sie an einem Punkt übereinstimmen. Grund: Ψ ist ein lokaler Homöomorphismus (Ψ′ ≠ 0). Die Menge der Übereinstimmung ist also offen und abgeschlossen.

**Satz Z.** Es gilt σ(sexp z) = e^{L(z+θ(z))} auf D und auf einem Streifen N = {x+iy : 0<x<1, 0<y<η(x)} über (0,1). Folgerungen:
- (i) θ(z) = 3 − z + h(u₃(z)) auf D, mit u_k = Log^k sexp z und dem **Hauptzweig** h = `lib3.Basis.h` (Zweigzahl 0).
- (ii) θ ist auf ℝ∖ℤ analytisch.
- (iii) 1 + θ′(z) = sexp′(z)·h′(sexp z) = sexp′(z)·h′(u₃)/(u₀u₁u₂).

*Beweis.* p₁ = e^{L(z+θ)} und p₂ = σ∘sexp sind beide stetige Hochhebungen von sexp|_D: Ψ∘p₁ = sexp nach (V1), Ψ∘p₂ = sexp nach U(c). Nach Lemma H genügt Übereinstimmung an einem Punkt (Zertifikat unten).

N ∪ D ist zusammenhängend, wenn η > 0 so klein ist, dass sexp(N) ⊂ ℍ (sexp′ > 0 auf (0,1)). Also gilt die Identität auch auf N.
- (ii): sexp bildet (0,1) auf (1,e) ab, und dieses Intervall trifft O nicht. Nach U(e) setzt sich log σ(sexp z) über (0,1) fort.
- (i): Zweimal U(d) anwenden. Die Zweigkonstante ist ganzzahlig, auf D konstant, und an vier Punkten numerisch 0 (|k| < 4·10⁻⁴⁹).
- (iii): Ableiten von e^{L(z+θ)} = σ(sexp z) und U(d). ∎

**Zertifikat** (`l1_schranken.py`, Teil B). P = e^{L(z₀+θ(z₀))} wird aus der θ-Reihe berechnet. Trennradius ρ = |σ′(w₀)|·Im w₀/4: Nach Koebe-¼ enthält σ(B(w₀, Im w₀)) die Scheibe B(σ(w₀), ρ). Ψ ist auf σ(ℍ) injektiv. Liegt also P im Abstand < ρ von σ(Ψ(P)), so ist P = σ(Ψ(P)).

| z₀ | \|σ(Ψ(P)) − P\| | ρ | \|Ψ(P) − sexp_Taylor(z₀)\| | Zweigzahl |
|---|---|---|---|---|
| 0.05i | 3.8·10⁻⁴⁸ | 0.121 | 3.4·10⁻⁵⁰ | 0 (±2·10⁻⁴⁹) |
| 0.07e^{iπ/3} | 6.9·10⁻⁴⁸ | 0.116 | 3.2·10⁻⁵⁰ | 0 |
| 0.08e^{3iπ/4} | 1.2·10⁻⁴⁸ | 0.096 | 4.1·10⁻⁵⁰ | 0 |
| 0.5+0.05i | 2.6·10⁻⁴⁸ | 0.038 | – | 0 |

Probe: Bei z = 0.06e^{0.6iπ} stimmen die θ-Reihe und 3 − z + h(u₃) auf 9·10⁻⁴⁹ überein. **Restlücke:** Die Rechnung ist Gleitkomma (dps 50), nicht Intervallarithmetik. Ein Punkt reicht. Ein Intervall-Zertifikat wäre kurz: θ-Reihe bei 0.05i, Ψ, σ. Nebenbei: In Knesers Konstruktion ist genau dieser Hauptzweig eingebaut.

**Fortsetzung in einen Sektor (für den Transfer).** Für α < π/2 sei S_α = {0 < |z| < r, −α < arg z < π+α}. S_α ist einfach zusammenhängend, und u₁ = Az(1+ε₁) ≠ 0. Also ist u₂ := log u₁ (stetiger Zweig) auf S_α definiert, mit Re u₂ < 0 und Im u₂ ∈ (−α−0.04, π+α+0.04). Weiter ist u₃ := log₊ u₂ mit Im u₃ ∈ (π/2, 3π/2). Die Funktion 3 − z + h(u₃) ist auf S_α holomorph und setzt θ fort. Am Rand von u₂ = negative Achse ist das die Fortsetzung aus U(e).

## 4. Satz M (Majorante)

**Werkzeuge aus der Theorie schlichter Funktionen** (f schlicht auf 𝔻, f(0) = 0, |ζ| = ρ):
- (Gr1) |log(ζf′(ζ)/f(ζ))| ≤ log((1+ρ)/(1−ρ)) (Grunsky; daraus kommt der Sternradius tanh(π/4)).
- (Gr2) |log(f(ζ)/(f′(0)ζ)) + log(1−ρ²)| ≤ log((1+ρ)/(1−ρ)) (Grunsky/Goluzin; der Betrag ist Koebes Wachstumssatz).

Literatur: Duren, *Univalent Functions* (1983), Kap. 2–3; Goluzin, *Geometric Theory of Functions*, Kap. IV. Die Logarithmen sind die auf 𝔻 stetigen Zweige mit Wert 0 bei ζ = 0.

**Lemma G (übertragen auf ℍ).** Sei d(w) = d_ℍ(L,w) = arcosh(1 + |w−L|²/(2 Im w·Im L)) und ζ(w) = (w−L)/(w−L̄). Für w ∈ ℍ gilt
- (G1) |log Q(w)| ≤ d(w), mit Q(w) := L·h′(w)·(w−L)(w−L̄)/(L−L̄);
- (G2) |Λ(w) + log(1−|ζ|²)| ≤ d(w), mit Λ(w) := log[σ(w)/(ζ(w)(L−L̄))] und 1−|ζ|² = 4 Im w Im L/|w−L̄|²;
- (G3) e^{−d}K(w) ≤ |h′(w)| ≤ e^{d}K(w), mit K = 2 Im L/(|L||w−L||w−L̄|). Daraus |h′(w)| ≤ 4 Im L/(|L||w−L|²) + 2/(|L| Im w).

Die Logarithmen in (G1) und (G2) sind die stetigen Zweige auf ℍ mit Wert 0 bei w = L.

*Beweis.* Sei m(ζ) = (L − L̄ζ)/(1−ζ). m bildet 𝔻 auf ℍ ab, mit m(0) = L und m′(ζ) = (L−L̄)/(1−ζ)². Die Funktion f = σ∘m ist schlicht (Lemma U), mit f(0) = 0 und f′(0) = L − L̄. m ist eine Isometrie, also d_𝔻(0,ζ) = log((1+|ζ|)/(1−|ζ|)) = d_ℍ(L,w). Aus σ′/σ = f′/(f m′) und ζm′(ζ) = (w−L)(w−L̄)/(L−L̄) folgt ζf′/f = Q. Ausserdem f(ζ)/(f′(0)ζ) = σ/(ζ(L−L̄)). Damit folgen (G1) und (G2) aus (Gr1) und (Gr2).

Für (G3): e^d ≤ 2 cosh d = 2 + |w−L|²/(Im w Im L), und |w−L̄| ≥ |w−L|. ∎

Probe (`l1_schranken.py`, Teil C): An 9 Punkten von w = 1+0.01i bis 60+1.6i sind die Verhältnisse |log Q|/d ≤ 0.38 und |Λ + log(1−|ζ|²)|/d ≤ 0.72. Beide sind ≤ 1, wie es sein muss.

**Die Kette auf D.** Sei u₀ = sexp z, u₁ = Log u₀ = Az(1+ε₁), ℓ₁ = ln(1/|u₁|) und u₂ = Log u₁ = −ℓ₁ + i arg u₁. Wegen |u₁| < 1 und Im u₁ > 0 gilt Re u₂ < 0 und 0 < Im u₂ < π, also arg u₂ ∈ (π/2, π). Für u₃ = Log u₂ folgt:

  **π/2 < Im u₃ < π,  Re u₃ = ln|u₂| ∈ [ln ℓ₁, ln(ℓ₁+π)],  Re u₃ ≥ X₀ = 0.780.**

Setze W(z) := π + |L| + ln(π + ln(1/(A(1−ε₁)|z|))). Dann gilt |u₃ − L| ≤ W und |u₃ − L̄| ≤ W, und mit (G3) und Im u₃ ≥ π/2 ist d(u₃) ≤ ln(2 + 2W²/(π Im L)) ≈ 2 ln W.

**Satz M** (Basis e, z ∈ D, r = 0.1):
- (M1) |θ(z)| ≤ C_θ + (4/|L|)·ln W(z), mit C_θ = 31.7 und 4/|L| = 2.910.
- (M2) |1+θ′(z)| ≤ C₊/(|z|·ln(1/(A(1+ε₁)|z|))). C₊ = 2.24 mit dem Gitter-Supremum H₁ = sup e^dK = 1.82 über dem Bereich von u₃. Rein analytisch gilt C₊ ≤ 19, aus (G3) mit |u₃ − L| ≥ 0.518.
- (M3) |1+θ′(z)| ≥ (1−0.056)·2 Im L / [|L| W²(2 + 2W²/(π Im L))·(1 + Ar(1.027))·(1.034)|z|·(π + ln(1/(A(0.966)|z|)))] ≍ 1/(|z| ln(1/|z|) (ln ln(1/|z|))⁴).
- (M4) Für den stetigen Zweig von log(1+θ′) auf D gilt |arg(1+θ′(z)) − arg(1+θ′(z₁))| ≤ 5π + 0.12 + d(u₃(z)) + d(u₃(z₁)), mit d(u₃(z)) ≤ ln(2 + 2W²/(π Im L)).

Alle Schranken fallen in |z|. Auf der Geraden z = x+iy sind sie also durch ihren Wert bei |z| = |x| beschränkt.

*Beweis.*
- (M1): Nach Satz Z(i) und U(d) ist Lθ(z) = L(3−z) + log(L−L̄) + Log ζ(u₃) + Λ(u₃) + 2πik. Log ζ ist auf u₃(D) stetig: Dort gilt |1−ζ| = 2 Im L/|u₃−L̄| ≤ 0.92, also |Log ζ| ≤ ln(1/0.08). Die Konstante k fällt weg, wenn man die Differenz zu z₁ = 0.05i bildet (|θ(z₁)| = 26.7 aus der Reihe). Nach (G2) gilt |Λ(u₃)| ≤ d + ln(|u₃−L̄|²/(4 Im u₃ Im L)) ≤ ln(2 + 2W²/(πIm L)) + ln(W²/(2π Im L)). Das ist ≤ 4 ln W + konst. C_θ ist das Maximum der Differenz über |z| ∈ (0, r) (`l1_schranken.py`, Teil D).
- (M2) und (M3): Satz Z(iii) mit |u₀| ∈ [1 − Ar(1.027), 1 + Ar(1.027)], |u₁| = A|z|(1±ε₁), ℓ₁ ≤ |u₂| ≤ ℓ₁ + π und (G3).
- (M4): Es gilt log(1+θ′) = log sexp′ + log Q(u₃) + log(L−L̄) − log L − log(u₃−L) − log(u₃−L̄) − Σ_{i<3} Log uᵢ + 2πik. Jeder Summand ausser log Q ist auf D stetig und hat ein Argument in einem Intervall der Länge ≤ π: u₃ − L und u₃ − L̄ liegen in ℍ, weil Im u₃ > π/2 > Im L, und uᵢ ∈ ℍ. Bei sexp′ ist die Länge ≤ 0.06. Für log Q gilt (G1). ∎

Vergleich mit dem wirklichen θ (`l1_schranken.py`, Teil D, längs z = t e^{iπφ}, φ = 0.02 / 0.5 / 0.98):

| t | \|θ\| | Schranke M1 | t·\|1+θ′\| | t·Schranke M2 | arg(1+θ′) | Im u₃ |
|---|---|---|---|---|---|---|
| 0.05 | 26.63–26.83 | 37.06 | 0.093–0.111 | 0.779 | −0.84…0.90 | 2.33–3.12 |
| 10⁻³ | 26.66–26.76 | 37.28 | 0.033–0.035 | 0.330 | −1.09…1.34 | 2.72–3.13 |
| 10⁻⁸ | 26.70–26.73 | 37.59 | 0.0092–0.0095 | 0.122 | −1.10…1.70 | 2.98–3.14 |
| 10⁻²⁰ | 26.74–26.75 | 37.89 | 0.0028 | 0.049 | −1.00…1.93 | 3.08–3.14 |
| 10⁻⁴⁰ | 26.77 | 38.11 | 0.00115 | 0.024 | −0.93…2.05 | 3.11–3.14 |

Die Schranken halten mit grossem Abstand. θ wächst in Wirklichkeit viel langsamer (wie ein Superlogarithmus), aber für die Integrierbarkeit genügt ln ln ln.

**Satz Z0′ (Mittelwertsatz, jetzt vollständig).** G sei holomorph und 1-periodisch in ℍ, stetig bis ℝ∖ℤ, mit G → G_∞ oben. Auf D gelte |G(z)| ≤ g(|z|) mit fallendem g ∈ L¹(0,r). Dann ist ∫₀¹ G(x+i0) dx = G_∞.

*Beweis.* Für y < r/2 ist |G(x+iy)| ≤ g(dist(x,ℤ)) auf dist(x,ℤ) < r/2. Auf [r/2, 1−r/2] × [0, r/2] ist G stetig, also beschränkt. Dominierte Konvergenz und d/dy∫₀¹G = 0 (Etappe 3, B-Z). ∎

Zulässig sind, nach (M1)–(M4) mit t = |z|, ℓ = ln(1/t) und W ≈ ln ℓ:
- **Φ(θ)** für ganze Φ mit log M_Φ(s) ≤ K e^{κs}. Dann ist g(t) ≤ exp(K′ W^{2.91κ}) ≤ t^{−1/2} für kleine t. Das deckt θ², θ³, e^{cθ}, e^{kLθ} (Z1, Z2) und e^{e^θ}. Für ganz beliebige Φ ist Z0 auf diesem Weg **nicht** gezeigt. Man könnte Φ entlang der Werte von θ beliebig schnell wachsen lassen.
- **(1+θ′)^s, Re s < 1** (Z3): |w^s| ≤ |w|^{Re s} e^{|Im s|·|arg w|}. Mit (M2)–(M4) ist g(t) ≲ t^{−Re s} ℓ^{−Re s} W^{2|Im s|} für 0 ≤ Re s < 1. Für Re s < 0 ist g(t) ≲ (tℓW⁴)^{|Re s|}·W^{2|Im s|}, beschränkt.
- **log(1+θ′)** (Z4): g(t) ≲ ln(1/t), integrierbar.
- **Fourier-Moden** G·e^{2πimz} (Z5): |e^{2πimz}| ≤ 1.
- **V3 (Binet):** E(z) = log sexp′(z) − Lz − Lθ(z). log sexp′ ist bei 0 holomorph, für θ gilt (M1).

## 5. Satz T (Transfer, S8/S9)

Basis e; q = e^{2πiz}, u = 1−q, c = ln(A/2π) + iπ/2, δ₁ = ½ + iκ/2π, H(q) = h(c + log u), bₙ = [qⁿ]H, b2ₙ = bₙ − δ₁[(n+1)b_{n+1} − 2n bₙ + (n−1)b_{n−1}] = [qⁿ](H − δ₁u²H′), vₙ = c − ln n.

**(T1) Δ-Analytizität.** θ (als Funktion von q) und H sind analytisch in Δ = {|q| < 1+η, q ≠ 1, |arg(q−1)| > φ₀}, für jedes φ₀ > 0 und ein η = η(φ₀) > 0.
- θ: Nach Satz Z(ii) ist θ analytisch auf ℝ∖ℤ, also auf einer Umgebung von |q| = 1 ohne die Stelle 1. Nahe q = 1 entspricht Δ dem Sektor S_α mit α = π/2 − φ₀ + o(1). Dort ist θ nach Abschnitt 3 holomorph.
- H: Für |q| ≤ 1, q ≠ 1, ist arg u ∈ (−π/2, π/2), also liegt c + log u in {0 < Im < π, Re ≤ ln(2A/2π) = −1.06}. Das ist in ℍ und fern von L (Re L = 0.32). Insbesondere ist H über |q| = 1 ohne die 1 analytisch. Nahe q = 1 in Δ gilt Re(c + log u) ≪ 0. Dort ist h(v) = 1 + h(log₊ v) holomorph (U(e)), auch wenn Im v < 0 oder > π.
- Wir benutzen für H denselben Zweig von h wie in θ (Zweigzahl 0). Ein anderer Zweig ändert nur b₀.

**(T2) aₙ − bₙ = O(1/(n² ln n)).**
Auf S_α gilt θ = 2 − z + h(u₂) (Satz Z(i), fortgesetzt), mit u₂ = c + log u + ε. Dabei ist ε = log(1+ε₁(z)) − log((1−q)/(−2πiz)) holomorph bei 0, mit ε(0) = 0, also |ε| ≤ C|u|. Die Zweige passen: Längs z = iy sind beide Seiten gleich ln(Ay) + iπ/2 + o(1).

Auf der Strecke von v = c + log u nach u₂ gilt h′(v) = h′(log₊v)/v. Weil Re v < 0, ist Im log₊ v ∈ (π/2, 3π/2). Nach (G3) ist dort |h′(log₊ v)| ≤ C, also |h′(v)| ≤ C/ln(1/|u|). Daraus folgt
  θ − H = 2 − z + E,  |E| ≤ C|u|/ln(1/|u|).
Ausserdem ist −z = (i/2π) log q = −(i/2π)u + O(u²). Für f := θ − H − 2 + (i/2π)u gilt also f = O(|u|/ln(1/|u|)) in Δ nahe 1. f ist Δ-analytisch, weil θ und H es sind und der Rest ein Polynom ist. Fern von 1 ist f auf dem Rand von Δ beschränkt.

Nach Flajolet–Odlyzko (Flajolet–Sedgewick, *Analytic Combinatorics*, Thm. VI.3; gilt für O-Terme mit beliebigen reellen α, β) folgt mit α = −1 und β = −1: [qⁿ]f = O(n⁻²(ln n)⁻¹). Für n ≥ 2 ist [qⁿ]f = aₙ − bₙ. ∎

Gebraucht wird nur eine Schranke, keine Entwicklung von h. Das war die eigentliche Lücke in S8.

**(T3) aₙ − b2ₙ = O(1/(n³ ln n)).** ε = δ₁u + O(u²) (Etappe 2, S7). h(v+ε) − h(v) = h′(v)ε + O(|ε|² sup|h″|), und es gilt
- h″(v) = (h″(w) − h′(w))/v² mit w = log₊ v,
- |h″(w)| ≤ C nach Cauchy, weil h′ auf {Im w ≥ π/4} durch (G3) beschränkt ist.

Wegen h′(v)·u = −u²H′(q) ist E = −δ₁u²H′ + O(|u|²/ln(1/|u|)). Mit −z = (i/2π)(−u − u²/2) + O(u³) ist f₂ := θ − H + δ₁u²H′ − 2 + (i/2π)(u + u²/2) = O(|u|²/ln(1/|u|)) und Δ-analytisch. Für n ≥ 3 ist [qⁿ]f₂ = aₙ − b2ₙ. Mit VI.3 (α = −2, β = −1) folgt die Behauptung. ∎

**(T4) Leitglied und Transfer-Reihe.** bₙ = −h′(vₙ)/n + O(1/(n ln² n)). Allgemein gilt bₙ = (1/n)Σ_{k=1}^{K} g_k(−1)^k h^{(k)}(vₙ) + O_K(1/(n ln^{K+1} n)), mit 1/Γ(s) = Σ g_k s^k. Dasselbe gilt für aₙ, weil aₙ − bₙ = O(1/(n² ln n)).

*Beweis* (Hankel-Kontur wie im Beweis von VI.3).
- Setze q = 1 + t/n. Dann ist c + log u = vₙ + log(−t). Zerlege H = h(vₙ) + h′(vₙ)(log u + ln n) + R.
- Die abgezogenen Teile haben exakt [qⁿ]1 = 0 und [qⁿ](log u + ln n) = −1/n für n ≥ 1.
- Auf dem inneren Teil |t| ≤ ln² n gilt |log(−t)| ≤ 2 ln ln n + π ≪ |vₙ|. Die Strecke von vₙ nach vₙ + log(−t) bleibt in {Re < 0, −π/2 < Im < 3π/2}, und dort ist |h″| ≤ C/|v|² (wie in T3). Also |R| ≤ C|log(−t)|²/ln² n, und dieser Teil trägt O((1/n)·ln⁻² n·∫|log(−t)|² e^{−Re t}|dt|) bei.
- Für ln² n ≤ |t| ≤ nδ ist |q|^{−n} ≤ e^{−c|t|} und |R| ≤ C ln n. Das gibt O(e^{−c ln² n}).
- Auf dem äusseren Kreis |q| = 1+η ist der Beitrag O((1+η)^{−n}).
- Für K Glieder: dieselbe Rechnung mit Taylor bis Ordnung K, |h^{(k)}(v)| ≤ C_k/|v|^k und [qⁿ](log u + ln n)^k = (1/n)·(d/ds)^k[1/Γ(−s)]_{s=0} + O(ln^k n/n²). ∎

**(T5) Grössenordnung.** vₙ liegt im zweiten Quadranten, also gilt Im Log vₙ ∈ (π/2, π) und |Log vₙ| ≈ ln ln n. Mit h′(vₙ) = h′(Log vₙ)/vₙ und (G3) folgt c/(ln n (ln ln n)⁴) ≤ |h′(vₙ)| ≤ C/ln n. Der Rest O(1/(n ln² n)) ist kleiner als die untere Schranke. Also gilt **c/(n ln n (ln ln n)⁴) ≤ |aₙ| ≤ C/(n ln n)** für n ≥ n₀, und ln|aₙ|/ln n → −1. ∎

Zahlen (`l2_transfer.py`, Daten aus Etappe 2; T1 = −h′(vₙ)/(n ln b)):

| n | n² ln n·\|a−b\| (e / 2) | n³ ln n·\|a−b2\| (e / 2) | ln n·\|a/T1 − 1\| (e) | n ln²n·\|b − T1\| (e / 2) |
|---|---|---|---|---|
| 16 | 0.0191 / 0.0190 | 0.0238 / 0.0238 | 0.695 | 0.115 / 0.115 |
| 256 | 0.0155 / 0.0155 | 0.0170 / 0.0170 | 0.795 | 0.142 / 0.141 |
| 1024 | 0.0134 / 0.0134 | 0.0142 / 0.0142 | 0.806 | 0.142 / 0.141 |
| 4000 | 0.0116 / 0.0116 | 0.0119 / 0.0119 | 0.810 | 0.141 / 0.139 |

Alle skalierten Reste sind beschränkt und fallen langsam. Das passt dazu, dass die Schranken um einen Faktor ~ln ln n zu grob sind. Der Leitterm allein hat den relativen Fehler ≈ 0.81/ln n. Die bewiesene Schranke dafür, O((ln ln n)⁴/ln n), ist gröber, weil die untere Schranke für |h′| grob ist.

## 6. Präzisierung zu Etappe 4

- **B-S1:** w = 1+θ′ ist nicht in L¹(0,1), denn ∫|θ′| ≥ |θ(½) − θ(ε)| → ∞. Damit sind F_a = 2gw, F_{V_k} = kg^{k−1}w und F_{e_m} = qᵐw nicht integrierbar. F_d = θ″/w = (log w)′ verhält sich wie −1/z (aus der Kette in Satz Z(iii)). F_c = −2θ″/w² ist mit (M3) integrierbar (≲ ℓW⁴).

  Richtig formuliert: Für p mit p(0) = 0 gilt p = Σ_k p̂(k)(e^{2πikx} − 1). Für k ≥ 1 erfüllt G = F_C·(q^k − 1) die Bedingungen von Z0′: |G| ≲ k|z||F_C| ≤ C k (ln W)^j/ℓ bzw. ≤ Ck. Ausserdem gilt G_∞ = −F_{C,∞} = 0. Also folgt
  **δC = Σ_{j≥1} F̃_C(j)·conj(p̂(j)),  F̃_C(j) := M[F_C·(e^{−2πijx} − 1)].**
  
  Das ist die Aussage von S1 mit regularisierten Koeffizienten. Die numerischen Prüfungen benutzen p ∝ (cos 2πkx − 1) bzw. sin 2πkx, sind also davon nicht betroffen. Die Richtung ⇒ in B-S2 (C¹ in endlich vielen Parametern) läuft mit denselben Majoranten, weil θ_T = θ∘φ + p und |φ(z)| ≍ |z| nahe 0.
- **B-S3 und Lemma R:** θ_T ∈ L¹ und lokal in H^p für alle p folgen aus (M1) für θ(φ(z)). Die Majoranten-Lücke ist dort geschlossen. Paulsen–Cowgill bleibt Literatur.

## 7. Allgemeine Basis b > η

Ersetze Log durch log_b = Log/ln b. Dann gilt log_b(ℍ) ⊂ {0 < Im < π/ln b} ⊂ ℍ, und L ist anziehend (|λ| > 1). Lemma U, G und H gelten wörtlich, mit h′ = σ′/(σ log λ) und K = 2 Im L/(|log λ||w−L||w−L̄|). In der Kette gilt Im u₃ = arg(u₂)/ln b ∈ (π/(2 ln b), π/ln b). Satz M gilt mit anderen Konstanten, Satz T wörtlich (c = ln(A/(2π ln b)) + iπ/2, vₙ = (c − ln n)/ln b). Numerische Konstanten sind nur für e ausgerechnet; der Transfer ist für e und 2 geprüft. Das Ein-Punkt-Zertifikat für S4 ist für b = 2 **nicht** gerechnet (θ-Reihe für Basis 2 liegt vor; ~1 min).

## 8. Grapheinträge: was ist jetzt geschlossen?

| Eintrag | Lücke vorher | jetzt |
|---|---|---|
| kn.satz.theta_singularitaet (id_kn_1) | S3 zitiert; S4 Hauptzweig nur numerisch; Lokalform (iii) numerisch | S3 **geschlossen** durch Lemma U (nur der benötigte Teil). S4 **geschlossen bis auf ein Gleitkomma-Zertifikat an einem Punkt** (Satz Z; Zweigzahl 0). (iii) ist jetzt **exakte Identität** θ = 2 − z + h(Log Log sexp z) auf D = {\|z\| < 0.1, Im z > 0}. Vorschlag: status skizze → bewiesen (rechnergestützt im Zweig). Basis 2: Zertifikat noch nicht gerechnet. |
| kn.beob.theta_koeffizienten (id_kn_1) | Transfer nicht streng | **Bewiesen** (Satz T), schärfer als angegeben: Rest O(n⁻²/ln n), zweite Ordnung O(n⁻³/ln n), Leitterm und K-Glieder-Reihe mit Rest O(1/(n ln^{K+1} n)), c/(n ln n (ln ln n)⁴) ≤ \|aₙ\| ≤ C/(n ln n), Exponent genau 1. Vorschlag: art beobachtung → satz. Offen: O-Konstanten nicht ausgeschrieben; „≍ 1/(n ln n ln ln n ⋯)“ bleibt Heuristik; Schwankung von arg aₙ. |
| kn.satz.theta_mittelwerte (id_kn_2) | Majorante abgelesen | **Geschlossen** (Satz M, Z0′). Präzisierung: Φ(θ) für ganze Φ mit log M_Φ(s) = O(e^{κs}), nicht für jede ganze Funktion. |
| kn.satz.binet_kappa (id_kn_2) | dieselbe Majorante | **Geschlossen** (E = log sexp′ − Lz − Lθ, Majorante aus M1). |
| br.satz.reelle_kennzeichnung (id_br_2) | θ_T ∈ L¹ und Majorante | **Geschlossen** (Abschnitt 6). Paulsen–Cowgill bleibt Literatur. |
| br.satz.endlich_viele_reichen_nicht (id_br_2) | Majorante für die Linearisierung | **Geschlossen, mit Korrektur der Formulierung:** F̂_C(j) durch die regularisierten F̃_C(j) ersetzen (F_C ist meist nicht in L¹). Gegenbeispiele unberührt. |
| br.frage.momentenfamilie (id_br_2) | benutzt S1 | Gleiche Präzisierung, sonst unberührt. |

Neue Quellen für die Einträge: dieser Bericht, `luecken/l1_schranken.py`, `luecken/l2_transfer.py`. Daten liegen in `luecken/daten/`. In `daten/` wurde nichts geändert.

## 9. Dateien und Aufrufe

Interpreter `C:/Users/janis/miniconda3/python.exe`, `PYTHONIOENCODING=utf-8`, kein Motor.

| Skript | Inhalt | Laufzeit | Ausgabe |
|---|---|---|---|
| `l1_schranken.py` | Taylor-Schranken, Univalenz r = 0.1, S4-Zertifikat, Grunsky-Probe, Konstanten C_θ, H₁, C₊, Vergleich mit θ | ~70 s | `daten/l1_schranken.json` |
| `l2_transfer.py` | skalierte Reste aₙ−bₙ, aₙ−b2ₙ, Leitterm (Basis e und 2) | ~1 s | `daten/l2_transfer.json` |
