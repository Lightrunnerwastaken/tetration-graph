
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
