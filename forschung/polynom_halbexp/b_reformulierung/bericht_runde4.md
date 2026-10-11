
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
b43_zweigschritt.py (b43.log).
