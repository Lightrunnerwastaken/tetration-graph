# Offenes Problem: Eckmoden der Abschnitte eines Kompositionsoperators (Hypothese (H₀))

Stand 11.10.2026, Strang W. Kennzeichnung: **bewiesen** / **numerisch** / **Vermutung**.

## 1. Rahmen
Sei L ≈ 0.3181 + 1.3372i der Fixpunkt e^L = L, ρ := |L|, D := {|x| < ρ}. Gesucht sind Polynome p_N vom Grad N mit T_N[p_N∘p_N − exp] = 0
(T_N = Abschneiden nach Grad N), die gegen die Kneser-Halbiteration f_K (f_K∘f_K = exp, singulär in L, L̄) konvergieren. Newton braucht die Abschnitte
J_N := T_N J|_{P_N} des linearisierten Operators **J δ := δ∘f + (f'∘f)·δ**, als Matrizen in ℓ²_ρ (Koeffizient k mit ρ^k gewichtet). Testmodell: ein reelles
Kubik f mit denselben Fixpunktdaten (f(L) = L, f'(L) = μ := √L); alle Zahlen unten für dieses Modell, Kneser gleich auf 2–3 Stellen.
Geometrie: Γin := {|x| = ρ, |f(x)| < ρ} (Bogen von L über −ρ nach L̄), Γout der Rest; h := der Zweig von f⁻¹ mit h(L) = L;
„A-Bogen“ := der Bogen des Kreises |x − f(0)| = |L − f(0)| ausserhalb D (er läuft von L nach L̄ an Γout vorbei).

## 2. Befund (numerisch, exakte Arithmetik bis N = 1000)
J_N⁻¹ hat r(N) ≈ N/250 „Eckmoden“ mit σ_i ~ e^{cN}, **c ≈ 0.0120** (kubisch und Kneser gleich). Die linken Singulärvektoren sitzen bei arg L
(Breite ~1/N, Höhe ~√N), die rechten (adjungierten) auf dem „Steg“ (Index k ≈ 0.95N, Winkel nahe π). Schreibe die adjungierten Moden
y_i := σ_iv_i als Laurent-Polynome Y_i(t) = Σ_{k≤N} y_{i,k}t^{−k}, t = x/ρ.

## 3. Hypothese (H₀)
Es gibt eine feste geschlossene Kurve γ* durch L, L̄ (auf der Γout-Seite in D durch h(A-Bogen), auf der Γin-Seite knapp ausserhalb von D), so dass
**|Y_i(t)| ≤ C·N^κ für t ∈ γ*, alle Eckmoden i, alle N ≥ N₀.**
**Bewiesen (Satz W1'')**: Aus (H₀) folgt für die Eckanteile der Pol-Antworten T_i(ζ) = ⟨u_i, J_N⁻¹T_NJ(k_ζ − T_Nk_ζ)⟩ die Schranke
|T_i(ζ)| ≤ C₁N^{κ'}e^{−(N+2)g_K(ζ)} ausserhalb eines Kompaktums K (g_K Green-Funktion), auch auf dem Schlitz nahe L. Weil die Singularität von f_K in L
sitzt und ihre Daten nur polynomial fallen, braucht der Konvergenzbeweis die **polynomiale** Form; e^{o(N)} statt C·N^κ genügt nicht.
**Gemessen bis N = 1000**: C = e^{−4.06}, κ = 0.88 auf γ*; Eckkonstanten C_E = 1, κ_E = 0.6; Reserve ε ≥ 0.60.

## 4. Exakte Ausgangsidentität (bewiesen)
Mit F̃(t) := f(ρt)/ρ, m̃ := (f'∘f)(ρ·), 𝒰(t) := Σ u_k t^{−k−1}, P Polynom:
**T(t) := m̃(t)Y(t)/t + Σ_{F̃(z)=t} Y(z)/(z F̃'(z)) − 𝒰(t) + P(t) = O(t^{−N−2}) bei ∞**,
T rational mit Polen nur in 0 und F̃(0) — äquivalent zu J_Nᵀy = u. Die Summe läuft über die drei Urbilder.

## 5. Grenzproblem (numerisch, kein Beweis)
(1/N)ln|Y| bzw. (1/N)ln|T| konvergieren gegen eine Funktion H, die auf der **geklebten Fläche ℳ** harmonisch ist: Man schneidet aus der Ebene das
Gebiet zwischen der „Bucht“-Kurve B (in D, Γout-Seite) und der Nullkurve Z_T von T (ausserhalb D) aus und klebt B an Z_T mit f (H = H∘h). ℳ hat
Geschlecht 0; bei 0 und ∞ ist H ~ −ln|t|, auf der Naht ist H = 0, B ist ein freier Rand. In L ist ℳ ein Kegel mit Deckverschiebung 2πi − ln μ:
H ≈ Re(c₁ζ^{κ₁}), **κ₁ = 2πi/(2πi − ln μ) = 1.118 − 0.032i** (Daten: 1.114 − 0.018i). Das freie Randwertproblem ohne freie Parameter liefert
**c = 0.0118 ± 0.0004** (gemessen 0.0117–0.0120), das Randprofil ab 8° auf 3e-4, ∂c/∂μ = 0.044 − 0.040i (gemessen 0.042 − 0.035i).
Im Inneren ist H die Green-Funktion g_{Ω_y}(·, 0) des Gebiets Ω_y zwischen B und der Nulllinie einer dünnen „Linse“ ausserhalb Γin, und die N Nullstellen
von Y_N verteilen sich nach dem harmonischen Mass von Ω_y (N = 1000: Anteil auf B 0.470 gegen 0.473; Winkel-KS 0.005; Abstände ~1.4 ln N/N).

## 6. Was nicht trägt (Befunde)
- **Lokal:** Die Bucht-Schranke ist über die Identität GLEICHWERTIG zu einer T-Schranke auf ∂K (Kreisschluss intrinsisch). Die Eckzone um L braucht
  am Aussenrand winkelaufgelöste, polynomial genaue Daten; die Innenscheibe |t − L| ≤ N^{−1/Re κ₁} folgt dann per Maximumprinzip.
- **Extremalität in ℓ²_ρ:** ‖J_N‖ = e^{0.47N}; ein rein aus dem Modell gebauter Testvektor hat Überlapp 0.997 mit v₁ (N = 800), aber
  Rayleigh-Quotient 10^{90}.
- **Gonchar–Rakhmanov in der Ebene:** Für die exakten Interpolanten Q_n (Q_n∘f + (f'∘f)Q_n = O(x^n)) gilt ∮Q_n(s)ψ_k(s)ds = 0 mit
  ψ_k = R_k + (f'∘f)(s)s^{−k−1}, R_k(s) = Σ_{f(x)=s}x^{−k−1}/f'(x). Jede Bedingung koppelt ein durch f transportiertes Moment in f(0) mit einem Moment
  in 0 bei gleichem Index — kein System {s^jW(s)}. Auf ℳ hätte das Problem GR-Form, aber ℳ existiert nur im Grenzwert (die Klebung T ≈ Y∘h ist eine
  Dominanzbilanz mit relativen Fehlern e^{−Nδ}, δ > 0 ausser in L, L̄).
- **Andere Skalenkreise:** in ℓ²_{ρ'} mit ρ' ≠ ρ gibt es kein ρ' ohne exponentielles Wachstum.

## 7. Fragen an Fachleute
1. **n-te Wurzel:** Grenzwerte der Raten erfüllen (bewiesen) eine „tropische Bilanz“: in der exakten Identität wird das Maximum der Termraten
   stets zweimal angenommen. Ist die Grenzkonfiguration (mit Normierung und den Singularitäten in 0, ∞) eindeutig, und ist c ihr Hauptwert
   (keine Konfiguration mit grösserer Rate)? Strukturell: ein Hindernisproblem (Kontaktmenge = Bucht) mit nichtlokaler Kopplung H = H∘h.
2. **Starke Asymptotik:** Gibt es einen Deift–Zhou-artigen Rahmen für eine „Schweiss-RH-Aufgabe“ (zwei Gebiete, längs einer Kurve mit einer
   stückweise analytischen Randabbildung verklebt) mit kleiner NICHTLOKALER Störung, mit lokaler Parametrix am Kegelpunkt L
   (Innenschicht der Breite 1/N, Vorfaktor N^{1/2})?
3. **Literatur:** Abschnitte gewichteter Kompositionsoperatoren mit Fixpunkt auf dem Abschneidekreis; Polynomlösungen approximativer
   Funktionalgleichungen mit Komposition; Verwandtschaft zu Muttalib–Borodin-Biorthogonalität (Komposition x ↦ x^θ) — dort fehlt aber der
   Summenterm (f'∘f)·δ, der hier die Ecke in L erzeugt.

Dateien und Belege: w_grenzschicht/BERICHT.md §15–§26 (Skripte w38–w53).
