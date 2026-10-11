# Strang X — Alternativen & Literatur (Phase 2)

Beginn 10.10.2026, 23:29 (Zeiten aus `date`). Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.
Quellen-Prüfstatus wie Strang C: **[O]** Originaltext (Seite/Satz) selbst gelesen, **[S]** sekundär (Abstract, Verlagsseite,
Zitat in einer selbst gelesenen Arbeit, ar5iv-Volltextauszug), **[U]** ungeprüft (Gedächtnis) — vor Weiterverwendung nachschlagen.
Literaturstand Phase 1: `c_literatur/BERICHT.md` §2 (hier nicht wiederholt, nur vertieft).
Bezeichnungen wie ARCHITEKTUR.md §1; zusätzlich μ = f'(L) = √L (|μ| = 1.1717, arg μ = 38.30°), h = f_K⁻¹, K̃δ := f'·(δ∘h),
T^{a₀}_N = Taylor-Abschnitt vom Grad N bei a₀.

Ziel: (1) Werkzeug-Literatur für W ((E_L), Grenzschicht der Eckmode) und O ((a′), Rest-Kontraktion), dazu M;
(2) Beweiswege, die (H) ganz umgehen — jeder numerisch so weit geprüft, dass klar ist, ob er trägt.

---

## 0e. Hauptbefund Meilenstein 6 (11.10. 04:10): explizite Konstanten — N₀ beweisbar ≈ 390, also Reserve-Zertifikat nötig (§3.9)

1. **Lemma X7 bewiesen** (bis auf zwei Gitter-Minima c_in = 0.674, c_out = 0.3456): |w_η|_{ℓ²_ρ} ≤ 6.3·√M auf ∂K' (gemessen 1.17√M) ⇒ C_R = 12.7, κ_R = 1.
2. **Lemma X8 (A_θ) bewiesen, elementar** — Koeffizienten von J_∞δ bis Zeile 5N statt Bernstein auf f(∂D): A_θ = 19.7 (c₀ = 0.54; gemessen 4.05).
   Kellogg–Warschawski (Lemma X6) wird für X4a nicht mehr gebraucht.
3. **Kontur d = 0.3:** Lage nachprüfbar (|x| ≤ 2.357, Konstante 17.1), braucht aber eine Sprungschranke am Schlitz bis u = 0.451 — dieselbe Lücke
   wie X5 mit R = 2; ohne sie Hybrid-Kontur, gleiches N₀ sobald A_θ ≥ 13.
4. **N₀ mit beweisbaren Konstanten ≈ 387** (Ecke C_E = 1, κ_E = 2), selbst mit perfekter Ecke ≈ 214; N₀ ∝ A_θ^{0.5}. ⇒ N₀ ≤ 193 nur mit besseren
   Konstanten; **Reserve: Zweig-Zertifikat bis 400 ≈ 12 h sequentiell / ½ h auf 32 Kernen, bis 1000 ≈ 32 Tage / 1–1.5 Tage** (Hochrechnung aus p24).
5. **(Fam-Eck) Bogenmodell bestätigt:** eine Glättung über einen Bogen fester Länge (Poisson-Tiefe ε ≈ 0.06, N-unabhängig) sagt ALLE Eckantworten
   (glatt, Zufall, p_N) auf ≈ 1–2 R voraus (Punktmodell: Faktor 10³–10⁴ daneben); ∂c/∂r ≈ 0.37.

## 0d. Hauptbefund Runde 4 (11.10. 03:25): Ecke reagiert polynomial auf Störungen; Bernstein selbsttragend; (a″) mit d ∝ N (§3.8)

1. **Lemma X6 (L²-Bernstein nahe f(∂D)), bewiesen** bis auf Kellogg–Warschawski (Pommerenke Thm 3.6): ‖δ'‖_{L²(Γ)} ≤ CN‖δ‖_{L²(f(∂D))} für Kurven Γ im
   Abstand ≤ c/N — damit ist Satz X4a ohne offene Literaturstelle. Eine fertige L²-Fassung in der Literatur habe ich nicht gefunden (Sup-Form: Kalmykov–Nagy–Totik 2017).
2. **(Fam-Eck) numerisch geklärt (N = 200/400/600):** Die Ecke reagiert auf JEDE glatte Sekanten-Störung nur polynomial (≈ N^{1.4–2}·R, am stärksten bei
   radialer Fixpunkt-Verschiebung: 32/131/230·R), nie wie σ₁ (1.84 → 8.9 → 68.4); das Matrixelement ⟨u₁, δAv₁⟩ ist selbst ∝ 1/σ₁. Die echte Lösung
   (kohärent: Buckel der Breite ~1/N bei L, L̄, oszillierend, lokale Mittel ≈ 0) gibt nur den gleichförmigen Bulk-Anteil −0.34·R. Verschärfte Vermutung:
   δσ_i/σ_i = O(R) + N·⟨∇c, geglättete Eckdaten⟩ — (Fam-Eck) ist die Glattheit der Ecken-Asymptotik c(L, μ, h|Γout) in den Eckdaten, d. h. sie fällt
   mit dem Beweis von (H₀) über Ws Modell mit ab, wenn dieser nur Eckdaten auf fester Skala benutzt. §3.7s „σ₁-Wachstum für Zufall“ ist widerlegt.
3. **(a″) mit kleinem d:** nötig ist d ≥ d₀(N) ≈ 0.02N − 2 (2/6/10/≳14 für N = 200/400/600/800), proportional zu N, nicht zur Eckmoden-Zahl.
   Grund: Überlaufprofile der Eckmoden trennen sich auf der Indexskala j ≈ 0.07N (Steg-Abstand). Ws Ratenmodell liefert (K_d) nicht direkt — es braucht
   die Vorfaktor-Ebene (Kegelharmonische); Anfrage an W.
4. **(H-i)-Text für ARCHITEKTUR §9** (Tangente nur auf 𝒵 verlangen) in §3.8.3, an P übergeben.

## 0c. Hauptbefund Runde 3 (11.10. 02:40): X3 abgeschlossen; Familie — eckfreier Teil bewiesen, Eckteil kohärent stabil (§3.5c, §3.7)

1. **Lemma X5 (Schwanz, bewiesen bis auf (K)-Fakten):** R = 1.649, Zeilen bis 5N, β = (M_R/ρ)(ρ/R)⁵ ≤ 0.798, Schlitzteil O(N^{−4.45}) ⇒
   (a″)_{4N} ⇐ (a″)_∞. Numerisch: Zeilen jenseits 2N ändern R̂, Rc um < 1e-11 (d = N und d = 4N gleich).
2. **Θ^c (d ≥ N, Kneser):** −0.039 / −0.012 / +0.020 / **+0.072** (N = 200/300/400/600), +0.028 je 100 Grade, getragen von der fernen Spitze.
   Grenzwert aus dem Profil 0.82 (0.99 mit Aufschlag) — **für d = 2N grenzwertig**; kleines d: 0.129 → 0.114 → 0.106 (N = 200/400/600), Θ_∞ ≈ 0.07–0.20.
3. **Satz X4a (eckfreier Teil der Familie, bewiesen modulo L²-Bernstein auf f(∂D)):** ‖(𝒥_P − J_∞)δ‖ ≤ θ_N‖J_∞δ‖, θ_N = O(N R_N), für ALLE P
   mit sup_{D̄}|P − f| ≤ R_N ⇒ (a″) und die eckfreie Pol-Schranke gleichmässig; Schwanz nur für f nötig. Bedingung κ < 3.69.
4. **Eckteil:** Satz X4b reduziert (Fam) auf die Eckdaten von 𝒥_P. Numerisch (x18): an der echten Lösung p_N ändern sich die Eckmoden relativ nur um
   0.34·R (Sekante) bzw. 0.68·R (Tangente), konstant von N = 400 bis 600, obwohl σ₁ von 8.9 auf 68.4 wächst; zufällige Sekanten-Störungen gleicher
   Grösse ≈ N^{1.5}·R, hochfrequente ≈ 0. **Keine exponentielle Eckverstärkung** — Vermutung (Fam-Eck) formuliert.
5. **Warnung:** Tangenten J_N(P) an beliebigen P der Kugel sind keine kleinen Störungen (η' wird auf f(D̄) ∋ 2.23 ausgewertet, Verstärkung 1.624^N);
   (H-i) für die Tangente gilt nur an Fast-Lösungen.

## 0a. Hauptbefund Runde 1 (11.10. 00:40): Eckmoden = Effekt der QUADRATISCHEN Abschneidung (§3)

- **Lemma 3.1 / Korollar 3.2 (bewiesen, elementar):** mit dem Rechteck-Abschnitt R = T_{N+d}J_∞|_{P_N} gilt J_N⁻¹ = R⁺_top + R⁺_bot S⁻¹BR⁺_top
  (S d×d), also σ_{d+1}(J_N⁻¹) ≤ 1/σ_min(R), und Π_N g = Π̂g + R⁺_bot S⁻¹ ρ̂(g) mit der eckfreien Kleinste-Quadrate-Projektion Π̂ und dem
  d-dimensionalen Überlauf-Residuum ρ̂(g).
- **Numerisch (kubisches Modell exakt, N ≤ 500; Kneser N = 40):** σ_min(R) = 0.7328 = min_{Γin}|m| = 1/‖σ‖_∞ (genau Ds/Os Rest-Konstante), sobald
  d ≥ 2, 2, 3, 5, 8 (N = 120 … 500); der quadratische Abschnitt hat dagegen 0.54 … 0.040. ⇒ (a′)/(E_R) ⇐ (a″) „Rechteck-Abschnitte gleichmässig
  injektiv“ — eine Aussage ohne Abschneidekante; lokal an L ist das Mellin-Symbol auf Re s = −½ durch |μ| = 1.17 von 0 weg.
- Konsequenz für die Architektur: die Eckproblematik reduziert sich auf eine d×d-Matrix S_N (d ≈ N/60 … N/100) und das Überlauf-Residuum.

## 0b. Runde 2 (11.10. 00:47–01:10): Satz X1/X2 im K'-Rahmen, Kneser-Daten (§3.5, §3.5a, §3.6)

- **Satz X1 (Reduktion bewiesen):** (F_{K'}) + (Univ) + (a″)_d + (Ẽ) ⇒ (LP). Der eckfreie Teil Π̂ bekommt (S) auf ganz ∂K' aus (a″) und Ps P-II-Argument
  (keine Eckmoden-Projektoren, kein r(N), kein (E_A), keine Bulk-Annahme). **Satz X2:** (Ẽ) ⇐ einteilige Randbedingung Θ^c_N ≤ 1 − ε für den
  kanonischen Eckteil c_N = Π_N − Π̂. Billige Formel: c_N(g) = J_N⁻¹Zᵀ(I + ZZᵀ)⁻¹o(g) (eine LU-Zerlegung).
- **Kneser (N = 200/300/400, Ps ∂K'-Gitter):** Θ^c = 0.129/0.129/0.114 (Ecken-Grenzwert −30, Schlitz ≤ −0.38), eckfreier Teil Θ = −0.24/−0.18/−0.16.
  Kubisch bis N = 600: eckfreier Teil auf ∂K' polynomial (fern 10.4 → 28.8), die ganze Spitze im Eckteil.
- **Ehrliche Einordnung (M 01:04):** (a″)_d mit kleinem d ist eine Kantenaussage vom Rang von (a′)_r. **Variante d = 2N bzw. d = N (Satz X3, §3.5b):**
  (a″) wird Os (a″)_∞ (unendlicher Operator, keine Kante; Kneser bis auf eine skizzierte Schwanzabschätzung mit f_K(2) = 3.43 < 8/ρ² = 4.23); der
  Eckteil ist dann gross, aber sein Randprofil flach: Kneser Θ^c = −0.039/−0.012/+0.020 (N = 200/300/400), Rc(L) ≈ 0.55·R(L).
  ⇒ Linearer Teil = (a″)_∞ [O] + EINE messbare Randungleichung Θ^c ≤ 1 − ε [W/X].
- **Eckoperator (§3.6):** 𝒞_N (Rang ≤ d) reproduziert σ₁, σ₂(J_N⁻¹) (N = 200…500) ohne SVD; Daten für W in x11_N*.npz.

## 0. Kurzfassung (Stand 10.10. 23:55)

- **Literatur:** Die drei passenden Theorien und je der Schritt, den sie liefern könnten:
  (i) *Splitting-Satz / Index-Kalkül endlicher Abschnitte* (Roch–Silbermann) = genau die Form von (a′) [O];
  (ii) *Rand- und Wellenpaket-Pseudomoden* (Böttcher–Grudsky: exakte Exponentialrate aus dem nullstellenfreien Ring [O];
  Trefethen–Chapman, Borthwick–Uribe, Dencker–Sjöstrand–Zworski: variable Koeffizienten, Rate = Wirkungsintegral) = Werkzeug für (E_L);
  (iii) *Kompositionsoperatoren mit Randfixpunkt* (Cowen, Bourdon [S]) — **unanwendbar** an L, und zwar beweisbar: Julia–Carathéodory
  verbietet einen drehenden Randfixpunkt für jede Selbstabbildung (Lemma 1.1, bewiesen).
- **Korrektur (numerisch, x02):** die von C (BRETT 23:05 (5), C §3.5) gemeldete Divergenz von Π_N f_K im Inneren des Halbmonds ab
  N ≈ 600–900 ist ein Rundungsartefakt der 60-stelligen Koeffizientendatei; echt ist nur das Wachstum am Scheitel a₀+R₀ (B b42).
- **Alternativwege:** Vorzeichen/Positivität/Monotonie in N tragen nicht (x01: Kneser-Koeffizienten wechseln das Vorzeichen mit Periode
  2π/arg L, Monotonie in N durch Phasendrehung von c_N widerlegt); Ersatzform-Kontraktion scheitert spektral (Bourdon-Scheibe ∋ −1).
  Bemerkenswert und ungenutzt: das Gewicht der negativen Kneser-Koeffizienten ist nur 0.0039 von 2.2405 (ℓ¹_ρ), und
  max_{|x|=r}|p_N| wird für alle N ≤ 1000 bei x = r > 0 angenommen.

---

## 1. Werkzeug-Literatur für O, W, M

### 1.1 Lemma (Julia-Hindernis; bewiesen, elementar)

*Aussage.* Sei U ein Gebiet mit L ∈ ∂U, das in L einen Tangentialkegel S (offener Sektor mit Scheitel L, Öffnung β ∈ (0, 2π))
besitzt, d. h. U ∩ B(L,r) ⊂ S_ε (Kegel S um ε aufgeweitet) und S_ε' ∩ B(L,r) ⊂ U für kleine r. Dann gilt weder f_K(U) ⊂ U noch
h(U) ⊂ U. Invariant können nur Gebiete sein, deren Rand ohne Tangentialkegel (z. B. logarithmisch-spiralig mit Steigung
arg μ/ln|μ| = 4.21) in L einläuft.
*Beweis.* Nach der Singularitätsformel (A1/B L8.2) ist f_K − L = μ t + O(|t|^{min(2, 1+Re γ)}) = μt + O(t²) für t = z − L im
Schlitzgebiet, also f_K(L + re^{iθ}) = L + |μ|r e^{i(θ + arg μ)}(1 + o(1)) gleichmässig für θ in kompakten Teilen von S. Da S − arg μ bzw.
S + arg μ für einen echten Sektor und arg μ ∉ 2πℤ nicht in S̄ liegt, gibt es eine Richtung θ ∈ S, deren Bild (für kleines r) ausserhalb
des aufgeweiteten Kegels, also ausserhalb U liegt. Für h ebenso mit 1/μ. ∎
*Bezug.* Das ist die lokale Form von Julia–Carathéodory: für eine holomorphe Selbstabbildung der Kreisscheibe ist die Winkelableitung in
einem Randfixpunkt reell und positiv [S: Encyclopedia of Mathematics „Julia–Wolff–Carathéodory theorem“, Cowen–Pommerenke, J. London
Math. Soc. 26 (1982) 271–289]. Für die Kreisscheibe D (Kegel = Halbebene, innere Normale bei 256.6°) heisst das: **Γout grenzt
zwangsläufig an L** — der Mechanismus hinter Ds Eckmode ist nicht ein Zufall der Kneser-Funktion, sondern erzwungen durch arg μ ≠ 0.
(Reelles μ > 1 erlaubt Γout = ∅ — passt zu Ds kubischem Test „μ reell ⇒ keine Eckmode“.)
*Geometrie (bewiesen, aus Lemma 1.1 mit den Normalen):* innere Normale von D in L: 256.6°, von D(a₀,R₀): 277.7° (Winkel 21.1°).
Das Bild von D(a₀,R₀) unter h (Drehung −38.3°) ragt deshalb in den Richtungen 187.7°–204.9° (ab L) aus D̄; das Bild von D unter f
(Drehung +38.3°) ergibt zusammen mit D den Kegel (166.6°, 384.9°), Innenwinkel 218.3°, **Aussenwinkel 141.7°** — exakt Strang Ps
Kompaktum K' = D̄ ∪ f_K(D̄) (P: 141.69°).

### 1.2 Für O — (a′) Rest-Kontraktion

| Quelle | Aussage | Technik | Welchen Schritt von (a′) |
|---|---|---|---|
| **Böttcher–Grudsky, Toeplitz band matrices with exponentially growing condition numbers, ELA 5 (1999) 104–125** [O, S. 104–109] | S. 105–106: Gohberg–Feldman: lim inf σ₁(T_n(a−λ)) > 0 ⇔ λ ∉ a(𝕋), wind(a,λ) = 0. Roch–Silbermann („splitting phenomenon“): bei wind = k ≠ 0 gehen genau \|k\| Singulärwerte gegen 0, σ_{\|k\|+1} bleibt von 0 weg. Thm 1.1: σ_{\|k\|} ≤ Ce^{−αn}; Gl. (3): α = min(log 1/r_λ, log R_λ), r_λ < 1 < R_λ der grösste Ring ohne Null- und Polstellen von a − λ. Thm 1.7: bei einer Nullstelle auf 𝕋 polynomial oder exponentiell je nach Windungszahl der Nachbarkomponenten. | Fredholm-Index der Grenzoperatoren T(a), T(ã); explizite Konstruktion fast-singulärer Vektoren aus Wurzeln von a − λ | (a′) ist genau eine **Splitting-Aussage** für I + K_N: r(N) Richtungen werden exponentiell klein, der Rest bleibt ≥ 1 − o(1). Ds Grenzwert 1.3663 = ‖σ‖_∞ ist der Beitrag des ersten Grenzoperators (Toeplitz mit Symbol σ, ‖T(σ)‖ = ‖σ‖_∞). |
| **Roch–Silbermann, Index calculus for approximation methods and singular value decomposition, J. Math. Anal. Appl. 225 (1998) 401–426** [S, Bibliographie TU Darmstadt + Zitat bei Böttcher–Grudsky] | Zahl der gegen 0 gehenden Singulärwerte = dim ker W(A) + dim ker W̃(A) (W, W̃ = starke Grenzwerte von A_N und der geflippten Folge R_N A_N R_N) | C*-Algebra der Folgen modulo Nullfolgen, „Fredholm-Folgen“ | Rahmen für „σ_{r+1} beschränkt“: (a′) ⇔ es gibt Rang-r-Korrekturen F_N mit sup‖(I + K_N + F_N)⁻¹‖ < ∞. **Hindernis:** W̃ existiert bei uns nicht in festen Flip-Koordinaten (D §11.2: Eckblock ändert sich um 40 % von N = 300 auf 400), und r(N) wächst (1 → 2 → 3 bei N ≈ 120, 300, 700). |
| **Silbermann, Modified finite sections for Toeplitz operators and their singular values, SIAM J. Matrix Anal. Appl. 24 (2003)** [S] und Rogozhin–Silbermann (Block-Toeplitz, ~2005) [S] | Splitting für modifizierte Abschnitte; lim σ_{k+1}(T_n(a)) = min(‖T(a)⁺‖⁻¹, ‖T(ã)⁺‖⁻¹) | Moore–Penrose-Inverse der Grenzoperatoren | gäbe den **Wert** der Rest-Schranke: (a′) mit Konstante 1 ⇔ ‖W̃⁺‖ ≤ ‖W⁻¹‖. |
| **Mascarenhas–Silbermann** (Toeplitz-Matrizen mit variablen Koeffizienten, Symbol a(x,θ) auf dem Zylinder [0,1]×𝕋, k-Splitting) [S, nur Seminarzusammenfassung IST Lissabon 04.07.2014; Originalarbeit nicht gefunden] | k-Splitting, k = Summe der Kerndimensionen zweier Toeplitz-Operatoren | C*-Algebra der Kac–Murdock–Szegő-Folgen | **der passende Rahmen**, wenn man J_N in k/N-Variablen schreibt (die Eckschicht skaliert mit N: Ds Tiefe d* ≈ 0.05N); unser Symbol ist allerdings an Γout nicht beschränkt (Komposition), also nicht in der Wiener-Algebra. |
| Garoni–Serra-Capizzano, Generalized Locally Toeplitz Sequences I/II (Springer 2017/2018) [S]; Tilli, Locally Toeplitz sequences, LAA 278 (1998) 91–120 [S]; Kac–Murdock–Szegő, J. Rational Mech. Anal. 2 (1953) 767–800 [S] | Verteilung der Singulärwerte variabler Toeplitz-Folgen nach |κ(x,θ)| | GLT-Algebra | nur Verteilungs-, keine Extremalaussagen; nützlich als Plausibilität „Bulk ≤ ‖σ‖_∞“. |
| Gohberg–Fel'dman, Convolution Equations and Projection Methods for their Solution, AMS Transl. Math. Monogr. 41 (1974) [S: Bibliographie Trefethen–Chapman]; Böttcher–Silbermann, Introduction to Large Truncated Toeplitz Matrices, Universitext 1999 [S]; Analysis of Toeplitz Operators, 2. Aufl. 2006 [U] | stückweise stetige Symbole: Abschnitte stabil ⇔ T(a), T(ã) invertierbar; Sprünge in ℓ² durch Strecken (mit Potenzgewicht durch Kreisbögen) überbrückt | lokales Prinzip (Gohberg–Krupnik, Douglas) | erklärt Ds Toeplitz-Test d20 exakt (instabil nur bei ψ = π: die Sprungstrecke von m(L) nach −Λ läuft durch/um 0). Für die echte Komposition auf Γout nicht anwendbar. |

Weitere über die Bibliographie von Trefethen–Chapman [O] belegte Quellen [S]: Reichel–Trefethen, Eigenvalues and pseudo-eigenvalues of Toeplitz matrices, LAA 162–164 (1992) 153–185; Böttcher, Pseudospectra and singular values of large convolution operators, J. Integral Equations Appl. 6 (1994) 267–301; Landau, On Szegő's eigenvalue distribution theory and non-Hermitian kernels, J. Analyse Math. 28 (1975) 335–357; Zworski, PAMS 129 (2001) 2955–2957; Davies, Proc. R. Soc. Lond. A 455 (1999) 585–599.

**Wichtigste Folgerung für O:** Die übliche Beweisform eines Splitting-Satzes braucht den *zweiten* Grenzoperator an der
Abschneidekante. Bei uns lebt die Kante semiklassisch (Skala k/N), deshalb ist das richtige Objekt ein Toeplitz-artiger Operator mit
variablen Koeffizienten in x = k/N nahe x = 1 — die Literatur dafür (Kac–Murdock–Szegő, Mascarenhas–Silbermann, GLT) setzt beschränkte
Symbole voraus; unser Symbol ist auf Γout eine Komposition (in Phasenraum-Sprache: komplexe Streckung k ↦ μk, §1.3).

### 1.3 Für W — Grenzschicht/WKB der Eckmode, (E_L)

| Quelle | Aussage | Technik | Welcher Schritt |
|---|---|---|---|
| **Böttcher–Grudsky ELA 1999, Thm 1.1 + Gl. (3)** [O] | exakte Exponentialrate α = min(log 1/r_λ, log R_λ): Abstand der nächsten Wurzel von a − λ zum Einheitskreis in log-Mass | Randmoden ∝ z₀^k mit Wurzel z₀ | Prototyp: **Rate = Im der komplexen Wellenzahl** der Wurzel des lokalen Symbols. |
| **Trefethen–Chapman, Wave packet pseudomodes of twisted Toeplitz matrices, CPAM 57 (2004) 1233–1264** [S, Abstract]; **Trefethen–Embree, Spectra and Pseudospectra (Princeton 2005), Kap. 7 „Toeplitz matrices and boundary pseudomodes“, Kap. 8 „Twisted Toeplitz matrices and wave packet pseudomodes“, Kap. 9, 10–11** [S, Inhaltsverzeichnis] | Bei variablen Koeffizienten a(x,k) existieren Pseudomoden mit exponentiell kleinem ε (Wellenpakete), wenn eine Windungs-/Twist-Bedingung (≈ Hörmanders Kommutatorbedingung) gilt; Symbol braucht in x nur an einer Stelle differenzierbar zu sein | WKB-Ansatz v_k = exp(N∫^{k/N} iξ(s)ds), ξ(s) Wurzel des lokalen Symbols | Ansatz für die adjungierte Eckmode: **c = ∫_{x*}^{1} log(1/\|z(x)\|) dx** über die Schicht zwischen Umkehrpunkt x* (\|z(x*)\| = 1) und Kante x = 1 (Vermutung; Ds Tiefe d* ≈ 0.05N ↔ x* ≈ 0.95). Das Randprofil V(arg L + τ) ist dann die Fourier-Transformierte eines Wellenpakets; analytische Phase S(τ) gibt λ(τ) = Re(S(τ) − S(0)) ~ τ² (p = 2), falls S in τ = 0 analytisch ist — Ds gemessenes p ≈ 1.5 deutet auf eine Nicht-Analytizität in L (Ecke, Exponent 1/α = 1.13) hin. |
| **Trefethen–Chapman, Volltext [O] (S. 1233–1258 gelesen)** | Def. 2.1/2.2: „twisted Toeplitz“ a_{jℓ} = c_{j−ℓ}(x_j), Symbol f(x,k) = Σ c_ℓ(x)e^{iℓk}. §3: WKB v_j = α(x_j)e^{−iNφ(x_j)}, Eikonal f(x, φ'(x)) = λ; Wellenpaket bei x* falls Im φ''(x*) < 0 ⇔ **Twist Im(∂_x f/∂_k f) > 0** (3.4) (Hörmander-Analogon). Thm 4.1: dann ‖(A−λ)v‖/‖v‖ ≤ M^{−N}, |v_j| ≤ C₁exp(−C₂N(x_j−x*)²), sofern f(x*,k) ≠ λ für k ≢ k*. Beweis über Transfermatrizen der Rekursion und exponentielle Dichotomie (Lemma 4.2/4.3, Kegelargument). §6: Twist ⇔ Windungszahl I(f,λ,x) fällt beim Durchgang durch x*. Thm 7.4: Symbol springt bei x* mit I⁻ > I⁺ ⇒ Pseudomode mit |v_j| ≤ C₁exp(−C₂N\|x_j − x*\|). Thm 7.5: bleibt richtig, wenn J feste Zeilen bei x* beliebig (auch voll) geändert werden. Kor. 7.6 = Reichel–Trefethen Thm 3.2 (Randmoden echter Toeplitz-Matrizen bei I(f,λ) < 0). | Transfermatrix + Kegel (konstruktiv); Raten nicht scharf (M ≈ exp(γΔx/4)), scharf nur per WKB (§8 dort) | Die Abschneidekante ist ein Sprung des Symbols (Matrix endet) ⇒ Thm 7.4/7.5-Situation: Randpseudomoden, Lokalisierung exp(−C₂N|x−1|) — passt zu Ds Eckvektoren (v₁ auf Zeilen N − d*, d* ∝ N). **Achtung:** alle Sätze dort setzen BANDMATRIZEN voraus; unser J̃_N ist voll (Komposition, f(0) ≠ 0) — in Ms radialem Modell ist die Spalte j eine Toeplitz-Spalte mit Symbol s(φ)^j = e^{Nx log s(φ)}, also exponentiell in N (kein festes Symbol). Ein „Spalten“-Analogon (Thm 7.3, Transposition) mit exponentiellem Symbol wäre zu entwickeln. |
| **Borthwick–Uribe, On the pseudospectra of Berezin–Toeplitz operators, Methods Appl. Anal. 10 (2003)** (arXiv math/0304182) [S] | deutet die Twist-Bedingung als Hörmanders Lösbarkeitsbedingung; Pseudomoden im Phasenraum lokalisiert | Berezin–Toeplitz-Quantisierung: Polynome vom Grad ≤ N = holomorphe Schnitte, h = 1/N | P_N ist genau eine Berezin–Toeplitz-Projektion (Hardy-Raum, h = 1/N) — der natürliche semiklassische Rahmen für J_N. |
| **Dencker–Sjöstrand–Zworski, Pseudospectra of semiclassical (pseudo-)differential operators, CPAM 57 (2004) 384–415** (arXiv math/0301242) [S] | Quasimoden dort, wo {Re p, Im p} ≠ 0 (Klammerbedingung); Resolventenschranken am Rand des Wertebereichs | Mikrolokalanalyse | Formale Begründung für „exponentiell kleine Quasimoden ⇔ Klammerbedingung am Phasenraumpunkt (x = arg L, k/N = 1)“. |
| **Flajolet–Sedgewick, Analytic Combinatorics (CUP 2009), §VIII.8 „Large powers“, S. 585ff.** [S, Inhaltsverzeichnis] | gleichmässige Sattelpunkt-Asymptotik von [z^k]A(z)^j für k/j in kompakten Bereichen | Sattel bei τA'(τ)/A(τ) = k/j | Einträge der Kompositionsmatrix: [x^k]f^j ≈ f(x_s)^j x_s^{−k}/√(2πj·…), Sattel j f'(x_s)/f(x_s) = k/x_s. **Sattel in L ⇔ k = μj** (komplex): die Komposition wirkt bei L als komplexe Streckung der Frequenz — Bs Streckband [N/\|μ\|, N] ist das Spaltenband, dessen Bild jenseits der Abschneidung liegt. Damit ist J_N ein Fourier-Integral-Operator mit Phase j log f(x) − k log x — W kann das lokale Symbol daraus ablesen. |
| Miller-Algorithmus: Gautschi SIAM Rev. 9 (1967) §3 [O nur S. 24–28 laut C]; **Olver, Error analysis of Miller's recurrence algorithm, Math. Comp. 18 (1964) 65–74** [S]; Olver, J. Res. NBS 71B (1967) 111–129 [S] (Randwertverfahren, automatische Fehlerkontrolle); Olver–Sookne, Math. Comp. 26 (1972) 941–947 [S] | Fehler des abgeschnittenen Problems = (f_N/g_N)·(dominante Lösung), exakt | Casoratian | eigene Herleitung (bewiesen, trivial): Lösung von y_0 = f_0, y_N = 0 ist y_n = f_0(f_n g_N − f_N g_n)/(f_0 g_N − f_N g_0) ⇒ relativer Fehler −(f_N/g_N)(g_n/f_n − g_0/f_0). Bs Ein-Moden-Befund c_N ≈ K·N·u_N/(δ_{−1})_N hat genau diese Form, bis auf den Faktor N (Grenzschicht). |

**Weitere Werkzeuge für W (kurz):**
- *Agmon-/Finsler-Abstand für Differenzenoperatoren mit h = 1/N:* Klein–Rosenberger, Agmon-type estimates for a class of difference operators,
  Ann. Henri Poincaré 9 (2008) 1177–1215 (arXiv 1706.06331) [S] — exponentieller Abfall von Eigenfunktionen mit Rate = Finsler-Abstand aus einer
  Eikonalgleichung; selbstadjungiert, also nur als Muster für „Rate c = Wirkungsabstand vom Umkehrpunkt zur Kante“.
- *q-Differenzengleichungen:* lokal bei L ist die Linearisierung δ(L + με) + μδ(L + ε) = g eine **Fuchssche** (reguläre) q-Differenzengleichung mit
  q = μ (|q| = 1.17, nicht reell): lokal trivial (Lösungen ε^s × q-elliptisch; die q-elliptischen Funktionen sind genau Knesers 1-periodische θ in
  Abel-Koordinaten). Nichttrivial ist nur das **Verbindungsproblem** zwischen L und L̄ — in der Möbius-Koordinate M = (x−L)/(x−L̄) zwischen 0 und ∞,
  Birkhoffs „generalized Riemann problem“ (Proc. Amer. Acad. 49 (1913)) [U]; moderne Klassifikation: Ramis–Sauloy–Zhang, Astérisque 355 (2013)
  (arXiv 0903.0853) [S] (behandelt v. a. irreguläre Gleichungen; für uns nur Begriffe). Für Os Riemann–Hilbert-Test heisst das: das RH-Problem lebt
  auf dem Ring ℂ*/μ^ℤ (Fundamentalbereich der Streckung), nicht auf dem Kreis.
- *Koopman/EDMD (Galerkin-Abschnitte von Kompositionsoperatoren):* Colbrook–Townsend, Rigorous data-driven computation of spectral properties of
  Koopman operators (arXiv 2111.14889, ResDMD) [S]: Pseudospektren von Kompositionsoperatoren ohne „spectral pollution“ über **rechteckige**
  Galerkin-Matrizen und Residuen — derselbe Kunstgriff wie §3 unten.

### 1.4 Kompositionsoperatoren (für O, M; Prüfung von Ds [U]-Zitat)

- **Bourdon, Spectra of some composition operators and associated weighted composition operators, J. Operator Theory 67 (2012) 537–560
  (arXiv 0908.1405)** [S, ar5iv-Volltextauszug]: Thm 4.1 (hyperbolisch, Denjoy–Wolff-Punkt 1, φ'(1) < 1): Sp(C_φ) = Sp_e(C_φ) = {|z| ≤ φ'(1)^{−1/2}};
  Thm 4.3 (parabolisch, a = φ''(1)): Spirale {e^{−at}: t ≥ 0} ∪ {0}; **Thm 5.5** (gewichtet): Sp(C_{g,φ}) = {|z| ≤ |g(1)|φ'(1)^{−1/2}}; zitiert
  Cowen, Composition operators on H², J. Operator Theory 9 (1983) 77–106 (Thm 2.1: Spektralradius φ'(ω)^{−1/2}) [S].
  → Ds Formel (D §9.3.5, dort [U]) ist damit belegt. Formal für K̃δ = f'·(δ∘h) (g = f', φ = h, φ'(L) = 1/μ): Radius |μ|·|μ|^{1/2} = |μ|^{3/2} = 1.268 > 1,
  also −1 im Inneren — das passt zu den unendlich vielen Kernfunktionen der unendlichen Linearisierung (Konkurrenten). **Aber** nach Lemma 1.1
  ist h keine Selbstabbildung (und hat zwei anziehende Randfixpunkte L, L̄ — für Selbstabbildungen unmöglich, nur ein Denjoy–Wolff-Punkt),
  der Satz ist also nur Heuristik.
- **Toeplitz-Kompositions-C*-Algebren:** Kriete–MacCluer–Moorhouse, J. Operator Theory 58 (2007) 135–156 [S] (C*(T_z, C_φ), φ linear-gebrochen mit
  Randpunkten); Quertermous, arXiv 1205.5786 [S, ar5iv-Auszug]: φ mit Randfixpunkt ζ, φ'(ζ) ≠ 1 (reell): C*(C_φ, 𝒦)/𝒦 ≅ Unitisierung von
  C₀([0,1]) ⋊_β ℤ, β_n f(x) = f(x^{φ'(ζ)^n}); Thm 5.5: Fredholm-Kriterium für Σ f_n U^n (f₀ ohne Nullstellen auf [0,1] und Σ f_n(1)z^n ohne
  Nullstellen in D̄), alle Fredholm-Indizes 0 (K₁ = 0); Thm 2.3: C_φT_w ≡ w(ζ)C_φ mod 𝒦. → **Für O das nächste fertige Symbolkalkül** einer
  Randfixpunkt-Komposition; der Fixpunkt liefert eine Intervall-(Mellin-)Variable. Bei uns ist φ'(ζ) = μ komplex und φ keine Selbstabbildung
  (Lemma 1.1) — die Theorie greift nicht direkt, zeigt aber die Form, die ein Symbol für J_∞ an L haben müsste.
- Glashoff–Ortlieb, arXiv 1705.00325 [S]: Finite-Section-Methode für Kompositionsoperatoren (Formabbildungen) — „not always successful“; keine
  Stabilitätskriterien, die übertragbar wären. Eine Theorie endlicher Abschnitte von C_φ auf H² mit Randfixpunkt habe ich **nicht** gefunden.

### 1.5 Für M — exaktes Modell mit drehendem abstossendem Randfixpunkt (bewiesen, Rechnung)

Loxodromisches Möbius-Modell f = M⁻¹(μM), M(x) = (x−L)/(x−L̄) (D d24: Eckmode, Rate 0.077 bei μ = √L; L̄ ist dort anziehend — anders
als bei exp, wo beide Fixpunkte abstossen; eine Möbius-Abbildung kann nicht zwei abstossende Fixpunkte haben). Exakt:
x − L̄ = (L−L̄)/(1−M), f(x) − L̄ = (L−L̄)/(1−μM). Jedes Polynom vom Grad ≤ N ist eindeutig (x−L̄)^N Q(M) mit Q Polynom vom Grad ≤ N, und
  C_f[(x−L̄)^N Q(M)] = (x−L̄)^N·((1−M)/(1−μM))^N·Q(μM),
das Abschneiden T_N bei x = 0 ist das Abschneiden bei M₀ = M(0) = L/L̄ (|M₀| = 1, arg M₀ = 2 arg L). In der M-Ebene ist D eine Halbebene
durch 0 (Rand durch L ↦ 0 und L̄ ↦ ∞), μ wirkt als loxodromische Streckung um 0. Die Eckmode sitzt also bei M = 0 (L) — der Fixpunkt der
Streckung liegt auf dem Rand der Halbebene, der Abschneidepunkt M₀ auf dem Einheitskreis.

**Satz 1.5 (exakte Faktorisierung für jede Möbius-Abbildung; bewiesen, geprüft x03 auf 1e-44).** Für f = (ax+b)/(cx+d), g = (a b; c d)
und P ∈ P_N ist P∘f = (cx+d)^{−N}·π_N(g)P mit (π_N(g)P)(x) := (cx+d)^N P(f(x)) = Σ_j p_j (ax+b)^j (cx+d)^{N−j} ∈ P_N — die
(N+1)-dimensionale irreduzible Darstellung von GL(2,ℂ) auf Polynomen vom Grad ≤ N (Matrixelemente in der Monombasis = Krawtchouk-Polynome,
Koornwinder SIAM J. Math. Anal. 13 (1982) 1011–1023 [S]). Daher
  T_N[P∘f] = L_N(r)·π_N(g)P,  r := (cx+d)^{−N},  L_N(r) = untere Dreiecks-Toeplitz-Matrix der Taylor-Koeffizienten von r bei 0,
und mit m = f'∘f (L_N bildet eine kommutative Algebra, L_N(r)⁻¹L_N(m) = L_N(m/r)):
  **J_N = L_N(r)·[π_N(g) + L_N(m/r)],  I + K_N = I + L_N(r/m)·π_N(g)**  (Ds Faktorisierung D8 in Möbius-Form).
π_N(g) ist diagonal in der Fixpunktbasis e_j = (x−L₁)^j(x−L₂)^{N−j}: π_N(g)e_j = (a−cL₁)^j(a−cL₂)^{N−j}e_j; mit a − cL₁ = cL₂ + d =: κ₁,
a − cL₂ = cL₁ + d =: κ₂ ist κ₁/κ₂ = f'(L₁) = μ, also Eigenwerte κ₂^N μ^j (geometrische Folge mit Quotient μ). *Beweis:* Fixpunktgleichung
cL² + dL = aL + b ⇒ (a − cL)(x − L) = ax + b − L(cx+d). ∎
Bedeutung: In Möbius-Modellen ist die ganze Nicht-Normalität in drei expliziten Bausteinen: Basiswechsel V (Monome ↔ e_j, Binomialmatrix
mit den Fixpunkten), Diagonale μ^j, und die Toeplitz-Faktoren L_N((cx+d)^{−N}) (Koeffizienten d^{−N}C(N+k−1,k)(−c/d)^k). Für
asymptotische Analysen (N → ∞, j/N, k/N fest) stehen die Krawtchouk-Asymptotiken zur Verfügung (Sattelpunkt bzw. Riemann–Hilbert:
Baik–Kriecherbauer–McLaughlin–Miller, Discrete Orthogonal Polynomials, Ann. Math. Stud. 164 (2007) [S]; Dai–Wong (Krawtchouk global,
Chin. Ann. Math. B) [S]) — mit „Band/Lücke/Sättigung“-Bereichen, die dem Streckband von Strang B entsprechen könnten (Vermutung).

---

## 2. Alternative Beweiswege — Prüfung

### 2.1 Positivität, Vorzeichen, Monotonie in N (x01, numerisch) — trägt nicht

| Prüfung | Ergebnis |
|---|---|
| Vorzeichen der Kneser-Koeffizienten | erste negative bei k = 4, 7, 10, 11, …; Anteil negativ → 0.498; Vorzeichenwechsel im Mittel alle 2.327 Grade (π/arg L = 2.349) |
| Gewicht der negativen Koeffizienten Σ_{a_k<0}\|a_k\|ρ^k | **0.00391** von Σ\|a_k\|ρ^k = 2.24047 (schon bei k ≤ 50 voll erreicht) |
| Vorzeichen von p_N gegen Kneser | gleich bis k/N ≈ 0.65–0.78 (N = 60…300), darüber Abweichungen (Kopf) |
| Vorzeichen von p_N − f_K und p_{N+1} − p_N | periodisch (Periode ≈ 4.7 = 2π/arg L): der Modus δ_{−1}, kein fester Kegel |
| Monotonie in N | a₀^{(N)} − a₀^K und p_N(1) − f_K(1) wechseln das Vorzeichen bei N = 2–3, 7, 14–16, 51–53, 215–217 (Phase von c_N dreht wie −2.116 log N; nächster Wechsel ≈ 216·e^{π/2.116} ≈ 950); \|a₀^{(N)} − a₀^K\| ist nicht monoton |
| max_{\|x\|=r}\|p_N\| | **für alle geprüften N (50…1000) und r = 0.5ρ, 0.9ρ, ρ bei x = r > 0 angenommen**, = p_N(r) (2.2326515 bei r = ρ) |
| Hankel-Determinanten det[a_{i+j}] | Vorzeichen −, −, +, +, −, −, + (n = 2…8) — keine Totalpositivität (Stieltjes) |

*Bewertung.* Eine Kegel-/Positivitätsstruktur, die eine a-priori-Schranke liefern könnte, gibt es in den Koeffizienten nicht — die
Oszillation mit arg L ist die Signatur der Randsingularität bei L (Pringsheim: positive Koeffizienten hätten die Singularität auf ℝ⁺).
Monotonie in N ist durch die logarithmische Phasendrehung von c_N widerlegt. **Rest-Idee:** „fast positiv“ (0.17 % negatives Gewicht) ⇒
sup_D|p| ≤ p(ρ) + 2·(negatives Gewicht) für jedes p mit Kneser-nahen Koeffizienten; das ersetzt aber nur sup_D durch einen reellen Wert und
braucht weiter Koeffizientenkontrolle.

### 2.2 Ersatzform bei a₀ (Taylor-Abschnitt bei a₀) — exakte Identität, aber keine Kontraktion

*Identität (bewiesen).* Für den linearen Kern: e_N := Π_N f_K − f_K erfüllt T_N[J_∞e_N] = 0; mit x = h(y) (h(a₀) = 0, lokal biholomorph) und
f'∘f∘h = f' folgt e_N(y) + f'(y)e_N(h(y)) = O((y−a₀)^{N+1}), also (da Π_N f_K ein Polynom vom Grad ≤ N ist)
  **Π_N f_K = T^{a₀}_N f_K − T^{a₀}_N[f'·(e_N∘h)].**
Nichtlinear ebenso p_N = T^{a₀}_N[exp∘p_N⁻¹] (C, README „Ersatzform“). Mit Cauchy auf D(a₀,R'):
|Π_N f_K(y) − f_K(y)| ≤ |f_K(y) − T^{a₀}_N f_K(y)| + sup_{D(a₀,R')}|f'|·sup_{h(D(a₀,R'))}|e_N| / (1 − |y−a₀|/R').
*Folgerung:* Konvergenz im Halbmond (Inneres von D(a₀,R₀)) folgt aus Kontrolle von e_N auf h(D(a₀,R')) — das liegt bis auf einen
schmalen Sektor bei L (Lemma 1.1: Richtungen 187.7°–204.9°) in D; dort braucht man e_N knapp ausserhalb W̄ (Bernstein–Walsh: Faktor
e^{N g}). Für festes y im Halbmond also Konvergenz, solange ε_N e^{c(y)N} klein ist — **für alle N offen**, durch Daten nicht widerlegt (2.3).
*Als Beweisweg (Kontraktion von p ↦ T^{a₀}_N[exp∘p⁻¹]):* Linearisierung −T^{a₀}_N K̃ mit K̃δ = f'·(δ∘h); formal Spektrum ⊃ Scheibe vom Radius
|μ|^{3/2} = 1.268 (Bourdon, §1.4) ∋ −1, auch nach Dämpfung p ↦ (1−θ)p + θ·Φ(p) bleibt ein Spektralpunkt vom Betrag 1 + 0.268θ > 1 —
**keine Kontraktion in irgendeinem Raum, der die Konkurrenz-Kernfunktionen enthält**; in Räumen mit Potenzgewicht bei L ist K̃ wegen der
Drehung unbeschränkt (Lemma 1.1). Trägt nicht.

### 2.3 Gegenprobe Halbmond (x02, numerisch) — Korrektur an Strang C

Rauschpegel ν_N(x) = Σ_k|b_k||x|^k·10^{−59} der gespeicherten Koeffizienten (60 Stellen) und empirisch (letzte Stelle zufällig gestört):

| N | x = 1.8: ν / empirisch / \|C-Wert − Grenzwert\| | x = 1.84: ν / emp. | x = 1.4+0.9i: ν / emp. |
|---|---|---|---|
| 500 | 3.8e-21 / 3.2e-23 / < 1e-11 | 1.3e-16 / 7.9e-19 | 3e-37 / 9e-39 |
| 600 | 1.7e-10 / 1.5e-11 / < 1e-11 | 4.9e-5 / 4.2e-6 (C: Abweichung 1.2e-6) | 6e-30 / 1e-30 |
| 650 | 8.4e-6 / 1.3e-6 (C: 4.6e-7) | 6.8 / 1.1 (C: 0.38) | 8e-27 / 8e-28 |
| 700 | 2.2 / 0.014 (C: 0.06) | 5.3e6 / 9.7e3 (C: 1.5e5) | 4e-23 / 4e-25 |
| 800 | 2.4e10 / 2.0e9 (C: 7.2e8) | 4.4e17 / 3.8e16 (C: 1.5e16) | 3e-16 / 2e-17 |
| 1000 | 3.2e31 / 5.2e29 (C: 2.7e29) | 4.2e40 / 8.7e38 (C: 4.5e38) | 0.105 / 0.0024 (C: 1.8e-3) |

Alle von C als Divergenz gemeldeten Werte liegen im Rauschen. Oberhalb des Rauschens konvergiert Π_N f_K im Inneren des Halbmonds bis
N ≈ 600 (1.8), 1000 (1.6, 1.4+0.9i). Echt ist dagegen Bs b42 am Scheitel a₀+R₀ = 1.848 (f_K dort über exp∘f_K∘Log genau; Fehler
1.3e-11 → 4.9e-11 → 1.4e-10 für N = 300, 400, 500 bei ν ≈ 1e-37 … 1e-15), Wachstumsrate ≈ 0.011/Grad.
Datei: `x02_halbmond_rauschen.py`, `x02.log`.

### 2.4 Lemma (Selektion nur am Rand; bewiesen unter (K)) — warum kein „weicher“ Weg tragen kann

*Aussage.* Sei U ⋐ D offen (Abschluss in D, also fern von L, L̄). Für jedes c ∈ ℂ, |c| < c₀(U), sei ψ_c(w) := c e^{2πiw} + c̄ e^{−2πiw}
(1-periodisch, reell auf ℝ, ψ_c(w + ½) = −ψ_c(w)) und g_c := T∘(id + ψ_c)∘A. Dann ist f_c := g_c⁻¹∘f_K∘g_c eine exakte, reelle Lösung von
f∘f = exp, holomorph auf einer Umgebung von Ū, f_c → f_K gleichmässig auf Ū für c → 0, und f_c ≠ f_K für c ≠ 0.
*Beweis.* (1) g_c vertauscht mit exp, wo definiert: exp∘g_c = T(A + ψ_c(A) + 1) = T((A+1) + ψ_c(A+1)) = g_c∘exp (T(w+1) = e^{T(w)},
A∘exp = A + 1). (2) Also f_c∘f_c = g_c⁻¹∘f_K∘f_K∘g_c = g_c⁻¹∘exp∘g_c = exp. (3) A ist holomorph auf D, A(Ū) kompakt, T ganz ⇒ g_c → id
gleichmässig auf Umgebungen von Ū und von f_K(Ū) (dort A über A∘f_K = A + ½ fortgesetzt; f_K(Ū) ∌ L, L̄), für kleines |c| injektiv ⇒
f_c holomorph um Ū und → f_K. (4) In Abel-Koordinaten ist f_c ↔ w ↦ φ_c⁻¹(φ_c(w) + ½), φ_c = id + ψ_c; das ist w + ½ genau dann, wenn
ψ_c(w + ½) = ψ_c(w), hier aber ψ_c(w + ½) = −ψ_c(w) ≢ ψ_c(w). (5) Realität: ψ_c reell auf ℝ, A, T reell. ∎
*Folgerung.* Jede Eigenschaft, die f_K mit allen f_c auf Ū teilt — lokale Konsistenz (p∘p = exp nahe 0), Schranken auf Kompakta,
Normalfamilien-Argumente auf Ū, Kontraktion auf Räumen über Ū —, kann die Grenzfunktion nicht festlegen. Die Auswahl von Kneser geschieht
ausschliesslich am Rand bei L, L̄ (dort hat f_c den starken Term ~c(x−L)^{1−γ}, unbeschränkt). Da die einzige Auswahlquelle das Abschneiden
ist und dieses nach B/C/D über die Eckschicht |x − L| ~ N^{−0.883} wirkt, **muss jeder Beweis die Eckschicht auflösen** — also eine
Aussage vom Typ (H₀)/(E_K) machen. Das ist kein Beweis, dass (H) nötig ist, aber ein Beweis, dass Eck-Information nötig ist.
*Schlupf (heuristisch, aus der gemessenen Ein-Moden-Struktur p_N − f_K ≈ c_Nδ_{−1} + konj. mit Abschneidung bei |x−L| ~ N^{−0.883}):*
sup_D|p_N − f_K| ≈ |c_N|·N^{0.883·3.447} = |c_N|N^{3.04}; für die README-Rückführung R1 (nur sup_D|p_N| < ∞, dann Montel + Eindeutigkeit)
genügte also c_N = O(N^{−3.04}) — wahr ist N^{−7.9}. Der Schlupf (Faktor N^{4.86}) ist derselbe wie in B-Max (κ < 4.69).

### 2.5 Direkter Weg in Kneser-Koordinaten (θ_N) — reduziert sich auf README-R1

p_N = T(A + ½ + θ_N(A)) definiert θ_N nur, solange p_N(T(w)) im Definitionsbereich von A liegt. Da p_N ein Polynom ist, ist p_N(L) ≠ L
(|p_N(L) − L| ~ N^{−4.85}), also θ_N(w) = A(p_N(T(w))) − w − ½ → A(p_N(L)) − w − ½ für Im w → ∞: θ_N ist **nicht** auf der ganzen
Halbebene klein, sondern nur bis zur Höhe Y_N ≈ (0.883/(2π·…))·log N (dort |T(w) − L| ~ N^{−0.883}; Ein-Moden-Grösse |c_N|e^{2πY_N} ~ N^{−4.86}
— konsistent). Die Gleichung θ(w) + θ(w + ½ + θ(w)) = ε_N(w) (ε_N aus p_N∘p_N − exp) bestimmt θ_N nur bis auf den Kern (antiperiodische
Moden, §2.4); die Auswahl kommt aus der Höhe Y_N, also aus der Eckschicht. Ausserdem ist ε_N klein nur, wo p_N(D_r) im Konvergenzbereich
liegt (|x| ≲ 1.1, da f_K(1.1) ≈ 1.8 schon nahe am Halbmond-Rand). Ergebnis: kein eigenständiger Weg; Kontrolle von θ_N bis zur Höhe Y_N
ist äquivalent zu sup_D-Schranken für p_N (R1) bzw. zu (H₀).

### 2.6 Weitere geprüfte Ideen (Kurzurteil)

- **Variationelle Charakterisierung / Galerkin-Koerzivität:** T_N ist in jedem diagonal gewichteten ℓ² ein Orthogonalprojektor. Wäre
  J_∞ in einer solchen Norm koerziv (Re⟨J_∞u,u⟩ ≥ c‖u‖²), so wäre sup_N‖J_N⁻¹‖ ≤ 1/c (Céa). Das widerspricht Ds exponentiell wachsender
  Eckmode (in ℓ²_ρ gemessen; D: auch alle getesteten Polynomgewichte (1+k)^s werden ab N ~ 240–500 erfasst). Für Normen, in denen T_N nicht
  orthogonal ist (z. B. H² auf W), ist T_N exponentiell unbeschränkt (Taylor-Reihe bei 0 divergiert auf D(a₀,R₀)∖D̄). **Trägt nicht.**
- **Hermite/Padé an 0 und a₀:** äquivalent zur Ersatzform (§2.2); die Zwei-Punkt-Hermite-Interpolation von f_K konvergiert nur auf der
  Lemniskate |x(x−a₀)| < |L||L−a₀|, p_N dagegen auf D (C §3.5) — p_N ist kein Hermite-Interpolant. **Kein neuer Hebel.**
- **Vitali in der Basis b:** p_N(b) und f_K(b) analytisch in b (Paulsen 2019 [S], Kneser für b > η); Konvergenz auf einer Menge mit
  Häufungspunkt + lokal gleichmässige Schranken in b ⇒ Konvergenz für alle b. Die Schranken sind aber genau (i) aus R1 (für komplexe b),
  und für b ↓ η entartet das Problem (L_b, L̄_b verschmelzen, parabolisch). Nahe η wird die Eckmode schwach (Ds c ≈ 0.074 ln|μ_b| → 0) und
  die Rate stark (Es Exponent 2s_b − 1: 27.5 bei b = 1.5) — ein Modell-Hinweis für M, kein Umweg. **Trägt nicht allein.**
- **Strukturierte Störungstheorie für (H₀) → (H):** unstrukturiert verstärkt J⁻¹EJ⁻¹ eine Störung der Grösse ε_N um σ₁² ~ e^{0.024N}
  (versagt ab N ~ 2000, passend zu ARCHITEKTUR 6.5); nötig ist ⟨v₁, E u₁⟩ ≪ 1/σ₁ für die kohärenten Störungen E = J_N(P) − J_N(f) —
  das misst Bs b46. Literatur: strukturierte Pseudospektren (Hinrichsen–Pritchard, Mathematical Systems Theory I, 2005, Kap. 5 [U];
  Rump, SIAM J. Matrix Anal. Appl. 28 (2006) [U]) — liefert nur Begriffe, keine Schranke.

---

## 3. Neuer Werkzeugbefund: die Eckmoden verschwinden in RECHTECKIGEN Abschnitten (für O, W)

### 3.1 Lemma 3.1 (Rechteck-Abschnitt ⇒ Splitting; bewiesen, elementar)

Sei R = [J_N; B] ∈ ℂ^{(N+1+d)×(N+1)} mit vollem Spaltenrang (B = Zeilen N+1…N+d von J_∞|_{P_N}), R⁺ = (R*R)⁻¹R* = [R⁺_top | R⁺_bot]
(Spalten N+1 | d) und J_N invertierbar. Dann ist S := I_d − B R⁺_bot invertierbar und
  **J_N⁻¹ = R⁺_top + R⁺_bot·S⁻¹·B R⁺_top,**  also  **σ_{d+1}(J_N⁻¹) ≤ ‖R⁺_top‖ ≤ 1/σ_min(R).**
*Beweis.* R⁺R = I gibt R⁺_top J_N = I − R⁺_bot B. Sylvester: det S = det(I − R⁺_bot B) = det(R⁺_top)·det(J_N) ≠ 0, da R⁺_top =
(R*R)⁻¹J_N* invertierbar. Probe: (R⁺_top + R⁺_bot S⁻¹BR⁺_top)J_N = I − R⁺_bot B + R⁺_bot S⁻¹(B − BR⁺_bot B) = I − R⁺_bot B + R⁺_bot S⁻¹SB = I.
Der zweite Summand hat Rang ≤ d; Weyl: σ_{d+1}(X + F) ≤ σ₁(X) für Rang F ≤ d. ∎
Dasselbe gilt wörtlich für I + K_N statt J_N (dann (a′) ⇐ σ_min(T_{N+d}(I + K)|_{P_N}) ≥ 1 − o(1)).

**Korollar 3.2 (Zerlegung der J-Projektion; bewiesen).** Mit der Kleinste-Quadrate-Projektion Π̂g := R⁺T_{N+d}J_∞g (reproduziert P_N) gilt
für jede Eingabe g
  **Π_N g = Π̂g + R⁺_bot·S⁻¹·ρ̂(g),  ρ̂(g) := [J_∞(Π̂g − g)]_{N+1…N+d}**  (Überlauf-Residuum der Kleinste-Quadrate-Lösung).
*Beweis.* w := T_N J_∞g. R⁺(w, [J_∞g]_over) = Π̂g ⇒ R⁺_top w = Π̂g − R⁺_bot[J_∞g]_over, also BR⁺_top w = BΠ̂g − (I − S)[J_∞g]_over und
S⁻¹BR⁺_top w = S⁻¹(BΠ̂g − [J_∞g]_over) + [J_∞g]_over. Einsetzen in Lemma 3.1. ∎
*Bedeutung:* Π̂ hat keine Eckmoden (‖R⁺‖ ≤ 1/σ_min(R)); Π̂k_ζ ist wie Π_Nk_ζ rational in ζ mit Polen nur in 0, f(0) und O(ζ^{−N−2}) (Cs Lemma 12.1
gilt wörtlich, R⁺ ist eine feste Matrix) — Bernstein–Walsh/Poisson–Jensen sind also auf Π̂ anwendbar. Der Eckteil ist **S⁻¹ angewandt auf ein
d-dimensionales Konsistenz-Residuum** — (E_K) wird zur Frage, ob ρ̂(k_ζ) bzw. ρ̂(f_K) die kleinen Singulärrichtungen von S meidet.

**Bedeutung.** Die gesamte Eck-Problematik steckt in der d×d-Matrix S_N (und den d Funktionalen B R⁺_top: „Überlauf der Kleinste-Quadrate-
Lösung in die Gleichungen N+1…N+d“). Der Rest R⁺_top ist durch den **Injektionsmodul des rechteckigen Abschnitts** beschränkt — eine
Aussage über J_∞ auf Polynomen ohne Abschneidekante. Toeplitz-Analogon: T_n(z − λ) (|λ| < 1) hat σ_min = O(|λ|^n), der Rechteck-Abschnitt
P_{n+1}T(z−λ)P_n ist gleichmässig injektiv — die kleinen Singulärwerte sind ein reiner Rand-/Index-Effekt der quadratischen Abschneidung
(Böttcher–Grudsky ELA 1999 [O]; Rechteck-/„uneven“ Abschnitte: Hansen, J. Funct. Anal. 254 (2008) [S], Colbrook–Roman–Hansen, Phys. Rev.
Lett. 122 (2019) [S], Lindner 2006 [S]).

### 3.2 Numerik (x04, x05; kubisches Modell, exakt, arb mit 2.6N + 400…500 Bit; ℓ²_ρ in Bild und Urbild)

σ_min(T_{N+d}J_∞|_{P_N}) (kleinste zwei Werte), d = Zahl der Zusatzzeilen:

| N | d = 0 (quadratisch) | d = 1 | d = 2 | d = 3 | d = 5 | d = 2N (= volle Abbildung) |
|---|---|---|---|---|---|---|
| 120 | 0.7198, 0.7329 | — | 0.7329 | — | 0.7329 | 0.7329 |
| 200 | 0.5399, 0.7328 | 0.7168, 0.7328 | **0.7328** | 0.7328 | 0.7328 | — |
| 300 | 0.2704, 0.5545 | 0.4747, 0.7063 | 0.6521, 0.7328 | **0.7328** | 0.7328 | 0.7328 |
| 400 | 0.1099, 0.2899 | — | 0.3829, 0.6034 | (d = 4: 0.6726, 0.7328) | **0.7328** (d = 5 … 13) | — |
| 500 | 0.0404 (= 1/σ₁, W) | — | — | (d = 4: 0.4170, 0.6036; d = 6: 0.6633, 0.7328) | **0.7328** (d = 8 … 20) | — |

Kneser selbst (x07, m = e^x/f_K', Koeffizienten aus der 4000er-Datei als exakte Daten): N = 40: 0.73195 für d = 0, 2, 5 (= 1/1.3663 = 1/‖σ‖_∞).
Benötigte Zusatzzeilen (bis der Bulk-Wert 0.7328 auf 4 Stellen erreicht ist) d(N) = 2, 2, 3, 5, 8 für N = 120, 200, 300, 400, 500 —
N/60 … N/100, wächst etwas schneller als linear (Vorsicht bei Extrapolation; für Lemma 3.1 genügt jede Schranke, z. B. N = 500, d = 6:
σ₇(J_N⁻¹) ≤ 1/0.6633 = 1.51). Zum Vergleich O: r(N) = 2 (N ≤ 400), 4 (600–800), 6 (1000) grosse Werte von (I+K_N)⁻¹ — d zählt auch die Moden
knapp über dem Bulk mit.

0.7328 = 1/1.36460 = min_{Γin}|m| = 1/‖σ‖_∞ (Ds Rest-Konstante im kubischen Modell). **Mit 2–3 Zusatzzeilen sind die Eckmoden vollständig
verschwunden**, und die volle Abbildung J_∞|_{P_N} (d = 2N) hat denselben Injektionsmodul. Kontrolle Lemma 3.1 (d = 6): σ(J_N⁻¹) per Formel
= 1.8524 / 1.3646 (N = 200), 3.6983 / 1.8035 / 1.3646 (N = 300) — exakt Ds Werte; ‖R⁺_top‖ = 1.3646 in beiden Fällen.
(Die Singulärwerte von S selbst sind bei d = 6 trivial winzig, 1e-24…1e-32, weil die überzähligen Zeilen fast linear abhängig sind; die
Faktoren R⁺_bot (1.7e-4) und BR⁺_top (2.9e-4) kompensieren das. Für eine saubere d×d-Analyse d minimal wählen bzw. orthonormieren.)

### 3.3 Folgerung für (a′)/(E_R) und Zusammenhang mit der lokalen Mellin-Analyse

Statt (a′) für den quadratischen Abschnitt genügt **(a″): inf_{δ∈P_N} ‖T_{N+d}J_∞δ‖_{ℓ²_ρ}/‖δ‖_{ℓ²_ρ} ≥ c > 0 gleichmässig in N, für ein
d = d(N) ≥ r(N)** (numerisch c = min_{Γin}|m|, d ≤ 3 bis N = 300). Das ist eine Ungleichung vom Gårding-Typ für J_∞ = C_f + M_m auf
Polynomen — ohne Abschneidekante, also für Lokalisierungsargumente zugänglich: auf Γin (Komposition glättet: f(Γin) ⊂ D) dominiert M_m mit
|m| ≥ 0.7328; auf Γout dominiert C_f; an L ist das lokale Mellin-Symbol μ^s + μ auf der H²-Geraden Re s = −½ von 0 weg:
**min_t |μ^{−1/2+it} + μ| = 1.1724 = |μ|** (numerisch, Minimum für t → ∞; die Nullstelle s = 1 + γ liegt bei Re s = 5.447, weit rechts).
Die Ecke ist für den Rechteck-Abschnitt also kein Hindernis; die exponentielle Eckmode ist ein Index-Effekt der QUADRATISCHEN Abschneidung.
(Lemma 3.1 verkraftet jedes d; wie schnell d(N) wachsen muss, siehe 3.4.)

*Symbolbild (Skizze, Ms radiales Modell f(e^{iφ}) = s(φ)e^{iφ}):* dort ist J_∞ = Op(σ) mit Kohn–Nirenberg-Symbol σ(φ,k) = s(φ)^k + m(φ).
Auf Γin (s < 1) gilt σ → m für k → ∞, |m| ≥ min_{Γin}|m|; auf Γout (s > 1) ist |σ| gross; am Fixpunkt ist s(φ)^k = e^{k s'(θ)(φ−θ)+…} > 0
reell, also |σ| ≥ min_{r>0}|r + μ| = |μ| (da Re μ > 0) — **das Symbol ist überall ≥ min_{Γin}|m| = 0.7328 im Betrag**, genau der gemessene
Injektionsmodul. Ein Beweis von (a″) bräuchte einen Gårding-/Calderón–Vaillancourt-Satz für dieses exotische Symbol (φ-Ableitungen kosten
einen Faktor k nahe θ, Klasse S⁰_{0,1}) — dort ersetzt man die Symbolrechnung durch die exakte Mellin-Analyse (Dilatation am Fixpunkt) und
klebt per lokalem Prinzip (Simonenko, Gohberg–Krupnik; Roch–Santos–Silbermann 2011 [S]). Literatur zu „nach unten beschränkt / abgeschlossenes
Bild“ für Kompositionsoperatoren (nur Selbstabbildungen): Cima–Thomson–Wogen, Indiana Univ. Math. J. 24 (1974/75) [S] (H²: abgeschlossenes Bild ⇔
Radon–Nikodym-Ableitung des Randbildmasses von 0 weg), Zorboska, Trans. AMS 344 (1994) 791–801 [S] (gewichtete Bergman-Räume), Galanopoulos–
Panteris (Hᵖ) [S].

### 3.4 Die eckfreie Projektion Π̂ hat KEINE Spitze und konvergiert wie Π_N (x07, x08; numerisch)

**Kneser selbst (x07, N = 300, Koeffizienten der 4000er-Datei als exakte Daten, 1280 Bit):** quadratisch σ_min = 0.27402, 0.55910 (kubisch
0.2704, 0.5545); **mit d = 3 Zusatzzeilen σ_min = 0.73191 (vierfach) = 1/1.36630** = Ds Kneser-Rest-Konstante. Projektionen auf f_K:

| | max_{|z|=0.98ρ}|P − f_K| | |P(L) − L| |
|---|---|---|
| Π_N f_K (quadratisch) | 8.30e-14 | 2.644e-12 (A: 2.6e-12 ✓) |
| Π̂ f_K, d = 3 | 7.97e-14 | 2.578e-12 |
| Π̂ f_K, d = 6 | 7.13e-14 | 2.413e-12 |
| Π̂ f_K, d = 10 | 6.23e-14 | 2.294e-12 |

⇒ Der Eckteil R⁺_bot S⁻¹ρ̂(f_K) aus Korollar 3.2 macht bei N = 300 nur ~10 % aus; Π̂ f_K konvergiert mit derselben Rate.

**Pol-Antworten R(η) = |Πk_η − T_Nk_η|_{ℓ¹_ρ} (x08, kubisch; Cs Grösse R_N, Werte quadratisch = Cs Kneser-Werte c23 auf 1 %):**

| η | N = 200: quadr. / rechteckig (d=3) | N = 300 (d=4) | N = 400 (d=6) | Faktor 300→400 quadr. / rechteckig |
|---|---|---|---|---|
| −ρ (ferner D-Bogen) | 17.98 / 10.44 | 45.17 / 15.86 | **124.6 / 20.57** | 2.76 / 1.30 |
| D-Bogen t = 0.5 | 22.78 / 14.21 | 53.18 / 21.35 | 131.8 / 25.18 | 2.48 / 1.18 |
| D-Bogen t = 0.1 | 40.0 / 33.8 | 68.6 / 54.0 | 109.1 / 71.5 | 1.59 / 1.32 |
| Ecke L | 67.12 / 66.59 | 105.8 / 105.1 | 145.8 / 144.7 | 1.38 / 1.38 |
| Schlitz L(1+N^{−0.883}) | 14.93 / 14.79 | 22.72 / 22.51 | 30.49 / 30.17 | 1.34 / 1.34 |
| A-Scheitel | 0.0111 / 0.0727 | 0.00113 / 0.0114 | 0.000117 / 0.00343 | fallend / fallend |

**Befund:** Die exponentielle Spitze auf dem fernen D-Bogen (Ds/Cs ≈ 1.0115^N, Kreuzung mit der Ecke bei N ≈ 420) ist **vollständig ein
Effekt der quadratischen Abschneidung**: Π̂ hat dort polynomiales Wachstum (≈ N^{0.9…1.0}), Ecke und Schlitz sind unverändert (N^{1.1}).
Für Π̂ scheint also Cs globale Bedingung (S) (polynomial auf ganz ∂W) zu gelten — damit trüge Cs Satz 12.3 (Bernstein–Walsh, bewiesen) für
Π̂ ohne jede Eck- oder Poisson–Jensen-Analyse (Π̂k_ζ ist rational mit Polen nur in 0, f(0), O(ζ^{−N−2}); Korollar 3.2).
**Vorschlag für die Architektur (Vermutung, numerisch gestützt bis N = 400):** Π_N f_K = Π̂ f_K + R⁺_bot S⁻¹ ρ̂(f_K). (i) Π̂ f_K → f_K über (S) für Π̂
(+ (F), Satz 12.3); (ii) Eckteil: |R⁺_bot S⁻¹ ρ̂(f_K)| → 0 — eine d-dimensionale Aussage (d ≈ N/60 … N/100). (S) für Π̂ verlangt eine
Pol-Antwort-Schranke für eine Matrix mit gleichmässig beschränkter Links-Inversen — das ist die Stelle, an der (a″) und Kohärenz zusammenkommen.

---

### 3.5 Satz X1 — Reduktion der linearen Vermutung auf (a″) und eine Eckbedingung für den kanonischen Eckteil (Beweis der Reduktion vollständig)

**Rahmen (P, ARCHITEKTUR §9.3):** K' = D̄ ∪ f_K(D̄), g' Green-Funktion von Ĉ∖K' (Pol ∞), Aussenwinkel bei L 141.69° = π − arg μ, α' = 0.787172,
g'(L(1+u)) ≥ c' u^{1/α'} (c' ≥ 0.754 für u ≤ 0.204), Schlitze S_L(δ), S_L̄(δ) bis zur Niveaulinie Γ'_δ = {g' = δ}, ω'_ζ harmonisches Mass.
(F_{K'}): f_K holomorph um K' bis auf die Schlitze, stetig bis an die Ufer, Sprung |Δf(L(1+u))| ≤ C_J u^{σ}, σ = Re(1+γ) = 5.447 [P-I, aus (K) +
Univalenz]. Alle Vektoren in ℓ²_ρ-Skalierung, |·|₁ := ℓ¹_ρ-Norm. J_∞ = J_∞(f_K), R_N := T_{N+d}J_∞|_{P_N} (d = d(N) ≤ N), Π̂ := R_N⁺T_{N+d}J_∞,
kanonischer Eckteil c_N(g) := Π_N g − Π̂g (= R⁺_bot S⁻¹ρ̂(g), Kor. 3.2).

**Annahmen.**
- **(a″)** σ_min(R_N) ≥ c₀ > 0 für alle N ≥ N₀ (numerisch c₀ = min_{Γin}|m| = 0.7319 bei Kneser, N = 40, 300; kubisch bis N = 500, §3.2).
- **(Univ)** f_K univalent und bi-Lipschitz auf einer Umgebung von D̄ ohne kleine Scheiben um L, L̄ (P: numerisch, min|f(x)−f(y)|/|x−y| = 0.3456).
- **(Ẽ)** Es gibt C, κ, ε > 0 mit |c_N(k_ζ)|₁ ≤ C N^κ e^{−εN g'(ζ)} für ζ ∈ S_L(δ) ∪ S_L̄(δ) ∪ Γ'_δ.

**Satz X1.** Unter (F_{K'}), (Univ), (a″), (Ẽ) gilt für N ≥ N₀
  α_N := |Π_N f_K − q_N|₁ ≤ C₁ N^{2}·N^{−α'(σ+1)} + C₂ N^κ (εN)^{−α'(σ+1)} + C₃ N^{2+κ} e^{−εNδ},   α'(σ+1) = 5.075,
also α_N → 0 für κ < 5.07, und sup_{D̄}|Π_N f_K − f_K| ≤ α_N + sup_{D̄}|f_K − q_N| → 0 (lineare Vermutung (LP)).
*Beweis.* (1) *Cauchy-Darstellung* (C Lemma 9.2/9.4, P-I): Π_N reproduziert P_N, also Π_N f_K − q_N = Π_N t_N mit dem Kneser-Schwanz t_N, und mit der
Darstellung t_N = (1/2πi)∫_{∂Ω'} μ(ζ)[k_ζ]_{>N}dζ (Ω' = {g' < δ} ohne Schlitze; μ = Δf auf den Schlitzufern, μ = f_K auf Γ'_δ; Fubini für |ζ| ≥ r₁ > |a₀|)
gilt Π_N f_K − q_N = (1/2πi)∫ μ(ζ) y_ζ dζ, y_ζ := Π_N k_ζ − T_N k_ζ.
(2) *Kanonische Zerlegung:* y_ζ = ŷ_ζ + c_N(k_ζ), ŷ_ζ := Π̂k_ζ − T_N k_ζ = R_N⁺T_{N+d}J_∞[k_ζ]_{>N} (wegen R_N⁺R_N = I).
(3) *Rationalität:* T_{N+d}J_∞k_ζ = T_{N+d}[1/(ζ − f_K)] + T_{N+d}[m/(ζ − x)]; die Taylor-Koeffizienten bei 0 sind rational in ζ mit Polen nur in
ζ = f_K(0) = a₀ bzw. ζ = 0; R_N⁺ ist eine feste Matrix. Also sind ŷ_ζ und c_N(k_ζ) = y_ζ − ŷ_ζ koordinatenweise rational mit Polen nur in
{0, a₀} ⊂ int K', und beide sind O(ζ^{−N−2}) bei ∞ (ŷ_ζ nach (2), y_ζ nach C Lemma 12.1). Also sind log|ŷ_ζ|₁ + (N+2)g'(ζ) und
log|c_N(k_ζ)|₁ + (N+2)g'(ζ) subharmonisch auf Ĉ∖K' und bei ∞ nach oben beschränkt.
(4) *Eckfreier Teil auf ∂K':* für η ∈ ∂K' = Γin ∪ f_K(Γout) ist |η| ≥ ρ und η ∉ f_K(D) (P-II). Dann |T_{N+d}J_∞k_η|_{ℓ²_ρ} ≤ C N^{3/2}: der Teil
m/(η − x) gibt ≤ ‖m‖_∞√(N+d+1)/ρ; für 1/(η − f_K) schätzt man die Koeffizienten über |x| = r_N = ρ(1 − 1/N) ab, |[x^k]|ρ^k ≤ (ρ/r_N)^k/dist(η,
f_K(D̄_{r_N})) ≤ e²·N/(c_f ρ) für k ≤ N + d ≤ 2N (bi-Lipschitz (Univ); nahe L trennt der Winkel 38.3°) — genau P-II, nur mit N+d statt N Zeilen.
Mit (a″): |ŷ_η|_{ℓ²_ρ} ≤ |T_{N+d}J_∞k_η|/c₀ + |T_Nk_η| ≤ C N^{3/2}, also |ŷ_η|₁ ≤ C N² auf ∂K'.
(5) *Bernstein–Walsh* (C Satz 12.2 in Ĉ∖K'): |ŷ_ζ|₁ ≤ C N² e^{−(N+2)g'(ζ)} für ζ ∉ K'.
(6) *Integrale:* Schlitz: 2·(|L|/2π)∫₀^{u₀}C_J u^σ·[C N² e^{−(N+2)c'u^{1/α'}} + C N^κ e^{−εN c'u^{1/α'}}]du ≤ C N² N^{−α'(σ+1)} + C N^κ(εN)^{−α'(σ+1)}
(∫₀^∞u^σe^{−λu^{1/α}}du = αΓ(α(σ+1))λ^{−α(σ+1)}); Niveaulinie Γ'_δ: |μ| ≤ max_{Γ'_δ}|f_K| und g' = δ. ∎
**Was damit gewonnen ist.** Gegenüber der SVD-Zerlegung (C 13.1, D8, P-II/P-IV) braucht man (i) keine Eckmoden-Projektoren Q_N, (ii) keine Aussage über
die Zahl r(N), (iii) keine (E_A)-Auslöschung auf dem A-Bogen (K' hat keinen A-Bogen in f_K(D)), (iv) keine Bulk-Annahme ‖Y_N P_{0.7N}‖: der ganze
eckfreie Teil folgt aus der EINEN Unterschranke (a″). Übrig bleibt (Ẽ) für den d-dimensionalen kanonischen Eckteil.
**(Ẽ) in Randform (C Lemma 13.2/Satz 13.3 + P-III, bewiesen):** (Ẽ) ⇐ ω'_ζ(λ^c_N) ≤ (1 − ε)g'(ζ) auf S_L ∪ Γ'_δ für die vorzeichenbehaftete Randrate
λ^c_N(η) := (1/N)log(|c_N(k_η)|₁/(CN^κ)) auf ∂K'. Nahe L genügt: Γin-Seite λ^c ≤ A|τ|^p (p ≥ 1/α' = 1.2704 ist grenzwertig, log-divergent) UND
f_K(Γout)-Seite λ^c ≤ −b|η − L|^{q} mit q·α' < 1 (dann ist ∫λ s^{−2}ds = −∞ und schlägt die log-Divergenz der Γin-Seite); fern von L: λ^c ≤ λ_max mit
λ_max·sup ω'/g' < 1.

**Wahl von d — Abwägung (numerisch, x09, N = 200):** Satz X1 gilt für jedes d ≤ 2N. Mit d = 2N (kubisch: volle Abbildung) ist (a″) genau Os
(a″)_∞ (unendlicher Operator, keine Kante) — ABER der Eckteil c^{(2N)} = Π_N − Π̂^{(2N)} ist dann nicht mehr klein bzw. abklingend: N = 200,
d = 400: f_K(Γout)-Seite τ = −0.06: Rc = 12.0 (gegen 0.056 bei d = 3), Schlitz v = 3: 0.62 (gegen 0.0075), Ecke L: 41 (gegen 2.1); Π̂^{(2N)}
weicht in breiter Front von Π_N ab (R̂(L) = 96.6 statt 66.6). Mit kleinem d ≈ d(N) ist der Eckteil „minimal“ (klein an L und auf dem Schlitz,
abklingend auf der f_K(Γout)-Seite) — dafür braucht (a″) dann die Kantenanalyse (wie viele Überlaufzeilen fangen die Eckmoden ein).
**Folgerung:** (a″)_∞ [O] + (Ẽ) für c^{(2N)} ist logisch hinreichend, aber (Ẽ) für c^{(2N)} sieht numerisch schlecht aus (kein Abfall bei
f_K(Γout)); der natürliche Weg ist d ≈ d(N) mit (a″)_{d(N)} [Kante] + (Ẽ) für den minimalen Eckteil.

**Satz X2 (einteilige Form; bewiesen).** Unter (F_{K'}), (Univ), (a″)_d und
  Rc_N(L) := |c_N(k_L)|₁ ≤ C N^κ (κ < 5.07),  Θ^c_N := sup_{ζ ∈ S_L(δ) ∪ S_L̄(δ) ∪ Γ'_δ} ω'_ζ(λ^c_N)/g'(ζ) ≤ 1 − ε,  λ^c_N(η) := (1/N)log(|c_N(k_η)|₁/Rc_N(L)) (η ∈ ∂K'),
gilt (Ẽ) mit derselben Konstante und damit (LP) nach Satz X1.
*Beweis.* c_N(k_ζ) ist koordinatenweise rational mit Polen nur in 0, a₀ ∈ int K' und O(ζ^{−N−2}) (Satz X1 (3)); log|c_N(k_ζ)|₁ + (N+2)g'(ζ) ist
subharmonisch auf Ĉ∖K', nach oben beschränkt, stetig bis ∂K'. Poisson–Jensen (C Lemma 13.2 mit K' statt W̄): log|c_N(k_ζ)|₁ ≤ log Rc_N(L) +
N·ω'_ζ(λ^c_N) − (N+2)g'(ζ) ≤ log(CN^κ) − (εN + 2)g'(ζ). ∎
(Satz X2 braucht keine Zerlegung in Moden; Θ^c_N und Rc_N(L) sind direkt aus Randdaten messbar — das ist die gesuchte „Eckbedingung in d Dimensionen“
in einteiliger Form. Die modale Form: 𝒞_N = Σ_{i≤d} s_i u_i v_iᵀ (§3.6), c_N(k_ζ) = Σ_i s_i u_i⟨v_i, T_{N+d}J_∞k_ζ⟩ — jede Komponente ist eine rationale
Funktion von ζ mit denselben Polen; Poisson–Jensen gilt für jede einzeln.)

**Numerik zu Satz X1 auf ∂K' (x09, kubisches Modell, d = 3/4/6/9/12 für N = 200/300/400/500/600, ℓ¹_ρ):**

| Punkt (τ in rad ab arg L) | R̂ (eckfreier Teil Π̂k − T_Nk), N = 200 … 600 | Rc (kanonischer Eckteil), N = 200 … 600 | R (quadratisch) |
|---|---|---|---|
| Ecke L | 66.6 / 105 / 145 / 185 / 226 | 2.12 / 2.33 / 3.41 / 5.39 / 7.48 | 67 … 228 |
| Γin τ = 0.01 | 62 / 97 / 132 / 167 / 202 | 2.50 / 3.04 / 4.92 / 8.51 / 12.9 | 63 … 210 |
| Γin τ = 0.1 | 42 / 67 / 91 / 110 / 127 | 6.8 / 11.8 / 25.7 / 55.5 / 102 | 46 … 210 |
| Γin τ = 0.4 | 20.6 / 31.1 / 39.1 / 55.7 / 58.7 | 16 / 39 / 97 / 210 / 474 | 31 … 512 |
| Γin −ρ | **10.4 / 15.9 / 20.6 / 24.1 / 28.8** | 15 / 41 / 116 / 362 / **1165** | 18 … 1179 |
| f(Γout) τ = −0.01 | 22.5 / 21.4 / 18.1 / 14.6 / 11.5 | 1.07 / 0.90 / 1.03 / 1.26 / 1.36 | 22 … 11 |
| f(Γout) τ = −0.03 | 3.1 / 1.31 / 0.55 / 0.24 / 0.11 | 0.34 / 0.16 / 0.10 / 0.069 / 0.040 | 3.0 … 0.08 |
| f(Γout) τ = −0.06 | 0.25 / 0.041 / 6.9e-3 / 1.3e-3 / 2.3e-4 | 0.056 / 0.010 / 2.4e-3 / 6.2e-4 / 1.4e-4 | |
| Schlitz L(1 + vN^{−α'}), v = 1 | 5.6 / 7.5 / 9.2 / 10.7 / 12.0 | 0.39 / 0.39 / 0.54 / 0.82 / 1.09 | 5.6 … 12.2 |
| Schlitz v = 3 | 0.046 (konstant) | 0.0075 / 0.0061 / 0.0072 / 0.0095 / 0.0113 | 0.046 |

⇒ (i) Der eckfreie Teil ist auf ganz ∂K' polynomial (Maximum in L, ~N^{1.1}; fern ~N^{0.93}; auf f_K(Γout) abklingend) — wie Satz X1 (4) aus (a″)
vorhersagt. (ii) Die gesamte exponentielle Spitze sitzt im Eckteil (fern Rate 0.0117/Grad; bei τ = 0.1 Rate 0.006–0.008, fallend); nahe L und
auf dem Schlitz ist der Eckteil polynomial (~N^{1.3…2}) und um 1–2 Grössenordnungen kleiner als der eckfreie Teil; auf der f_K(Γout)-Seite fällt
er exponentiell. Genau dieses Profil braucht (Ẽ) in Randform (Poisson–Jensen auf ∂K'); Θ für Kneser folgt (x12/x12b).

**Kneser selbst: einteilige Poisson–Jensen-Bedingung im K'-Rahmen für die drei Teile (x12 + x12b; Ps Randgitter p14, Ps kprime-Modul für
τ', g', ω'; Θ = sup_{Schlitz ∪ Γ'_δ} ω'(λ)/g', λ = (1/N)log(Q(η)/Q(L))):**

| N (d) | Θ für R (quadr.) | Θ für R̂ (eckfrei) | Θ für Rc (kanon. Eckteil) | Rc(L) | Ecken-Grenzwert (1/κ₀')∫λ/s² für Rc | Schlitz max ω'/g' (Rc) |
|---|---|---|---|---|---|---|
| 200 (3) | −0.189 | −0.235 | **+0.129** (Γ'_0.05) | 2.13 | −32.7 | −0.383 |
| 300 (4) | −0.117 | −0.183 | **+0.129** | 2.34 | −30.9 | −0.386 |
| 400 (6) | −0.070 | −0.155 | **+0.114** | 3.43 | −29.7 | −0.388 |
| 600 (12) | −0.009 | −0.125 | **+0.106** | 7.53 | −28.3 | −0.388 |

(Θ für R stimmt mit Ps p15 bis auf die Γ'_δ-Diskretisierung überein: P −0.13/−0.09/−0.06.)
⇒ **(Ẽ) gilt numerisch in Randform mit ε ≥ 0.87** (N ≤ 400): der Eckteil ist am Schlitz negativ, sein Ecken-Grenzwert ist stark negativ (−30, weil
er auf der f_K(Γout)-Seite abklingt — K' ist für DIESEN Eckteil also nicht grenzwertig, anders als für Ps volle Eck-Projektion), positiv nur auf der
Niveaulinie durch die ferne Spitze (≤ 0.13). Der eckfreie Teil hat Θ < 0 (wie Satz X1 (4) es aus (a″) erzwingt: polynomial auf ∂K').

### 3.5b Variante d = 2N: (a″) wird zur Aussage über den UNENDLICHEN Operator (Os Teil), Θ^c bleibt klein (Kneser, x12)

| N (d) | Rc(L) | Θ für R̂ | **Θ^c** | λ^c max Γin / min f(Γout) | Ecken-Grenzwert Rc |
|---|---|---|---|---|---|
| 200 (400) | 41.1 | −0.095 | **−0.039** | +0.0007 / −0.016 | −27.1 |
| 300 (600) | 60.6 | −0.069 | **−0.012** | +0.0006 / −0.011 | −25.0 |
| 400 (400 = N) | 79.8 | −0.055 | **+0.020** | +0.0012 / −0.009 | −23.7 |

(N = 400: die ferne Spitze überholt Rc(L) (Rc(−ρ) = 122.7 > 79.8), daher jetzt Θ^c > 0 auf Γ'_0.05; Schlitz −0.03. Erwartet: Θ^c steigt mit N gegen die
Spitzenschranke ≈ 0.1–0.2, wie in allen Varianten.)

Mit allen Zeilen bis 3N ist der Eckteil gross (Rc(L) ≈ 0.6·R(L), wächst ~N^{0.96}) und klingt auf f_K(Γout) nur schwach ab — aber sein Randprofil ist
relativ zu Rc(L) flach, und Θ^c ist bisher NEGATIV. Für grosse N trägt die ferne Spitze (Rate 0.0117 auf dem D-Bogen) auch hier positiv bei (erwartet
Θ^c ≈ 0.1–0.2 wie bei allen Varianten); Läufe N = 400, 600 mit d = N laufen.
**Satz X3 (Variante, Reduktion bewiesen bis auf die markierte Schwanzabschätzung):** (F_{K'}) + (Univ) + **(a″)_∞:** inf_{δ∈P_N}‖J_∞δ‖_{H²(D_ρ)}/‖δ‖ ≥ c_∞ > 0
+ **(Schwanz)** ‖T_{>3N}J_∞δ‖_{ℓ²_ρ} ≤ CN^{−1}‖δ‖ für δ ∈ P_N (aus der Singularitätsentwicklung (K): Koeffizienten von f_K^j und m jenseits 3N sind
O(N·k^{−5.4}ρ^{−k}) — Skizze, nicht ausgeführt) + **Θ^c_{2N} ≤ 1 − ε, Rc_{2N}(L) ≤ CN^κ** ⇒ (LP). Denn (a″)_∞ + (Schwanz) ⇒ (a″)_{2N} mit c₀ = c_∞ − o(1),
dann Satz X1/X2 mit d = 2N (der P-II-Schritt (4) braucht nur k ≤ 3N). Beim kubischen Modell entfällt (Schwanz) (J_∞P_N ⊂ P_{3N}).
**Aufteilung damit:** O: (a″)_∞ (Satz O7 + Injektivität; für Kneser D4|(G2)); X/W: Θ^c (gemessen, Randform) — keine Kanten- und keine Rest-Abschätzung
für endliche Abschnitte mehr nötig.
**(Schwanz) — Skizze mit Zahlen.** δ ∈ P_N, ‖δ‖_{ℓ²_ρ} = 1 ⇒ |δ(y)| ≤ √(N+1)(|y|/ρ)^N für |y| ≥ ρ. Cauchy für F = δ∘f_K über den Kreis |x| = R mit
Schlüsselloch um die radialen Schlitze von L, L̄ (f_K holomorph in |x| < 7.6 ausser auf den Schlitzen, (K)): (i) Kreisteil |F_k|ρ^k ≤ √(N+1)(M_R/ρ)^N(ρ/R)^k;
mit Zeilen bis (1+q)N braucht man (M_R/ρ)(ρ/R)^{1+q} < 1. Gemessen: M_2 := max_{|x|=2}|f_K| = f_K(2) = exp(f_K(log 2)) = **3.4313** (Maximum auf der positiven
Achse, |f_K(2e^{iφ})| fällt für |φ| ≤ 1.18: 3.43, 3.34, 3.08, 2.70, 2.40, 2.27), also für q = 2 (Zeilen bis 3N): (3.4313/ρ)(ρ/2)³ = **0.810 < 1**;
für q = 3 genügt R = 1.65 (u ≤ 0.2, Ps geprüfter (F)-Bereich): 1.98·(ρ/1.65)⁴ = 0.954 < 1. (ii) Schlitzteil: Sprung ΔF ≈ δ'(f)Δf, |Δf| ≤ C_J u^{5.447},
|f_K(L(1+u))|/ρ ≈ 1 + 0.92u ⇒ Beitrag ≲ N^{3/2}∫u^{5.447}e^{−(k−0.92N)u}du ≲ N^{3/2}(k − 0.92N)^{−6.447}; Summe der Quadrate über k > (1+q)N:
O(N^{−4.4}). m·δ analog (m = eˣ/f_K' hat bei L die Singularität t^{γ}). ⇒ ‖T_{>(1+q)N}J_∞δ‖_{ℓ²_ρ} ≤ C N^{−4.4} + exp. klein. (Skizze; offen nur die
Buchhaltung der Schlitzintegrale und M_R für |φ| > 1.18 — dort ist f_K(2e^{iφ}) über exp∘f_K∘Log zu rechnen.)

### 3.5c X3 abgeschlossen (Runde 3, 11.10. 01:35–02:00): Schwanz-Lemma X5, Θ^c-Verlauf, Grenzwert — mit ehrlicher Einschränkung

**Lemma X5 (Schwanz; bewiesen unter den genannten (K)-Fakten).** Sei R := 1.649, u₁ := R/ρ − 1 = 0.1997, β := (M_R/ρ)(ρ/R)^5. Voraussetzungen:
(K1) f_K holomorph auf {|x| < R + 0.01} ohne die radialen Schlitze S_L, S_L̄, stetig bis an die Ufer (B L8.4 / A belegen |x| < |L| + 0.15 = 1.52;
für 1.52 ≤ |x| ≤ 1.66 [Skizze]: f_K = sexp(slog x + ½) ist singulär nur, wo slog singulär ist (L, L̄) oder slog x + ½ ∈ {−2, −3, …}, d. h.
x ∈ {sexp(−2.5), sexp(−3.5), …}; das nächste ist sexp(−2.5) = log(log a₀) = −0.362 ± πi mit |x| = 3.16 — im Kreis |x| < 3.16 also nur L, L̄;
dasselbe trägt die Fassung R = 2, nicht von P geprüft); (K2) Sprung |Δf(L(1+u))| ≤ C_J u^σ, σ = 5.447, C_J ≤ 0.95 für u ≤ 0.2 (A k25, P); (K3) M_R := max_{|x|=R}|f_K| = f_K(R) ≤ f_K(1.65) = 2.7252
(x14: Maximum auf der positiven Achse, |f_K(1.65e^{iφ})| fällt monoton bis 0.171 bei φ = π) ⇒ **β ≤ 0.798** (ρ = 1.37456); (K4) auf den
Schlitzufern |f_K(L(1+u))| ≤ ρ(1 + c₁u), c₁ ≤ 1.2 (lokale Formel f_K = L + μLu + O(u²), |μ| = 1.172); (K5) m = eˣ/f_K' holomorph auf demselben Gebiet
(f_K'(x) = e^{h(x)}/f_K'(h(x)) ≠ 0), max_{|x|=R}|m| < ∞, Sprung |Δm| ≤ C u^{Re γ}.
*Aussage:* für δ ∈ P_N gilt ‖T_{>5N}J_∞δ‖_{ℓ²_ρ} ≤ τ_N‖δ‖_{ℓ²_ρ}, τ_N = C₁√N β^N + C₂N^{−4.4}.
*Beweis.* |δ(y)| ≤ √(N+1)‖δ‖max(1,|y|/ρ)^N, |δ'(y)| ≤ N√(N+1)‖δ‖ρ⁻¹max(1,|y|/ρ)^{N−1}. Cauchy für F = δ∘f_K über den Kreis |x| = R mit Schlüssellöchern um
S_L, S_L̄ (kleine Kreise um L, L̄ verschwinden, F ist dort beschränkt). Kreisteil: |F_k|ρ^k ≤ √(N+1)‖δ‖(M_R/ρ)^N(ρ/R)^k, Quadratsumme über k > 5N
≤ (N+1)‖δ‖²β^{2N}/(1 − (ρ/R)²). Schlitzteil: |ΔF| ≤ sup|δ'|·|Δf| ≤ N√(N+1)‖δ‖ρ⁻¹(1 + c₁u + C_J u^σ/ρ)^N C_J u^σ und |x|^{−k−1}ρ^k = ρ⁻¹(1+u)^{−k−1}
≤ ρ⁻¹e^{−0.833ku} (log(1+u) ≥ u/1.2 für u ≤ 0.2); für k > 5N ist 0.833k − 1.2N ≥ 0.58k, also |F_k^{Schlitz}|ρ^k ≤ C N^{3/2}‖δ‖Γ(σ+1)(0.58k)^{−σ−1},
Quadratsumme über k > 5N: O(N³·N^{−2σ−1}) = O(N^{−8.9}) ⇒ Wurzel O(N^{−4.45}). m·δ: ‖T_{>5N}(mδ)‖ ≤ ‖δ‖Σ_{i≥4N}|m_i|ρ^i (Young), |m_i|ρ^i ≤ C(ρ/R)^i + C i^{−1−Re γ}
(gleiche Kontur) ⇒ O(N^{−4.45}). ∎
*Robustheit:* Der Schlitzteil braucht nur c₁ < 0.833·5 = 4.1 (dann ist 0.833k − c₁N ≥ c k für k > 5N); |μ| = 1.172 lässt viel Reserve, der O(u²)-Rest bei u ≤ 0.2
ist nicht ausgewertet (x14-Rekursion scheitert am Schlitz, weil log L = L) — (K4) also [Skizze, Reserve Faktor 3.5].
*Die Fassung der Leitung (R = 2, Zeilen bis 3N):* f_K(2) = 3.4313 < R³/ρ² = 8/ρ² = 4.2342 ⇔ β₂ = (M₂/ρ)(ρ/2)³ = 0.810 < 1 — das ist die richtige Bedingung,
aber der Beweis braucht dann (K1) bis |x| = 2 (Skizze oben über slog; numerisch: x14 wertet f_K = exp∘f_K∘Log auf dem ganzen Kreis |x| = 2
ohne Störung aus) und eine Sprungschranke |Δf(L(1+u))| bis u = 0.455 (nur für u ≤ 0.2 belegt). Deshalb hier die Fassung
R = 1.649 mit Zeilen bis 5N, d = 4N (Θ^c-Werte für d = 4N und d = 2N sind identisch, s. u.); mit (K1) bis |x| = 2 gilt dasselbe für d = 2N.
*Folgerung:* σ_min(T_{5N}J_∞|_{P_N}) ≥ c_∞ − τ_N, also **(a″)_{4N} ⇐ (a″)_∞** (Os unendlicher Operator). (K3) ist numerisch (eine Maximumsuche auf
einem Kreis), (K1), (K2), (K5) sind die bekannten (K)-Folgerungen; (K4) aus der lokalen Formel.
*Numerische Bestätigung, dass Zeilen jenseits 3N nichts ändern:* Kneser N = 200, d = 400 und d = 800: R̂(L) = 96.6214 und Rc(L) = 41.0750 identisch;
N = 300, d = 600 und d = 1200: 150.8484 / 60.5691 identisch; σ_min = 0.73191 in allen Fällen (x13).

**Θ^c-Verlauf (alle Zeilen; Kneser, Ps ∂K'-Gitter; x13 über Normalgleichungen, x12b):**

| N | d | σ_min(R) | Rc(L) | Rc(−ρ) | λ max Γin | λ min f(Γout) | Schlitz max | Γ'_0.05 | **Θ^c** |
|---|---|---|---|---|---|---|---|---|---|
| 100 | 400 | 0.73191 | 21.01 | 7.43 | — | — | — | — | −0.071 |
| 200 | 200 / 800 | 0.73191 | 41.08 | 18.87 | +0.00074 | −0.0158 | −0.058 | −0.039 | **−0.039** |
| 300 | 1200 | 0.73191 | 60.57 | 45.60 | +0.00064 | −0.0112 | −0.041 | −0.012 | **−0.012** |
| 400 | 400 | — | 79.76 | 122.7 | +0.00120 | −0.0088 | −0.030 | +0.020 | **+0.020** |
| 600 | 600 | 0.73190 | 117.34 | 1139.6 | +0.00379 | −0.0061 | −0.018 | +0.072 | **+0.072** |

**d = N genügt:** N = 200 mit d = 200 und d = 800 gibt dieselben R̂, Rc auf 1e-11 (alle 112 Punkte) — Zeilen jenseits 2N tragen nichts
(Lemma X5 ist hier sehr vorsichtig). Der Lauf N = 600, d = 600 (x13, 461 s) ersetzt den abgebrochenen x12-Lauf.
**Trend:** +0.027/+0.032/+0.026 je 100 Grade (linear ≈ +0.028/100); getragen von der fernen Spitze auf der Γin-Seite (Rc(−ρ)/Rc(L) = 0.46 → 9.7
für N = 200 → 600; Γ'_0.05 ist das Maximum, Schlitz und Ecke bleiben tief negativ). f(Γout)-Abfall wird schwächer (−0.016 → −0.006).
**Grenzwert — ehrliche Einschränkung (x15b):** Mit d ≥ 2N ist c^{(2N)} = c^{(klein)} + D, D = Π̂^{(klein)} − Π̂^{(2N)} (Differenz zweier eckfreier Projektionen,
polynomial auf ∂K', aber OHNE Abfall auf f_K(Γout)). Asymptotisch hat c^{(2N)} also auf der Γin-Seite das Spitzen-/Eckprofil (Ws Messung: λ = 0.0009/0.0024/
0.0047/0.0090/0.0115 bei τ = 2/4/8/30/180°), auf der f_K(Γout)-Seite aber λ = 0. Für dieses Profil gibt Poisson–Jensen in K': **Θ_∞ ≈ 0.82** (mit 20 %
Aufschlag auf das W-Profil 0.99), Ecken-Grenzwert endlich (+3.6, schon durch κ₀' geteilt), weil der Γin-Exponent p ≈ 1.4–1.5 nur knapp über 1/α' = 1.27 liegt — **K' ist für die
Variante d = 2N asymptotisch grenzwertig** (Θ < 1, aber ohne Reserve). Für den MINIMALEN Eckteil (kleines d), der auf f_K(Γout) abklingt (gefittet
λ = −0.0047 bei 1°, −0.053 bei 10°), gibt dasselbe Γin-Profil **Θ_∞ ≈ 0.20 (δ = 0.05) bzw. 0.07 (δ = 0.1)** (inkl. 20 % Aufschlag) mit Ecken-Grenzwert −27.
**Reserve gegen 1 (Antwort auf die Leitung):** linear extrapoliert erreichte Θ^c die Grenze 0.82 bei N ≈ 3300, 1 erst bei N ≈ 3900; der
Grenzwert aus dem Profil ist 0.82 (ohne Aufschlag) bzw. 0.99 (mit 20 % Aufschlag auf Ws Profil) — **Reserve 0.18 bzw. 0.01, also knapp**.
⇒ Zwei ehrliche Varianten: (X3) (a″)_∞ [O, ohne Kante] + Θ^c_{2N} — numerisch Θ ≤ 0.072 bis N = 600, asymptotisch ≈ 0.82 (grenzwertig); (X2) (a″)_{d(N)}
[Kantenaussage] + Θ^c_{klein} — asymptotisch ≈ 0.07–0.20 (Reserve 5–14). Im Hybrid K* (W̄ an der Ecke) wäre die Ecke unkritisch, aber dort liegt der
A-Bogen in f_K(D) und die eckfreie Schranke aus (a″) gilt dort nicht (Ps (E_A)).

### 3.5a Literatur zu Rechteck-Abschnitten — was davon greift (Auftrag der Leitung, Punkt 2)

- **Hansen, On the approximation of spectra of linear operators on Hilbert spaces, J. Funct. Anal. 254 (2008)** [S] (n-Pseudospektren) und
  Colbrook–Roman–Hansen, Phys. Rev. Lett. 122 (2019) [S]: für beschränktes A wächst σ_min(P_m A P_n) monoton in m gegen σ_min(AP_n) ≥ γ(A)
  (Injektionsmodul); Rechteck-Abschnitte mit genügend vielen Zeilen sind also genau dann stabil, wenn A nach unten beschränkt ist. **Greift**
  (sinngemäss, J_∞ unbeschränkt aber auf P_N wohldefiniert) für d = 2N bzw. d = ∞ ⇒ (a″)_∞ = Os Satz-O7-Aussage. Für kleine d (≈ N/60) liefert
  die Theorie nichts — dort ist es Kantenanalyse.
- **Toeplitz:** T(a), a in der Wiener-Algebra, T(a) injektiv (Windungszahl ≥ 0 in der passenden Konvention): P_{n+l}T(a)P_n mit l ≥ Bandbreite
  ist exakt T(a)P_n, also stabil; im Wiener-Fall gehen die Schwänze ‖(I−P_{n+l})T(a)P_n‖ ≤ Σ_{j>l}|a_j| → 0 gleichmässig in n. Quadratische Abschnitte
  haben dagegen |wind| exponentiell kleine Singulärwerte (Splitting; Roch–Silbermann 1996/98, Böttcher–Grudsky 1999 [O]). Moore–Penrose-Inverse
  von Abschnitten: Heinig–Hellinger, IEOT 19 (1994) 419–446 [S]; Silbermann, SIAM J. Matrix Anal. Appl. 24 (2003) („modified finite sections“,
  stückweise stetige Symbole) [S]. Pseudospektren rechteckiger Matrizen: Wright–Trefethen, IMA J. Numer. Anal. 22 (2002) 501–519 [S];
  stückweise stetige Symbole konvergieren nur algebraisch: Böttcher–Embree–Trefethen, SIAM J. Matrix Anal. Appl. (2002) [S].
  **Greift nicht direkt:** J̃_∞ ist in ℓ²_ρ nicht bandbeherrscht (Spalte j hat Masse bis Zeile ~1.04j und darüber, |J_N| ~ 1.63^N), und die „Windung“
  ist unendlich (Spirale des Mellin-Symbols, O3), daher wächst die benötigte Zeilenzahl d(N) — kein fester Index.
- **Mechanismus-Analogie (Skizze):** J_∞ auf H²(D_ρ) ist (numerisch) injektiv mit Modul min_{Γin}|m|; die quadratische Abschneidung „erzwingt
  Surjektivität“ in endlicher Dimension und erzeugt die Fast-Kerne — wie T_n(z−λ). Die Zahl der dabei sichtbaren Kokern-Richtungen wächst mit N
  (Os Paare bei N ≈ 60, 450, 950; mein d(N) = 2, 2, 3, 5, 8 für N = 120…500).

### 3.6 Der kanonische Eckoperator 𝒞_N (Rang ≤ d) — billige Formel und Gestalt (bewiesen + numerisch)

**Formel (bewiesen).** Mit Z := BJ_N⁻¹ (d×(N+1)) ist die Kleinste-Quadrate-Lösung von [J_N; B]δ ≈ (w_top, w_over) gleich
δ = J_N⁻¹(w_top − e), e = Zᵀ(I_d + ZZᵀ)⁻¹(Zw_top − w_over) (setze δ = J_N⁻¹(w_top − e), minimiere |e|² + |Zw_top − w_over − Ze|²; Push-through).
Also **Π̂g = Π_N g − 𝒞_N w,  𝒞_N := J_N⁻¹Zᵀ(I_d + ZZᵀ)⁻¹[Z | −I_d],  w = T_{N+d}J_∞g**, und c_N(g) = 𝒞_N w = J_N⁻¹Zᵀ(I+ZZᵀ)⁻¹·o(g) mit
o(g) = [J_∞(Π_N g − g)]_{N+1…N+d} (Überlauf-Residuum der QUADRATISCHEN Lösung). Rechnerisch braucht man nur eine LU-Zerlegung von J_N
(x09, x11; N = 1000 in Minuten statt Stunden). (I + ZZᵀ) ≥ I ist gutartig; die exponentielle Verstärkung steckt in J_N⁻¹Zᵀ.

**Gestalt (x11, kubisch, N = 200, d = 3):** σ(𝒞_N) = 1.809, 0.751, 0.062 (σ₁(J_N⁻¹) = 1.852). Mode 1: Eingang (rechter Singulärvektor) zu
90 % in den Zeilen k/N ∈ [0.895, 0.960], Median 0.945, Anteil der Überlaufzeilen ≈ 0 — genau Ws „Steg“ (x ≈ 0.95); Ausgang |u(ρe^{iφ})| =
1.6e-8 (0°), 1e-5 (40°), 0.022 (70°), **4.2 (76.6° = arg L)**, 2.9 (80°), 0.7–0.9 (120°–180°) — Ds/Ws Eckmode u₁. Der kanonische Eckoperator
IST also die Eckmode, aber ohne SVD von J_N⁻¹ und mit der Eingabe-Normierung des Rechteck-Abschnitts.

| N (d) | σ(𝒞_N) (robust: Eigenwerte von (YᵀY)(FFᵀ) in arb/mpmath) | σ₁, σ₂ von J_N⁻¹ (W) | Eingang Mode 1/2/3: Median k/N | Ausgang Mode 1 bei arg L |
|---|---|---|---|---|
| 200 (3) | 1.809, 0.751, 0.062 | 1.852, 1.365 | 0.945 / 0.940 / 0.965 | 4.2 |
| 300 (4) | 3.624, 1.655, 0.257, 0.0133 | 3.698, 1.803 | 0.950 / 0.930 / 0.953 | 6.3 |
| 400 (6) | 9.072, 3.362, 0.811, 0.117, 0.0085, 3e-4 | 9.100, 3.45 | 0.953 / 0.938 / 0.912 | 7.9 |
| 500 (9) | 24.75, 7.968, 1.214, 0.561, 0.103, 0.0102, 6e-4, … | 24.77, 8.00 | 0.952 / 0.938 / — | 9.2 |

(Bei N = 600 war die float64-Rang-d-SVD unbrauchbar (Dynamik von Z), daher die robuste Variante; N = 600 läuft neu.) Die beiden grössten Werte sind
die exponentiell wachsenden Eckmoden (Rate (ln 24.75 − ln 9.07)/100 = 0.0100, wie Ds lokale Rate bei N ≈ 450); der Rest des Ranges d ist klein.
Eingänge auf dem „Steg“ k/N ≈ 0.93–0.95 (W), Ausgänge bei arg L. Dateien x11_N{N}_d{d}.npz (Uc, Vc, sc) für W.

## 3.7 Familien-Stabilität (H₀) → (H) im X1/X3-Rahmen (Auftrag Runde 3)

**Ziel (ARCHITEKTUR §5.1, §9.3 (Fam)):** (H-i)–(H-iii) für die Sekanten-Operatoren 𝒥_P v = v∘P + M_P v, M_P = (f∘P − f∘f)/(P − f), an allen Fast-Lösungen
P ∈ 𝒮_N ∩ B_N, gleichmässig. Im X1-Rahmen heisst (H-ii): (a″)^P + Eingangsschranke auf ∂K' (P-II) + (Ẽ)^P für den Eckteil c^P = Π^P_N − Π̂^P.

**Satz X4a (eckfreier Teil der Familie; bewiesen modulo (a″)_∞ [O], Lemma X5, einer L²-Bernstein-Ungleichung auf glatten Jordankurven [Lit., s. u.]
und (Univ)).** Es gibt θ_N = O(N·R_N) (R_N = K₀C'N^{κ−4.69} Radius von B_N, ARCHITEKTUR §4.2), so dass für alle P mit sup_{D̄}|P − f_K| ≤ R_N gilt:
  σ_min(T_{5N}𝒥_P|_{P_N}) ≥ (1 − θ_N)(c_∞ − τ_N) − θ_Nτ_N,
und |Π̂^P k_ζ − T_Nk_ζ|₁ ≤ C N² e^{−(N+2)g'(ζ)} ausserhalb K', gleichmässig in P.
**Bedingung:** θ_N → 0 verlangt κ < 3.69 (strenger als κ < 4.69 in (H)); für κ = 3.69 genügt eine kleine Konstante K₀C'. Gemessen (x17) liegt die
echte Lösung bei sup_{|x|=ρ}|p_N − f_K| = 1.95e-11 / 6.53e-13 / 9.14e-14 (N = 200/400/600, Steigung −4.90/−4.85), also N·sup ≈ N^{−3.9} — auf dem
kohärenten Pfad ist θ_N winzig; der Bernstein-Faktor ist echt: sup|p_N' − f_K'|/sup|p_N − f_K| = 0.54N / 0.51N / 0.46N.
**Normen:** alle Normen in (1)–(2) sind L²(∂D_ρ)-Normen und T_{5N} die Fourier-Abschneidung (Kontraktion) — so stört es nicht, dass M_P an
Schlitzkreuzungen (P(x) auf der anderen Seite von S_L) nicht holomorph sein muss. **Nur Sekante:** für die Tangente J_N(P) (in (H-i) ebenfalls
verlangt) braucht (1) zusätzlich |P' − f'| klein auf P(D̄) ≈ f(D̄) ⊄ D̄ — das folgt nicht aus P ∈ B_N, sondern erst aus der Fast-Lösungs-
Eigenschaft (P − f)∘P = [P∘P − f∘f] − [f∘P − f∘f] [Skizze; für p_N numerisch gleich gut, s. Tabelle].
*Beweis.* (1) *Relative Störung des unendlichen Operators:* 𝒥_Pδ − J_∞δ = (P − f)∫₀¹δ'(f + s(P−f))ds + (M_P − m)δ. Zweiter Term: M_P − m ist der Differenzenquotient
von f' (f'' beschränkt nahe f(D̄), da σ − 2 = 3.45 > 0; Schlitzkreuzungen geben nur ~R_N^{4.4}) ⇒ |M_P − m| ≤ C R_N ⇒ ‖·‖ ≤ C R_N‖δ‖ ≤ (C R_N/c_∞)‖J_∞δ‖.
Erster Term: E := f(∂D) ist eine glatte Jordankurve (f univalent auf D̄, C^{5,0.44} an L), und für Polynome vom Grad ≤ N gilt die L²-Bernstein-Ungleichung
‖p'‖_{L²(E)} ≤ C_E N‖p‖_{L²(E)} (Beweis: p/Φ^N ist im Äusseren von E holomorph, also ‖p‖_{L²(E_r)} ≤ C‖p‖_{L²(E)} auf den Niveaukurven E_r = {g_E = 1/N};
Cauchy-Formel für p' über E_r mit Abstand ~1/N, Schur-Test). Dasselbe auf den um ≤ R_N ≪ 1/N verschobenen Kurven Γ_s = (f + s(P−f))(∂D). Mit |f'| auf ∂D
nach oben und unten beschränkt (Univ): ‖δ'(f + s(P−f))‖_{L²(∂D)} ≤ C N‖δ∘f‖_{L²(∂D)} ≤ C N(1 + ‖m‖_∞/c_∞)‖J_∞δ‖. Zusammen ‖(𝒥_P − J_∞)δ‖_{L²(∂D_ρ)} ≤ θ_N‖J_∞δ‖_{L²(∂D_ρ)}.
(2) *Übergang zum Rechteck-Abschnitt OHNE Schwanz für P:* ‖T_{5N}𝒥_Pδ‖ ≥ ‖T_{5N}J_∞δ‖ − ‖(𝒥_P − J_∞)δ‖_{L²} ≥ ‖T_{5N}J_∞δ‖ − θ_N(‖T_{5N}J_∞δ‖ + τ_N‖δ‖) (Lemma X5 nur für f!).
(3) *Eingänge:* 𝒥_Pk_η = 1/(η − P) + M_P/(η − x); für η ∈ ∂K' liegt P⁻¹(η) innerhalb R_N/c_f von h(η) ∉ D, also gilt Ps P-II-Abschätzung mit r_N = ρ(1 − 1/N)
unverändert (R_N ≪ 1/N). (4) Rest wie Satz X1 (4)–(5). ∎
*Bemerkung:* Ein Schwanzargument für P selbst wäre falsch — P = (1−λ)f + Q ist ausserhalb W̄ exponentiell gross; der Trick (2) braucht den Schwanz nur für f.
**Bernstein-Literatur:** L^p-Bernstein-/Markov-Ungleichungen auf glatten Jordankurven (Totik; Nagy–Totik) [U, genaue Stelle nachschlagen]; der Beweis oben ist
selbsttragend bis auf die Standardschätzung ‖p‖_{L²(E_r)} ≤ C‖p‖_{L²(E)} (Smirnov-Raum des Äusseren, Kellogg–Warschawski für glatte E).

**Satz X4b (Reduktion der Familien-Annahme auf die Eckteile; Beweis = X1 für 𝒥_P).** (F_{K'}) + (Univ) + (a″)_∞ + Lemma X5 + [für alle P ∈ 𝒮_N ∩ B_N:
S^P_N invertierbar (⇔ (H-i), da det S = det R⁺_top·det J_N) und Θ^{c,P} ≤ 1 − ε, Rc^P(L) ≤ CN^κ] ⇒ (H-ii) gleichmässig ⇒ mit (H-iii) Satz H1/H2.
Damit ist die Familien-Frage auf die d-dimensionalen Eckteile der Sekanten-Operatoren reduziert; der ganze eckfreie Rest ist gleichmässig bewiesen.
*Exakte Störungsformel für die Eckteile (bewiesen):* Π^P_N g − Π_N g = −(T_N𝒥_P|_{P_N})⁻¹T_N[(𝒥_P − J_∞)(Π_N g − g)] (aus T_NJ_∞(Π_N g − g) = 0); mit Bs B-Id
fallen darin alle Terme mit J_∞(Π_Ng − g) = x^{N+1}ψ bis auf Kopfterme (N+1)[J_∞E]_{N+1}…·x^N weg — die Familienstörung wird also von Pol-Fehlern auf D̄
(eckfreie Eingaben) und den ersten Überlauf-Residuen o_i(g) = [J_∞(Π_Ng − g)]_{N+i} getragen — genau den Grössen, die in Kor. 3.2/§3.6 den Eckteil steuern.

**Numerik (x13, x16; Tangente J_∞(p_N) an der echten Lösung p_N, gemeinsam/pN_N100_bis_1000.json, gegen J_∞(f_K); Ps ∂K'-Gitter):**

| N, d | σ_min (f_K / p_N) | max rel. Abweichung R / R̂ / Rc über alle 112 Randpunkte | Θ^c (f_K / p_N) |
|---|---|---|---|
| 100, 3 | 0.73191 / 0.73191 | identisch auf 5 Stellen | 0.0715 / 0.0715 |
| 200, 3 | 0.73191 / 0.73191 | 1.0e-11 / 1.1e-11 / 2.8e-11 (sup|p_N − f_K| ≈ 1.9e-11) | 0.1288 / 0.1288 |
| 400, 6 | 0.73191 / 0.73191 | 5.2e-13 / 4.5e-13 / 1.5e-12 (sup|p_N − f_K| ≈ 6.5e-13) | 0.1141 / 0.1141 |
| 600, 12 | — / 0.73190 | 7.0e-14 / 5.2e-14 / 2.4e-13 (sup|p_N − f_K| ≈ 9.1e-14) | 0.1062 / 0.1062 |
| 200, 800 | 0.73191 / 0.73191 | R 1e-11; R̂, Rc: Median 0.9 %, max 87 % (ferne f_K(Γout)-Punkte, Werte dort ~1) | −0.039 / −0.044 |

⇒ Mit kleinem d ändern sich ALLE Grössen (auch der Eckteil) relativ nur um O(sup|P − f|) — die Familie ist numerisch so stabil wie überhaupt messbar.
Mit d ≥ 2N sehen die Zusatzzeilen P ausserhalb D̄ (wo p_N ≠ f_K), daher ~1 % Abweichung in R̂/Rc; die Unterschranke σ_min bleibt identisch (wie Satz X4a sagt).

**Eckmoden unter Störung (x18, x18b; 11.10. 02:05–02:30).** Gemessen werden die grössten Singulärwerte σ₁ ≥ σ₂ ≥ … von A_P⁻¹ (A_P = T_N𝒥_P|_{P_N}
bzw. Tangente J_N(P)); σ₁ (und ab N ≈ 400 σ₂) sind die Eckmoden, der Rest Bulk (1.3646…). Störungen gleicher Grösse R := sup_{|x|=ρ}|p_N − T_Nf_K|:
pS/pN = echte Lösung p_N (kohärent), st/sh/sv = f + η, η zufällig reell mit Koeffizienten nur in k ≤ N/2 / k ≥ 0.9N / alle k. Eintrag = (δσ₁/σ₁)/R; bei Zufall je Ziehung (Seeds 1/2/3).

| N | σ₁(J_N⁻¹) | Sekante pS (kohärent) | Sekante st (tief) | Sekante sv (alle) | Sekante sh (hoch) | Tangente pN | Tangente zh / zv / ft |
|---|---|---|---|---|---|---|---|
| 200 | 1.8410 | −0.13 | +2.6 | +1.6 | 3e-7 | −0.26 | σ₁ → 102 / 107 / 13.1 (Faktor 55 / 58 / 7) |
| 400 | 8.9154 | −0.34 (σ₂: −0.30) | +9.2 / +3.2 / +7.2 | +6.2 / +1.8 / +4.8 | 3e-6 | −0.68 (σ₂: −0.61) | — |
| 600 | 68.410 | −0.34 (σ₂: −0.32) | +16.1 / +1.1 / +10.6 | +0.8 | 2e-6 | −0.68 (σ₂: −0.65) | — |

*Befunde.* (1) **Kohärente Familie (p_N) — keine Eckverstärkung:** δσ₁/σ₁ = −0.13R / −0.34R / −0.34R (Sekante) bzw. −0.26R / −0.68R / −0.68R
(Tangente) für N = 200/400/600, während σ₁ um ×4.8 und ×7.7 wächst. Ab N = 400 ist der Quotient **konstant** — die Eckmoden der Sekanten- und
Tangenten-Operatoren an der echten Lösung sind relativ O(R) stabil, ohne Faktor σ₁ und ohne Bernstein-Faktor N.
(2) **Zufällige Sekanten-Störungen gleicher Grösse:** tieffrequent +2.6R / +9.2R / +16.1R — wächst polynomial (≈ N^{1.5}), NICHT wie σ₁
(drei Zufallsziehungen: N = 400 st 9.2/3.2/7.2, sv 6.2/1.8/4.8; N = 600 st 16.1/1.1/10.6, sv 0.8 — Mittel st 2.6 → 6.5 → 9.3,
also ×1.4 von 400 auf 600 bei ×7.7 in σ₁); hochfrequent (k ≥ 0.9N) praktisch 0 (1e-6·R). Der naive Störungssatz (|δσ₁/σ₁| ≤ σ₁·θ_N mit θ_N ~ N R) überschätzt um
den Faktor σ₁·N/16 ≈ 2600 (N = 600). Strukturell heisst das: das Eck-Matrixelement ⟨u₁, T_N(𝒥_P − J_∞)v₁⟩ ist selbst ∝ s₁ = 1/σ₁ — die Störung
erhält die Fast-Singularität der Eckmode (Vermutung (Fam-Eck): |⟨u_i, T_N(𝒥_P − J_∞)v_i⟩| ≤ C N^a s_i·sup_{D̄}|P − f|; numerisch a ≈ 1.5 für
Zufall, a ≈ 0 kohärent).
(3) **Tangenten an beliebigen P der Kugel sind KEINE kleinen Störungen:** in P'∘P wird η' auf f(D̄) ausgewertet, und f(D̄) reicht bis f(ρ) = 2.2327
(ausserhalb D̄ und sogar ausserhalb W̄: |f(ρ) − a₀| = 1.734 > R₀ = 1.349). x18b, N = 200: ‖T_N[η'∘f]‖_{ℓ²_ρ} = 1.2e10 (k ≤ N/2) bzw. 2.3e31 (k ≥ 0.9N)
bei ‖η‖ = 1e-11 — Verstärkung (f(ρ)/ρ)^N = 1.624^N (= 1.3e42 bei N = 200, gemessen 1.9e42). Schon P = T_Nf (nur Abschneiden) verschiebt σ₁
von 1.84 auf 13.1, Zufalls-η auf 102. Für p_N dagegen ‖T_N[p_N'∘p_N] − m‖ = 3.8e-5 und δσ₁/σ₁ = −0.26R. ⇒ (H-i) für die Tangente J_N(P) ist
**keine Kugel-Eigenschaft**; sie braucht P ≈ f in C¹ auf P(D̄) ≈ f(D̄), was nur Fast-Lösungen haben (As „Zufallsrichtungen ≳ 1.17^N“ in anderer
Gestalt). Die Sekante 𝒥_P (nur f wird bei P ausgewertet) hat dieses Problem nicht.
(4) **Folgerung für den Beweis.** Der eckfreie Teil ist gleichmässig bewiesen (X4a). Für den Eckteil reicht eine Störungsrechnung um f genau dann,
wenn (Fam-Eck) gilt — dann ist s_i(A_P)/s_i(J_N) = 1 + O(N^a R_N) und mit S^P, Θ^{c,P} ebenso (Kor. 3.2 ist eine rationale Funktion der Eckdaten);
Bedingung κ < 4.69 − a. Ohne (Fam-Eck) bliebe nur die naive Schranke σ₁N R_N (versagt ab N ≈ 3000–3500, ARCHITEKTUR §6.5) oder ein
struktureller Beweis von (H₀) für jedes 𝒥_P aus P-stabilen Grössen (Fixpunkt L_P, μ_P, (Univ), Randkurve P(∂D)).

## 3.8 Runde 4 (11.10. 02:50–03:30): Bernstein-Lemma, Eckempfindlichkeit der Familie, Tangente in (H-i)

### 3.8.1 Lemma X6 (L²-Bernstein nahe f(∂D)) — bewiesen bis auf Kellogg–Warschawski

**Aussage.** Sei E := f(∂D) (Jordankurve wegen (Univ); C^{5,0.44}, denn f(L + ζ) = L + μζ + … + Cζ^{5.447}…; f' ≠ 0 auf ∂D, P: min|f'| = 0.3456).
Es gibt C_E und c > 0 mit: für jedes Polynom δ vom Grad ≤ N und jede rektifizierbare Kurve Γ im Abstand ≤ c/(2N) von E (insbesondere
Γ_s = (f + s(P − f))(∂D) für sup_{D̄}|P − f| ≤ c/(2N)) gilt ‖δ'‖_{L²(Γ, |dz|)} ≤ C_E N‖δ‖_{L²(E, |dz|)}; ebenso ‖δ‖_{L²(Γ)} ≤ C_E‖δ‖_{L²(E)}.
**Beweis.** (i) *Aussen:* Φ: ext E → {|w| > 1} konform, Φ(∞) = ∞; F := δ/Φ^N ist holomorph in ext E ∪ {∞}. Für Ψ = Φ⁻¹ ist F∘Ψ holomorph in
|w| > 1 inkl. ∞, also M₂(r) = (1/2π)∫|F∘Ψ(re^{iθ})|²dθ = Σ_k|b_k|²r^{−2k} fallend in r. Nach Kellogg–Warschawski (E ∈ C^{1,α}) ist Ψ' stetig
und ≠ 0 auf |w| ≥ 1, also |dz| ≍ dθ auf allen Niveaukurven E_r = Ψ(|w| = r), 1 ≤ r ≤ 2. Damit ‖δ‖_{L²(E_r)} ≤ C r^N‖δ‖_{L²(E)}, für r ≤ 1 + 2c'/N
also ≤ C e^{2c'}‖δ‖_{L²(E)}. (ii) *Innen:* φ: D(0,1) → int E konform; δ∘φ ∈ H^∞, M₂(r) wachsend in r, φ' stetig ≠ 0 auf dem Rand (K–W) ⇒
‖δ‖_{L²(E^{in}_r)} ≤ C‖δ‖_{L²(E)} auf den inneren Niveaukurven. (iii) *Cauchy, gemittelt:* für analytisches δ ist |δ'(z)|² ≤ C t^{−4}∬_{D(z,t)}|δ|²dA
(Cauchy-Formel auf Kreisen der Radien t/2…t mitteln, Cauchy–Schwarz). Mit t = c/N und z ∈ Γ: ∫_Γ|δ'|²|dz| ≤ C t^{−4}·t·∬_U|δ|²dA, U = 2t-Umgebung
von E (jeder Punkt w liegt in D(z,t) nur für z auf einem Γ-Bogen der Länge ≤ Ct). U ist durch innere und äussere Niveaukurven geblättert
(Ko-Fläche, |∇| ≍ 1 nach K–W), also ∬_U|δ|² ≤ C t·sup_{Niveau}‖δ‖²_{L²} ≤ C t‖δ‖²_{L²(E)} nach (i), (ii). Zusammen ≤ C t^{−2}‖δ‖² = C N²‖δ‖². ∎
**Verwendung in X4a:** ‖δ'(f + s(P−f))‖_{L²(∂D)} ≤ C‖δ'‖_{L²(Γ_s)} (|(f + s(P−f))'| ≥ c auf ∂D) ≤ C N‖δ‖_{L²(E)} ≤ C N‖δ∘f‖_{L²(∂D)}. Bedingung
R_N ≤ c/(2N): gilt, da R_N = K₀C'N^{κ−4.69} und κ < 3.69.
**Literatur (Prüfstand):** Kellogg–Warschawski in der C^{m,α}-Form = Ch. Pommerenke, *Boundary Behaviour of Conformal Maps*, Springer 1992, Thm 3.6
[S: so zitiert in arXiv 1803.09458; Buchtext nicht selbst gelesen]. Eine fertige **L²**-Bernstein-Ungleichung auf glatten Jordankurven habe ich
nicht gefunden; die scharfe **Supremums**-Form auf C²-Kurven (Konstante über Normalableitungen der Green-Funktion) ist Kalmykov–Nagy–Totik,
*Bernstein- and Markov-type inequalities for rational functions*, Acta Math. 219 (2017), arXiv 1610.06706 [S]; Übersicht über L^p-Fassungen:
Queffélec–Zarouf, arXiv 1903.10801 [S]. Der obige Beweis ist selbsttragend (nur K–W und Monotonie der H²-Mittel) — eine eigene Literaturstelle
für die L²-Bernstein-Ungleichung ist damit **nicht mehr nötig**, nur K–W.

### 3.8.2 Eckempfindlichkeit: die Ecke reagiert POLYNOMIAL auf jede glatte Störung, auf kohärente fast gar nicht (x18 Formen, x19, x19b)

Störungen P = T_Nf + η der Sekante 𝒥_P (erste Ordnung: M_P = m + ½(f''∘f)η), alle auf sup_{|x|=ρ}|η| = R normiert (R = sup|p_N − T_Nf|).
Eintrag (δσ₁/σ₁)/R; in Klammern die Grösse des Matrixelements |⟨u₁, δA v₁⟩|/R = Eintrag/σ₁.

| η | lokale Daten bei L (x19) | N = 200 (σ₁ = 1.84) | 400 (8.92) | 600 (68.4) | Wachstum 400→600 |
|---|---|---|---|---|---|
| e0 = konst. | Fixpunkt radial +1.36R | +32 (17) | +131 (15) | +230 (3.4) | ×1.76 ≈ N^{1.4} |
| e1 = x | Fixpunkt tangential 1.36R | +4.3 | +15.9 | +27.5 | ×1.73 |
| eL = q, q(L) = 0 | nur μ ändert sich (0.57iR) | +2.0 | +7.3 | +12.1 | ×1.66 |
| eL2 = q² | Fixpunkt und μ fest | −0.25 | −0.76 | — | |
| eQ = q·x^{N/2} | mittelfrequent | 0.001 | −0.01 | — | |
| Zufall tief (k ≤ N/2) | | +2.6 | +3…+9 | +1…+16 | |
| Zufall hoch (k ≥ 0.9N) | | 3e-7 | 3e-6 | 2e-6 | |
| **p_N (kohärent)** | Fixpunkt radial −1.4R, Betrag δμ ≈ 185R | **−0.13** | **−0.34** | **−0.34** | **×1.0** |
| Bulk σ₃, σ₄ (alle Formen) | | ±0.25…0.30 | ±0.29 | ±0.30 | |

**Befunde.** (1) **Keine Form wächst wie σ₁** (×7.7 von 400 auf 600): die Antwort der Ecke auf jede glatte Störung ist polynomial, ≈ N^{1.4–2}·R,
am grössten für eine radiale Verschiebung des Fixpunkts (e0). Das Matrixelement ⟨u₁, δAv₁⟩ ist also selbst ∝ s₁ = 1/σ₁ (e0: 17 → 15 → 3.4 R).
(2) Die echte Lösung verschiebt den Fixpunkt des Polynoms p_N genauso stark (−1.4R radial, x19) und ändert P'(L) um 185R — und trotzdem antwortet
die Ecke nur mit −0.34R, d. h. **nur mit dem gleichförmigen Bulk-Anteil** (≈ −0.3R bei allen Formen, auch bei σ₃, σ₄). Grund (x19b): η_p = p_N − T_Nf
ist ein **Buckel der Breite ~1/N bei L und L̄** (Anteil mit |η| ≥ R/2: 0.024 → 0.012, L¹-Masse/R: 0.036 → 0.017 für N = 200 → 400) und **oszilliert**:
sein Mittel über ±3° um L ist 0.022R bzw. 0.004R, über ±10° 0.001R bzw. < 0.001R. Die Ecke sieht die Störung nur geglättet über Bögen fester
Länge — dort ist die kohärente Störung praktisch null; ihr Fixpunkt-Versatz ist ein 1/N-Effekt, den die Ecke nicht auflöst.
(3) Damit ist die Lage klarer als in §3.7: dort hatte ich aus zwei N-Werten auf σ₁-Wachstum für Zufallsstörungen geschlossen — **durch N = 600 widerlegt**.

**Vermutung (Fam-Eck), verschärft:** Für P mit sup_{D̄}|P − f| ≤ R gilt
  δσ_i/σ_i = β_N(P) + N·⟨∇c, Δ_ε(P)⟩ + (höhere Ordnung), |β_N| ≤ C R (gleichförmiger Bulk-Anteil),
wobei Δ_ε(P) die über Bögen fester Länge ε um L, L̄ geglätteten Eckdaten von P sind (radiale Lage des Fixpunkts, Multiplikator, h auf Γout) und c die
Ecken-Rate aus Ws Leitordnungsproblem (w §20: c hängt nur von h auf Γout, L, μ ab). Folgerungen: allgemein |δσ_i/σ_i| ≤ C N R (gemessen lokal N^{1.4};
e0: 230R/600 ≈ 0.38·NR); für kohärente P (Δ_ε(P) ≪ R/N) |δσ_i/σ_i| ≤ C R.
**Beweisweg:** (Fam-Eck) ist keine eigene Störungsrechnung, sondern die **Glattheit der Ecken-Asymptotik in den Eckdaten**: Wer (H₀) über Ws Modell
beweist (s₁(J_N) = e^{−cN}N^{κ}(1 + o(1)) mit c = c(L, μ, h auf Γout)), beweist mit denselben Schritten s₁(A_P) für jedes 𝒥_P mit c(Δ_ε(P)) — sofern
der Beweis nur Eckdaten auf Skala ε benutzt. Für (H) reicht schon die schwache Form |δσ_i/σ_i| ≤ C N^a R_N mit a < 4.69 − κ.
**Konsequenz für X4b:** mit (Fam-Eck) für a = 2 folgt (H) aus (H₀) für alle P ∈ B_N (nicht nur kohärente), solange κ < 2.69; der eckfreie Teil (X4a)
braucht κ < 3.69.

### 3.8.3 (H-i): Tangente nur an Fast-Lösungen — Vorschlag für ARCHITEKTUR §9 (an P über das BRETT)

Vorschlagstext: (H-i) verlangt die Invertierbarkeit der **Tangente** J_N(P) = T_N[δ∘P + (P'∘P)δ] und der **Sekante** T_N𝒥_P. Für die Sekante ist (H-i)
eine Kugel-Eigenschaft (eckfrei X4a; Ecke polynomial stabil, §3.8.2). Für die Tangente NICHT: P'∘P wertet P' auf P(D̄) ≈ f(D̄) aus, und f(D̄) reicht bis
f(ρ) = 2.2327 (ausserhalb W̄: |f(ρ) − a₀| = 1.734 > R₀ = 1.349); eine Störung η mit sup_{D̄}|η| = R ändert T_N[P'∘P] um bis zu R·(f(ρ)/ρ)^N = R·1.624^N
(x18b: 2.3e31 bei N = 200, R = 1.8e-11); σ₁(J_N(P)⁻¹) springt dann von 1.84 auf ~100 (x18 zh/zv), schon P = T_Nf allein auf 13.1. An Fast-Lösungen
ist die Tangente harmlos (p_N: ‖T_N[p_N'∘p_N] − m‖ = 3.8e-5, δσ₁/σ₁ = −0.68R bei N = 400 und 600). Satz H1 braucht die Tangente nur in Schritt (2)
(Satz über implizite Funktionen) und nur **an den Lösungen (λ, Q(λ)) ∈ 𝒵** — dort ist P = P_λ eine Fast-Lösung; (H-i) für J_N(P) ist also nur auf 𝒵 zu
verlangen, nicht auf ganz 𝒮_N ∩ B_N. Dass auf 𝒵 P ≈ f in C¹ auf P(D̄) gilt, folgt aus (P − f)∘P = [Φ(P) − Φ(f)] − [f∘P − f∘f] [Skizze].

### 3.8.4 Kanten-Aussage für (a″) mit kleinem d (Auftrag 1; x20, Kneser; Abstimmung mit W über das BRETT)

**Messung (x20, genau: (R_dᵀR_d)⁻¹ = J⁻¹J⁻ᵀ − J⁻¹Z_dᵀ(I + Z_dZ_dᵀ)⁻¹Z_dJ⁻ᵀ in arb, dann float-eigvalsh).** σ_min(R_d), R_d = T_{N+d}J_∞|_{P_N}:

| N | r (Eckmoden) | σ₁, σ₂(, σ₃) von J_N⁻¹ | d = 0 | 1 | 2 | 3 | 4 | 6 | 8 | 12 | d₀ (σ_min = 0.7319) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 200 | 1 | 1.84 | 0.543 | 0.719 | 0.7319 | 0.7319 | | | | | 2 |
| 400 | 2 | 8.92, 3.39 | 0.112 | 0.236 | 0.387 | 0.543 | 0.676 | 0.7319 | | | 5–6 |
| 600 | 2 | 68.4, 19.8 | 0.0146 | 0.038 | 0.077 | 0.136 | 0.215 | 0.424 | 0.645 | 0.7319 | 9–12 |
| 800 | 3 | 592, 147, 1.96 | 0.0017 | 0.005 | 0.011 | 0.023 | 0.041 | 0.106 | 0.223 | 0.573 | > 12 (d = 16: Präzision reicht nicht, Wert unbrauchbar) |

⇒ **d₀(N) ≈ 0.02N − 2** (2 / 6 / 10 / ≳ 14), also **proportional zu N**, nicht zur Zahl r(N) der Eckmoden (r = 2 bei N = 400 und 600, d₀ 6 → 10).
Unterhalb von d₀ steigt σ_min(R_d) gleichmässig (Faktor 2–3 je Zeile), oberhalb ist es exakt die Bulk-Konstante 0.7319 (= x13, gleich für d = N, 4N).
Die Läufe von §3.5/§3.6 mit „kleinem d“ (3/4/6/12 für N = 200/300/400/600) lagen also genau auf d ≈ d₀(N).

**Mechanismus.** J_∞ ist auf H²(D_ρ) unbeschränkt (f(D_ρ) ⊄ D_ρ); die Überlaufzeilen sind auf JEDEM Vektor riesig (Eckmode 1, Zeile N+1:
1.5e6 / 1.5e12 / 8.6e17 für N = 200/400/600; Bulk-Mode bei N = 600: 3.8e51). „Einfangen“ heisst also nicht „genug Überlauf“, sondern: die
Überlaufprofile j ↦ [J_∞v_i]_{N+j} der Eckmoden müssen in den ersten d Zeilen von einander und von den Bulk-Überläufen **linear unabhängig** sein.
Sie sind glatt in j/N: das Verhältnis Mode 1/Mode 2 ist linear in j mit relativer Steigung 0.035 (N = 400) bzw. 0.024 (N = 600) je Zeile, d. h. die Profile
trennen sich auf der Indexskala j ≈ 0.07N — genau der Abstand des Stegs (Eingänge der Eckmoden bei k/N ≈ 0.93) von der Kante. Daher d₀ ∝ N.

**Kanten-Aussage (K_d), genau formuliert.** Es gibt ε₀ > 0 und c₀ > 0, so dass für d = ⌈ε₀N⌉ (gemessen genügt ε₀ = 0.02 bis N = 600; bei N = 800 ist d₀ > 12, die Gerade 0.02N − 2 sagt 14 voraus — d = 16 dort wegen Präzision nicht auswertbar):
  σ_min(T_{N+d}J_∞|_{P_N}) ≥ c₀.
Äquivalent (Lemma 3.1): ‖J_N⁻¹(I + Z_dᵀZ_d)^{−1/2}‖ ≤ 1/c₀, Z_d = T_{(N,N+d]}J_∞J_N⁻¹. Mit Satz X1 (Variante „kleines d“) ergibt das zusammen mit
Θ^c_{d} ≤ 1 − ε (gemessen 0.129/0.114/0.106 für N = 200/400/600 mit d = 3/6/12, Grenzwert 0.07–0.20) die lokale Pol-Schranke (LP).
**Liefert Ws Eckmodell (K_d)?** Nicht direkt. Ws Modell (w §19–§20) gibt die Eckmoden in der Ebene auf exponentieller Ordnung (Rate H, c = 0.0118) und
den Kegel-Exponenten κ₁; in Koeffizienten übersetzt (Sattelpunkt) bestimmt die Rate das Überlaufprofil nur bis auf e^{o(N)}. Alle Eckmoden haben
dieselbe führende Rate; sie unterscheiden sich in den Kegelharmonischen Re(c_n w^n), n = 1, 2, … (w §19), also auf Vorfaktor-Ebene. (K_d) braucht
genau diese Ebene: die Überlaufprofile der Eckmoden auf j ≤ εN bis auf relative Fehler o(1) — oder wenigstens die gleichmässige lineare Unabhängigkeit
der n-ten Kegelharmonischen auf der Indexskala j/N. **Anfrage an W** (BRETT): Gibt das Kegelmodell die Fortsetzung der Eckmoden über die Kante k > N
(Überlaufprofil als Funktion von j/N) mit Vorfaktor? Dann ist (K_d) eine endliche Determinanten-Bedingung (Gram-Matrix der r Profile auf [0, ε₀]).
Gegenprobe ohne Modell: (K_d) folgt auch aus (a″)_∞ [O] + Lemma X5, wenn man d = 4N nimmt — das ist die grenzwertige Variante (Θ^c → 0.82).

### 3.8.5 (Fam-Eck) quantitativ: die Eckantwort ist eine lineare Funktion der lokalen Fixpunktdaten (x18 se2–se4, x19, N = 400)

Glatte Störungen η (sup = R) von P = T_Nf + η, Sekante; lokale Daten des Fixpunkts L_P von P: radialer Versatz δr = |L_P| − ρ, tangentialer δt,
Multiplikator δμ = P'(L_P) − μ (alle in Einheiten R, relativ zu T_Nf):

| η | δr | δt | δμ | (δσ₁/σ₁)/R gemessen | Modell 90.7·δr + 16.9·Re δμ + 14.1·Im δμ | ohne diesen Punkt angepasst |
|---|---|---|---|---|---|---|
| 1 | +1.357 | +0.169 | −0.153 + 0.766i | +131 | +131.4 | +133.1 |
| x | +0.150 | +1.359 | −0.053 + 0.028i | +15.9 | +13.1 | +13.1 |
| q (Null in L) | 0 | 0 | 0.575i | +7.3 | +8.1 | +8.2 |
| x² | −1.288 | +0.460 | +0.128 + 0.663i | −106 | −105.3 | −105.0 |
| x³ | −0.746 | −1.146 | −1.265 + 0.606i | −80 | −80.4 | −83.1 |
| x⁴ | +0.943 | −0.991 | −1.351 − 1.646i | +39 | +39.5 | +42.5 |

⇒ **Für glatte Störungen ist die Eckantwort auf 2 % (Kreuzvalidierung ≤ 8 %) eine lineare Funktion von (δr, δμ)**; der tangentiale Versatz
zählt nicht (Koeffizient 1.6 bei freier Anpassung). Der Radialkoeffizient wächst linear in N: a(N) ≈ 22 / 90 / 158 für N = 200/400/600
(aus e0, nach Abzug des μ-Anteils), d. h. a(N) ≈ 0.34·N − 46 — das ist N·∂c/∂r mit **∂c/∂r ≈ 0.34**: Rückt der Fixpunkt um δr vom
Abschneidekreis nach aussen, steigt die Eckrate c um 0.34·δr. Der μ-Teil ≈ N·Re((0.042 − 0.035i)·δμ) (N = 400).
**p_N passt NICHT in dieses lokale Modell** (Punktdaten δr = −1.36, δμ = −21 + 183i würden +2200 vorhersagen, gemessen −0.34): seine Punktdaten
bei L stammen vom 1/N-Buckel; die Ecke sieht nur geglättete Daten (eQ = q·x^{N/2} mit denselben Punktdaten wie q antwortet mit 0.001 statt 7.3).
**Präzisierte Vermutung (Fam-Eck):** δσ₁/σ₁ = β(P) + N·[∂_r c·δr_ε(P) + Re(∂_μ c·δμ_ε(P))] + …, mit geglätteten Fixpunktdaten (Skala ε fest),
∂_r c ≈ 0.34, ∂_μ c ≈ 0.042 − 0.035i, |β| ≤ 0.3R. **Prüfstein für W:** Ws Kegelmodell muss ∂_r c ≈ 0.34 und ∂_μ c (über κ₁(μ) = 2πi/(2πi − log μ))
wiedergeben; dann ist (Fam-Eck) die Ableitung von c(L, μ, h) nach den Daten — genau die Form, in der die Leitung sie führen will.

## 3.9 Meilenstein 6 (11.10. 03:55–04:10): Lemma X7 und X8 mit expliziten Konstanten, Kontur d = 0.3, Zertifikats-Reserve, Bogenmodell

### 3.9.1 Lemma X7 (Eingangsnorm des eckfreien Teils) — bewiesen bis auf zwei Gitter-Minima [num.]

**Aussage.** Für η ∈ ∂K' (= Γin ∪ f_K(Γout), beide Hälften) und M ≥ 1 gilt mit w_η := T_M[J_∞k_η] = T_M[1/(η − f) + m·k_η]:
  |w_η|_{ℓ²_ρ} ≤ C_X7·√M,  C_X7 = e^{1/2 + 1/(4M)}(1 − 1/(4M))^{−1/2}·(1/c_K + max_{∂D}|m|)/ρ,  c_K = min(c_in, c_out),
mit c_in := min_{η∈Γin, x∈∂D}|η − f(x)|/|η − x| = **0.6741** und c_out := min_{ξ∈Γout, x∈∂D}|f(ξ) − f(x)|/|ξ − x| = **0.3456** (= min|f'|), max|m| = 2.3376
(x26, Gitter 6000 Punkte) ⇒ **C_X7 = 6.30 (M = 100) … 6.28 (M → ∞)**; gemessen 1.17 (x22).
**Beweis.** (i) Für g holomorph auf D̄_r, r < ρ: |g_k|ρ^k = (ρ/r)^k|g_k|r^k, also ‖T_Mg‖_{ℓ²_ρ} ≤ (ρ/r)^M‖g‖_{L²(∂D_r, dφ/2π)}.
(ii) r := ρ(1 − 1/(2M)): (ρ/r)^M ≤ e^{1/2 + 1/(4M)} (−log(1−y) ≤ y + y² für y ≤ ½), und ρ² − r² = (ρ²/M)(1 − 1/(4M)).
(iii) Für |ξ| ≥ ρ ist (1/2π)∫|ξ − re^{iφ}|^{−2}dφ = 1/(|ξ|² − r²) ≤ M/(ρ²(1 − 1/(4M))).
(iv) *Untere Schranken:* Für η ∈ Γin ist q(x) := (η − f(x))/(η − x) holomorph in D, und 1/q ist stetig auf D̄ (1/q(η) = 0), denn η ∉ f_K(D̄) für
η ∈ Γin ∖ {L, L̄}, bei η = L ist q(x) → f'(L) = μ ≠ 0. Maximumprinzip für 1/q ⇒ |η − f(x)| ≥ c_in|η − x| auf D̄. Für η = f(ξ), ξ ∈ Γout ⊂ ∂D, ist
q_ξ(x) := (f(ξ) − f(x))/(ξ − x) holomorph auf D, stetig und nullstellenfrei auf D̄ (Univalenz; q_ξ(ξ) = f'(ξ) ≠ 0) ⇒ |η − f(x)| ≥ c_out|ξ − x|.
In beiden Fällen ist der Vergleichspol (η bzw. ξ) vom Betrag ≥ ρ.
(v) Also ‖1/(η − f)‖_{L²(∂D_r)} ≤ (1/c_K)·√M/(ρ√(1 − 1/(4M))) und ‖m k_η‖_{L²(∂D_r)} ≤ max_{D̄}|m|·√M/(ρ√(1 − 1/(4M))) (|η| ≥ ρ auf ∂K';
max_{D̄}|m| = max_{∂D}|m| nach Maximumprinzip, m = eˣ/f' holomorph auf D, stetig auf D̄). Mit (i), (ii) folgt die Aussage. ∎
**Folge (Satz P-II scharf, X1):** ŷ_η = Π̂k_η − T_Nk_η = R⁺w_η − T_Nk_η, |R⁺| ≤ 1/c₀, |T_Nk_η|_{ℓ²_ρ} ≤ √(N+1)/ρ ⇒
|ŷ_η|_{ℓ¹_ρ} ≤ √(N+1)·[C_X7√(1+α)/c₀ + 1/ρ]·√(N+1) =: C_R(N+1), **C_R = 12.7** (c₀ = 0.54, α = d/N = 0.05) statt bisher C·N² (P-II).
Zu zertifizieren [P]: die beiden Minima c_in, c_out (Intervallrechnung auf ∂D × Bogen, braucht zertifizierte f_K-Werte auf ∂D) und max|m|.

### 3.9.2 Lemma X8 (Familien-Konstante A_θ) — elementar, OHNE Bernstein auf f(∂D) und ohne Kellogg–Warschawski

**Aussage.** Für δ ∈ P_N und P mit sup_{D̄}|P − f| ≤ R, N·R ≤ 1/(4A'):
  ‖(𝒥_P − J_∞)δ‖_{L²(∂D)} ≤ [A_θ·N + ½M2/c₀]·R·(1 + 2A'NR)·‖J_∞δ‖,
  A_θ := (1/min_{∂D}|f'|)·[K/ρ + (max|m|/ρ + max|m'|/N)/c₀] + τ'_N,  K = 5 (Zeilen bis 5N, Lemma X5), τ'_N ≤ 1e-8 (N ≥ 100).
Zahlen (x26): min|f'| = 0.3456, max|m| = 2.3376, max|m'| = 1.3839, M2 = 0.8183 ⇒ **A_θ = 19.7 (c₀ = 0.54), 17.3 (c₀ = 0.7319), 22.2 (c₀ = 0.4249)**;
mit Zeilen bis 3N (braucht X5 in der Fassung R = 2) 15.5 / 13.1 / 18.0. Gemessene Bernstein-Konstante zum Vergleich: K_B·(1 + max|m|/c₀) = 4.05.
**Beweis.** (𝒥_P − J_∞)δ = (P − f)∫₀¹δ'(f + s(P − f))ds + (M_P − m)δ. *Hauptterm:* F := δ∘f = J_∞δ − mδ ist auf D holomorph mit
F' = f'·δ'∘f, also ‖δ'∘f‖ ≤ ‖F'‖/min|f'|. Koeffizientenweise ist ‖(J_∞δ)'‖_{L²(∂D_ρ)}² = Σk²|c_k|²ρ^{2k−2} ≤ (5N/ρ)²‖J_∞δ‖² + Σ_{k>5N}k²|c_k|²ρ^{2k−2},
und der Rest ist nach den Koeffizientenschranken von Lemma X5 (Kreisteil ~ β^N, Schlitzteil ~ k^{−σ−1}N^{3/2}) ≤ (τ'_N N)²‖δ‖². Ferner
‖(mδ)'‖ ≤ max|m'|‖δ‖ + max|m|‖δ'‖ ≤ (max|m'| + max|m|N/ρ)‖δ‖ (Bernstein auf dem Kreis in ℓ²_ρ: ‖δ'‖ ≤ (N/ρ)‖δ‖) und ‖δ‖ ≤ ‖J_∞δ‖/c₀.
*Verschiebung s(P − f):* δ'(f + t) = Σ_{j<N}δ^{(j+1)}(f)t^j/j! (δ Polynom, endliche Summe); δ^{(j)}∘f rekursiv über Faà di Bruno aus F^{(j)}, und
‖F^{(j)}‖ ≤ (5N/ρ)^j‖J_∞δ‖ + … wie oben ⇒ ‖δ^{(j)}∘f‖ ≤ (A'N)^j‖J_∞δ‖ mit A' ≈ A_θ, also Faktor e^{A'NR} ≤ 1 + 2A'NR. *Multiplikator:* M_P − m =
½(f''∘f)(P − f) + O(R²) (Taylor; Schlitzkreuzungen nahe L tragen ≤ C_J R^{4.447}), ‖·‖ ≤ ½M2·R/c₀·‖J_∞δ‖. ∎
**Bemerkung:** Lemma X8 ersetzt Lemma X6 in Satz X4a vollständig; Kellogg–Warschawski wird für X4a nicht mehr gebraucht (die Literaturfrage
§3.8.1 ist damit für den Beweis gegenstandslos; Original weiterhin nicht gelesen). Preis: Faktor ≈ 4.9 in A_θ gegenüber der gemessenen
Bernstein-Konstante (19.7 statt 4.05).

### 3.9.3 Satz P-I mit Kontur d = 0.3 — was nötig ist und was fehlt

Satz P-I: sup_{D̄}|Π^P_N f − f| ≤ (Schlitze L, L̄ bis u₀(d)) + (1/2π)∮_{Γ_d}|f_K||E_N||dζ| mit |E_N(ζ,·)| ≤ C N^κ e^{−εNg(ζ)}/dist(ζ, K).
**Nachprüfbare Kurvenlage (x24, x26; Green-Funktion von W̄ explizit: T(z) = (z−L)/(z−L̄) bildet ext W̄ auf einen Sektor der Öffnung 158.934°
ab, dann Potenz 1/α und Halbebenen-Green-Funktion):** Γ_{0.3}: Länge 12.871, δ_d = dist(Γ_d, W̄) = 0.5095, max|z| = 2.357, Abstand zu den nächsten
Singularitäten −0.362 ± πi ≥ 1.231, M_d = max_{Γ_d}|f_K| = 4.251 (180 Punkte, f_K = exp∘f_K∘Log, keine Auswertung gescheitert) ⇒ Konturkonstante
M_d|Γ_d|/(2πδ_d) = **17.1** (d = 0.1: 37.9; 0.2: 22.0; 0.5: 14.2). Schlitzlänge (Green-Funktion auf L(1+u)): u₀(0.1) = 0.151, u₀(0.2) = 0.297,
**u₀(0.3) = 0.451** (|L(1+u₀)| = 1.995), u₀(0.5) = 0.799.
**Was fehlt:** (i) Holomorphie von f_K innerhalb Γ_{0.3} bis auf die Schlitze — [Skizze] (slog-Argument §3.5c; numerisch keine Störung); (ii) eine
Sprungschranke |Δf(L(1+u))| für **u ≤ 0.451** — heute nur bis u = 0.2 (C_J = 0.95). Das ist **dieselbe** Lücke wie in Lemma X5, Fassung R = 2
(u ≤ 0.455): mit einer Sprungschranke bis u ≈ 0.46 [P, Kneser-Motor] werden beide Reparaturen auf einmal sauber.
**Ohne (ii) — Hybrid-Kontur K3:** Niveaulinie d = 0.3 fern von L, L̄, Schlitz nur bis u = 0.2 (g = 0.134), auf [0.2, 0.451] grob |Δf| ≤ 2M_d und die
Schlitz-Rate ε_S ≈ 0.89 (W1: g_{K_T}/g_W̄ = 0.89–0.97). Zusatzterm (2M_d/π)C N^κe^{−0.12N}log(2.26)/c_d — ab A_θ ≥ 13 ändert er N₀ nicht (x23b).

### 3.9.4 N₀ mit beweisbaren Konstanten (x23b)

C_R = 12.7, κ_R = 1 (Lemma X7), A_θ beweisbar, Kontur K1 (d = 0.3); N₀ = kleinstes N mit θ_N ≤ ½:

| A_θ | 4.05 (gemessen) | 13.1 | 15.5 (3N-Zeilen) | **19.7 (bewiesen)** | 22.2 (c₀ = 0.4249) |
|---|---|---|---|---|---|
| K1 d = 0.3 | 172 | 314 | 342 | **387** | 412 |
| K2 d = 0.1 | 370 | 401 | 408 | 422 | 433 |
| K3 Hybrid | 188 | 314 | 342 | 387 | 412 |

(Ecke C_E = 1, κ_E = 2, ε = 0.69.) **Empfindlichkeit:** N₀ ≈ 121/172/244/314/387/482/629 für A_θ = 2/4/8/13/20/30/50 — etwa N₀ ∝ A_θ^{0.5}.
Eck-Konstanten bei A_θ = 19.7 (bzw. 15.5): κ_E = 1: N₀ ≈ 214 (193) unabhängig von C_E ≤ 1 — **auch mit perfekter Ecke liegt N₀ mit beweisbarem A_θ
bei ≈ 214**; κ_E = 2: 216/234/387/1186 für C_E = 0.01/0.1/1/10; κ_E = 2.5: 245/565/3009 für C_E = 0.01/0.1/1.
⇒ **Mit beweisbaren Konstanten ist N₀ ≤ 193 nicht erreichbar** (der Rest allein gibt ≈ 214); N₀ ≈ 390 braucht C_E ≤ 1, κ_E ≤ 2. Das macht die
Zertifikats-Reserve (§3.9.5) zum eigentlichen Hebel.

### 3.9.5 Reserve über 193: Kosten des Zweig-Zertifikats bis N = 400 bzw. 1000 (Schätzung aus Ps Protokollen p24_lauf1/2.log)

Gemessen (p24, N = 49…193, NC = 7 Tschebyschow-Knoten, Präzision 2M + 320 Bit): Zeit je Übergang **t ≈ 42 s·(M/193)^{2.8}** (Fit M ≥ 100), Summe bis 193:
2228 s. Wrapping-Term t3 (relativ, Reserve-Verbrauch) **∝ 1.88^M·r** (Fit: 1e-42 bei M = 49 → 3.3e-8 bei M = 193; r = Röhrenradius ≈ 1e-55…1e-50,
begrenzt durch den Interpolationsrest R^{−NC}, R ≈ 4e11…6e12).
**Hochrechnung:** aus t3 = 3.3e-8 bei r = 8.5e-51 (M = 193) folgt t3 ≈ 3.7e-11·1.88^M·r; t3 ≤ 1e-8 verlangt r ≤ 1e-8/(3.7e-11·1.88^M) ⇒ **r ≈ 5e-108 bei M = 400,
≈ 2e-272 bei M = 1000** (P schätzte 1.62^M — dann weniger). Heute r ≈ 1e-55 bei NC = 7 und R ≈ 1e11–1e12 (je Knoten ≈ 11.5 Dekaden) ⇒ **NC ≈ 12 (M = 400)
bzw. ≈ 26 (M = 1000)**. Die Präzision 2M + 320 Bit (1120 bzw. 2320 Bit ≙ 1e-337 bzw. 1e-698) reicht dafür; die
1.17^N-Verstärkung zufälliger Fehler (A) betrifft nur den nicht-strengen Prädiktor und ist durch 0.23N + 100 Bit gedeckt.
| bis | Zeit je Übergang (Ende) | Summe sequentiell | parallel (32 Kerne) | Speicher je Prozess |
|---|---|---|---|---|
| 400 | ≈ 42·(400/193)^{2.8}·(12/7) ≈ 550 s | ≈ 12 h | ≈ ½ h | ≈ 0.2 GB (arb ≈ 190 B/Eintrag, 4–6 Matrizen 401²) |
| 1000 | ≈ 42·(1000/193)^{2.8}·(26/7) ≈ 4.3 h | ≈ 32 Tage | ≈ 1–1.5 Tage | ≈ 2 GB (≈ 340 B/Eintrag, 1001²) ⇒ 32 Prozesse ≈ 64 GB |
**Parallelisierbar**, weil jeder Übergang N → N+1 mit einer vorab (nicht streng) berechneten Startbox geprüft werden kann; die Kette ist gültig, wenn
anschliessend jede End-Box in der nächsten Start-Box liegt (billige Prüfung). Gehört auf aicgolling (geteilter Server: Grafikkarte/Kerne vorher prüfen).
**Nutzen:** Zertifikat bis 400 deckt das Szenario „beweisbare Konstanten, Ecke C_E ≤ 1, κ_E ≤ 2“ (N₀ ≈ 387); bis 1000 zusätzlich κ_E = 2 mit C_E ≲ 5
oder κ_E = 2.5 mit C_E ≲ 0.1. Dazu bleibt (A) — zertifizierte f_K auf D̄ auf ≈ R_{N₀}/10 (≈ 1e-6 bei N₀ ≈ 400) — Ps Satz P-V.

### 3.9.6 (Fam-Eck): Bogenmodell — die Ecke sieht die Störung über einen Bogen fester Länge geglättet (x25, x25c)

Modell: Antwort = a·δr_ε + Re((c − id)·δμ_ε) mit den linearen Fixpunktdaten δL = η(L_ε)/(1 − μ), δr = Re(δL·L̄/ρ), δμ = η'(L_ε) + f''(L)δL,
ausgewertet am Poisson-Punkt L_ε = L(1 − ε) statt in L (= Mittel von η über einen Bogen der Länge ~ερ). (a, c, d) aus den sechs glatten Formen
(1, x, q, x², x³, x⁴), dann **Vorhersage** für p_N und die Zufallsstörungen (je Seed; st: k ≤ N/2, sv: alle k), Einheit R:

| N | bestes ε | RMS-Fehler der Vorhersagen | p_N vorhergesagt / gemessen | Zufall vorhergesagt / gemessen | Punktmodell (ε = 0) RMS |
|---|---|---|---|---|---|
| 200 | 0.07–0.08 | 0.31–0.40 | −0.07…−0.15 / −0.13 | 2.9/2.6, 3.3/2.7, 2.4/2.7, 1.7/1.6, 2.2/1.8, 1.5/1.7 (ε = 0.08) | 109 |
| 400 | 0.06 | 0.92 | +0.25 / −0.34 | 8.7/9.2, 1.1/3.1, 7.7/7.2, 5.9/6.2, 0.6/1.8, 5.2/4.8 | 1700 |
| 600 | 0.06 | 2.2 | +0.15 / −0.34 | 11.6/16.1, 1.0/1.1, 12.2/10.6, 0.7/0.8 | 3700 |

⇒ **Ein einziger Glättungsparameter ε ≈ 0.06, nahezu unabhängig von N (0.07–0.08 / 0.06 / 0.06), erklärt alle gemessenen Eckantworten** —
glatte Störungen, Zufallsstörungen und die echte Lösung — auf ≈ 1–2 R, wo das Punktmodell um den Faktor 10³–10⁴ danebenliegt (p_N: 8082 statt −0.34).
Die Glättung über einen Bogen **fester** Länge (≈ ερ ≈ 0.08, d. h. ±3.5° um L) ist damit numerisch belegt. Der Radialkoeffizient a(N) ≈ 23/96/170
(N = 200/400/600) wächst linear (≈ 0.37·N − 50) ⇒ **∂c/∂r ≈ 0.37**.
**Präzisierte Vermutung (Fam-Eck), Bogenform:** δσ₁/σ₁ = N·[∂_r c·δr(P_ε) + Re(∂_μ c·δμ(P_ε))] + O(R), P_ε(x) := P(x(1 − ε))-Daten bei L, ε ≈ 0.06.
Für die echte Lösung sind die geglätteten Daten klein (1/N-Buckel bei L, geglättet; Modellvorhersage −0.15…+0.25 R), also ist ihre Eckantwort O(R) —
in Übereinstimmung mit den Messungen. **Beweisweg mit W:** Ws Kegelmodell muss liefern, dass c nur von den auf Skala ε geglätteten Daten (L, μ, h) abhängt, mit ∂_r c ≈ 0.37
und ∂_μ c ≈ (0.042 − 0.035i) (N = 400: (c − id)/N); dann ist (Fam-Eck) die Ableitung von c, und der Beweis von (H₀) über das Modell trägt die Familie mit.

## 3.10 Meilenstein 7 (11.10. 04:15–): Zweig-Zertifikat 193 → 400 auf aicgolling; A_θ geschärft

### 3.10.1 A_θ — wo Lemma X8 Faktoren verschenkt (x27, x26)

**Die beste Konstante** des Hauptterms, A*(N) := max_{δ∈P_N} ‖δ'∘f‖_{L²(∂D)}/(N‖J_∞δ‖_{ℓ²_ρ}) (J_∞δ bis Zeile 6N; verallgemeinertes
Eigenwertproblem, Cholesky von T_{6N}J_∞ in arb, Werte δ_j'(f(x_i)) aus EINER hochgenauen Auswertung von f an 8N Punkten):

| N | 50 | 100 | 150 | 193 |
|---|---|---|---|---|
| A*(N) | 0.6099 | 0.6205 | 0.6246 | 0.6266 |
| Präzisionsprobe | PRECF 12/20 gleich | 8/12 gleich | 6/16 gleich | (8) |

(Achtung: mit zu wenig Präzision — < 12N + 600 Bit bei N = 50, < 8N + 600 bei N = 100 — kommen Werte bis 1e13 heraus; mit nur 3N Zeilen ist das
Problem schlecht gestellt.) **Lemma X8 gab 19.7, die gemessene Bernstein-Rechnung 4.05, das Optimum ist ≈ 0.63.**
**Aufschlüsselung von 19.7:** (i) Ableitung von J_∞δ über den höchsten Zeilenindex: (5/ρ)/min|f'| = 10.5 — das Extremal hat seine Koeffizienten
bei Index ≈ N, nicht 5N (Faktor ≈ 5); (ii) 1/min|f'| = 2.89 statt eines gewichteten Werts; (iii) der Term mδ' mit max|m|/(min|f'|ρc₀) = 9.1.
**Beweisbare Verbesserung (Lemma X8′, sofort):** δ'∘f = [(J_∞δ)' − m'δ − mδ']/f' termweise geteilt, also max|m/f'| = 2.1176 statt max|m|/min|f'|
= 6.76 und max|m'/f'| = max|f''∘f| = 0.8183:
  A_θ = (K/ρ)/min|f'| + max|m/f'|/(ρc₀) + M2/(c₀N) = **13.4** (K = 5, c₀ = 0.54); 12.6 (c₀ = 0.7319); **9.2** mit K = 3 (braucht X5 mit R = 2, d. h.
die Sprungschranke bis u = 0.455, dieselbe wie für die Kontur d = 0.3).
**Was darüber hinaus nötig wäre:** Der Faktor 5 aus (i) lässt sich mit Koeffizientenschranken allein nicht holen (zwischen 1.2N und 3N gibt es keine
Abklingschranke, β_α > 1). Der Weg zu A ≈ 1 ist eine L²-Bernstein-Ungleichung auf E = f(∂D) mit expliziter Konstante (Lemma X6 mit zertifizierter
konformer Abbildung) — rechnergestützt machbar, aber ein eigener Auftrag.

**N₀ (Kontur d = 0.3, C_R = 12.7, κ_R = 1) für die A_θ-Stufen:**

| A_θ | κ_E = 2, C_E = 1 | κ_E = 2, C_E = 0.1 | κ_E = 2.5, C_E = 0.1 | κ_E = 2.5, C_E = 1 | κ_E = 1, C_E = 1 |
|---|---|---|---|---|---|
| 0.63 (optimal) | 93 | 78 | 93 | 222 | 68 |
| 2 | 121 | 92 | 126 | 529 | 83 |
| 9.2 (X8′, 3N) | 262 | 166 | 334 | 1685 | 154 |
| **13.4 (X8′, bewiesen)** | **318** | 196 | 432 | 2243 | 182 |
| 19.7 (X8) | 387 | 234 | 565 | 3009 | 215 |

⇒ Mit X8′ (bewiesen) und einem Zertifikat bis 400 ist N₀ ≤ 400 für κ_E = 2, C_E ≤ 1 gesichert; mit dem Optimum 0.63 läge N₀ unter 100.


### 3.10.2 Zweig-Zertifikat 193 → 400 — vorbereitet und getestet, NICHT gelaufen (Team angehalten, 04:48)

Skripte in `zert/` (Kopie von Ps p24 mit NC-Parameter und vollen Newton-Schritten an den Knoten; Start-Boxen 2^{−(N+150)}ρ^{−k}; Kette per
zert_kette.py), Startanleitung `zert/STARTANLEITUNG.md`. Auf aicgolling getestet: 193 → 194 → 195 ok und Kette ok (je ≈ 27 s), 300 → 301 ok
(NC = 10, t1 = 2.8e-72, t3 = 9.6e-7, 309 s). Zwei Fehlstarts lehrten: (1) Ps Sehnen-Newton an den Knoten begrenzt r auf ~1e-51 ⇒ ab M ≈ 230
scheitert jeder Übergang (t3 bis 1e43) — volle Newton-Schritte beheben das; (2) die Startbox muss grösser sein als die Ausgabebox des Vorgängers
(2^{−(M+160)}), sonst reisst die Kette. Erwartete Dauer des vollen Laufs ≈ 2.5 h mit 12 Arbeitern.

## 4a. Wiederaufnahme (Stand 11.10. 04:48, Team angehalten)

**Stand der Buchhaltung (KONSTANTEN.md §00):** N₀ heute bewiesen ≈ 318 (Lemma X7: C_R = 12.7; Lemma X8′: A_θ = 13.4; Kontur d = 0.3 mit
Hybrid-Ausweg; Ecke C_E = 1, κ_E = 2 angenommen), repariert/gemessen bis ≈ 47–150; ohne Ecken-Annahme ist N₀ nicht bestimmt (W). Optimum
A*(N) ≈ 0.63 ⇒ N₀ < 100 möglich. **Nächste Schritte:** (1) Zertifikat bis 400 laufen lassen (zert/STARTANLEITUNG.md), P prüft die Kette;
(2) P: c_in, c_out zertifizieren; Sprungschranke am Schlitz bis u ≈ 0.46 (Kneser-Motor) — macht Kontur d = 0.3 und X5 mit R = 2 sauber
(A_θ dann 9.2); (3) W: Eck-Konstanten C_E, κ_E explizit; (Fam-Eck) als Ableitung der Ecken-Rate c (Bogenmodell §3.9.6 als Prüfstein);
(4) X: A_θ → A* über eine L²-Bernstein-Ungleichung auf f(∂D) mit expliziter Konstante. **Dateien:** BERICHT §0e, §3.7–§3.10; KONSTANTEN.md;
x17–x27 (Skripte/Logs); zert/ (Zertifikat).

## 4. Offene Punkte / nächste Schritte (Stand 11.10. 03:25)

1. **(Fam-Eck) beweisen:** Glattheit der Ecken-Asymptotik in den Eckdaten (mit W/O). Prüfbar jetzt: Hängt δσ₁/σ₁ für glatte η nur von geglätteten
   Daten ab (radialer Fixpunkt-Versatz, μ)? Test: zwei Formen mit gleichen geglätteten Daten, verschiedener Feinstruktur ⇒ gleiche Antwort.
2. **(K_d):** Überlaufprofile der Eckmoden aus Ws Kegelmodell (Vorfaktor-Ebene) oder direkt aus Os T_L; Gram-Bedingung auf [0, ε₀].
   N = 800 mit höherer Präzision nachmessen (d = 16 war unbrauchbar).
3. **Tangente in (H-i):** nur auf 𝒵 verlangen (Text §3.8.3) — P trägt in ARCHITEKTUR §9 ein.
4. Lemma X5: Konstante c₁ (Schlitzufer) am Ufer selbst auswerten; Fassung R = 2 braucht die Sprungschranke bis u = 0.455.

## 5. Dateien (x_alternativen/)

x01–x03 (Vorzeichen, Halbmond-Rauschen, Möbius-SL(2)) · x04, x05 (Rechteck-Abschnitte, Lemma 3.1; Logs x04_N300, x05_*) · x06 (Modelle kubisch/radial) ·
x07 (Kneser Rechteck, Π̂f_K) · x08 (Pol-Antworten quadratisch/rechteckig) · x09 (kubisch auf ∂K': R, R̂, Rc; json/npz/log je N) ·
x10 (Überlauf der Eckmoden; Befund: Überlauf ist für alle Vektoren riesig, kein brauchbares Mass) · x11 (Eckoperator 𝒞_N; npz für W) ·
x12 (Kneser auf Ps ∂K'-Gitter; json je N,d) · x12b (Θ im K'-Rahmen mit Ps kprime) · x13 (Normalgleichungen, grosse d, Kneser oder p_N; σ_min) ·
x14 (M_R auf |x| = 1.65, 2) · x15, x15b (Θ_∞-Profil) · x16 (Familienvergleich zweier json) · x17 (sup|p_N − f_K|, Bernstein-Faktor) ·
x18 (Eckmoden σ₁…σ₄ von A_P⁻¹ unter kohärenter/zufälliger Störung, Sekante und Tangente; Logs x18_tangente.log, x18_sekante.log) ·
x18b (Diagnose: Grösse der Tangenten-Störung) · x19 (Fixpunkt-Versatz und Multiplikator der Störungen) · x19b (Gestalt von p_N − T_Nf auf dem Kreis) ·
x20 (σ_min(R_d) für kleine d, Überlaufprofile der Eckmoden; Logs x20.log, x20_N800.log) · x21–x26 (Konstanten, Eingangsnormen, N₀-Bilanz gemessen/beweisbar, Kontur, Bogenmodell; siehe KONSTANTEN.md) · Logs x18_formen*.log, x18_seeds.log · brett_eintrag.py, brett1–15.txt.
