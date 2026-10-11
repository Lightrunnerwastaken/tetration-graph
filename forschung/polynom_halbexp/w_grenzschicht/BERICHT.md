# Strang W — Grenzschicht/WKB der Eckmoden → (E_L)

Stand: 11.10.2026, 00:22, ergänzt 00:55 (§7–§11) (Zeiten aus `date`). Kennzeichnung: **bewiesen** / **Skizze** / **numerisch** / **Vermutung**.
Testbett: reelles kubisches Modell von D (§9.3), f(x) = x + (x² − 2Re(L)x + |L|²)(α + βx), f'(L) = √L; umskaliert t = x/ρ, ρ = |L|.
Bezeichnungen wie D: J̃⁻¹ = Σσ_i u_i v_iᵀ (ℓ²_ρ), u_i Ausgabe, v_i Eingabe (Gleichungen); **y := J̃^{−T}u₁ = σ₁v₁** (y-Normierung:
J̃ᵀy = u₁, also |y| = σ₁); adjungierte Mode als Laurent-Polynom y(t) = Σ_k y_k t^{−k} (= σ₁V₁(ρt)). Γin = {|f| < ρ} (Bogen von L über −ρ
nach L̄), Γout der Rest. Phasenraum = (Winkel φ auf |t| = 1, Index x = k/N).

## 00000000. Kurzfassung Weg (A), n-te-Wurzel-Stufe (Stand 04:32)

**Klebungslücken (§26.1):** bindend ist die Naht Γin mit δ = Λ(τ) (0.001 bei 1°, 0.012 bei 103°); beide Lücken > 0 ausser in L, L̄ — auf n-ter-Wurzel-Ebene
genügt das. **Satz W11 (bewiesen):** Kompaktheit + „tropische Bilanz“ der Grenzraten ohne jede Orthogonalität; die Fehler sind darin enthalten. Offen ist
**Vermutung W12** (Eindeutigkeit/Extremalität des Grenzproblems, ein Hindernisproblem mit nichtlokaler Kopplung) ⇒ c und (A1) in n-ter Wurzel.
**Entscheidend (§26.3):** e^{o(N)} reicht an der Ecke NICHT, auch nicht für Stufe (i) — f_Ks Singularität sitzt in L mit nur polynomialem Abfall;
(i) braucht polynomiales (E_L) und damit starke Asymptotik in einer festen Umgebung von L, die nur global (Deift–Zhou-artig) zu haben ist.
Problemblatt: OFFENES_PROBLEM.md.

## 0000000. Kurzfassung „(H₀) global“, Schritt 1–2 (Stand 04:22)

**Schritt 1 entschieden: keine Gonchar–Rakhmanov-Orthogonalität (§25.2).** Die exakten Relationen koppeln in jeder Bedingung ein durch f
transportiertes Moment in a₀ = f(0) mit einem Moment in 0 (gleicher Index) — nicht von der Form {s^jW}. Auf der Grenzfläche ℳ wäre es GR-Form,
aber ℳ existiert nur im Grenzwert (Klebung = Dominanzbilanz, keine Identität). **Schritt 2:** klassische S-Eigenschaft nicht anwendbar (geschlossener
Träger, einseitiges Green-Potential); auf ℳ gilt das Analogon (Glattheit über die Naht, numerisch 5e-4). **Neu:** J_N⁻¹ = Σ Q_nA_nᵀ/h_n mit exakten
Interpolanten; ‖Q_N‖‖A_N‖/|h_N| trägt die Rate c. **Vorschlag:** GR mit Stabilität unter exponentiell kleinen Klebungsfehlern auf approximativen
Flächen ℳ_N (+ Os Innenschicht); Vorbild Muttalib–Borodin für die Kompositionsstruktur.

## 000000. Kurzfassung Runde 7 (Stand 04:07)

**Wachstumsschranke (I) nicht bewiesen (§24.1):** σ₁(ℓ²_ρ) ≤ (ρ'/ρ)^N‖J_N⁻¹‖_{ℓ²_{ρ'}} ist bewiesen, aber es gibt keinen Skalenkreis ρ' > |L| ohne Wachstum
(Eckmode nur verzögert, ab 1.03|L| ein zweiter Modus). **Eckzone (§24.2):** Ratenebene summierbar, Zone nur 5.6·ln N h-Schritte, für y ausserhalb D eine
kontrahierende Rekursion (Faktor 1/|L| = 0.727, Lemma W9) — aber das Maximumprinzip braucht am Aussenrand winkelaufgelöste, polynomial genaue Daten
(globale Gestalt); die Eckzone ist kein lokales Störungsproblem. **Vorschlag:** (H₀) als benannte Hypothese (gemessen bis 1000: C_E = 1, κ_E = 0.6,
ε_eff ≥ 0.60; an X gemeldet), Beweis als eigenes langes Vorhaben.

## 00000. Kurzfassung Runde 6 (Stand 03:48, Nachtrag 03:58: §23.6 — Modell trifft Xs ∂c/∂μ, W1'' ohne N^{κ+2})

**Extremalität trägt nicht (§23.1):** ‖J̃‖ = e^{0.47N}, daher ist der ℓ²-Rayleigh-Quotient für Testvektoren unbrauchbar — obwohl ein rein aus dem
Modell gebauter Testvektor die Eckmode in ℓ² auf 1–5 % trifft (Überlapp 0.985/0.997 bei N = 400/800). **Reduktion (§23.2):** fern von L genügt eine
Wachstumsschranke σ₁ ≤ e^{c'N}; an der Ecke scheitert Ps K'-Rahmen asymptotisch (Kegel-Exponent 1.118 < 1.27), dort braucht es (H*)_loc im W̄-Rahmen;
die Innenscheibe |t − L| ≤ r_N = A·N^{−0.894} folgt per Maximumprinzip aus einer Kreisbedingung (KB) — numerisch polynomial (~N^{0.9}), und die
Kreisprofile kollabieren genau auf dieser Skala. **Konstanten für X (§23.3):** |Y| ≤ e^{−4.06}N^{0.88} auf γ*, |T| ≤ e^{−4.60}N^{1.03} auf ∂K_T,
ε_eff ≥ 0.60/0.80 (d = 0.05/0.1) für 150 ≤ N ≤ 1000 (gemessen). **Schnittstelle O (§23.4):** O liefert (KB), ich den Aussenteil.

## 0000. Kurzfassung Runde 5 (Stand 03:31)

**(H-V1)/(H-V3) für alle N: nicht bewiesen, aber auf EINEN Standardsatz reduziert (§22).** Lemma W3: die Buchtschranke ist dieselbe Aussage wie
die T-Schranke auf dem A-Bogen (Kreisschluss intrinsisch). Satz W1'': eine Hypothese (H*) auf einer Kurve γ* ersetzt (H-I), (H-V1), (H-V3).
Das Leitordnungsmodell ist die Green-Funktion des Gebiets Ω_y zwischen Buchtrand und Linsen-Nulllinie; die Nullstellen der Eckmode verteilen sich
nach dessen harmonischem Mass (N = 1000: Anteil innen 0.470 gegen 0.473, Winkel-KS 0.005). Satz W5: Nullstellenverteilung + Normierung ⇒ (H*) bis
e^{o(N)}. Fehlend: diese Nullstellenverteilung für alle N (Extremalität, Gonchar–Rakhmanov-Typ) und die Ecke (O, Doppel-Skalenmodell).
Korrektur: c_Modell = 0.0118 ± 0.0004 (Kondition der Ladungssimulation geprüft, §22.7).

## 000. Kurzfassung Runde 4 (Stand 02:40)

**Lokale Ecktheorie (§19, Skizze + numerisch):** nahe L lebt die Ratenfunktion auf einem Kegel mit Deckverschiebung 2πi − log μ; H ≈ Re(c₁ζ^{κ₁}),
**κ₁ = 2πi/(2πi − log μ) = 1.118 − 0.032i** (Daten: κ = 1.114 − 0.018i). Re κ₁ > 1 ist strukturell (⇔ arg μ(2π − arg μ) > ln²|μ|) — die Bedingung p > 1
von Satz W1. **c hergeleitet (§20, numerisch):** Leitordnungsproblem §15 + Eckfunktionen ergibt ohne freie Parameter **c = 0.0118** (gemessen 0.0117/0.0120),
Profil ab 8° auf ≤ 0.0003, Buchtrand 0.814 (0.815). Innenschicht trägt kein c (nur den Eckvorfaktor). Offen: Profil nahe L (τ ≤ 4°, Daten unter dem
Modell), Beweis von (H-V1)/(H-V3) für alle N (Plan §21).

## 00. Kurzfassung Runde 2 (Stand 01:29)

**Entscheidung (§14, numerisch, kubisch UND Kneser, N ≤ 1000, genaue Daten):** der (E_L)-Weg trägt in der W̄-Ecke. Poisson–Jensen-Kriterium
Θ_∞ = 0.22 (d = 0.05) / 0.10 (d = 0.1); auf dem Schlitz negativ; sogar termweise ohne A-Bogen-Auslöschung 0.23/0.12. Eckteil-Exponent auf dem
D-Bogen nahe L 1.5–1.6 (> 1.13); A-Bogen fällt exponentiell; Engpass ist der ferne D-Bogen auf Γ_d (Reserve 4.5 bzw. 10). Meine Vermutung W1
(1.27 = 1/α') zurückgezogen. Liste der durch die acb_series-Kappung verfälschten D-Zahlen an P (BRETT 00:56).
**O(c)-Theorie (§15, Skizze):** geschlossenes Leitordnungsproblem (H harmonisch auf Ĉ∖Bucht, H = H∘h auf dem Halbmond, freier Buchtrand) trifft
Buchtrand und Stegfrequenzen, aber nicht c (0.018 statt 0.012) — offen ist die Ecke L.

## 0. Kurzfassung (Runde 1)

1. **Phasenraum-Gestalt (numerisch, N ≤ 1000):** beide Eckmoden leben auf demselben „Bogen“ im Phasenraum (Bild w06_phasenraum_N1000.png):
   zwei **Beine** bei φ ≈ arg L und φ ≈ −arg L (dicht innerhalb Γin, alle Indizes 0 < x < 0.85) und ein **Steg** bei x ≈ 0.90–0.97 quer
   über ganz Γin. u₁ sitzt auf den Beinen (Masse bei x ≈ 0.4–0.75, Winkel L), v₁ auf dem Steg (x ≈ 0.95, Maximum bei 180°). Der Steg ist
   zur Kante x = 1 hin gaussisch abgeschnitten (|v_k| ≈ e^{−30N(x−0.953)²}). Die Eckmode ist also eine **Randmode des Phasenraumrechtecks
   Γin × [0, 1]** (Sprünge des Toeplitz-Symbols bei L, L̄ = Beine; Abschneidekante = Steg).
2. **Leiter von Modenpaaren (numerisch):** σ₃, σ₄ folgen σ₂ mit Verzögerung ~500–600 in N (σ₃(900) = 3.54 ≈ σ₂(400) = 3.45), gleiche
   Endrate; Fit ln σ = A + cN − a ln N über N = 850–950: **c = 0.0120 (σ₁ und σ₂), a₁ = 0.81, a₂ = 1.29** (kubisch = Kneser-Wert 0.0120).
   Passt zu O (00:14: r(N) wächst ~N/250).
3. **Asymptotisches Randprofil (numerisch, Trennung κ log N + λN über N = 600…1000):** λ(103°) = 0.01198 = c. Γin-Seite:
   λ(τ) − λ(0) = 0.0010 / 0.0014 / 0.0024 / 0.0032 / 0.0045 / 0.0054 / 0.0067 / 0.0089 / 0.0110 bei τ = 1.5 / 2 / 3 / 4 / 6 / 8 / 12 / 24 / 48°.
   **Lokaler Exponent 1.26–1.29 auf 1.5°–3°**, danach fallend (konkav: 1.06, 0.88, 0.75, 0.62, … bei 3–4°, 4–5°, 5–6°, 6–8°). λ/τ^{1.27}
   ist auf [1.5°, 3°] konstant (0.000595 ± 1 %); τ^{1.13} und τ^{1.5} passen schlechter. Unterhalb ~1° Innenschicht (polynomial).
   **Ds „p ≈ 1.44–1.52“ enthält den polynomialen Vorfaktor** (Raten aus zwei N statt Trennung); die Rohraten fallen mit N (600→700: 1.51,
   900→1000: 1.32 auf 2°–4°).
4. **Γout-Seite: λ = 0, nicht negativ (numerisch, wichtig für Ps Formulierung):** dort ist |y| auf dem Sockel e^{−17…−20}·N^{−2.4}
   (Fit: κ = −2.4, λ = −0.0002 ± 0.0004). Die bei festem N gemessenen negativen „Raten“ (P 00:04: −0.0055 … −0.0157) sind
   (1/N)ln(Sockel/Eckwert), nicht exponentiell; für N → ∞ geht das gegen 0. **Asymptotisch gilt also λ ≤ 0 auf Γout (genauer λ = 0 mit
   polynomial kleinerem Vorfaktor), aber NICHT λ ≤ −B|τ|^q.**
5. **n-te-Wurzel-Struktur = gespiegelte Green-Funktion (numerisch):** Y(w) := Σ_k y_k w^k (w = 1/t) ist asymptotisch ein
   **Bernstein–Walsh-extremales Polynom** einer Menge K_y ⊃ D̄ (im Inneren: (1/N)ln|y(t)| ≈ g_{Ĉ∖K_y}(1/t̄)); Nullstellen von y auf Γin
   knapp ausserhalb |t| = 1 und auf einer „Bucht“-Kurve im Inneren auf der Γout-Seite (r ≈ 0.815 bei 0°, stabil N = 400…800).
   Lokale Frequenz auf dem Steg = Normalableitung der Green-Funktion: x(φ)/∂_n g_{Ĉ∖W̄} = 1.011 (135°–180°), N-unabhängig (N = 500, 700,
   1000). Mit Zwei-Kreis-Mengen passt K_y mit rechtem Scheitel b ≈ 1.24 (statt W̄: 1.346) auf ±0.01; W̄ selbst weicht um +0.025 ab.
   Die Rate c ≈ 0.012 und das Profil λ sind eine O(c)-Korrektur zu dieser Leitordnung.
6. **Mechanismus der Bilanz (numerisch, radiales Modell von M, N = 150):** die Stegbeiträge des Multiplikator-Terms (m·Y) werden in
   (Jᵀy)_k = 0 vor allem vom Kompositionsterm AUS DEN BEINEN (|φ ∓ θ| < 6°) ausgeglichen, nicht von Γout — die Beine koppeln über die
   komplexe Streckung x ↦ μx alle Frequenzen.
7. **Vermutung W1 (λ ≍ τ^{1/α'}, α' = 1 − arg μ/π) — durch den arg-μ-Test NICHT bestätigt (§11):** das Profil nahe L ist fast universell
   (Exponent 1.1–1.3 auf 2°–4° für arg μ = 0.3/0.669/1.0). Neu (§7–§10): zwei Rechenfallen (float innen, acb_series-Kappung bei 10 Termen in
   Ds d22/d28/d29), kopffreier Eckanteil auf ∂W̄ genau (A-Bogen fällt exponentiell; aussen Θ ≈ 0), Leitordnungsmodell (Dipol-Green-Funktion),
   und die drei Strukturaussagen (V1)–(V3), auf die sich (E_K) zurückführen lässt.

## 1. Daten und Methode

- w_kern.py: Modelle, J_N in arb (wie D: J[k][j] = [x^k]f^j + [x^{k−j}]f'∘f), volle Inverse umskaliert (Präzision 0.9N+300 Bit), SVD in float64.
- w01_eckmode_daten.py: N = 150 … 1000 (N = 1000 in 9.5 min). σ (kubisch): 1.4904 (150), 1.8524 (200), 3.6983 (300), 9.1002 (400), 24.766 (500),
  70.710 (600), 207.39 (700), 619.08 (800), 1871.7 (900), 5713.7 (1000); σ₂: 1.80 (300), 3.45, 8.00, 20.48, 55.28, 153.9, 438.1, 1267.3 (1000);
  σ₃ = 1.42, 2.00, 3.54, 7.32 und σ₄ = 1.40, 1.95, 3.42 (N = 700 … 1000).
- y = J̃^{−T}u₁ hochgenau (arb, 1.2N+400 bzw. 1.5N+400 Bit) für Auswertungen innerhalb |t| < 1 (float-Horner ist dort wertlos: Summanden
  bis 10^{80}). Kontrolle: 𝒴 = Σ_m c_m (t−ã₀)^{−m−1} (c_m = ℓ_y((F̃−ã₀)^m)) stimmt mit Σ_{F̃(z)=t} y(z)/(zF̃'(z)) überein (3 Punkte, 5 Stellen).

## 2. Phasenraum (w02, w06)

Gabor-Transformation im Index (Fensterbreite 0.6√N), Bild w06_phasenraum_N1000.png. Kennwerte (N = 600, v₁): Fensterzentrum k/N → Winkel
des Maximums: 0.6 → 77.4°, 0.7 → 78.8°, 0.8 → 83.5°, 0.85 → 90.5°, 0.88 → 99.8°, 0.9 → 111.5°, 0.92 → 131.5°, ≥ 0.94 → 180°.
Indexprofil v₁: log10|v_k| = −4.1 (0.5), −3.0 (0.7), −2.5 (0.8), −1.6 (0.9), −0.7 (0.95, Max), −2.9 (0.97), −6.3 (0.98), −11.4 (0.99), −19.8 (1.0).
u₁: Masse 5/50/95 % bei k/N = 0.39/0.59/0.75, Plateau 10^{−2.2} auf [0.8, 0.95], dann −15 (0.99).

## 3. Profil nahe L (w03, w08)

Λ_N(φ) := (Σ_{i≤2}σ_i²|V_i(ρe^{iφ})|²)^{1/2} (invariant gegen die Mischung der beiden reellen Moden). Fit log Λ_N = a + κ log N + λN über
N = 600, 700, 800, 900, 1000 (w08): siehe Kurzfassung 3. κ(0) = 0.63 (Ecke ~N^{0.6}); κ(τ) fällt nach aussen auf −0.5 (Steg).
Rohraten zwischen benachbarten N (w03, 900→1000): 0.00081/0.00143/0.00202/0.00306/0.00390/0.00521/0.00759/0.01091 bei τ = 2/3/4/6/8/12/24/103°.
Innenschicht τ ≲ 1°: Rohrate ≈ +0.00025 = (κ(τ) − κ(0))/N, kein exponentieller Anteil.

## 4. Gestalt in der Ebene (w04, w05, w07, w10, w11, w13, w15)

- Nullstellen (float, nur nahe |t| = 1 verlässlich): Dichte im Winkel gleichverteilt (≈ N/2π); u₁: Γin innen (N(r−1) ≈ −12…−16), Γout aussen
  (+23…+30); v₁ spiegelbildlich. Der Abstand der Stegnullstellen von |t| = 1 ist (λN + ln(Steg/Sockel))/(xN) — für N → ∞ endlich (≈ λ/x).
- Genaue Werte (arb, N = 400/600/800): auf der Γout-Seite bleibt y bis tief ins Innere auf dem Sockel (Bucht bis r ≈ 0.815 bei 0°, 0.84 bei
  30°, 0.89 bei 50°; N-unabhängig). y(h(A-Bogen)) ist daher klein (10^{−7.3} am Scheitel, N = 400), ebenso 𝒴 auf dem A-Bogen.
- Leitordnung: (1/N)ln|y(t)| ≈ g_{Ĉ∖K_y}(1/t̄) + O(c), K_y ⊃ D̄ (Zwei-Kreis-Fit: Scheitel b ≈ 1.24; Tabelle w15). Grund (Skizze): die
  Padé-Restfunktion Q := m̃y + t𝒴 − (Polynom) − t𝒰 ist O(t^{−N−1}) bei ∞ mit Polen nur in 0, ã₀ (exakt, aus Jᵀy = u), also Bernstein–Walsh-
  extremal im Äusseren; auf dem Steg ist Q ≈ m̃y, und die harmonische Fortsetzung von −g über den D-Bogen ist +g(1/t̄) (Schwarz-Spiegelung) —
  das erklärt die gespiegelte Green-Funktion im Inneren und die Frequenz = Gleichgewichtsdichte auf dem Steg.

## 5. Exakte Umformung der adjungierten Gleichung (bewiesen)

Für ein Polynom f vom Grad d (hier d = 3) und y ∈ ℂ^{N+1} sei y(t) = Σ_{k≤N} y_k t^{−k}, 𝒴(t) := Σ_{F̃(z)=t} y(z)/(z F̃'(z)). Dann
(i) 𝒴(t) = Σ_{m=0}^{N} c_m (t − ã₀)^{−m−1} mit c_m = Σ_k y_k [t^k](F̃ − ã₀)^m (rational, einziger Pol in ã₀ = F̃(0) der Ordnung N+1;
    die Pole an kritischen Werten heben sich in der Summe über die Äste weg);
(ii) Σ_{k≥0} (J̃ᵀy)_k t^{−k−1} = m̃(t)y(t)/t + 𝒴(t) + P(t) als Laurent-Reihe bei ∞ (P Polynom), d. h. Jᵀy = u ⇔
    m̃y/t + 𝒴 − 𝒰 = T − P mit 𝒰 = Σ_{k≤N}u_k t^{−k−1} und T = O(t^{−N−2}) analytisch auf Ĉ∖{0, ã₀}.
*Beweis.* (J̃ᵀy)_k = (1/2πi)∮ y(z)[F̃(z)^k + z^k m̃(z)]dz/z; Erzeugende Σ_k(·)w^{−k−1} = (1/2πi)∮ y(z)[1/(w−F̃(z)) + m̃(z)/(w−z)]dz/z
für grosses |w|; Residuen ausserhalb des kleinen Kreises (Urbilder von w, z = w, ∞) geben (ii); (i) durch Entwicklung von 1/(w − F̃(z)) um ã₀,
(F̃ − ã₀)^m = O(z^m) bricht bei m = N ab; Symmetrie in den Ästen beseitigt die Verzweigungspunkte. ∎ (Numerisch geprüft, w07.)
Das ist Ds transponierte Gleichung (D §11.3) in Laurent-Form, mit Cs Rationalitätsstruktur (Pole 0, ã₀).

## 6. Offene Punkte / nächste Schritte

1. Test von Vermutung W1 (Exponent folgt 1/α' = 1/(1 − arg μ/π)?) — läuft (w16, arg μ = 0.3 und 1.0).
2. Grössere N (≈ 1300) für die Kleine-τ-Asymptotik (Innenschicht von der äusseren Potenz trennen).
3. Theorie der O(c)-Korrektur: freies Randwertproblem (Stegnullstellenkurve, Bucht) mit Quellen bei L, L̄; daraus c und λ.
4. Kneser (A_f mit Kneser-Koeffizienten) zum Abgleich.

---

# Runde 1b (00:25–00:55): Eckanteil auf ∂W̄ genau, Leitordnungsmodell, zwei Rechenfehler im Team

## 7. Zwei numerische Fallen (wichtig für alle)

1. **float-Horner innerhalb |t| < 1 ist wertlos** (Summanden bis 10^{30…80}); y = σ₁v₁ muss dort aus J̃^{−T}u₁ in arb ausgewertet werden
   (w10, w17). Ps p10-Werte „Urbildterm auf dem A-Bogen +0.031“ waren dieses Rauschen (P hat zurückgezogen, BRETT 00:37). Genau:
   (1/N)ln|y(h(η))| auf dem A-Bogen = +0.0012/−0.0020/−0.0051/−0.0089/−0.0145/−0.0188/−0.0241 bei |η−L| = 0.02/0.05/0.1/0.19/0.38/0.6/1.0 (N = 600).
2. **python-flint `acb_series(...).inv()` liefert nur ctx.cap = 10 Terme** (Rest Nullen). Betroffen: Ds d22, d28, d29 (kein ctx.cap) ⇒ Ds
   Schlitz- und Eckwerte (§9.4, §10.4, §11.3) sind falsch; D-Bogen-Werte kaum (Urbild dort ausserhalb D). Gefunden, weil der kopffreie Eckanteil
   die Bernstein–Walsh-Schranke verletzte. Neumessung (w23, mit Fix; Eckanteil |⟨u₁,·⟩| mit Kopf / kopffrei):

| N | L | u = 0.001 | u = 0.003 | u = 0.01 | u = 0.03 | u = 0.1 | τ = 2° | τ = 8° | τ = 24° |
|---|---|---|---|---|---|---|---|---|---|
| 400 | 2.47 / 3.35 | 2.23 / 2.75 | 1.81 / 1.85 | 0.826 / 0.457 | 0.0894 / 7.8e-3 | 1.5e-3 / 3.1e-9 | 2.84 | 3.14 | 4.97 |
| 600 | 3.14 / 4.28 | 2.71 / 3.22 | 2.00 / 1.81 | 0.643 / 0.239 | 0.0351 / 6.2e-4 | 5.5e-4 / 2.1e-13 | 3.87 | 6.69 | 22.2 |
| 800 | 3.65 / 5.00 | 3.03 / 3.46 | 2.03 / 1.64 | 0.467 / 0.117 | 0.0168 / 4.8e-5 | 2.7e-4 / 1.5e-17 | 5.06 | 15.3 | 108 |
| 1000 | 4.10 / 5.62 | 3.25 / 3.58 | 2.00 / 1.44 | 0.334 / 0.0562 | 0.0095 / 3.7e-6 | 1.5e-4 / 1.1e-21 | 6.61 | 37.0 | 547 |

   (D, N = 400, cap 10: L 1.81, u = 0.003: 1.19, u = 0.03: 0.041.) Qualitativ bleibt Ds Bild (Ecke ~N^{0.6}, Schlitz fallend); kopffrei fällt der
   Eckanteil auf dem Schlitz wie e^{−εNg} mit **ε ≈ 0.73 (u = 0.01, 0.03, 0.1; Raten 800→1000)** (g(L(1+u)) aus C §12.2).

## 8. Kopffreier Eckanteil T(η) = ⟨u₁, y_η⟩ auf ∂W̄ und aussen (w20, w21; genau, mit Fix)

T ist rational (Pole 0, a₀), O(η^{−N−2}) ⇒ Poisson–Jensen auf Ĉ∖W̄ (Kontrolle: Schranke an ζ = 3, −3, 2+2i, 5 erfüllt, fast Gleichheit).
Raten (1/N)ln|T|, N = 300…800 (u-Normierung, Ecke L ≈ N^{0.6}):
- **D-Bogen:** +0.0014 (τ = 2°) … +0.011 (180°) — wie das Profil von V₁ (§3).
- **A-Bogen: exponentieller ABFALL** −0.005 (|η−L| = 0.1), −0.008 (0.19), −0.0124 (0.38), −0.020 (1.0) — nicht nur Sockel.
- **Aussen auf der Normalen aus L:** (1/N)ln|T| = −0.0122/−0.048/−0.100 bei r = 0.03/0.1/0.2, für N = 400 und 800 gleich, ≈ −g(η) = −0.012/−0.045/−0.098:
  Θ ≈ 0 … −0.06. Der A-Bogen-Abfall gleicht den D-Bogen-Anstieg nahe L mehr als aus (vorzeichenwechselndes Profil im Sinne von C §13.4 (a)).
**Folgerung für die Zielform:** im W̄-Rahmen ist die richtige zweiseitige Form von (E_L) die auf ∂W̄ (D-Bogen UND A-Bogen), nicht auf Γin/Γout.
Termweise (Eckformel D6) erhält man auf dem A-Bogen nur λ_A ≤ 0 (beide D6-Terme auf dem Sockel, siehe §10 (V3)); dann muss der D-Bogen allein
p > 1/α = 1.13 liefern. Der wahre Abfall auf dem A-Bogen entsteht durch Auslöschung der D6-Terme (Sockel gegen Sockel).

## 9. Leitordnungsmodell (n-te Wurzel): Dipol-Green-Funktion der Linse Λ_h (Skizze + numerischer Test, w18, w19)

Bausteine (Skizze): (i) T = m̃y/t + 𝒴 − 𝒰 + P exakt (§5); (ii) T nullstellenfrei aussen ⇒ Φ_T = ω(λ) − g; auf dem Steg T ≈ m̃y/t, also
setzen sich Φ_y (innen) und Φ_T (aussen) harmonisch über Γin fort; (iii) y bleibt auf dem Sockel in der Bucht, deren Rand = h(A-Bogen)
(aussen jenseits des A-Bogens ist T klein, also y(h(t)) klein); (iv) Pole: Φ_y ~ −log|t| bei 0, Φ_T ~ −log|t| bei ∞.
⇒ **Φ_B = g_𝒟(t,0) − g_𝒟(t,∞)** auf 𝒟 = Ĉ∖Λ_h, Λ_h = Linse zwischen A-Bogen und h(A-Bogen) (Öffnungswinkel arg μ = 38.3° bei L); eindeutig, keine
freien Parameter. Vorhersage gegen Messung (kubisch):
- **Stegfrequenz x = −∂_rΦ_B:** 0.9563 (180°) / 0.939 (125°) / 0.879 (93°) / 0.775 (81°) gegen gemessen 0.957 / 0.936 / 0.88 / 0.787 — sehr gut.
- **Inneres:** Φ_y − Φ_B ≈ konstant −0.018 ± 0.004 auf r = 0.3…0.9 (Γin-Seite) — Form richtig.
- **Rate und Profil falsch:** c_pred = Φ_B(−1) = 0.025 gegen 0.0120; nahe L Φ_B ~ r^{0.56} (erste Eckmode des 321.7°-Winkels) gegen gemessen ~τ^{1.1…1.3}.
  ⇒ die Leitordnung erklärt die grosse Struktur (Frequenzen, Inneres), aber c und λ(τ) sind eine Korrektur derselben Grössenordnung, die
  Eckstruktur bei L (Nullstellen von T/y nahe L, Innenschicht) fehlt im Modell. Bucht gemessen bei r = 0.815 (0°) gegen h(A-Scheitel) 0.827.

## 10. Was ein Beweis von (E_K) über die Eckformel konkret braucht (W̄-Rahmen; Skizze)

Mit D6 und Poisson–Jensen auf Ĉ∖W̄ genügen drei Eigenschaften der adjungierten Eckmoden (für alle r(N) Moden gemeinsam, vgl. O 00:33):
- **(V1) Profil auf dem D-Bogen:** σ|V(ρe^{i(arg L+τ)})| ≤ CN^κ e^{λ(τ)N}, λ ≤ Aτ^p mit p > 1/α = 1.1326 nahe L, λ < 0.068 sonst.
- **(V2) Linsenbegrenzung:** σ|V(ζ)| ≤ CN^κ für ζ ausserhalb D̄ ausserhalb einer dünnen Linse (Breite ≈ λ/x) längs Γin — dann ist der D6-Urbildterm
  auf dem D-Bogen polynomial (Ps formale Lücke (2) verschwindet; die äussere Poisson-Schranke überschätzt, weil V auf Z_bar Nullstellen hat).
- **(V3) Bucht:** σ|V(ζ)| ≤ CN^κ auf h(A-Bogen) (⊂ Bucht) — dann sind die D6-Terme auf dem A-Bogen polynomial (Ps Lücke (1) verschwindet).
Numerisch: (V1) mit p ≈ 1.1…1.3 auf 2°–4° (N ≤ 1000, siehe unten), (V2) und (V3) erfüllt (§4, §7). Alle drei sind Aussagen über die
n-te-Wurzel-Struktur der adjungierten Mode — genau das, was die Grenzschicht-Theorie liefern muss.

## 11. Test Vermutung W1 (Exponent = 1/α' = 1/(1 − arg μ/π)) — nicht bestätigt

Kubisch mit |μ| = 1.1724, arg μ = 0.3 / 0.669 / 1.0 rad (w16, N = 500…900, Fit wie w08): c = 0.0120 / 0.0121 / 0.0135; λ(4°) = 0.0036 / 0.0032 /
0.0030; λ(24°) = 0.0092 / 0.0089 / 0.0095; lokale Exponenten auf 2°–4°: 1.28/0.93/0.93 · 1.17/1.27/1.15 · 1.38/1.08/1.24 (Rauschen ±0.2).
Vorhersage W1: 1.106 / 1.270 / 1.467 — **kein Trend sichtbar**; das Profil nahe L ist fast universell. W1 zurückgestellt. Die Übereinstimmung
p ≈ 1.27 mit 1/α' (K') ist vermutlich Zufall; im W̄-Rahmen (1/α = 1.13) reicht das gemessene Profil (P: Θ konvergiert, 0.22).


## 12. Kantenschwanz (Schnittstelle zu O, Lemma O6)

‖σ₁v₁ auf k < cN‖₂ (kubisch, w01): c = 0.5: 0.124 → 0.193, c = 0.7: 0.357 → 0.623 (N = 400 → 1000, ~N^{0.6}, lokaler Exponent steigt 0.63 → 0.74);
c = 0.8: 0.62 → 2.17, c = 0.85: 0.91 → 9.5, c = 0.9: 1.76 → 88.6 (exponentiell). **Warnung (Vermutung aus §9):** die Stegfrequenz fällt zur Ecke
hin gegen 0 (x ~ ∂_n g ~ τ^{1/α−1}); dann trägt jede Indexschwelle c eine (sehr kleine) Rate λ(τ(c)) > 0, z. B. c = 0.7 ↔ τ ≈ 2° ↔ 0.0014.
„Polynomial für c ≤ 0.7“ wäre dann nur bis N ~ 10³ richtig; asymptotisch nur mit Rate λ_c → 0 für c → 0.

## 13. Nächste Schritte (Stand 00:55; überholt, siehe §16)
1. O(c)-Theorie: A-Bogen-Daten λ_A und Buchtverschiebung selbstkonsistent bestimmen (Leitordnung §9 liefert Frequenzen, aber c um Faktor 2 zu gross).
2. (V1)–(V3) als Beweisziel mit P abstimmen; Hybrid: Rest über K' (P), Eckteil über W̄ (gemessen Θ ≤ 0 nahe L).
3. Kante als d×d-Problem (X: S_N) — WKB-Beschreibung von S_N.


---

# Runde 2 (ab 00:56, Auftrag der Leitung): Entscheidungsrechnung — trägt der (E_L)-Weg?

## 14. Genaue Randdaten und Poisson–Jensen-Kriterium (w24, w25; numerisch)

**Methode.** Kubisches Modell, N = 400…1000. Für 85 Punkte (L; D-Bogen τ = 0.05°…103.38°; A-Bogen |η−L| = 0.002…1.9; Schlitz L(1+u),
u = 1e-3…0.151; Niveaulinien Γ_d, d = 0.05, 0.1) in arb (0.9N+400 Bit, ctx.cap gesetzt): y_η = J⁻¹T_N J∞k_η − T_N k_η, R_N = |y_η|_{ℓ¹_ρ},
Eckteil K_N = |Q_N y_η| (Q_N = Projektion auf die r(N) Eck-Ausgänge u_i mit σ_i > 1.40; r = 2…4) und K₁ = |⟨u₁, y_η⟩|. Asymptotische Rate je Punkt:
Fit ln X = a + κ ln N + λN über N = 500…1000; Fehler = Abweichung der Fits über 500–900 bzw. 600–1000. Poisson–Jensen in Ĉ∖W̄ über die
explizite Abbildung auf die obere Halbebene (Randdaten in s = τ(η) interpoliert, Spiegelpunkte direkt abgebildet), Θ = sup ω_ζ(λ)/g(ζ) auf dem
Schlitz und auf Γ_d (C Lemma 13.2/Satz 13.4, P-I(ii)). Kontrolle: Bernstein–Walsh/Poisson–Jensen an ζ = 3, −3, 2+2i, 5 erfüllt (fast Gleichheit).

**Asymptotisches Profil λ − λ(L) (Eckteil K_N, ± Teilfenster):**

| D-Bogen τ | 0.5° | 1° | 2° | 3° | 4° | 8° | 15° | 30° | 60° | 103.4° |
|---|---|---|---|---|---|---|---|---|---|---|
| λ | +0.0001 | +0.0003 | +0.0009 | +0.0016 | +0.0024 | +0.0047 | +0.0067 | +0.0090 | +0.0108 | +0.0115 ± 0.0003 |

| A-Bogen |η−L| | 0.002 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 | 0.3 | 0.6 | 1.0 | 1.9 |
|---|---|---|---|---|---|---|---|---|---|---|
| λ | −0.0003 | −0.0013 | −0.0020 | −0.0033 | −0.0053 | −0.0084 | −0.0108 | −0.0157 | −0.0195 | −0.0228 ± 0.0004 |

Lokale Exponenten D-Bogen (K_N): 1.54/1.61/1.59/1.54/1.50/1.50/1.55/1.58/1.52/1.34 auf 0.05→0.1→0.2→0.35→0.5→0.75→1→1.5→2→3→4°, danach
konkav. (K₁: 1.48 → 1.10 auf 1°–4°; R_N: 1.6 auf 0.75°–4°.) A-Bogen nahe L ~ −dist^{0.6…0.9}. Schlitz direkt: λ_slit ≈ −1.05…−1.2·g.

**Kriterium Θ (W̄):**

| Grösse | Variante | Schlitz u = 1e-3 / 0.01 / 0.151 | Γ_0.05 | Γ_0.1 | Θ (d = 0.05 / 0.1) |
|---|---|---|---|---|---|
| K_N | endlich N = 1000 | −0.76 / −0.38 / −0.07 | 0.122 | 0.056 | 0.122 / 0.056 |
| K_N | asymptotisch | −0.55 / −0.26 / −0.04 | 0.219 | 0.104 | 0.219 / 0.104 |
| K_N | asymptotisch + Fehler | −0.50 / −0.23 / −0.04 | 0.225 | 0.107 | 0.225 / 0.107 |
| K_N | termweise: λ_A := 0, D + Fehler | +0.063 / +0.049 / +0.021 | 0.231 | 0.113 | 0.231 / 0.113 |
| R_N | asymptotisch (einteilig) | −0.51 / −0.23 / −0.03 | 0.239 | 0.114 | 0.239 / 0.114 |
| R_N | termweise | +0.056 / +0.048 / +0.023 | 0.248 | 0.121 | 0.248 / 0.121 |
| K₁ | asymptotisch | −0.42 / −0.20 / −0.03 | 0.232 | 0.110 | 0.232 / 0.110 |

(Schlitzwerte für u → 0 laufen in der vorzeichenbehafteten Variante gegen −∞, termweise gegen ≈ 0.10 — endlich, weil p·α ≈ 1.5·0.883 > 1.)

**Ergebnis.** Der (E_L)-Weg trägt in der W̄-Ecke mit Reserve: Θ ≈ 0.22 (d = 0.05) bzw. 0.10 (d = 0.1), also ε ≈ 0.78 bzw. 0.90, und zwar
(a) mit den echten zweiseitigen Daten (A-Bogen fällt) und (b) sogar termweise ohne Auslöschung (nur D-Bogen-Profil mit p ≈ 1.5 + Bucht/Sockel auf dem
A-Bogen). Engpass ist der ferne D-Bogen auf Γ_d (Rate 0.0115 gegen λ* ≈ 1.03d), nicht die Ecke; Reserve dort Faktor ≈ 4.5 (d = 0.05) bzw. 10 (d = 0.1).
Ds „Reserve 11“ (∫λdν = 0.086) war als Eck-Grenzwert richtig (termweise hier 0.06–0.10). **Korrektur zu W 00:24:** der Exponent 1.26–1.29 galt für das
Kreisprofil σ|V| zweier Moden (N = 600–1000) und ist nicht robust (mit allen Moden 1.3–1.8, verrauscht); für die massgebliche Grösse (Eckteil der
kopffreien Pol-Antwort) ist p = 1.5–1.6. Vermutung W1 zurückgezogen. Damit ist auch Ps K'-Ecke (Schwelle 1.2704) für diese Grösse nicht grenzwertig.
**(E_A)** ist erfüllt und schärfer als verlangt: R_N fällt auf dem A-Bogen exponentiell (λ = −0.0014/−0.0046/−0.010/−0.019 bei |η−L| = 0.01/0.1/0.3/1.0,
κ ≈ 0…1); R_N(L) ~ N^{1.12}, λ(L) = 0.
**Kneser (w26 + w25b, A_f mit den 50-stelligen Koeffizienten aus gemeinsam/ als exakte Daten, m = eˣ/f_K', N = 500…1000; σ₁ = 5400.7 bei
N = 1000 = D):** praktisch identisch: Eckteil K_N asymptotisch Θ = 0.215 (d = 0.05) / 0.103 (d = 0.1), mit Fehler 0.221/0.106, termweise 0.226/0.122;
R_N einteilig 0.238/0.114; D-Bogen-Exponent nahe L 1.48–1.59 (K_N) bzw. 1.6 (R_N); A-Bogen −0.0013 (0.01) … −0.0222 (Scheitel); Schlitz
λ ≈ −1.1·g. **Status:** numerisch (N ≤ 1000, kubisch und Kneser).


## 15. O(c)-Theorie, Stand 01:20 — ein geschlossenes Leitordnungsproblem (Skizze + numerische Tests)

**Neue exakte Beobachtung (numerisch auf 1e-3 genau, w29–w31).** Der kopffreie Eckanteil T ist ausserhalb D̄ analytisch und (numerisch)
nullstellenfrei; im Halbmond (zwischen Γout und dem A-Bogen) gilt **(1/N)ln|T(t)| = (1/N)ln|y(h(t))|** (T ≈ 𝒴, Urbildterm dominiert), z. B. N = 600
auf der reellen Achse t = 1.05/1.1/1.2/1.25/1.3/1.33ρ: 0.3005/0.2352/0.1166/0.0621/0.0099/−0.0235 gegen 0.3010/0.2357/0.1172/0.0628/0.0105/−0.0228;
ebenso auf Strahlen in den Halbmond hinein bis 0.01 an L heran (Abweichung ≤ 5e-4). T fällt glatt (ohne Knick) durch den A-Bogen; Nulldurchgang
bei 1.31ρ, und h(1.31ρ) = 0.815ρ = gemessener Buchtrand.
**Damit (Skizze):** H := Φ_y auf Ω₁ (Inneres ohne Bucht), := Φ_T auf Ĉ∖D̄, ist EINE harmonische Funktion auf 𝒟' = Ĉ ∖ Bucht̄ (die Bucht = Gebiet
zwischen Γout und B), mit H ~ −log|t| bei 0 und bei ∞, H = 0 auf B, und **H = H∘h auf dem Halbmond** (offene Menge ⇒ Cauchy-Daten auf Γout).
Der Buchtrand B ist frei; konsistent ist automatisch B = h({H = 0}). Keine freien Parameter mehr: c = H(−1), λ = H|_{Γin}.
**Numerische Lösung (w33/w33b, Ladungssimulation, B mit 3 Formparametern, Kleinste Quadrate):**

| | Buchtrand (0°) | Stegfrequenz 180° / 120° / 90° | c = H(−1) | Profil 1° / 4° / 15° / 60° |
|---|---|---|---|---|
| gemessen (kubisch) | 0.815 | 0.957 / 0.936 / 0.88 | 0.0120 | 0.0003 / 0.0024 / 0.0067 / 0.0108 |
| Dipol-Modell (w19, H = 0 auf A-Bogen und h(A-Bogen)) | 0.827 (fest) | 0.956 / 0.939 / 0.879 | 0.0252 | 0.0053 / 0.0103 / 0.0174 / 0.0242 |
| freier Rand, gleichmässiges Gitter (w33) | 0.8196 | 0.958 / 0.939 / 0.898 | 0.0246 | (Ecke nicht aufgelöst) |
| freier Rand, Ecken-gradiert (w33b) | **0.8156** | **0.9565 / 0.935 / 0.867** | 0.0176 | 0.0037 / 0.0071 / 0.0120 / 0.0169 |

Buchtrand und Stegfrequenzen werden vom Modell getroffen; die Rate c und das Profil nahe L nicht (c um Faktor 1.5–2 zu gross, nahe L ~τ^{0.5}
statt τ^{1.5}). Der Engpass ist die **Ecke L** des freien Randwertproblems: lokal wirkt h als komplexe Streckung (log-polar: Verschiebung um
−log μ), und eine harmonische Funktion mit H = H∘h auf einem Keil ist dort eine Summe von Re(A_n ζ^{2nγ}) (ζ = t − L, Re 2γ = 8.89) mit
log-spiralförmigen Nullkurven — ein solches lokales Verhalten ist viel flacher als das gemessene τ^{1.5}. Das Leitordnungsproblem ist also bis auf
die Ecke richtig gestellt; die Ecke verlangt eine eigene (vermutlich von der Innenschicht 1/N beeinflusste) Analyse. **Status: Vermutung/Skizze**
(Herleitung der Gleichungen aus den exakten Beziehungen §5 + Nullstellenfreiheit, die nur numerisch belegt ist).
**Nebenbefund:** die Menge, die den Eckteil trägt, ist keine Zwei-Kreis-Menge (w27: A-Bogen-Daten nahe L verlangen Scheitel b ≈ 1.15–1.2ρ, fern ≈ 1.3ρ).


## 16. Nächste Schritte (Stand 01:29)
1. Lokale Analyse der Ecke L im freien Randwertproblem §15 (Mellin in log ζ mit der komplexen Verschiebung −log μ, Kopplung Keil ↔ Bild-Keil);
   Ziel: Eckexponent des Profils (gemessen 1.5–1.6) und die Konstante, dann c aus dem korrigierten Modell.
2. Verbindung zur Innenschicht (Skala 1/N) — vermutlich trägt sie die Normierung (Eckwert polynomial) und damit c.
3. Beweisrichtung: termweise genügen (V1) [D-Bogen-Profil, p > 1.13, sup < 0.068] und (V3) [Sockel in der Bucht]; beide sind Aussagen über das
   freie Randwertproblem §15 — dessen Lösungstheorie wäre der Weg zu einem Beweis.


---

# Runde 3 (ab 01:30): Satz W1 — (E_K) aus polynomialen Randschranken auf einem vergrösserten Kompaktum (termweise, OHNE Exponentenbedingung p > 1.13)

## 17. Satz W1 (bewiesen; Hypothesen numerisch)

**Bezeichnungen.** f Polynom (kubisches Modell; für Kneser siehe Bemerkung 3), ρ = |L|, J_N invertierbar, (σ_i, u_i, v_i), i ≤ r, die Eckmoden der
ρ-umskalierten Inversen, y_i := σ_i v_i ∈ ℂ^{N+1}, Y_i(x) := Σ_{k≤N} y_{i,k}(ρ/x)^k (= σ_iV_i, D §10.1), T_i(ζ) := ⟨u_i, y_ζ⟩ mit Cs kopffreier
Pol-Antwort y_ζ = Π_Nk_ζ − T_Nk_ζ (Eckteil, C 13.0). W̄ = D̄ ∪ D̄(a₀,R₀); Γin, Γout wie oben; A-Bogen := ∂D(a₀,R₀)∖D̄; h := der Zweig von f⁻¹ mit
h(η) ∈ D für η auf dem A-Bogen. P_ext[Λ] := Poisson-Integral des Äusseren der Kreisscheibe |x| < ρ (Randdaten Λ auf |x| = ρ).

**Hypothesen** (Konstanten C ≥ 1, κ ≥ 0, x₀ ∈ (0, 1], Λ: ∂D → [0, ∞) stetig, Λ = 0 auf Γout):
- **(H-I) Indexschwanz:** |y_{i,k}| ≤ C N^κ für 0 ≤ k < x₀N (und x₀N > deg m).
- **(H-V1) Kreisprofil:** |Y_i(ρe^{iφ})| ≤ C N^κ e^{NΛ(φ)} für alle φ.
- **(H-V3) Bucht:** |Y_i(h(η))| ≤ C N^κ für η auf dem A-Bogen.
- **(G) Geometrie** (endlich prüfbar): mit ℬ := {t : |t| ≥ ρ, P_ext[Λ](t) > x₀ log(|t|/ρ)} („Linse“) und K_T := W̄ ∪ ℬ̄ gilt: (G1) für η ∈ ∂K_T
  liegen alle Urbilder x_j ≠ h(η) von η ausserhalb D ∪ ℬ, für η ∈ ∂ℬ auch h(η); (G2) |x_j| ≥ c₁, |f'(x_j)| ≥ c₂ auf diesen Urbildern; Ĉ∖K_T ist
  zusammenhängend und regulär, 0, a₀ ∈ int K_T.

**Aussage.** Für alle i ≤ r und ζ ∉ K_T:  **|T_i(ζ)| ≤ C₁ N^{κ+2} e^{−(N+2) g_{K_T}(ζ)}**, C₁ = C₁(f, C, c₁, c₂) explizit; und
g_{K_T} ≥ g_W̄ − M_ℬ mit M_ℬ := max_{ℬ̄} g_W̄. Insbesondere gilt (E_K) auf Γ_d mit ε = 1 − M_ℬ/d, auf dem Schlitz mit ε = inf g_{K_T}/g_W̄.

**Beweis.** (1) *Aussenschranke (Lemma W1a).* Y_i = Y_low + Y_high, Y_low := Σ_{k<m}y_{i,k}(ρ/x)^k, m := ⌈x₀N⌉. Für |x| ≥ ρ ist |Y_low| ≤ mCN^κ (H-I).
G(x) := (x/ρ)^m Y_high(x) ist holomorph in |x| > ρ einschliesslich ∞; auf |x| = ρ ist |G| = |Y_high| ≤ |Y_i| + |Y_low| ≤ (m+1)CN^κ e^{NΛ} (H-V1).
log|G| ist subharmonisch und bei ∞ beschränkt ⇒ log|G(x)| ≤ log((m+1)CN^κ) + N·P_ext[Λ](x). Also |Y_i(x)| ≤ mCN^κ + (m+1)CN^κ·
exp(N[P_ext[Λ](x) − x₀ log(|x|/ρ)]) ≤ (2m+1)CN^κ für |x| ≥ ρ, x ∉ ℬ.
(2) *Eckformel (D6) für T_i:* T_i(ζ) = Y_i(ζ)m(ζ)/ζ + Σ_j Y_i(x_j)/(x_jf'(x_j)) + Σ_{k<deg m} y_{i,k}ρ^k[Q_ζ]_k − ⟨u_i, T_Nk_ζ⟩, x_j = Urbilder von ζ,
Q_ζ = (m(x) − m(ζ))/(ζ − x). Für ζ ∈ ∂K_T: erster Term nach (1) (ζ ∉ ℬ, |ζ| ≥ ρ; auf Γin ∩ ∂K_T ist Λ = 0 und (H-V1) direkt); Urbildterme nach (G1)+(1)
bzw. (H-V3) für h(η) auf dem A-Bogen; dritter Term nach (H-I) (|[Q_ζ]_k| ≤ c(f) auf dem kompakten ∂K_T); Kopfterm ≤ ‖T_Nk_ζ‖_{ℓ²_ρ} ≤ √(N+1)/|ζ|.
Also max_{∂K_T}|T_i| ≤ C₁N^{κ+2}.
(3) T_i ist rational in ζ mit Polen nur in 0, a₀ ∈ int K_T und T_i(ζ) = O(ζ^{−N−2}) (C 12.1/13.0) ⇒ Bernstein–Walsh (C Satz 12.2) in Ĉ∖K_T.
(4) g_W̄ − g_{K_T} ist harmonisch und beschränkt auf Ĉ∖K_T und ≤ M_ℬ auf ∂K_T ⇒ ≤ M_ℬ überall. ∎

**Bemerkungen.** (1) Es wird KEINE Exponentenbedingung an Λ nahe L gebraucht, nur dass ℬ den Schlitz nicht trifft und K_T in L eine Ecke mit
positivem Aussenwinkel behält: für Λ ≤ Aτ^p mit p > 1 ist die Linsenbreite O(τ^p) und g_{K_T}/g_W̄ → 1 am Schlitz (Grenzwert (1/κ₀)∫ψs^{−2}ds mit
ψ ~ s^{1−α+αp}, endlich für p > 1); auch p = 1 mit kleinem A genügt (dann nur ein etwas kleinerer Eckwinkel). Damit ist Ps grenzwertige Frage
„p > 1.13 bzw. 1.27?“ für den termweisen Weg erledigt. (2) Ps Lücken (1) und (2) (A-Bogen-Urbilder, äusseres Poisson-Mittel) verschwinden: (1) ist
(H-V3), (2) ist Lemma W1a (Zerlegung Y_low/Y_high statt Poisson-Mittel der vollen Funktion). (3) Kneser: D6 durch O6 + (E_sing) ersetzen (P-IV);
sonst wörtlich. (4) Für alle r(N) Moden gemeinsam: |Q_Ny_ζ| ≤ Σ_i|T_i| kostet einen Faktor r ≤ N.

## 18. Numerische Prüfung der Hypothesen (kubisch, N = 400…1000, w35–w37; Kneser-Daten w26 zeigen dasselbe Profil)

- **(H-V1):** asymptotisches Kreisprofil (alle Eckmoden σ_i > 1.40, Fit a + κ ln N + λN, + Fehler als Umhüllende): Λ = 0.0009 (2°), 0.0025 (4°),
  0.0049 (8°), 0.0069 (15°), 0.0091 (30°), 0.0109 (60°), **λ_max = 0.0116** (180°); nahe L (≤ 1°) 0 ± 0.0005 (Innenschicht ~1/N, Kollaps in τN, w34);
  Exponent 1°–4.5°: 1.42. Γout: λ = 0.
- **(H-I):** max_{k<xN}|y_k| = 0.031 → 0.035 (x = 0.5), 0.081 → 0.113 (x = 0.7, ~N^{0.48}), 0.18 → 0.63 (x = 0.8, ~N^{1.8}) für N = 500 → 1000 ⇒ x₀ = 0.7.
- **Linse (x₀ = 0.7):** Breite s = λ/x₀-artig: 0.0018 (2°), 0.0037 (4°), 0.0097 (15°), 0.0164 (180°); **M_ℬ = 0.0156**; C₁ := sup_{Γout} ∂_s P_ext[Λ] = 0.20 < x₀
  (Γout-Seite frei); Schlitz ausserhalb ℬ (u = 1e-3: P_ext = 1.6e-4 < 7.0e-4).
- **(G1):** Urbilder der Linsenrandpunkte: |h|/ρ − 1 = 0.0028 (0.25°) … 0.040 (4°) … 1.45 (180°) gegen Linsenbreite dort 0.0006 … 0.0023 … 0.015;
  weitere Urbilder bei |x| = 2.4–5.5ρ.
- **(H-V3):** |Y(h(η))| auf dem A-Bogen: nahe L polynomial (ln = 1.11 → 1.71, ~N^{0.66} bei |η−L| = 0.005), weiter weg fallend bis auf den
  Sockel (ln −15 … −19, N ≥ 800).
- **Direkte Probe der Aussage:** T auf ∂K_T (genau): max ln|T| auf dem Linsenrand 1.30 → 1.79, in L 1.52 → 2.18, auf dem A-Bogen ≤ 1.44 (N = 400 → 1000):
  polynomial (Fits λ ≤ 3e-4, auf dem A-Bogen negativ).
- **ε (w37, Ladungssimulation für g_{K_T}):** Schlitz g_{K_T}/g_W̄ = 0.885 (u = 0.003), 0.92 (0.01), 0.94 (0.03), 0.97 (0.15); **Γ_0.05: ε_Γ = 0.694, Γ_0.1: 0.851**.
**Ergebnis:** Satz W1 liefert (E_K) für das kubische Modell mit ε ≥ 0.69 (d = 0.05) bzw. 0.85 (d = 0.1) aus Hypothesen, die bis N = 1000 erfüllt sind.
Reserve am fernen D-Bogen: ε_Γ > 0 ⇔ M_ℬ < d; M_ℬ ≈ 1.35·λ_max (x₀ = 0.7) ⇒ trägt, solange λ_max < 0.037 (d = 0.05) bzw. 0.074 (d = 0.1): Faktor 3.2 bzw. 6.4
über dem gemessenen λ_max = 0.0116.
**Offen (für alle N):** (H-I), (H-V1), (H-V3) — Aussagen über die adjungierten Eckmoden; Status numerisch (N ≤ 1000).


---

# Runde 4 (ab 02:00): lokale Ecktheorie bei L und Herleitung von c

## 19. Lokale Ecktheorie bei L: ein Kegel mit Deckverschiebung 2πi − log μ (Skizze; numerisch geprüft, w38/w39)

**Herleitung (Skizze, aus dem Leitordnungsproblem §15).** Nahe L ist h(L + ζ) = L + ζ/μ + O(ζ²). Schneidet man den Keil zwischen dem Buchtrand B
und der Nullkurve Z_T von T im Halbmond aus (Öffnung arg μ, denn B = h(Z_T)) und klebt B an Z_T mit h, so entsteht ein Kegel. In s = log ζ ist
die Deckverschiebung **2πi − log μ** (eine volle Drehung um L minus die komplexe Streckung). Einwertige Koordinate w = ζ^{κ₁},
**κ₁ = 2πi/(2πi − log μ)**. Beschränkte harmonische Funktionen mit H(L) = 0 sind dort H = Re Σ_{n≥1} c_n w^n; jede Potenz ζ^{nκ₁} erfüllt
H = H∘h samt Normalableitung auf Γout exakt (für das linearisierte h). Folgerungen:
- (a) Lokaler Exponent Re κ₁ = 2π(2π − β)/((2π − β)² + ℓ²) mit β = arg μ, ℓ = ln|μ|; für μ = √L: **κ₁ = 1.1182 − 0.0317i**.
- (b) **Re κ₁ > 1 ⇔ β(2π − β) > ℓ²**; das gilt für 0 < β < π und kleines |ℓ| immer (hier 3.75 gegen 0.025). Das ist genau die Bedingung p > 1 aus
  Satz W1 (§17, Bem. 1), jetzt strukturell begründet (Skizze).
- (c) Nullkurven durch L: die Naht (Z_T aussen, B = h(Z_T) innen) und eine zweite im Winkelabstand π/Re κ₁ = 161.0° (Log-Spiralen mit
  Steigung Im κ₁/Re κ₁).
- (d) Vergleich mit dem Keilmodell §15 (Re 2γ = 8.89): das ist Os Mellin-Bild (Nullstellen von 1 + μ^{s−1}: s = 1 + (2k+1)iπ/ln μ,
  Re s − 1 = 4.44·(2k+1)). Der Keil sieht nur den Halbmond, nicht die volle Drehung um L; in 2D läuft die Ratenfunktion aber um L herum (Φ_y innen,
  Φ_T aussen, harmonisch über Γin), daher die 2πi. Formal: ersetzt man in Os Multiplikator ln μ durch ln μ − 2πi (eine Drehung um L), so liegen die
  Nullstellen bei **s = 1 ± κ₁/2** (Re = 1.559 bzw. 0.441). 1.559 liegt auffällig nahe am Eckexponenten p ≈ 1.5–1.6 von K_N (D, §14). **Vermutung**,
  nicht erklärt.

**Numerische Prüfung (w38: genaue arb-Werte ln|T| aussen, ln|Y| innen auf Kreisen um L; w39: Raten aus a + κ ln N + λN, N = 400…1000,
Fehler ≤ 2e-4).** ψ = Winkel relativ zur äusseren Normalen in L (+90° = Γin-Tangente, −90° = Γout-Tangente, A-Bogen-Tangente −68.9°).

| r = \|t − L\| | H(ψ = 0) | H(ψ = 180°) | Nullstellen ψ | Abstand |
|---|---|---|---|---|
| 0.03 | −0.01307 | +0.01310 | −79.6°, +89.0° | 168.6° |
| 0.06 | −0.02744 | +0.02974 | −77.2°, +89.1° | 166.3° |
| 0.12 | −0.05756 | +0.06708 | −76.7°, +90.0° | 166.7° |

- Fit H = Σ_{n≤2} Re(c_n w^n) mit festem κ₁: rms 8e-4 (Werte bis 0.067). **κ frei: κ = 1.114 − 0.018i** (Theorie 1.118 − 0.032i).
  Radialexponenten je Strahl 0.99–1.19 (Mittel ≈ 1.1). Nahtabstand 166–169° statt 161° (die n = 2-Korrektur erklärt nur ~1–2°).
- Naht Z_T bei ψ = −77.8° (im Halbmond zwischen Γout und A-Bogen); B = h(Z_T) bei ψ ≈ −116°. Probe: bei ψ = −105° (in der Bucht) fällt |Y|
  exponentiell mit Rate −0.0033/−0.0058 (r = 0.03/0.06) = harmonische Fortsetzung von H über B (vorhergesagt −0.0026/−0.0057), danach Sockel.
- **Rolle der Innenschicht (w34 neu ausgewertet, numerisch + Skizze):** Kollaps in ξ = τN (τ in Grad) mit β = 1 nur für ξ ≲ 100, dort
  ln(Λ_N(τ)/Λ_N(0)) ≈ 0.0018ξ (N = 400…1000), also τ ≲ 0.1–0.25°. In inneren Variablen ist der äussere Eckterm N·H = N^{1−κ₁}|c₁|(Nr)^{κ₁}·(Winkel)
  = O(N^{−0.118}) bei festem ξ: die Innenschicht trägt **keinen exponentiellen Anteil**, nur den polynomialen Eckvorfaktor N^{κ(0)} (κ(0) ≈ 0.5–0.6).
  **c wird NICHT von der Innenschicht bestimmt** (meine Vermutung §16.2 ist damit überholt), sondern vom äusseren Problem — siehe §20.

## 20. Freies Randwertproblem MIT Eckfunktionen: c = 0.0118 hergeleitet (Skizze des Modells + numerisch, w40)

**Verfahren (w40_freier_rand_ecke.py).** Wie w33b (H = −log|t| + Σq log|t − s| + C, H = 0 auf B, Cauchy-Daten H = H∘h auf Γout, Σq = 0), dazu
(a) **Eckfunktionen** Re(c_n Z^{nκ₁}), n ≤ 3, Z = e^{iγ}(t − L̃)/(t − P), Schnitt = Kreisbogen von L̃ nach P = 0.913 IN der Bucht (geprüft),
gespiegelt bei L̄; (b) **Eckwinkel der Bucht frei** (B nicht mehr an h(A-Bogen) geheftet). Eingaben nur f (h auf Γout), L, μ — **kein Fit an Daten**.

| Variante | Rest | c = H(−1) | Bucht (0°) | Profil 1° / 2° / 4° / 8° / 15° / 30° / 60° / 103° |
|---|---|---|---|---|
| w33b (ohne Eckfunktionen, Ecke geheftet) | — | 0.0176 | 0.8156 | 0.0037 / — / 0.0071 / — / 0.0120 / — / 0.0169 / — |
| ohne Eckfunktionen, Eckwinkel frei, K = 8 | 5.6e-3 | 0.0153 | 0.814 | 0.0027 / 0.0039 / 0.0055 / 0.0076 / 0.0099 / 0.0125 / 0.0146 / 0.0153 |
| mit Eckfunktionen, K = 4 / 6 / 8 / 10 / 12 / 16 | 1.4e-3 … 3.3e-4 | 0.0131 / 0.0126 / 0.0123 / 0.0121 / 0.0119 / 0.0119 | 0.78–0.88 | |
| **mit Eckfunktionen, K = 20** | **2.6e-4** | **0.01184** | **0.814** | **0.0012 / 0.0020 / 0.0032 / 0.0050 / 0.0070 / 0.0093 / 0.0112 / 0.0118** |
| gemessen (kubisch, w35/w08) | | 0.0117 (λ(103°) − λ(L)); σ-Rate 0.0120 | 0.815 | 0.0003 / 0.0009 / 0.0025 / 0.0049 / 0.0069 / 0.0091 / 0.0109 / 0.0116 |

- Robust: Schnittrichtung −100°/−104° liefert dasselbe (0.01184); verlässt der Schnitt die Bucht (−108°), scheitert die Lösung (Rest 7.7e-3).
  Konvergenz in K: 0.01193 → 0.01187 → 0.01184 (K = 12/16/20); c_Modell = 0.0118 ± 0.0002 (korrigiert: ± 0.0004, Kondition geprüft, §22.7).
- Lokal (r = 0.12): H(0) = −0.0572 / H(180°) = +0.0674 / Nullstellen −76.1°, +89.4° gegen gemessen −0.0576 / +0.0671 / −76.7°, +90.0°;
  r = 0.06 auf 2 %, r = 0.03 auf 5 % (Nullstellen −76.0°, +85.3° gegen −79.6°, +89.0°).
- **Ergebnis:** Das Leitordnungsproblem §15 mit der Eckstruktur §19 trifft **c (0.0118 gegen 0.0117/0.0120)**, das Profil ab 8° (auf ≤ 0.0003),
  den Buchtrand und die lokalen Werte um L für r ≥ 0.06 — ohne freie Parameter. Der Unterschied zu §15 (0.0176) lag allein an der Ecke.
- **Nicht hergeleitet: das Profil nahe L (τ ≤ 4°).** Dort liegen die Daten UNTER dem Modell (der zweite Knoten liegt in den Daten fast genau auf der
  Γin-Tangente, im Modell 4° davor). Ursache offen: Modell für r < 0.03 unteraufgelöst (Restmaximum in der Ecke) und/oder Daten vor-asymptotisch
  (r = 0.03 ⇔ rN ≤ 30 bei N ≤ 1000). Das gemessene „τ^{1.5}“ auf 1°–4° ist nach §19 ein Übergangsexponent; asymptotisch (τ → 0) sagt die Theorie
  Re κ₁ = 1.118 voraus (2.118, falls der Knoten exakt tangential an Γin läge). Für Satz W1 ist das gleichgültig (p > 1 genügt, und das Modell liegt
  nahe L ÜBER den Daten, ist also eine verträgliche obere Hülle).
- Reserve am fernen D-Bogen: mit dem Modellwert λ_max = c = 0.0118 ist M_ℬ ≈ 1.35·0.0118 = 0.016 < d = 0.05 (ε_Γ ≈ 0.68) — die Entscheidung §14/§18
  hängt also nicht an Messfeinheiten, sondern ist eine Eigenschaft des Leitordnungsproblems.

## 21. Stand und nächste Schritte (02:40)
1. Prüfung durch P: §19 (Kegel, κ₁, Bedingung (b)) und §20 (Modell = Leitordnung §15 + Eckfunktionen; Zahlen).
2. Beweisweg für (H-V1)/(H-V3) (Plan, Skizze): subharmonisches Kleben auf einer festen Vergleichsfläche ℳ'. Aus der exakten Identität §5 folgt im
   Halbmond ausserhalb der Linse T = Y∘h/(h·F'∘h) + O(N^κ) (die übrigen Terme sind nach Lemma W1a polynomial); u_T := (1/N)ln max(|T|, N^K) und
   u_Y∘h stimmen dort bis ε_N = O(ln N/N) überein und lassen sich mit einer harmonischen Kragenfunktion zu einer subharmonischen Funktion auf ℳ'
   kleben ⇒ u ≤ H_ℳ' + const. Offen: (i) die Bucht-Schranke (V3) als Eingang ist zirkulär — braucht ein eigenes Argument (Tropik der vier Terme
   m̃y/t, 𝒴, T, 𝒰 in der Bucht); (ii) die Konstante (Normierung über den polynomialen Eckwert); (iii) eine rechnergestützte Kontrolle von H_ℳ'
   (Ladungssimulation mit Fehlerschranke, z. B. über das Maximumprinzip für den Kollokationsrest).
3. Nahes Profil (τ ≤ 4°): feinere Eckauflösung im Modell (mehr Eckterme mit Regularisierung, Ladungen gradiert) gegen Daten bei grösserem N.


---

# Runde 5 (ab 02:50): Beweisweg für (H-V1)/(H-V3) für alle N — Befund, Vereinfachung, Programm

## 22. Die zirkuläre Stelle ist intrinsisch; eine einzige Hypothese genügt

### 22.1 Lemma W3 (bewiesen): Bucht- und Γout-Schranke sind Aussagen über T

Bezeichnungen wie §17; D6 für die Eckmode i (Index weggelassen):
T(η) = Y(η)m(η)/η + Σ_j Y(x_j)/(x_jf'(x_j)) + Σ_{k<deg m} y_kρ^k[Q_η]_k − ⟨u, T_Nk_η⟩, x_j = Urbilder von η.
**(a)** Für η auf dem A-Bogen ist h(η) ∈ Bucht eines der Urbilder; also exakt
  Y(h(η)) = h(η)f'(h(η))·[T(η) − Y(η)m(η)/η − Σ_{x_j ≠ h(η)} Y(x_j)/(x_jf'(x_j)) − Σ_k y_kρ^k[Q_η]_k + ⟨u, T_Nk_η⟩].
  Sind Y auf dem A-Bogen und an den fernen Urbildern sowie y_k (k < deg m) polynomial beschränkt (das liefert Lemma W1a aus (H-I), (H-V1)),
  so gilt **(H-V3) ⇔ |T| ≤ CN^κ auf dem A-Bogen**.
**(b)** Ebenso für t ∈ Γout mit η = f(t) ∈ f(Γout) (t ist dort das Urbild h(η)): **Λ = 0 auf Γout ⇔ |T| ≤ CN^κ auf f(Γout)**.
*Beweis:* D6 nach dem jeweiligen Urbildterm aufgelöst; |h f'(h)| und die übrigen Faktoren sind auf den kompakten Bögen beschränkt. ∎
**Folgerung.** Die „zirkuläre Stelle“ ist keine Schwäche einer Technik: (H-V3) ist bis auf polynomiale Faktoren DIESELBE Aussage wie die
Schlussfolgerung von Satz W1 auf dem A-Bogen-Teil von ∂K_T. Inhaltlich sagt sie: die Nullkurve Z_T von T trennt Γout vom A-Bogen und von f(Γout)
— eine Aussage über die LAGE des freien Randes. T auf dem A-Bogen hängt (Maximumprinzip im Halbmond) von T auf Γout ab, und dort ist T ≈ Y∘h mit
h(Γout) ⊂ Ω₁ (Innenbereich, exponentiell gross). Keine Umformung, die nur Aussenwerte benutzt (D6, Konturdarstellung unten, Kleben auf einer festen
Fläche), kann sie liefern; man braucht eine globale Vergleichsaussage (§22.6, §22.8).
**Nebenbei (bewiesen): Konturdarstellung des Eckanteils.** Mit dem Abschneiderest ρ_ζ = T_N[(f/ζ)^{N+1}/(1 − f/ζ)] ist
  **T(ζ) = ⟨y, ρ_ζ⟩ = (1/2πi)∮_γ y(x)·(f(x)/ζ)^{N+1}/(1 − f(x)/ζ)·dx/x**
für JEDE Kontur γ, die 0 umläuft, aber kein Urbild von ζ (Integrand analytisch bis auf 0 und f⁻¹(ζ)). D6 ist der Grenzfall γ → ∞ (Residuen an den
Urbildern und bei ∞, dort mit Jᵀy = u). Damit (1/N)ln|T(ζ)| ≤ max_γ[(1/N)ln|y| + ln|f|] − ln|ζ| + o(1) — gibt allein nur Kreisschranken (|ζ| > max_γ|f| ≈ 1.3), nicht die Ecke.

### 22.2 Satz W1'' (bewiesen): EINE Hypothese auf EINER Kurve statt (H-I), (H-V1), (H-V3)

Sei ℬ irgendeine Linse längs Γin wie in §17 (Rand ∂ℬ ausserhalb D̄, von L über −ρ(1+s) nach L̄, Breite → 0 in L, L̄) und
**γ* := ∂ℬ ∪ h(A-Bogen) ∪ Spiegelbild** (geschlossene Jordankurve durch L, L̄, umläuft 0 und ã₀; auf der Γout-Seite IN der Bucht).
**(H*)** |Y_i(x)| ≤ CN^κ für x ∈ γ*.   **(G*)** = (G1), (G2) aus §17 (Urbilder von ∂ℬ ausserhalb von ℬ ∪ D̄; |x_j|, |f'(x_j)| nach unten beschränkt).
**Aussage:** wie Satz W1: |T_i(ζ)| ≤ C₁N^{κ+2}e^{−(N+2)g_{K_T}(ζ)} ausserhalb K_T = W̄ ∪ ℬ̄.
*Beweis.* Y_i ist analytisch in Ext(γ*) ∪ {∞} (Polynom in 1/x) ⇒ |Y_i| ≤ CN^κ auf ganz Ext(γ*) (Maximumprinzip). Cauchy auf |x| = 2:
|y_{i,k}|ρ^k ≤ 2^k CN^κ für k < deg m. Für ζ ∈ ∂K_T liegen ζ, alle fernen Urbilder und (für ζ ∈ ∂ℬ nach (G1)) h(ζ) in Ext(γ*); für ζ auf dem
A-Bogen ist h(ζ) ∈ γ*. Also sind alle D6-Terme ≤ C(f)N^{κ+1}; Kopfterm ≤ √(N+1)/|ζ|. Bernstein–Walsh wie in §17. ∎
**Gewinn:** (H-I) — die unhandliche Phasenraum-Aussage über die niedrigen Indizes — und das Kreisprofil (H-V1) entfallen ganz; die Linse ist
ein freier Konstruktionsparameter (nur (G*) und der Schlitz zählen). **Numerik** (genau, arb, N = 400…1000): auf h(A-Bogen) polynomial nahe L
(κ = 0.59), weiter weg fallend (w36); auf ∂ℬ der Mess-Linse (x₀ = 0.7) polynomial nahe L (κ ≤ 0.73, λ ≤ 2e-4) und fallend weiter weg (λ bis −0.0036).
**Linse aus dem MODELL (w41, ohne Daten definiert):** Λ := Γin-Profil des Leitordnungsmodells (w40), ℬ_x = {P_ext[Λ] > x·ln(|t|/ρ)}:
x = 0.7: M_ℬ = 0.0159, ε_Γ ≥ 0.682 (d = 0.05) / 0.841 (d = 0.1); x = 0.9: M_ℬ = 0.0124, ε_Γ ≥ 0.752 / 0.876; Schlitz ausserhalb der Linse.
(H*) auf dem Rand der Modell-Linse (w41, genau, N = 400…1000, alle Eckmoden σ_i > 1.40): **x = 0.7: alle Raten ≤ +0.0002** (nahe L
polynomial, κ ≤ 0.75; ab 2° fallend bis −0.004); x = 0.9: bei 8°–30° Raten +0.0003…+0.0006 (Linse dort zu dünn) ⇒ x = 0.7. (w41 benutzte noch die
schlecht konditionierte §20-Lösung; deren Γin-Profil weicht vom Referenzmodell §22.7 um ≤ 2e-4 ab.)

### 22.3 Warum Kleben auf einer festen Vergleichsfläche allein nicht schliesst (Analyse, Plan §21 geprüft)

Auf ℳ₀ = (Ĉ∖D̄) ∪ D' mit Γout ∼ h(Γout) (feste Fläche, kein freier Rand) wäre U := (1/N)ln|T| aussen, (1/N)ln|y| innen subharmonisch und über
Γin, Γout bis auf die Restterme der Identität (m̃y/t auf Γout, 𝒴 auf Γin, 𝒰, P) verklebt. Drei Hindernisse:
(i) Die Restterme enthalten y auf Γout und in der Bucht — genau die Grössen aus (H-V1)|Γout und (H-V3) (Lemma W3).
(ii) Kleben mit Böden (max(|·|, Boden)) braucht Bodenfunktionen, die dort aus- und dort eingeschaltet werden, wo U den Boden kreuzt — das ist der
     freie Rand B ∪ Z_T selbst.
(iii) Auf der kompakten Fläche ist eine subharmonische Funktion mit festen Polen bei 0, ∞ bis auf eine KONSTANTE bestimmt; die Konstante (= die
     Normierung c bzw. H(L) = 0) kommt aus der Quelle u, die nur auf Bodenniveau sichtbar ist (längs B).
Also: das Kleben liefert die Gestalt (H bis auf eine Konstante), sobald der freie Rand bekannt ist; Lage des Randes und Konstante sind global.

### 22.4 Bootstrap in N (Option der Leitung) — geprüft, ungeeignet für (V3)

Rändern: J̃_{N+1}ᵀy' = u' ⇒ ȳ' := y'|_{≤N} = J̃_N^{-T}(ū' − ηc) (η = y'_{N+1}, c = Zeile N+1 von J̃_∞). Der Eckanteil von ȳ' ist parallel zu den
Eckmoden von N (das vererbt die Gestalt), aber der Volumenanteil Σ_{i>r}σ_iv_i⟨u_i, ·⟩ ist IN der Bucht exponentiell gross (|v_i(t)| bis
(ρ/|t|)^N), ebenso der Randterm η t^{−N−1}. Die polynomiale Grösse von y' in der Bucht ist eine Auslöschung zwischen diesen beiden (stationäre Phase
im Index: der Steg schwingt mit der Frequenz von 180°, die Bucht liegt bei 0°). Eine Koeffizienten-Induktion sieht diese Auslöschung nicht.

### 22.5 Xs Form (Lemma 3.1) — dieselbe Aussage, homogen formuliert (bewiesen)

Sei R = [J_N; B] (σ_min(R) ≥ c₀ nach O + X5) und Q eine Orthonormalbasis von ker Rᵀ (dim d), q = (q_top, q_bot). Dann ist (I − RR⁺) = QQᵀ,
S = Q_botQ_botᵀ und der kanonische Eckanteil
  **c_N(k_ζ) = R⁺_bot (Q_botᵀ)⁻¹ [⟨q_top, ρ_ζ⟩ + ⟨q_bot, (J_∞ tail_ζ)_over⟩]_{q ∈ Q}.**
Die Spalten q sind NORMIERTE annullierende Funktionale (Rᵀq = 0: J_Nᵀq_top = −Bᵀq_bot) — homogenes Problem ohne Quelle; ‖q_bot‖ ~ e^{−cN},
q_top ≈ normierte Eckmode. Also: (Ẽ) ⇔ T-Schranken für die normierten Funktionale × ‖(Q_botᵀ)⁻¹‖ ~ e^{cN}. Xs Θ^c (Satz X2) ist eine
Randdaten-Bedingung derselben Art, keine a-priori-Schranke.

### 22.6 Programm für alle N (drei Skalen) — Stand

**(A) Aussen: Vergleichsprinzip mit dem Modell als Barriere.** Zu zeigen: limsup (1/N)ln|y_N| ≤ H_Modell (+ o(1)) auf Ext(γ*) — damit (H*).
Nach Satz W5 (§22.8) ist das GLEICHWERTIG mit der Nullstellenverteilung (A1) samt Normierung (A2), d. h. mit asymptotischer Extremalität.
Werkzeug: Gonchar–Rakhmanov–Stahl (n-te-Wurzel-Asymptotik nicht-hermitescher Orthogonalität). Entsprechungen: Orthogonalitätsrelationen
∮y(x)[F̃(x)^j + x^j m̃(x)]dx/x = 2πi u_j (j ≤ N), Rest T = O(t^{−N−2}) (Interpolation bei ∞), S-Eigenschaft ↔ Cauchy-Daten-Bedingung H = H∘h
mit gleicher Normalableitung auf Γout (im Modell eingebaut), Riemannsche Fläche ℳ vom Geschlecht 0 (H = −ln|φ| + const). Neu gegenüber Stahl:
die Komposition (Spur über die Äste von F̃⁻¹) und die Quelle u. **Offen — das ist der eine fehlende globale Satz.**
**(B) Ecke:** auf dem Kegel (§19) sind beschränkte harmonische Funktionen mit H(L) = 0 genau Re Σc_nw^n ⇒ H = O(r^{Re κ₁}), Re κ₁ > 1 (bewiesen
für das Modell); für y Teil von (A).
**(C) Innenschicht (Skala 1/N):** Amplitudenebene; Vorfaktor N^{κ(0)} (κ(0) ≈ 0.51) — hier gehört Os T_L/Mellin hin (BRETT 02:59 an O).
**Agmon-Variante (mit O):** Konjugation mit analytischem Gewicht g = e^{Nφ} erhält die Operatorform: g⁻¹J(gδ) = e^{N(φ∘f−φ)}δ∘f + mδ. Eine
(FA)-Abschätzung für diese Familie (Os Konturverschiebung §30) gäbe gewichtete Schranken für alle N — Hindernis: Abschneidung T_N vertauscht nur
mit radialen Gewichten.
**Ehrlicher Stand:** siehe §22.9.

### 22.7 Korrektur zum Modell §20: Kondition der Ladungssimulation geprüft (numerisch)
Die Lösung K = 20 aus §20 war schlecht konditioniert: Ladungen bis |q| = 1.1·10⁴, und ZWISCHEN den Kollokationspunkten auf B Fehler bis 1.5
(nahe 0°, wo der Schnittendpunkt P und die Ladungen dicht an B liegen). Neu mit Prüfung zwischen den Kollokationspunkten (w40, Argumente QF/QS/BS):

| Diskretisierung | max\|q\| | Rest (Kollok.) | Fehler zwischen Punkten: B / Γout-Wert / Γout-Ableitung | c | Profil 1° / 4° / 15° |
|---|---|---|---|---|---|
| §20 (Ladungen 0.2/0.5/0.8, jede 2.) | 1.1e4 | 2.6e-4 | 1.5 / 1.8e-3 / 0.32 | 0.01184 | 0.0012 / 0.0032 / 0.0070 |
| Ladungen 0.4/0.6, jede 2., B dicht | 10.7 | 4.9e-3 | 3.8e-4 / 6.7e-4 / 3.6e-3 | **0.01176** | 0.0010 / 0.0030 / 0.0069 |
| Ladungen 0.35/0.5/0.65, jede 3., B dicht | 3.0 | 5.9e-3 | 6.1e-4 / 6.7e-4 / 4.3e-3 | **0.01181** | 0.0010 / 0.0031 / 0.0069 |
| zu grob (0.4/0.6, jede 3.) | 2.4–2.7 | 2.9e-2 … 7.5e-2 | — | 0.0110 / 0.0099 (verworfen) | |

**Ergebnis bleibt:** gut konditionierte UND aufgelöste Varianten geben c = 0.0118; Fehlerbalken ehrlich **±0.0004** (statt ±0.0002 in §20/BRETT 02:42).
Das Profil nahe L bleibt in allen brauchbaren Lösungen über den Daten. Referenzmodell jetzt w40_modell_ref.npz (Ladungen 0.4/0.6).

### 22.8 Das Modell als Green-Funktion — und Satz W5: von der Nullstellenverteilung zur Schranke (H*)
**Beobachtung (numerisch, Referenzmodell):** Sei Ω_y das Gebiet um 0, das vom Buchtrand B und der Nulllinie von H ausserhalb Γin (Linsenrand)
begrenzt wird. Dann ist H|_{Ω_y} = g_{Ω_y}(·, 0) (Green-Funktion mit Pol in 0): H = 0 auf ∂Ω_y, H ~ −ln|t| bei 0. Der Fluss von ∇H durch ∂Ω_y ist
2π·0.9999; das Potential −ln|t| + ∫ln|t − z|dν(z) der harmonischen Massverteilung ν (= Mass |∇H|ds/2π auf ∂Ω_y) stimmt mit H auf Ω_y und mit 0
aussen auf 1e-4 überein. Vorhersage: **alle N Nullstellen der Eckmode liegen asymptotisch auf ∂Ω_y, Anteil auf B = 0.4728**, und
**(1/N)ln|y_N| (höchster Koeffizient) → ∫ln|z|dν(z) = −0.0650**.
**Satz W5 (bewiesen, Standard-Potentialtheorie).** Sei P_N(w) := Σ_k y_k w^k (w = 1/t), ν̃_N sein normiertes Nullstellenmass, ν̃ das Bild von ν.
Gelten (A1) ν̃_N → ν̃ schwach-* mit Trägern in einem festen Kompaktum, und (A2) limsup (1/N)ln max_{|t|=2}|y_N(t)| ≤ 0 (dort ist G = 0), so ist für jedes
Kompaktum E: limsup_N max_E[(1/N)ln|y_N| − G] ≤ 0 mit G := g_{Ω_y}(·, 0) auf Ω_y, 0 sonst. Insbesondere |y_N| ≤ e^{o(N)} auf γ* (⊂ Ĉ∖Ω_y) und mit
Satz W1'' |T_i(ζ)| ≤ e^{o(N)}e^{−(N+2)g_{K_T}(ζ)} ausserhalb K_T.
*Beweis.* (1/N)ln|P_N(w)| = (1/N)ln|lc(P_N)| + ∫ln|w − ω|dν̃_N. Abstiegsprinzip (Saff–Totik I.6.8): limsup ∫ln|w_N − ω|dν̃_N ≤ ∫ln|w − ω|dν̃ für w_N → w,
also gleichmässig auf Kompakta. Mittel über |w| = 1/2 (Jensen: ∫ln|w − ω|dθ/2π = max(ln½, ln|ω|), keine Singularität, stetig in ω) liefert mit (A2)
limsup (1/N)ln|lc| ≤ −∫max(ln½, ln|ω|)dν̃(ω). Einsetzen und die Potentialidentität oben (G = 0 auf |t| = 2). ∎
**Umkehrung (Blatt–Saff–Simkani):** n-te-Wurzel-Extremalität von P_N auf K̃ := 1/(Ĉ∖Ω_y) ⇒ (A1). Also: **(H*) auf n-te-Wurzel-Ebene ⇔ (A1)+(A2)
⇔ „die Eckmode ist als Polynom in 1/t asymptotisch extremal auf K̃“.** Das ist der eine fehlende globale Satz, jetzt in Standardform.
**Test (w42/w43: Nullstellen von P_N mit arb/flint, alle N; Vergleich mit dem Referenzmodell, kein Fit):**

| N | Anteil innen (Modell 0.4728) | Abstand zu B (Median) | Abstand zur Linsen-Nulllinie (Median) | Winkel-KS innen / aussen | Ausreisser \|t\| > 1.1 |
|---|---|---|---|---|---|
| 400 | 0.4650 | 0.0209 | 0.0426 | 0.009 / 0.007 | 4 |
| 600 | 0.4683 | 0.0149 | 0.0289 | 0.007 / 0.003 | 4 |
| 800 | 0.4688 | 0.0115 | 0.0230 | 0.006 / 0.004 | 4 |
| 1000 | **0.4700** | **0.0096** | **0.0185** | **0.005 / 0.002** | 4 |

Normierung (A2'): (1/N)ln|y_N| = −0.0703/−0.0688/−0.0680/−0.0674 (N = 400/600/800/1000), extrapoliert −0.0647 (a + b·lnN/N + c/N) bzw. −0.0655
(a + c/N) gegen ∫ln|z|dν = −0.0650. Alle Nullstellen ausser 4 (beschränkt, N-unabhängige Zahl, |t| ≤ 5) liegen in 0.818 ≤ |t| ≤ 1.1; die Abstände zur
Modellkurve fallen etwa wie ln N/N. **⇒ (A1) und (A2) sind numerisch sehr genau erfüllt**; die Nullstellen sind ein feines, vom Fit unabhängiges
Prüfobjekt — das Modell ist damit über Raten hinaus bis zur Nullstellenverteilung bestätigt.

### 22.9 Stand Runde 5 (ehrlich)
- **Bewiesen:** Lemma W3 (Bucht/Γout-Schranke ⇔ T-Schranken; Kreisschluss intrinsisch), Konturdarstellung von T, Satz W1'' (eine Hypothese (H*) auf
  γ*), Identität 22.5 (Xs Eckanteil über ker Rᵀ), Satz W5 (Nullstellenverteilung + Normierung ⇒ (H*) und (E_K) bis e^{o(N)}).
- **Numerisch (genau, N ≤ 1000):** (H*) auf der Modell-Linse (x = 0.7) und auf h(A-Bogen); (A1), (A2) wie oben; c_Modell = 0.0118 ± 0.0004.
- **NICHT bewiesen:** (H-V1)/(H-V3) bzw. (H*) für alle N. Fehlend sind genau zwei Bausteine: **(A1) für alle N** — asymptotische Extremalität der
  Eckmode auf K̃ (Satz vom Gonchar–Rakhmanov–Stahl/Blatt–Saff–Simkani-Typ für die Orthogonalität ∮y[F̃^j + x^jm̃]dx/x = 2πiu_j; neu ist die
  Komposition) — und **die Ecke** (e^{o(N)} → polynomial nahe L: Innenschicht, Amplitudenebene). Den zweiten Baustein übernimmt O mit seinem
  Doppel-Skalenmodell x = Le^{−w/N} (BRETT 03:05); Anschluss an den Kegel bei r ~ N^{−1/Re κ₁} = N^{−0.894} (Os Grenze N^{−0.9}).
- Bootstrap in N: ungeeignet (§22.4). Agmon-Gewichte: O rechnet nach, dass die Restverschiebung exp(0.68·N r^{κ₁}) kostet — nur bis r ~ N^{−0.9}
  gleichmässig; danach braucht es die Eikonal-Bedingung, und T_N vertauscht nicht mit Winkelgewichten (BRETT 03:05). Nicht weiter verfolgt.
- Os Mellin-Nullstellen 1 + μ^{s−1}: O und ich sind uns einig, dass 1.559 Zufall ist (Amplitude des Vorwärts-δ ohne Drehung um L).


---

# Runde 6 (ab 03:35): Extremalität, explizite Konstanten, Schnittstelle zu O

## 23. Befunde Runde 6

### 23.1 Auftrag 1: Folgt (H*) aus Extremalität + Remez/Bernstein–Walsh? — Nein (w44; Befund + Begründung)
Variationell ist 1/σ₁ = min_{v≠0} ‖J̃ᵀv‖/‖v‖ (ℓ²_ρ), die Eckmode v₁ ist der Minimierer. **Test (w44, arb):** Testvektor P(w) = Π(w − 1/z_j) mit
N Modellpunkten z_j auf ∂Ω_y (Quantile des harmonischen Masses ν, Anteil 0.4728 auf B — KEIN Datum benutzt), normiert.

| N | 200 | 300 | 400 | 600 | 800 |
|---|---|---|---|---|---|
| Überlapp \|⟨v, v₁⟩\| (\|⟨v, v₂⟩\|) | 0.919 (0.110) | 0.948 (0.198) | 0.985 (0.036) | 0.954 (0.233) | 0.997 (0.025) |
| Rayleigh-Quotient R = ‖J̃ᵀv‖ | 5.0e20 | 8.4e32 | 5.5e44 | 1.2e67 | 6.1e90 |
| wahres 1/σ₁ | 0.54 | 0.27 | 0.11 | 0.014 | 0.0016 |

Das Modell trifft den Eckunterraum in ℓ² auf 1–5 % (stärkster bisheriger Test, ohne Fit) — aber der Rayleigh-Quotient ist unbrauchbar, weil
**‖J̃‖ = e^{0.47N}** (gemessen σ_max(J̃) = 1.6e20 / 1.2e41 / 9.8e61 bei N = 100/200/300; Komposition mit |f| bis 1.35ρ auf ∂D).
**Begründung, warum der Weg nicht trägt:** (a) „Obere Schranke“ (σ₁ ≥ e^{(c−ε)N} per Testvektor) braucht Defektgenauigkeit e^{−cN} gegen
‖J̃‖ = e^{0.47N}, also einen ℓ²-Fehler ≤ e^{−0.48N} — aus einem asymptotischen Modell unmöglich (und für (H*) ohnehin die falsche Richtung).
(b) „Untere Schranke“ (‖J̃ᵀv‖ ≥ e^{−(c+ε)N}‖v‖_gewichtet für ALLE v) ist genau die gewichtete Koerzivität = (A1)/Agmon. (c) Die r(N) ~ N/250
Eckmoden haben dieselbe Rate (σ₁/σ₂ polynomial): keine exponentielle Lücke isoliert v₁; ℓ²-Nähe e^{−εN} wird durch Bernstein–Walsh auf der
Buchtseite (|w| bis 1.23) um 1.23^N verstärkt — Remez/BW können sie nicht in sup-Schranken auf K̃ verwandeln. **Ergebnis:** Extremalität in ℓ²_ρ
plus Remez/BW liefert (H*) nicht. In Os gewichteten Normen (E²(p), J_∞ beschränkt und koerziv) entfiele (a); (b) bleibt der globale Kern.

### 23.2 Minimale Eingänge für (E_K) — Reduktion (Skizze; (iii) bewiesen) und warum die Ecke den W̄-Rahmen braucht
**(i) Fern von L genügt eine Wachstumsschranke (I): σ_i ≤ CN^κe^{c'N}** für die Eckmoden. Auf ∂K' (K' = D̄ ∪ f(D̄)) liegen alle Urbilder auf oder
ausserhalb ∂D, also gibt D6 grob |T_i| ≤ C(f)√N·σ_i; mit Poisson–Jensen (P-III) trägt das, solange c'·ω_ζ(fern)/g(ζ) < 1 (nicht neu: Ps K'-Rahmen, dort Θ' ≤ 0.11 mit dem echten Profil).
**(ii) An der Ecke scheitert K' asymptotisch — jetzt erklärt:** der K'-Eckexponent ist 1/α' = 1.27, das Γin-Profil hat nach dem Kegel §19 den
Exponenten Re κ₁ = 1.118 < 1.27 (Ps „an der Grenze“ wird für N → ∞ zu „verletzt“). Im W̄-Rahmen genügt p > 1 (Satz W1/W1''). Also braucht die Ecke
**(H*)_loc: |Y_i| ≤ CN^κ auf γ* ∩ {|t − L| ≤ r₀}** (Linsen- und Buchtteil), r₀ ≈ 0.02–0.03 (Rand der Eckzone, wo c'·ω/g < 1 wird).
**(iii) Innenschicht durch Maximumprinzip (bewiesen):** Y_i und T_i sind in einer Umgebung von L analytisch (Pole nur in 0 bzw. 0, ã₀). Aus der
**Kreisbedingung (KB): max_{|t−L| = r_N} (|Y_i|, |T_i|) ≤ C_E N^{κ_E}, r_N := A·N^{−1/Re κ₁} = A·N^{−0.894}**, folgt dieselbe Schranke auf der
ganzen Scheibe |t − L| ≤ r_N — keine eigene Innenschicht-Analyse nötig. (KB) passt zum Kegel: dort ist N·G ≤ |c₁|·N·r_N^{1.118} = |c₁|A^{1.118} = O(1).
**Zusammen (Skizze): (E_K) ⇐ (I) + (H*)_loc auf r_N ≤ |t − L| ≤ r₀ + (KB).** Zu prüfen bleibt das Kleben: im Maximumprinzip für Y nahe L geht die
Grobrate des fernen Teils von γ* mit dem harmonischen Mass ein (Beitrag ~ N·c'·ω_ζ(fern)), also wieder ein Poisson–Jensen-Kriterium c'·sup ω/g < 1
in der Eckzone (dafür den fernen Teil von γ* auf ∂D legen, wo die Grobschranke σ₁ gilt). Gegenüber Runde 5 wäre der globale Teil damit auf eine
Wachstumsschranke (I) reduziert; die Gestalt (A1) würde nur noch in der Eckzone gebraucht.
**Test von (KB) (w46, arb, alle Eckmoden):** max über den Kreis |t − L| = A·N^{−0.894} von ln|Y| (innen) bzw. ln|T| (aussen):

| A \ N | 200 | 300 | 400 | 600 | 800 | 1000 | Fit |
|---|---|---|---|---|---|---|---|
| 0.5 | 0.69 | 1.34 | 1.58 | 1.83 | 2.03 | 2.19 | ~N^{0.87}, λ = −0.0014 |
| 1 | 1.22 | 1.86 | 2.10 | 2.35 | 2.55 | 2.70 | ~N^{0.86} |
| 2 | 2.33 | 2.96 | 3.19 | 3.43 | 3.62 | 3.77 | ~N^{0.83} |
| 4 | 4.73 | 5.34 | 5.54 | 5.77 | 5.94 | 6.07 | ~N^{0.78} |
| (Wert in L: ln\|T(L)\|) | 0.57 | 1.25 | 1.52 | 1.79 | 2.01 | 2.18 | ~N^{0.94} |

(KB) gilt numerisch mit C_E N^{κ_E}, κ_E ≈ 0.9 (keine exponentielle Komponente: λ < 0). Und: der Abstand zum Eckwert hängt nur von A ab
(A = 1: 0.65 → 0.52; A = 4: 4.16 → 3.89 über N = 200…1000; die Zunahme von A = 2 auf 4 ist 2.15-mal die von 1 auf 2 ≈ 2^{1.1}) — **die Kreisprofile
kollabieren auf der Skala N^{−1/Re κ₁}**: direkte Bestätigung der Anpassungsskala aus dem Kegel (bei Skala 1/N müsste der A = 4-Wert um ≈ 0.7 fallen).

### 23.3 Auftrag 2: explizite Konstanten (gemessen, kubisch, alle Eckmoden σ_i > 1.40; für Xs KONSTANTEN.md)

| N | 150 | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900 | 1000 |
|---|---|---|---|---|---|---|---|---|---|---|
| H*(N) = max_{γ*} ln\|Y\| | −0.02 | 0.33 | 0.98 | 1.15 | 1.29 | 1.39 | 1.48 | 1.58 | 1.67 | 1.74 |
| M_T(N) = max_{∂K_T} ln\|T\| | 0.15 | 0.57 | 1.25 | 1.52 | 1.68 | 1.79 | 1.90 | 2.01 | 2.10 | 2.18 |
| ε_eff (d = 0.05) | 0.663 | 0.625 | 0.599 | 0.606 | 0.615 | 0.622 | 0.628 | 0.632 | 0.635 | 0.639 |
| ε_eff (d = 0.1) | 0.831 | 0.813 | 0.800 | 0.803 | 0.807 | 0.811 | 0.814 | 0.816 | 0.818 | 0.819 |
| ln σ₁ − 0.0120N | −1.40 | −1.78 | −2.29 | −2.59 | −2.79 | −2.94 | −3.07 | −3.17 | −3.27 | −3.35 |

(γ*, K_T mit der Modell-Linse x = 0.7; ε_eff := ε_Γ − M_T/((N+2)d), ε_Γ ≥ 1 − M_ℬ/d = 0.682 / 0.841.) **Hüllen über alle gemessenen N:**
|Y_i| ≤ e^{−4.06}N^{0.883} auf γ*, |T_i| ≤ e^{−4.60}N^{1.026} auf ∂K_T, σ₁ ≤ e^{0.0120N − 1.40} (genauer σ₁ ≈ e^{0.0118N + 2.4}N^{−0.81}).
**Damit für (E_K) auf Γ_d: |T_i(ζ)| ≤ e^{−(N+2)·d·ε_eff}, ε_eff ≥ 0.60 (d = 0.05) bzw. 0.80 (d = 0.1) für 150 ≤ N ≤ 1000** — von ε_Γ ≥ 0.68/0.84
bleiben nach dem Vorfaktor 0.60/0.80 als Reserve; der Vorfaktor wächst nur wie N^{1.03}, ε_eff steigt ab N = 300 wieder.
**Quantitative Form von (A1) (gemessen, w43):** Nullstellen im Abstand η_B = 1.39·ln N/N von B und η_Linse = 2.7·ln N/N von der Linsen-Nulllinie,
Winkel-Diskrepanz ≈ 0.7·ln N/N. **Lemma W5q (bewiesen, elementar):** liegen alle Nullstellen innerhalb η von Modellpunkten mit Quantil-Zuordnung,
so ist (1/N)ln|y(t)| ≤ G(t) + 2η/dist(t, ∂Ω_y) + ε_disc(t) + ε_lc (|ln|w−ω| − ln|w−ω̃|| ≤ 2η/d für d ≥ 2η; ε_disc = Quadraturfehler der Quantilpunkte,
O(ln N/N) abseits der Kurve; ε_lc aus der Normierung auf |t| = 2). Mit η = a·ln N/N wird der Faktor e^{N·2η/d} = N^{2a/d}: **polynomial** auf jedem
festen Stück von γ* abseits von L; zur Ecke hin wächst der Exponent wie 1/d ~ r^{−1.118} — dort übernimmt (KB). Die ln N/N-Verschiebung ist genau die
Spur der polynomialen Vorfaktoren (Verschiebung = Δκ·ln N/(N|∇G|)).
**Ehrlich:** alle Zahlen sind für N ≤ 1000 gemessen (arb, Eckmoden aus float-SVD), NICHT für alle N ≥ N₀ bewiesen.

### 23.4 Auftrag 3: Schnittstelle zu O (Ecke)

| | liefert | Form | Stand |
|---|---|---|---|
| **O** (Innenmodell x = Le^{−w/N}, Wiener–Hopf) | (KB) auf dem Kreis r_N = A·N^{−0.894}, d. h. \|w\| ≈ A·N^{0.106}/\|L\| | max(\|Y_i\|, \|T_i\|) ≤ C_E N^{κ_E}, κ_E, C_E explizit | offen (O) |
| **W** (Kegel §19, Modell §22) | Aussenteil (H*)_loc auf r_N ≤ r ≤ r₀ und Anschluss an (I) | \|Y_i\| ≤ CN^κ auf γ* ∩ {r_N ≤ r ≤ r₀}, r₀ ≈ 0.03 | numerisch (w38, w41, w46) |
| **beide** | Anschlussbedingung im Überlapp r ~ r_N | O: Wachstum des Fast-Kerns für \|w\| → ∞ höchstens exp(\|c₁\|·N^{1−κ₁}\|Lw\|^{Re κ₁}(1+o(1))) mit Naht + 2. Knoten im Abstand 161°; W: Kegelentwicklung mit polynomialem Vorfaktor bis r_N | |

Daraus übernehme ich: (KB) ⇒ Innenscheibe per Maximumprinzip (23.2(iii)) ⇒ (H*) auf γ* ∩ {r ≤ r_N}. Prüfgrössen für O: Eckwert N^{κ(0)}, κ(0) ≈ 0.51
(L_T- und L_Y-Fit); Kollaps ln(Λ_N(τ)/Λ_N(0)) ≈ 0.0018·τN (τ in Grad, τN ≲ 100) (w34); Kreiswerte w46.

### 23.5 Stand Runde 6
- Auftrag 1: Extremalität in ℓ²_ρ trägt NICHT (‖J̃‖ = e^{0.47N}); aber das Modell reproduziert den Eckunterraum in ℓ² auf 1–5 % (w44).
- Neu bewiesen: Lemma W5q, Innenscheibe per Maximumprinzip aus (KB) (23.2(iii)). Skizze: Reduktion (E_K) ⇐ (I) + (H*)_loc + (KB) (23.2).
- Erklärt: warum Ps K'-Rahmen an der Ecke asymptotisch scheitert (1.118 < 1.27).
- Offen für alle N: (I) (Wachstumsrate, global), (H*)_loc auf r_N ≤ r ≤ r₀ (lokale Gestalt), (KB) (O).

### 23.6 Nachträge (Antworten an X und O, 03:57)
**(a) Konstanten C_E, κ_E auf dem Weg W1'' (w47).** Ohne Lemma W1a fällt der Verlust N^{κ+2} weg: D6 auf ∂K_T gibt
max_{∂K_T}|T_i| ≤ C(f)·sup_{Ext γ*}|Y_i| + max(1/|ζ|)·√(N+1), mit **C(f) = 4.95** (|m/ζ| 1.66, Urbilder 0.77, niedrige Koeffizienten über Cauchy auf
|x| = 2ρ 2.52) und max 1/|ζ| = 0.728. Mit der gemessenen (H*)-Hülle e^{−4.06}N^{0.883}: max|T| ≤ 16.1 / 19.0 / 25.8 / 42.1 / 61.0 (N = 150/193/300/600/1000),
also **C_E = 1, κ_E = 0.6** für 150 ≤ N ≤ 1000 (der Kopfterm √N dominiert bis N ≈ 300). Direkt gemessen ist M_T ≤ e^{−4.60}N^{1.03} (kleiner).
**(b) Ableitungen der Eckrate aus dem Modell (w48, w49; Xs Prüfstein 03:49).** Kubische Familie mit festem L und μ = √L + δμ:
c(δμ = ±0.02) = 0.01264 / 0.01089, c(δμ = ±0.02i) = 0.01256 / 0.01097 ⇒ **∂c/∂μ = 0.044 − 0.040i** (in Xs Schreibweise δc = 0.044 Re δμ + 0.040 Im δμ);
X misst an p_N **0.042 − 0.035i** — Übereinstimmung auf 5 % bzw. 13 %, ohne Fit. Lokalität: eine Störung ε(x² − 2Re L·x + |L|²)² (verschwindet doppelt
in L, L̄; ε = ±0.003 ändert f auf Γout um ≤ 1.1 %) verschiebt c nur um ∓0.00013 (1.1 %) — pro relativer Änderung von f auf Γout etwa 2.5-mal schwächer
als die μ-Störung: c hängt vorwiegend, aber nicht ausschliesslich von (L, μ) ab. ∂c/∂r (Fixpunkt neben dem Kreis) kann das Modell in seiner jetzigen Form
nicht (es setzt L auf den Kreis); dafür wäre die Ecke mit Versatz neu aufzusetzen.
**(c) Os Vorfaktor N^{1/2} (BRETT 03:52) und mein κ(0) ≈ 0.51.** O misst die normierte VORWÄRTS-Eckmode u₁ (Breite 1/N um arg L, daher Höhe √N durch die
Normierung). Die exakte Identität bei t = L (D6/§5 im Fixpunkt) koppelt die Eckwerte: y(L)·(m̃(L)/L + 1/(L·f'(L))) − T(L) = 𝒰(L) − P(L) − (ferne Urbilder),
und 𝒰(L) ist die Erzeugende von u₁ in L ~ √N. Damit ist der adjungierte Eckwert |y(L)| ~ N^{1/2} (gemessen κ = 0.50 für Y, 0.51 für T in L) bis auf eine
Auslöschung zwischen y(L) und T(L) durch u₁ festgelegt — die Innenschicht von O liefert also genau den Vorfaktor, den (KB) braucht.


---

# Runde 7 (ab 03:59): Wachstumsschranke, Eckzone in der Streifenkoordinate

## 24. Befunde Runde 7

### 24.1 Auftrag 1 — Wachstumsschranke σ₁ ≤ e^{c'N}: Versuch über einen grösseren Skalenkreis (scheitert, w50)
**Bewiesen (Normäquivalenz auf P_N):** für ρ' ≥ ρ ist ‖δ‖_ρ ≤ ‖δ‖_{ρ'} ≤ (ρ'/ρ)^N‖δ‖_ρ, also **σ₁(ℓ²_ρ) ≤ (ρ'/ρ)^N·‖J_N⁻¹‖_{ℓ²_{ρ'}}**. Gäbe es ein
ρ' > |L| (Fixpunkt im Innern der Scheibe), in dem die Abschnitte gleichmässig stabil sind, so folgte (I) mit c' = ln(ρ'/ρ) aus einer Stabilitätsaussage
OHNE Ecke. **Gemessen (arb-Inverse, float-SVD; ln σ₁ in ℓ²_{ρ'}, ρ' = |L|(1 + s)):**

| s \ N | 100 | 200 | 300 | 400 | 500 | 600 | 800 |
|---|---|---|---|---|---|---|---|
| 0 (ℓ²_ρ) | 0.31 | 0.62 | 1.31 | 2.21 | 3.21 | 4.26 | 6.43 |
| 0.005 | 0.31 | 0.49 | 0.98 | 1.65 | 2.43 | 3.27 | 5.01 |
| 0.01 | 0.31 | 0.40 | 0.73 | 1.21 | 1.80 | 2.45 | 3.84 |
| 0.02 | 0.32 | 0.33 | 0.44 | 0.65 | 0.93 | 1.28 | 2.08 |
| 0.025 | 0.32 | 0.32 | 0.37 | 0.49 | 0.68 | 0.91 | 1.55 |
| 0.03 | 0.32 | 0.32 | 0.34 | 0.40 | 0.54 | 0.90 | 2.67 |
| 0.05 | 0.33 | 0.33 | 0.49 | 1.86 | 3.90 | 6.26 | 11.54 |
| 0.1 | 0.35 | 2.67 | 7.91 | 13.99 | 20.49 | 27.26 | 41.31 |

(s < 0, also ρ' < |L|, wächst schneller.) Für 1.005 ≤ ρ'/|L| ≤ 1.025 ist die Eckmode nur verzögert (Rate bei N = 600–800 noch 0.003–0.008, steigend);
ab ≈ 1.03 setzt ein zweiter Modus ein (Rate 0.026 bei 1.05, 0.07 bei 1.1; nicht von Nullstellen von m = f'∘f, die liegen bei ≥ 3|L|). **Es gibt keinen
Skalenkreis ohne Wachstum**; der Weg liefert bestenfalls c' ≈ ln 1.02 + 0.004 ≈ 0.024 und auch das nur numerisch. **(I) bleibt unbewiesen.** In Xs
Form ist (I) eine untere Schranke für σ_min(S_N) = σ_min(Q_bot)² (§22.5), also eine Schranke für J_N^{−T} auf den Überlauf-Funktionalen — dieselbe
globale Art wie (A1).

### 24.2 Auftrag 2 — Eckzone r_N ≤ |t − L| ≤ r₀ in der Streifenkoordinate (Analyse)
- **Ratenebene (Skizze, für das Modell):** in ζ = log((x−L)/(x−L̄)) wirkt h als ζ ↦ ζ − ℓ + O(e^{Re ζ}); die Störung ist relativ O(r) und über die Zone
  summierbar (∫ r d ln r ≤ r₀). Die Kegelentwicklung H = Σ Re(c_n w^n) + O(r^{Re κ₁ + 1}) gilt gleichmässig in der Zone.
- **Länge der Zone:** ln(r₀/r_N)/ln|μ| ≈ 0.894·ln N/0.159 = 5.6·ln N h-Schritte. Selbst ein konstanter Verlustfaktor C je Schritt gäbe nur N^{5.6·ln C}.
- **Lemma W9 (bewiesen bis auf die Restterme): kontrahierende Rekursion für y ausserhalb D.** Aus der Identität:
  y(t) = (t/m̃(t))·[T(t) − a(t)·y(h(t)) − R(t)], a = 1/(h·F̃'∘h), R = ferne Urbilder + 𝒰 − P. In L ist |a·t/m̃| = (1/|μ|)/|μ| = 1/|μ|² = 1/|L| = 0.727 < 1.
  Längs einer h-Bahn (Drehung um −arg μ = −38.3°, Schrumpfung um 1/|μ|) ist also |y(t)| ≤ C·Σ_j 0.727^j·|T(h^j t)| + (Bucht-Werte nach Eintritt in D)
  + Σ_j 0.727^j|R| — geometrisch, kein Verlust über die Zonenlänge. Die fernen Urbilder in R brauchen |y(x)| für |x| ≥ 2.4ρ, also bei der Grobschranke
  σ₁ ≤ e^{0.012N} nur die niedrigen Koeffizienten |y_k| ≤ poly für k ≤ 0.012N/ln 2.4 ≈ 0.014N (gemessen: max_{k<0.3N}|y_k| ≈ 0.008).
- **Satz W10 (bedingt, Skizze):** auf der festen Vergleichsfläche (Naht Γout ∼ h(Γout), dieselbe Deckverschiebung 2πi − log μ, also dieselben
  Exponenten nκ₁) ist Ψ = Re(c₁w) + Re(c₂w²) harmonisch; Maximumprinzip auf dem Kegelring r₁ ≤ |t − L| ≤ r₀, dann Lemma W9 und Lemma W3, ergibt
  (H*)_loc, FALLS (a) (KB) bei r₁ = r_N·(C ln N)^{0.894}, (b) U ≤ Ψ + C·ln N/N auf |t − L| = r₀ (winkelaufgelöst), (c) die Klebung an den Nähten im
  Rahmen bleibt, (d) niedrige Koeffizienten wie oben.
- **Warum das nicht lokal schliesst:** (b) ist die globale Gestalt am Aussenrand — die Grobschranke (I) gibt dort nur FLACHE Daten U ≤ c', und deren
  harmonische Fortsetzung (Hadamard, log-Interpolation zwischen r₁ und r₀) ist auf γ* positiv. (c) braucht entweder untere Schranken an den Nähten
  oder Böden; Böden vertragen sich nicht mit den negativen Randwerten von Ψ. Beides ist wieder Information über den positiven Lappen, also globale Gestalt.
  **Die Eckzone ist kein in sich geschlossenes Störungsproblem;** sie erbt am Aussenrand den globalen Eingang.

### 24.3 Ehrliche Lage nach Runde 7
Bewiesen: Normäquivalenz-Schranke (24.1), Lemma W9 (Kontraktion 1/|L| ausserhalb D), dazu alles aus §22–§23. Nicht bewiesen und nach meiner Analyse
nicht mit lokalen Mitteln beweisbar: (I) und (A1) auf einem Kreis um L (polynomial genau, winkelaufgelöst) samt niedrigen Koeffizienten. Diese
Aussagen sind Teile EINER globalen Aussage über das freie Randwertproblem (Lage von Buchtrand und Linse bzw. die Normierung c). Vorschlag: im Hauptsatz
als benannte Hypothese (H₀) führen — gemessen bis N = 1000 mit C_E = 1, κ_E = 0.6, ε_eff ≥ 0.60 (BRETT 04:06 an X) — und den Beweis als eigenes,
langes Vorhaben ansetzen (Gonchar–Rakhmanov/Riemann–Hilbert mit Verschiebung auf der Fläche ℳ).


---

# Langes Vorhaben „(H₀) global“ (ab 04:08)

## 25. Schritt 1 und 2 entschieden: keine Gonchar–Rakhmanov-Orthogonalität; S-Eigenschaft nur auf der Grenzfläche ℳ

### 25.1 Die exakten Relationen (bewiesen)
**(a) Adjungiert (Eckmode y, ‖u‖ = 1):** für alle q ∈ P_N und JEDE Kontur Γ um 0 (Integrand analytisch bis auf 0):
  ∮_Γ y(z)·[q(F̃(z)) + m̃(z)q(z)] dz/z = 2πi⟨u, q⟩.
**(b) Vorwärts, exakte Interpolanten:** Q_n monisch vom Grad n mit Q_n∘f + m·Q_n = O(x^n) bei 0 (Kern des Rechteckabschnitts T_{n−1}J_∞|P_n,
existiert stets). Äquivalent, auf jeder grossen Kontur Γ um 0 und a₀ = f(0):
  ∮_Γ Q_n(s)·ψ_k(s) ds = 0 (k < n),  ψ_k(s) = R_k(s) + m(s)s^{−k−1},  R_k(s) := Σ_{f(x)=s} x^{−k−1}/f'(x) = −[x^k](f(x) − s)^{−1}
  (R_k rational, einziger Pol in a₀ der Ordnung k+1; Erzeugende Σ_k ψ_k(s)z^k = 1/(s − f(z)) + m(s)/(s − z)).
**(c) Biorthogonalität:** J_N ist die Matrix der Bilinearform ⟨⟨A, Q⟩⟩ = (1/2πi)∮ A(1/x)·[Q(f(x)) + m(x)Q(x)] dx/x; mit den biorthogonalen Paaren
(Q_n, A_n) und h_n = ⟨⟨A_n, Q_n⟩⟩ gilt J_N⁻¹ = Σ_{n≤N} Q_nA_nᵀ/h_n (falls alle Hauptminoren ≠ 0). **Gemessen (w51, w52):** das letzte Paar liegt im
Eckunterraum (Q_N in span{u₁, u₂} zu 0.99; |⟨A_N, v₁⟩| = 0.92–0.97) und ‖Q_N‖‖A_N‖/|h_N| = 1.2 / 1.7 / 4.1 / 11.1 gegen σ₁ = 1.5 / 1.9 / 3.7 / 9.1
(N = 150/200/300/400), Rate von ln(‖Q_n‖‖A_n‖/|h_n|) bei n ≈ 500: 0.0105 → c. Die Nullstellen von A_N liegen aber anders als die der Eckmode
(Median |t| = 1.06 statt Linse bei 1.00–1.02, w53) — die Eckmode ist eine Mischung vieler letzter Paare, nicht A_N allein.

### 25.2 Entscheidung zu Schritt 1: es gibt KEINE Orthogonalität vom Gonchar–Rakhmanov-Typ — woran es scheitert
GR (und Stahl) brauchen: die Unbekannte ist ein Polynom P_n in EINER Koordinate, orthogonal zu {s^j·W_n(s)}_{j<n} mit EINEM (ggf. von n abhängigen)
Gewicht auf einer Kontur. Hier gilt:
1. **Zwei Koordinaten mit gleichem Index.** Jede Bedingung ist eine Summe aus einem Moment in a₀, durch f transportiert (R_k = Potenzen der
   Umkehrzweige), und einem Moment in 0 (Multiplikation mit m). span{ψ_k : k < n} ist ein n-dimensionaler Teilraum von
   span{(s−a₀)^{−l−1}, m(s)s^{−k−1}} (2n-dimensional), aber NICHT von der Form {s^j W}. Getrennte Bedingungen in 0 und a₀ wären Mehrpunkt-Padé
   (das hat GR-Form); die Kopplung derselben Ordnung k über f ist neu.
2. **Kontur auf ∂Ω_y verschoben (dominante Terme):** auf B dominiert der Kompositionsterm (|F̃| > 1), auf der Linse der Multiplikationsterm. In der
   t-Koordinate der T-Seite lautet (a) dann ∫_{Z_T} q(t)·y(h(t))·ω₁ dt + ∫_{Linse} q(t)·y(t)·ω₂ dt ≈ ⟨u,q⟩: Testfunktionen Polynome in t, die Unbekannte
   aber als y∘h auf Z_T und als y auf der Linse — eine Riemann–Hilbert-Aufgabe mit nicht-Carlemanscher Verschiebung h (Fixpunkte L, L̄).
3. **Der eigentliche Grund:** Auf der geklebten Fläche ℳ (B ∼ Z_T) würden beide Teile EINE Orthogonalität (Unbekannte in L(N·[0]), Testfunktionen in
   L(N·[∞]); bei Geschlecht 0 über die Uniformisierende φ Polynome in 1/φ bzw. φ — GR-Form!). Aber ℳ existiert nur im GRENZWERT: die Klebung T ≈ y∘h
   ist eine Dominanzbilanz, keine Identität (die übrigen Terme der Identität sind nur relativ klein, wo y∘h exponentiell gross ist). Bei endlichem N
   ist y keine meromorphe Funktion auf einer festen Riemannschen Fläche mit nur einem Pol. Damit fehlt die exakte Orthogonalität, von der GR ausgeht.
4. **Auch keine klassische RH-Aufgabe:** die Daten sind polynomial (f, m Polynome), es gibt kein Gewicht mit Sprüngen; die Sprungkurven (B, Linse)
   entstehen erst aus der N-abhängigen Bilanz (wie Szegő-Kurven), sind also frei.

### 25.3 Schritt 2: S-Eigenschaft
- **In der Ebene: nicht anwendbar bzw. verletzt.** Der Träger des Nullstellenmasses ist die GESCHLOSSENE Kurve ∂Ω_y (trennt 0 von ∞); das
  Green-Potential ist einseitig (g_{Ω_y} innen, 0 aussen). Die S-Eigenschaft (gleiche Normalableitungen beidseits) gehört zu Bögen mit
  zusammenhängendem Komplement; hier ist die Normalableitung aussen 0. Das Nullstellenmass ist das harmonische Mass von Ω_y (Balayage) — die
  Lage asymptotisch extremaler Polynome auf einem Kompaktum mit Innerem (Blatt–Saff–Simkani), nicht nicht-hermitescher Orthogonalität auf S-Bögen.
- **Auf ℳ: erfüllt (numerisch).** Das richtige Analogon ist die Glattheit von H über die Naht B ∼ Z_T (gleiche Cauchy-Daten); w29–w31: H = H∘h
  auf einer offenen Menge auf 5e-4; im Modell eingebaut (§15, §20). Das ist ein starkes Indiz für ein GR-artiges Prinzip AUF der Grenzfläche —
  aber ohne exakte Orthogonalität (25.2.3) fehlt der Hebel.

### 25.4 Nächstbeste globale Struktur (Vorschlag)
**(A) „Approximative Riemannsche Flächen“ ℳ_N.** Für jedes N eine Fläche ℳ_N (Klebung entlang einer festen Kurve nahe B ∼ Z_T) und darauf die
geklebte Funktion; die Klebungsfehler sind relativ e^{−N·H} (H > 0 an der Naht), also exponentiell klein ausser bei L, L̄ (dort Os Innenschicht).
Ziel: ein GR-Argument „bis auf exponentiell kleine Klebungsfehler“ — man bräuchte eine Stabilitätsversion von GR (Nullstellenverteilung stabil
unter relativ kleinen Störungen der Orthogonalitätsrelation). Offen; mir ist keine Literatur dazu bekannt.
**(B) Exakte Interpolanten statt Fast-Kern (25.1(b,c)).** Die Eckrate c steckt in ‖Q_N‖‖A_N‖/|h_N| — exakte Objekte ohne Fast-Singularität.
Für Q_N gilt die Funktionalgleichung Q_N∘f = −m·Q_N + x^N S_N; im Grenzwert liefert sie eine Gleichung für das logarithmische Potential
U∘f − U = (Bilanz) — wieder ein freies Randwertproblem, aber für exakte Grössen. Analytisches Vorbild: Muttalib–Borodin-Biorthogonalität
(Komposition x ↦ x^θ; für rationales θ Mehrfach-Orthogonalität mit (r+1)×(r+1)-RH, Deift–Zhou, lokale Parametrix mit Meijer-G). Unser Fall: Komposition
mit dem Kubik f (d = 3 Zweige, also Struktur der Grösse d + 1 = 4) PLUS Toeplitz-Term m — die Summe zweier Formen hat dort kein Vorbild.
**(C) RH mit Verschiebung (Litvinchuk; Karlovich–Kravchenko).** Fredholm-Theorie für Operatoren mit nicht-Carlemanscher Verschiebung mit Fixpunkten
existiert; ein nichtlinearer steilster Abstieg für solche Aufgaben ist mir nicht bekannt.
**Einschätzung:** Keiner der drei Wege ist ein Standardwerkzeug; (A) ist dem Problem am nächsten und am ehesten mit Os Innenschicht verträglich.

### 25.5 Literatur (Kennzeichnung)
- Kuijlaars, Molag: The local universality of Muttalib–Borodin biorthogonal ensembles with parameter θ = 1/2, Nonlinearity 32 (2019), arXiv:1810.00741 —
  **nur zitiert** (Zusammenfassung gesehen: θ = 1/2 → Mehrfach-Orthogonalität Typ II, 3×3-RH, Deift–Zhou, Parametrix in 0 mit Meijer-G).
- Wang, Zhang: A vector Riemann–Hilbert approach to the Muttalib–Borodin ensembles, arXiv:2103.10327 — **nur zitiert** (Abstract: zwei 1×2-Vektor-RH,
  ganzzahliges θ, Parametrix (θ+1)×(θ+1)).
- Claeys, Romano: Biorthogonal ensembles with two-particle interactions, Nonlinearity 27 (2014) 2419–2444, arXiv:1312.2892 — **nur zitiert** (Abstract:
  Rekursion/Christoffel–Darboux für rationales θ, Vektor-RH, Gleichgewichtsmass).
- Gonchar–Rakhmanov (Equilibrium distributions and degree of rational approximation, 1987), Stahl (Padé-Theorie, 1985–86), Blatt–Saff–Simkani (1988),
  Saff–Totik (Logarithmic Potentials with External Fields, 1997), Fokas–Its–Kitaev, Deift–Zhou, Litvinchuk (BVP with shift) —
  **nur zitiert, aus dem Gedächtnis** (nicht nachgelesen in dieser Sitzung).


---

# Weg (A), Stufe n-te Wurzel (ab 04:24)

## 26. Fehlertoleranz, Klebungslücken, und warum e^{o(N)} für Stufe (i) nicht reicht

### 26.1 Klebungslücken (gemessen am Referenzmodell)
Auf einer festen Vergleichsfläche ℳ₀ werden y (innen) und T (aussen) an zwei Nähten verklebt; der relative Klebungsfehler an einem Nahtpunkt p ist
e^{−N·δ(p)}, δ = (Rate des dominanten Terms) − (Rate der übrigen Terme der Identität):
- **Naht Γout ∼ h(Γout)** (T ≈ y∘h/(hF̃'∘h)): δ_out(φ) = H(h(e^{iφ})) = 0.0010 / 0.0024 / 0.0058 / 0.0139 / 0.033 / 0.072 / 0.167 / 0.38 bei
  0.5 / 1 / 2 / 4 / 8 / 15 / 30 / 76.6° vor L.
- **Naht Γin** (T ≈ (m̃/t)·y): δ_in(τ) = Λ(τ) = 0.0005 / 0.0010 / 0.0018 / 0.0030 / 0.0048 / 0.0069 / 0.0092 / 0.0111 / 0.0118 bei
  0.5 / 1 / 2 / 4 / 8 / 15 / 30 / 60 / 103°. Das ist die BINDENDE Lücke. Legt man die Naht um η ins Innere, wächst sie auf ≈ Λ(τ) + 0.2·η
  (Steigung des Stegs 0.95 gegen Quelle ≤ 0.75), am Rand von D muss sie aber durch L.
- **Beide Lücken verschwinden genau in L, L̄** (wie r^{1.1}); sonst sind sie positiv.
**Die Ungleichung, nach der die Leitung fragt:** Auf n-ter-Wurzel-Ebene genügt δ(p) > 0 für alle Nahtpunkte p ≠ L, L̄ — das ist erfüllt
(min δ ausserhalb U_r(L) = δ_in: ≈ 0.0018 bei r = 0.03, ≈ 0.004 bei r = 0.1, ≈ 0.007 bei r = 0.25). Bei endlichem N ist der Fehler
erst klein, wenn N·δ ≳ 3: bei N = 1000 ab τ ≈ 4°, bei N = 400 ab τ ≈ 15°.

### 26.2 Fehlertoleranz — die eigentliche Frage ist Eindeutigkeit, nicht Toleranz
In der n-ten Wurzel verschwinden alle Fehler mit strikt kleinerer Rate. Das sieht man am saubersten in einer Kompaktheitsform, die GANZ ohne
Orthogonalität auskommt:
**Satz W11 (bewiesen, Standard-Potentialtheorie).** Sei v_N der normierte Eckvektor (‖v_N‖ = 1, J̃ᵀv_N = ε_N u_N, ‖u_N‖ = 1) und V_N(t) = Σ v_{N,k}t^{−k}.
(a) u_N := (1/N)ln|V_N| ist lokal gleichmässig nach oben beschränkt (≤ max(0, ln(1/|t|)) + o(1)); jede Teilfolge hat eine L¹_loc-konvergente
Teilfolge mit Grenze u ≢ −∞ und max_{|t|=1} u = 0 (Hartogs). Ebenso für w_N := (1/N)ln|T_N| und für die Raten der Quelle.
(b) Die Grenzraten der Terme der exakten Identität T = m̃V/t + 𝒱 − ε𝒰 + P erfüllen überall die **tropische Bilanz**: jede Rate ist ≤ dem Maximum
der übrigen (das Maximum wird mindestens zweimal angenommen). Beweis: |a_i| ≤ Σ_{j≠i}|a_j| ⇒ (1/N)ln|a_i| ≤ max_{j≠i}(1/N)ln|a_j| + ln 4/N; Grenzübergang
in L¹_loc (max ist 1-Lipschitz), dann f. ü. ⇒ überall (subharmonisch).
**Folgerung:** Die Klebungsfehler sind in (b) automatisch enthalten; es bleibt nur zu zeigen, dass das Grenzsystem (b) samt Normierung und
Singularitäten (u ~ ln(1/|t|) in 0, w ~ −ln|t| in ∞) **eindeutig** ist bzw. dass die Quellrate −c* nicht unter −c_Modell liegen kann:
**Vermutung W12 (Extremalität des Grenzproblems).** (a) Jede Grenzkonfiguration hat c* ≤ c_Modell; (b) Gleichheit nur für die Modellkonfiguration
(Green-Funktion von Ω_y bzw. H auf ℳ). Aus W11 + W12 folgen lim (1/N)ln σ₁ = c_Modell (das wäre (I) scharf!) und (A1) auf n-ter-Wurzel-Ebene.
Struktur von W12: ein **Hindernisproblem mit nichtlokaler (Verschiebungs-)Bedingung** — die Bucht ist die Kontaktmenge (Rate = Boden der Quelle),
B ihr freier Rand, die Klebung H = H∘h die nichtlokale Kopplung; c ist sein Hauptwert. Für lokale Hindernisprobleme gibt es Eindeutigkeit über
Vergleichsprinzipien; für die nichtlokale Kopplung ist mir nichts bekannt.

### 26.3 Wichtig für die Planung: e^{o(N)} reicht an der Ecke NICHT, auch nicht für Stufe (i)
Die Kneser-Singularität sitzt IN L; ihre Daten fallen nur polynomial (Koeffizienten ~ ρ^{−k}k^{−σ−1}). Der Eckanteil des Residuums ist ein Integral
über den Schlitz nahe L: ∫₀^{u₀} |T(ζ(u))|·u^σ du. Mit |T| ≤ e^{ε_N N}·e^{−N a u^b} (n-te Wurzel, ε_N → 0) ist der Beitrag von u ≲ (ε_N/a)^{1/b}
von der Grösse e^{ε_N N}·(ε_N/a)^{(σ+1)/b} — das wächst, sobald ε_N·N ≫ ln N. **Also braucht schon Stufe (i) (Existenz für grosse N, p_N → f_K)
eine POLYNOMIALE Eckschranke (E_L)**, nicht nur e^{o(N)}. Und die Eckzone lässt sich nach §24.2 (Satz W10(b)) nur mit polynomial genauen
Randdaten am Aussenrand r₀ schliessen — n-te-Wurzel-Asymptotik liefert dort nur o(1)-Genauigkeit. **Folgerung:** Für (i) wie für (ii) braucht man
starke (polynomial genaue) Asymptotik in einer festen Umgebung von L, und die ist nur global zu haben (Deift–Zhou-artig: äussere Parametrix aus H auf ℳ,
lokale Parametrix in L von O, Kleinheitsargument für die nichtlokale Störung). Die n-te-Wurzel-Stufe (W11 + W12) bleibt als Teilergebnis wertvoll:
sie liefert c und die Gestalt, also die g-Funktion für den steilsten Abstieg, und den fernen Teil von (E_K).
**Was zwischen (i) und (ii) liegt:** (ii) braucht zusätzlich explizite Konstanten in jedem Schritt der starken Asymptotik (Kleinheitsnorm der
Störung, Anschluss der Parametrizes), (i) nur deren Existenz.


---

# W12 — Eindeutigkeit des Grenzproblems (lange Strecke ab 04:35), Zwischenstand

## 27. Formulierungen und erste Befunde zu W12

### 27.1 Das Grenzproblem als Schweissproblem (Umformulierung, bewiesen als Äquivalenz der Bedingungen)
{H = 0} ist auf ℳ eine geschlossene Kurve (Naht B ∼ Z_T plus Linsen-Nulllinie λ), die 0 von ∞ trennt. Mit Ω_y (Gebiet um 0, berandet von B und λ) und
E (Aussengebiet, berandet von Z_T = f(B) und λ) ist ℳ = Ω_y ∪_g E, geklebt längs ∂Ω_y ↔ ∂E durch g = id auf λ, g = f auf B. H = −ln|φ| mit der
Uniformisierenden φ: ℳ → Ĉ (φ(0) = 0, φ(∞) = ∞). **Die freie Randbedingung ist: φ bildet die Klebekurve auf einen KREIS ab**, d. h. mit den
Riemann-Abbildungen φ_in: Ω_y → 𝔻 (0 ↦ 0) und φ_out: E → Ĉ∖𝔻̄ (∞ ↦ ∞) gilt φ_out∘g = φ_in auf ∂Ω_y (konforme Schweissung mit Kreisbedingung).
Unbekannte: zwei Kurven (B, λ); Bedingung: eine reelle Funktion auf B ∪ λ — Anzahl passt. c = −ln|φ(−1)|. In L trifft id auf f: dort ist die
Schweissung nicht analytisch (komplexe Streckung μ), das ist der Kegelpunkt mit κ₁.

### 27.2 Lokale Eindeutigkeit des Modells (numerisch)
Jacobi-Matrix der Kleinste-Quadrate-Lösung (Referenzmodell, 20 Formparameter des Buchtrands): Singulärwerte 3.43 … 0.009 (voller Rang,
Kondition ≈ 400). Schwächste Richtung = Gesamttiefe der Bucht: Schritte ±0.05 / ±0.1 / ±0.2 längs ihr erhöhen den Rest 4.9e-3 → 5.2e-3 / 8.5e-3 /
2.8e-2 und ändern c nur in 2. Ordnung (0.01176 → 0.01171 / 0.01155 / 0.0110). **Das Modell ist eine (flache) isolierte Lösung**; c ist robust,
die Buchttiefe bei 0° nur auf ≈ ±0.01 bestimmt.

### 27.3 Verbandsstruktur der Unterlösungen (bewiesen) — und wo die Normierung herkommt
Nenne (u, w) eine tropische Unterlösung zur Quellrate s, wenn u, w subharmonisch sind (mit den Singularitäten wie in W11) und für die Terme der
Identität gilt: jede Rate ≤ Maximum der übrigen (Quelle mit Rate s). **(a)** Das Maximum zweier Unterlösungen ist eine Unterlösung (alle Terme sind
monoton in (u, w): Identität und Komposition); die Klasse ist also ein ⊕-Verband — Perron-Konstruktion möglich. **(b)** Die Hüllkonfiguration
lim (1/N)ln‖J̃_N⁻¹e_t‖ (Supremum über alle Quellen u) ist eine Unterlösung zur expliziten Quellrate s_max = max(0, ln(1/|t|)).
**(c) Befund zur Normierung (wichtig):** In der Modellkonfiguration steht die Quelle fast nirgends im Maximum — in Ω₁ balancieren y und T, im
Halbmond und in der Bucht T und y∘h, weit draussen y gegen seine eigenen Urbild-Terme (alle auf dem Sockel). **Tropisch unsichtbar ist die Quelle
also ausser an der Ecke**: innen nahe L wächst die Rate von y wie x_eff·(Tiefe) mit x_eff = 0.59–0.73 (w38, ψ = 180°), die der Quelle u₁ (sitzt bei arg L,
Indizes 0.4–0.75N) wie ≈ 0.75·(Tiefe) — dort balancieren y und Quelle. Folge: das Grenzsystem ist bis auf die Ecke homogen; ohne den Eckanteil
bestimmt es nur den ANSTIEG (Maximum minus Sockelniveau), nicht das absolute Niveau. Die absolute Wachstumsrate c = Anstieg + (Sockel − Quelle)
hängt an der Ecke: y(L)·(μ + 1/μ)/L̃ = T(L) + 𝒰(L) − P(L) − (ferne Urbilder) — die beiden y-Terme im Fixpunkt löschen sich NICHT aus (μ + 1/μ ≠ 0).
**W12 ist deshalb mit der Normierung in der Ecke zu formulieren** (u(L) = Quellniveau), nicht mit max_{|t|=1} u = 0 allein; ob das tropische System
u(L) = 0 erzwingt (T kann im Fixpunkt ausgleichen), ist offen.

### 27.4 Nächste Schritte
1. W12 präzise mit Eck-Normierung; Prüfen, ob die Perron-Hülle der Unterlösungen (mit Eck-Normierung und Singularitätsstärke 1 in 0, ∞) endlich ist.
2. Diskrete Perron-Iteration (Gitter, Mittelwert-Bedingung für Subharmonie, Bilanz mit den drei Urbildern, Ecke als Randbedingung) als
   unabhängiger Test: reproduziert die maximale Unterlösung c ≈ 0.012?
3. Vergleichsprinzip für monotone, additiv homogene Systeme (Shapley/Collatz–Wielandt-Typ): u ≤ S(u) ⊕ b mit S aus Mittelung (Subharmonie) und
   Komposition/Maximum — der einfache Collatz–Wielandt-Schluss braucht eine Super-Lösung mit beschränkter Differenz; die Ecke ist die Stelle, wo
   die Differenz kontrolliert werden muss.

### 27.5 Entscheidung: W12 in tropischer Form trägt die Wachstumsschranke NICHT (Begründung)
**Satz W13 (bewiesen, elementar).** Sei (u, w) eine tropische Unterlösung (Raten von y und T samt Quelle s und Polynomteil p wie in W11), und gebe es
keinen Punkt, an dem das Maximum NUR von einem y- oder T-Term und der Quelle angenommen wird (die Quelle wird nie allein zum Ausgleich gebraucht).
Dann ist (u + δ, w + δ) für jedes δ ≥ 0 wieder eine Unterlösung.
*Beweis.* Wo das Maximum von zwei Nicht-Quell-Termen angenommen wird, bleibt das nach der Verschiebung so; T ist die Summe der übrigen Terme, seine
Rate ist (ohne Auslöschung) gleich deren Maximum, steht also im Maximum mit; wo T durch Auslöschung klein ist, müssen zwei y-Terme gleichauf liegen. ∎
**Befund (numerisch am Modell):** Die Voraussetzung von W13 ist erfüllt — in Ω₁ gleichen y und T aus, im Halbmond und in der Bucht T und y∘h, weit
draussen y gegen seine eigenen Urbild-Terme (alle auf dem Sockel); im Fixpunkt L sind die beiden y-Terme zwar dieselbe Zahl, aber T(L) gleicht aus.
**Folgerung:** Das Grenzsystem aus W11 ist (bis auf Eck-Effekte, die es nicht sieht) VERSCHIEBUNGSINVARIANT. Es kann höchstens die GESTALT und den
ANSTIEG relativ zum Sockel festlegen (W12' = Eindeutigkeit bis auf Konstanten, weiterhin offen), aber nicht die absolute Rate c = lim (1/N)ln σ₁.
Die absolute Normierung — wo der Sockel relativ zur Quelle liegt — ist eine Anregungsaussage: |y_k| ≤ poly für kleine k ((H-I)), gleichwertig mit
„die Eckmoden werden von glatten Daten x^k nicht exponentiell angeregt“, und J_∞⁻¹x^k ist in L schwach singulär (Mellin-Exponent Re s = 5.44), also
ist (H-I) selbst eine (E_L)-Aussage. Ebenso liefert die Identität im Fixpunkt y(L) ~ 𝒰(L) ~ N^{1/2} nur, wenn T(L) polynomial ist — wieder (E_L).
**Damit ist der Teilweg „n-te-Wurzel-Eindeutigkeit ⇒ scharfe Wachstumsschranke“ gescheitert:** Normierung und Ecke lassen sich nicht von der
Gestalt abkoppeln. Die Normierung entsteht in L, wo die Quelle (vorwärts-Eckmode, Breite 1/N, Höhe √N) die adjungierte Mode trägt — genau dort, wo
Os Parametrix sitzt. Folgerung für die Planung: direkt auf der starken Ebene arbeiten (Schweiss-RH-Aufgabe mit Quelle als Inhomogenität, äussere
Parametrix aus H auf ℳ, lokale Parametrix in L mit der Quelle); dort ist die Normierung Teil der Lösung und nicht ein zusätzlicher Eingang.
W11 bleibt als Werkzeug (Gestalt der g-Funktion), W12' (Gestalt bis auf Konstanten) als Nebenziel.

### 27.6 Stand W12 beim Anhalten (11.10., Team wird angehalten)
- **Bewiesen:** W11 (Kompaktheit der Raten, tropische Bilanz der Grenzwerte); 27.3(a,b) (Unterlösungen bilden einen ⊕-Verband; die Hüllkonfiguration
  ist Unterlösung zur Quellrate s_max); W13 (Verschiebungsinvarianz, wenn die Quelle nie allein ausgleicht); 27.1 Umformulierung als konforme
  Schweissung mit Kreisbedingung.
- **Numerisch:** Modell lokal eindeutig (Jacobi-Singulärwerte 3.4 … 0.009, c robust in 2. Ordnung); Voraussetzung von W13 am Modell erfüllt.
- **Vermutet / offen:** W12' (Gestalt und Anstieg bis auf Konstanten eindeutig); ob die Ecke die Normierung tropisch erzwingt (eher nein, T(L) gleicht aus).
- **Gescheitert:** „W12 ⇒ scharfe Wachstumsschranke“ (Normierung = (E_L)-Aussage, §27.5).
- **Offene Ansätze:** (1) starke Ebene direkt: Schweiss-RH-Aufgabe mit Quelle, Parametrix in L zusammen mit O (Anpassungsform BRETT 04:35);
  (2) W12' über ein Vergleichsprinzip für monotone, additiv homogene Systeme (Collatz–Wielandt) bzw. diskrete Perron-Iteration;
  (3) hochgenaue Schweiss-Rechnung (Fourier/Wegmann) statt Ladungssimulation, mit Linearisierung für lokale Eindeutigkeit.
