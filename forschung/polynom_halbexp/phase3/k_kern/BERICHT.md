# Strang K (Kern) — Bericht, Phase 3

Auftrag (Leitung, 11.10.2026 19:50): analytischer Beweisweg für (H₀) bzw. seine Umgehung, neue Mathematik, Spielzeugmodelle mit Randfixpunkt.
Kennzeichnung: **bewiesen** / **Skizze** (Lücken benannt) / **numerisch** / **Vermutung**. Zeiten aus `date`. Bezeichnungen wie ARCHITEKTUR §1;
zusätzlich r_t := J_N⁻¹T_Nk_t ∈ P_N (Antwort des Abschnitts auf den Cauchy-Kern bei t) und **G_N(t) := |t|·‖r_t‖_{ℓ²_ρ} = ‖J̃_N⁻¹ẽ_t‖**
(Os SVD-freie Majorante, O §33.1; ẽ_t = ((ρ/t)^k)_{k≤N}). Es gilt G_N(t) = sup_{‖u‖=1}|Y_u(t)| (alle adjungierten Vektoren, Ecke und Bulk).

## 0. Kurzfassung (Stand 11.10. 20:30; Nachtrag 23:39)
### 0.A ABSCHLUSSNOTIZ PHASE 3 (Strang K, 12.10. 14:10; Konsolidierung bis zur Entscheidung des Autors, keine neuen Versuche)

**Erreicht, bewiesen UND von V gegengelesen** (Liste mit Fundstellen: HAUPTSATZ_KR0.md §7):
- Kurvenkriterium K1/K1′/K1″, Lemma K2, Sätze K3–K5 (Strukturidentität, Kopf-Zerlegung), Lemma K6.
- Explizite Eckprofile K10/K11 mit Integraldarstellung, Lemma K13 (Kopfentwicklung mit beliebigen Konturen), Rang-1-Zuwachs.
- Bootstrap K16 (bedingt), **K19** (Datenseite = klassische Approximation), **K20** (primale Rate ≡ 0 auf D̄), O1, Spiralfenster.
- **K22 = (S_E)** auf dem gelesenen Satz 2 von Andrievskii–Pritsker–Varga, Geometrie von E numerisch (k54, v55).
- **(CD1)**, die exakte Christoffel-Darboux-Identität (§40, V 13:35, auch für Kneser).
- **Satz A** (Logik V) und der bedingte Hauptsatz in der K_{r₀}-Fassung (HAUPTSATZ_KR0.md Fassung 3, Budget γ = 0, Rate N^{−3.96}).

**Numerisch gesichert** (kubisch, meist = Kneser):
- (EP): G_n ~ n^{−3.2}, R_A = 1.008 bis n = 2000.
- (FA): Reserven auf Linse und Γout.
- Am Steg ist das Produkt Kopf × Adjungierte flach und explizit (§39–§40).
- Vs CF-Reserve 0.0024 für die Familie, Linsenfenster S ≈ 0.013–0.015.
- Vs T1: gleichmässige Konvergenz auf D̄ ist wohl falsch, daher K_{r₀}.

**Widerlegt oder zurückgezogen:**
- K17; (EP⁺); F3 in der alten Form; die Steg-Lücke „+0.003“ von (KK) (ein Messartefakt).
- Lemma ρ mit expliziten Parametrizes (§37); der Abschneide-Rest als Hebel (§38).
- Dyson-Kontraktion aus CD1 (§41); lokal integrabler Eckkern (§42).

**Offen (benannte Hypothesen):**
- (EP) = (G-P) + (A-lok); (FA) (⇐ W12′); für die Familie (BU) bzw. (Fam-koh)^{r₀}.
- Gemeinsamer Kern: Kohärenz an exakt ausgeglichenen Stellen. Am Steg ist das Produkt flach (Marge 0). Die Ecke ist echt zweikanalig (§42). Bei (BU) entscheiden die Phasen der Paarungen (§43).

**Nächste Kandidaten (für die Entscheidung des Autors):**
1. K3-4: RH-Problem mit Verschiebung bzw. konformer Verschweissung, mit einer lokalen Parametrix an L, die BEIDE Kanäle enthält (mindestens 3×3). Langfristig.
2. Ein Beweis von (BU) über die Phasen der Paarungen π_n. Vorher misst V die R_P und R_A auf der kritischen Linse S₀ ≈ 0.0125 bis N = 3000 und das gewichtete CF.
3. (FA) aus Ws Grenzproblem W12′.

**Für Fachleute:** phase3/OFFENES_PROBLEM_2.md (deutsch) und OPEN_PROBLEM_2_en.md (englisch), je etwa 2–3 Seiten, in sich geschlossen.

### 0.0 ÜBERGABE an einen frischen K-Lauf (Stand 12.10. 12:50, nachgeführt 13:07 nach §39; ersetzt keine der Notizen darunter, geht ihnen aber vor)

**Ziel.** Beweis des Hauptsatzes in der K_{r₀}-Fassung: p_N → f_K gleichmässig auf K_{r₀} = D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄)), mit Rate. Das massgebliche Dokument ist
phase3/HAUPTSATZ_KR0.md, **Fassung 3** (§1–§7, ca. 200 Zeilen; die Anhänge A/B sind alte Fassungen). Satz A (linear, bei f_K) ist bewiesen bis auf drei Annahmen.
Der Hauptsatz braucht dazu die Familie.

**In 30 min einsteigen (Lesereihenfolge):**
1. HAUPTSATZ_KR0.md §1–§2 (15 min): Bezeichnungen, Kontur 𝒞, E, Satz A, Annahmen (EP) = (G-P) absolut + (A-lok), (FA), (S_E), Familie F1/F2/F3′/F4, Kern (KK) in §2.6.
2. KERNPROBLEM.md §4 (5 min): der gemeinsame Kern (KK) mit Zahlen.
3. Dieser Bericht: §37.4 und §38.4 (Ergebnisse), dann nur die Überschriften §26–§38 (10 min). Ältere Teile nur bei Bedarf.
4. ../BRETT.md, oberste 5 Einträge; PLAN.md, die letzten 3 „Anpassung K“.

**Stand der Annahmen.** (S_E): Satz K22 auf gelesenem Satz (Andrievskii–Pritsker–Varga 2) + numerisch geprüfte Geometrie von E (k54, Vs v55), s_E = 3.96.
(FA): offen, folgt aus W12′. (EP): numerisch klar (V bis n = 2000: G_n ~ n^{−3.2}, R_A = 1.008), offen. Familie: (F3′) = (Fam-koh)^{r₀} offen.
Der gemeinsame Kern ist **(KK′)** (KERNPROBLEM §4, Nachtrag 2; HAUPTSATZ §2.6). Die frühere Fassung (KK) mit „Steg-Lücke +0.003“ war ein Messartefakt (13:07, §39.2).
Am Steg ist das Produkt primal × adjungiert polynomial und explizit: ΔK_n(t,t)·m(t) → κ(t), K_n(t,t)·m(t) = t(log δ_n)′ + Δ(t) + o(1) (§39.3).
Es bleiben: die Komposition der Familien-Daten (P̃ punktweise; Linsenfenster laut V), die Eck-Wachstumszone ausserhalb D nahe L, und für den Parametrix-Weg
zu (EP) die Primale an f(Γout) (≈ ½·ln(|w|/ρ) je Grad, §38).

**Gescheitert, und warum (nicht wiederholen):**
- H-A, Wronski-Paarung auf dem Torus AN DER ECKE (§16): Die Köpfe sind in komplementären Sektoren homogen (primal innen, adjungiert aussen), nirgends beide. Q weicht um 10^{±3} ab.
- Betragsabschätzung auf einer Kontur, auch mit gemischten Datendarstellungen (§34.3): Die Bedingungen an ∇U^ν(L) widersprechen sich in erster Ordnung.
- (FA) per superharmonischer Vergleichsfunktion (§33): O1 (Laplace-Masse der Köpfe), O2 (λ-Naht ⇔ W12′).
- F3 nahe L lokal über (EP⁺) (§37.1): Fenster 1 (Σ₁) ist primal ein Wachstumsfenster.
- Lemma ρ / explizite Parametrizes (§37): Taylor gibt Kompositionsüberlauf, der E-Approximant gibt ρ = O(1). Im Wachstumssektor gilt das innere Profil nur bis |ζ| ≈ 2–3.
- Identität über den Abschneide-Rest R_n (§38): R_n = u_n∘f + m u_n − T_nk_L, auf Γout exponentiell (bis +0.252 je Grad). Sie ist tautologisch.
- Ältere Wege vor Phase 3: KERNPROBLEM.md §2 (Tabelle).

**Nächste Kandidaten.**
1. Erledigt in §39–§40: Das Produkt am Steg ist polynomial und explizit, und es gibt eine exakte Christoffel-Darboux-Identität (CD1), bewiesen in §40.1.
   Offen ist der Dyson-Term Ξ_n, also der Kompositionskanal. Die Dyson-Kontraktion ist geprüft und gescheitert (§41: Spektralradius > 1 ab n ≈ 58).
   Offen als Kandidat: B(t) exakt bestimmen.
   Ideenrunde 3 (IDEEN.md): K3-2 integrabler Eckkern gescheitert (§42: Ecke echt zweikanalig). K3-3 Gerüst für (F3′) steht (§43); offen ist EINE
   Ungleichung (BU), die termweise nur mit Marge 0 auf der kritischen Linse S₀ ≈ 0.0125 schliesst. Robust nur über die Phasen-Kohärenz (Vs 0.0024).
2. (FA) aus W12′ (Ws Randwertproblem).
3. (Fam-koh)^{r₀}: Vs CF-Test (V §3.24): Repräsentant + Adjungierte ≈ −0.0024 je Grad, linsenunabhängig ⇒ Linsenfenster S ∈ (0.0126, 0.0147). Nächster V-Test:
   CF und Adjungierte bei S ≈ 0.0137 auf der ganzen Linse; danach das gewichtete CF. Je Kopf genügt der Keim mod x^{n+1} (§39.4).

**Werkzeuge und Fallen.**
- Python: C:/Users/janis/miniconda3/python.exe. Skripte mit dem Write-Werkzeug anlegen, weil Bash-Heredocs Backslashes verlieren.
- k_kern/k_kern.py: J_tilde, modell, kubik_koeff.
- Exakte Eck-Pol-Antwort u_n wie in k55c/k56: J̃ solve mit bits = 1.3n + 400. Bei arb → mpmath immer .mid().str(nd, radius=False) verwenden.
- k54: kontur(u1, u2, s, r1, m) liefert die Teile von 𝒞.
- flint: ctx.cap = 10 für Reihen, Präzision mit N wachsen lassen.
- Läufe auf dem Laptop ≤ 20 min, längere über V auf dem Server. Nur V benutzt den Kneser-Motor.
- Kein git. Schreiben nur in phase3/ (gemeinsame Dateien, k_kern/) und auf dem BRETT.
- BRETT-Kopf `## [Strang K] dd.mm. HH:MM — Titel`, neueste oben, direkt vor dem ersten `## [Strang` nach „kennzeichnen.“.
- Zeiten immer aus `date`. Kennzeichnung bewiesen / Skizze / numerisch / Vermutung. Ein Satz zählt als Grundlage erst, wenn V ihn gelesen hat.

**Korrekturen, die man kennen muss.**
- K17 ist zurückgezogen (§29).
- Die Primale ist auf ∂D nicht exponentiell klein (K20, §31).
- Fenster 1 ist kein primales Ruhefenster (§37.1).
- E ⊃ D̄ ist falsch, richtig ist E ⊃ K_{r₀}. Die Öffnung ist 131.7°, nicht 99°.
- F3 in der Form „P nahe f_K auf 𝒞“ ist falsch (V §3.20).
- Die Schlitzkonstante in Fassung 1 ist 8.2, nicht 16.3.

**Für den Autor (12.10. 09:10) — welche Fassung wir jetzt anstreben und warum:** Wir wollen beweisen, dass p_N gleichmässig auf der Kreisscheibe |x| ≤ |L| ohne
zwei kleine, feste Kreise um die Fixpunkte L und L̄ gegen Knesers halbe e-Funktion f_K konvergiert (K_{r₀}-Fassung, z. B. r₀ = 10⁻³). Daraus folgt alles, was die
ursprüngliche Vermutung verlangt: jeder Koeffizient a_k^{(N)} geht gegen Knesers Koeffizienten, und p_N(1) → T_e(½). Die schärfere Fassung „gleichmässig bis an L
heran“ verfolgen wir nicht weiter, weil sich das Eckmuster bei L mit wachsendem N langsam dreht (komplexer Exponent κ₁) und eine gleichmässige Randaussage
dadurch für sehr grosse N unbeweisbar und womöglich falsch ist (Prüfung läuft, Test T1); die K_{r₀}-Fassung ist davon nicht betroffen und braucht von der Ecke nur
noch eine Zahl je Grad (das Eckprodukt) plus lokale, explizite Eckprofile.
**Für den Autor (12.10. 13:07):** Eine Korrektur zu meiner Notiz von 12:45. Ich hatte zwei Zahlen verglichen, die auf verschiedene Art gemessen waren.
Richtig gemessen wird knapp ausserhalb der Scheibe die eine Hälfte der Rechnung genau so schnell gross, wie die andere klein wird. Das Produkt bleibt harmlos und
folgt sogar einer einfachen Formel. V hat dazu passend gefunden, dass die Nachbarpolynome dort mit etwas Spielraum gut genug sind, wenn man die „Linse“ richtig dick wählt.
Offen bleiben vor allem zwei Dinge: Die Nachbarpolynome müssen dort punktweise nahe an Knesers Funktion bleiben, und nahe beim Punkt L gibt es eine Zone, in der wirklich etwas wächst.

**Für den Autor (12.10. 12:45; Korrektur 13:07: der zweite Teil stimmt so nicht mehr, siehe die Notiz darüber):** Wir haben vier Wege probiert, die erste Annahme (EP) zu beweisen. Alle vier scheitern an derselben Stelle: Knapp ausserhalb der
Kreisscheibe wird die eine Hälfte der Rechnung exponentiell gross, und die andere wird fast, aber nicht ganz so schnell klein. Das Produkt bleibt in Wirklichkeit
harmlos, aber nur weil sich Vorzeichen gegenseitig aufheben, nicht weil einzelne Teile klein sind. Ein Beweis muss dieses Aufheben direkt zeigen.
Das steht jetzt kurz in phase3/KERNPROBLEM.md §4. (EP) bleibt vorerst eine ausdrücklich genannte Annahme.

**Für den Autor (12.10. 11:05; ersetzt die Notizen von 10:01 und 10:23):** Das Ergebnis steht in phase3/HAUPTSATZ_KR0.md (Fassung 2): Wenn drei Aussagen über das lineare
Problem stimmen — zwei Antwortfolgen an der Ecke wachsen höchstens wie eine Potenz (EP), die adjungierten Köpfe sind fern der Ecke exponentiell kleiner als in der Ecke (FA),
und Knesers Funktion lässt sich auf einem etwas grösseren Gebiet mit Polynomen annähern (S_E) — und dasselbe für die Nachbarpolynome im Newton-Schritt gilt, dann konvergieren
Deine Polynome auf der Scheibe bis |L| ausser in zwei kleinen Kreisen um L, L̄. Alle drei Aussagen sind gemessen bzw. in Messung, bewiesen sind sie noch nicht.
**Nachtrag zweiter K-Lauf (12.10. 08:44; §18–§21):** Die Grenzaufgabe K → ∞ ist gelöst: adjungierter Kopf → F_{β_A}(λ_Aζ), primaler → F_{β_P}(λ_Pζ),
F_β(u) = Σ Γ(β − iαk)(−u)^k/(k!Γ(β)), α = log μ/2π, β_A = ½ − 2iα, β_P = ½ + iα, λ_P = −μ^{−1/2}λ_A (numerisch auf 1e-7 bzw. 6 Stellen; Defekt-Tabellen und
Konvergenz n^{−a} auch im Anpassungsbereich); Lemma K10.1, Satz K10 (Existenz), K11/K11′ (Eindeutigkeit unter Mellin-Randbedingung) bewiesen. Die explizite
Lösung zeigt: komplexes κ₁ ⇒ Eckmuster dreht sich um 1.45°·ln N gegen die feste Geometrie ⇒ tangentiale Linse versagt (von V bestätigt: τ_N steigt ab 1000),
Keil trägt bis N_×(φ), mitgeführter Keil bis ~10¹⁶; linsenunabhängig wird die primale Kreisnorm ab ~10⁸ nicht-lokal. (H1) in Keil-Geometrie reduziert sich auf
B_k beschränkt + Anpassungsfehler + Dämpfung des Stegs; Bessel nur polynomial (primale Köpfe kohärent).
0. **Nachtrag 23:39 (§9):** A2a entschieden — Eck-Korrektur D_N universell (hängt nur von μ und der Linse bei L ab; kubisch = Kneser auf 4e-4 im führenden Paar,
   flache Störung dort ≤ 2e-5), ein eigenständiges Eckmodell mit Gradabschnitt gibt es nicht (affines Modell: J nicht injektiv, Rate 0.70 exakt erklärt).
   A2b entschieden — Satz K3: Ws Tropik allein liefert die Schleifenbedingung in keiner Koordinate; (SB) ⇐ W12′ + Eck-Normierung u(L) = 0.
   Vermutung S′: (B) + (SB) ⇒ (S_γ). Familie: Satz K4 (P'∘P = (eˣ + R_P')/P') — Newton-Matrix an p_N und Π_N f_K so stabil wie J_N, exponentiell instabil bei
   T_N f_K (‖Ĵ⁻¹‖ = 2.3e25 bei N = 600); Literaturabgleich §9.4.
   **Nachtrag 23:59 (§10):** Linsen-Spannung entschieden — Ĵ_N^P hängt nur vom Keim von P ab, K1″(b, c) gelten wörtlich für die Familie; die
   Abweichung von P auf der Linse spielt für (S_γ) keine Rolle (S = 0.1, N = 600: 8.1e9 auf der Linse, ‖Ĵ⁻¹‖ = 0.965153 wie bei f_K). P auf γ
   nur in Daten-Schlüssen (H-ii)/H1, dort Kosten e^{NF₀}·Kopf (Carathéodory–Fejér), Schwelle N ≈ 6300–7800 inkohärent; Rest: (Fam-koh).
   **Nachtrag 12.10. 00:25 (§11):** A2c Schritt 1 „ohne Ecke“ endgültig gescheitert — L ∈ γ, also enthält n-te-Wurzel-(S_γ) die linsenfreie
   Aussage G_N(L) = e^{o(N)}; fern von L reduziert Lemma K6 auf das Grad-Antwortprofil (Rate 0 für ξ ≤ 0.8, exponentiell nur im Kantenband 0.85–0.97).
   Satz K5 (exakte Kopf-/LU-Zerlegung): Ĵ_N⁻¹ = Σ_{n≤N} Köpfe, U unabhängig von k, D_N = projizierte künftige Köpfe; n-te Wurzel ⇔ Köpfe
   subexponentiell (gemessen: w-Köpfe sättigen bei 0.015 bis N = 1000, x-Köpfe ~ e^{0.011n}). Wendung: A2c um die Ecke organisieren (E₀/E₁ in Kopfform).
1. **Satz K1 (bewiesen, elementar): (H₀) folgt aus einer polynomialen Schranke für G_N auf EINER festen Kurve γ = Γout ∪ Linse(Γin) ∪ konj.**
   Der Beweis ist eine Cauchy-Darstellung der Abschnittslösung — er braucht **keine Eckformel D6, kein (E_sing), kein (a′)_r, keine Eckprojektoren
   Q_N, keinen Kantenschwanz** und gilt wörtlich für Knesers f_K (nur Analytizität von f_K, m auf D̄ ∪ Linse und Stetigkeit bis L). Rahmen ist
   K_γ = K' ∪ Linse; jedes κ < 2.79 in G_N ≤ CN^κ genügt. Mit Gs Summenkriterium G10 reicht sogar **sup_γ|𝒟_n| ≤ C für alle n** (κ = 1).
   Damit ist Vs Lücke B6.2 (W1'' nur für polynomiales f) geschlossen, und der Engpass hat eine einzige, messbare Form: **(H*_G) G_N ≤ CN^κ auf γ.**
2. **Numerik k01 (kubisch, arb):** G_N ist auf Γout konstant (N = 100…1200, keine Rate), auf der Linse konstant fern von L und höchstens ~√N nahe L;
   auf dem Γin-Kreis selbst wächst G exponentiell (Kontrolle). Siehe §1.5.
3. Kompositions-Vorkonditionierung (Lemma K2, bewiesen): Linksmultiplikation mit T_N C_k (k(0) = 0) vertauscht mit dem Abschneiden, lässt Π_N
   unverändert und wirkt auf Grad-n-Anteile wie das Phasenraumgewicht e^{n·ε Re φ}. (H*_G) auf k(∂D) ⇔ polynomiale Stabilität von T_N C_k J_∞|_{P_N}.

## 1. Satz K1 — das Kurvenkriterium

### 1.1 Lemma K1.0 (Cauchy-Darstellung der Abschnittslösung) — bewiesen
Sei γ eine rektifizierbare Jordankurve mit 0 ∈ Int γ und W analytisch in Int γ, stetig auf Int γ‾. Dann
  **J_N⁻¹ T_N W = (1/2πi) ∮_γ W(t) r_t dt**,  also  **‖J_N⁻¹T_NW‖_{ℓ²_ρ} ≤ A_γ · max_γ G_N · max_γ |W|**,  A_γ := (1/2π)∮_γ|dt|/|t|.
*Beweis.* Cauchy-Formel W(x) = (1/2πi)∮_γ W(t)/(t − x) dt für x ∈ Int γ (für W ∈ A(Int γ) ∩ C(Int γ‾) und rektifizierbares γ gültig). Die
Taylor-Koeffizienten bei 0 sind [x^k]W = (1/2πi)∮_γ W(t) t^{−k−1} dt, also T_N W = (1/2πi)∮_γ W(t)·T_Nk_t dt (endlichdimensionales
vektorwertiges Integral). Lineare Abbildung J_N⁻¹ anwenden; Dreiecksungleichung für Integrale mit ‖r_t‖ = G_N(t)/|t|. ∎

### 1.2 Satz K1 — Aussage
**Voraussetzungen.**
- (γ) γ = Γout ∪ ∂ℬ ∪ (Spiegelbild), wobei ℬ eine „Linse“ über Γin ist: ℬ = {ρe^{s}e^{iφ} : φ ∈ Γin, 0 ≤ s ≤ s(φ)}, s stetig, s > 0 im Inneren
  von Γin, s(φ) = O(|φ − arg L|^p) mit p > 1 an L (ebenso an L̄). Also Int γ‾ = D̄ ∪ ℬ̄, und γ ⊂ {|t| ≥ ρ}.
- (A) f, m = f'∘f sind analytisch auf einer Umgebung von D̄ ∪ ℬ̄ ohne {L, L̄} und stetig bis L, L̄ (für f_K: (K1)/(F) + (K3) f_K(L) = L,
  f_K'(L) = √L; die Schlitze S_L, S_L̄ treffen D̄ ∪ ℬ̄ nicht und auch f(D̄ ∪ ℬ̄) nicht, siehe 1.4).
- (Geo) K_γ := (D̄ ∪ ℬ̄) ∪ f(D̄ ∪ ℬ̄); Ĉ∖K_γ zusammenhängend, 0 und a₀ = f(0) in int K_γ; g_γ := Green-Funktion von Ĉ∖K_γ (Pol ∞) ist
  Lipschitz bis ∂K_γ (gilt für stückweise Dini-glatten Rand mit konvexen Ecken); J_N invertierbar.
- (H*_G) sup_{t∈γ} G_N(t) ≤ C_G N^κ.

**Aussage.** Für η ∉ K_γ und z ∈ D̄:
  (i) |Π_N k_η(z)| ≤ √(N+1)·A_γ·(1 + M_γ)·C_G N^κ / dist(η, K_γ),  M_γ := max_{D̄∪ℬ̄}|m|;
  (ii) für η ∉ K_γ^{1/N} (1/N-Umgebung): |E_N(η, z)| ≤ C₁ N^{κ+3/2} e^{−(N+2) g_γ(η)} + e^{−(N+1)g_γ(η)}/(|η| − ρ),
also (H₀^{K_γ}) in der Form von Satz P-I mit ε = 1 und Exponent κ + 3/2.

### 1.3 Beweis
(i) W := J_∞k_η, W(x) = 1/(η − f(x)) + m(x)/(η − x). Für η ∉ K_γ ist W analytisch auf Int γ und stetig auf Int γ‾ (Nenner ≠ 0, da f(x) ∈ K_γ
und x ∈ K_γ; (A)). Auf γ ist |W| ≤ (1 + M_γ)/dist(η, K_γ). Lemma K1.0 mit Π_Nk_η = J_N⁻¹T_NW, dann sup_{D̄}|P| ≤ √(N+1)‖P‖_{ℓ²_ρ} für P ∈ P_N.
(ii) y_η := Π_Nk_η − T_Nk_η ist koordinatenweise rational in η mit Polen nur in 0 und a₀ und y_η = O(η^{−N−2}) bei ∞ (C Lemma 12.1; gilt für jedes
analytische f, auch f_K: [x^k]1/(η − f(x)) hat nur den Pol a₀, [x^k] m(x)/(η − x) nur den Pol 0). Auf ∂K_γ^{1/N} ist nach (i) sup_{D̄}|y_η| ≤
√(N+1)A_γ(1+M_γ)C_G N^{κ+1} + (N+1)/ρ =: C₂N^{κ+3/2}. Bernstein–Walsh (C Satz 12.2) in Ĉ∖K_γ^{1/N}: sup_{D̄}|y_η| ≤ C₂N^{κ+3/2}e^{−(N+2)g_{1/N}(η)};
wegen g_{1/N} ≥ g_γ − sup_{∂K_γ^{1/N}} g_γ ≥ g_γ − Lip/N folgt (ii) für y_η. Schliesslich E_N = y_η − [k_η]_{>N} und |[k_η]_{>N}(z)| ≤
(ρ/|η|)^{N+1}/(|η| − ρ) ≤ e^{−(N+1)g_γ(η)}/(|η| − ρ), da K_γ ⊃ D̄ ⇒ g_γ(η) ≤ log(|η|/ρ). ∎

**Folgerung (LP).** Gilt zusätzlich (F_{K_γ}) (f_K holomorph um K_γ bis auf S_L, S_L̄, Sprung ≤ C_J u^{5.447}) und g_γ(L(1+u)) ≥ c u^{1/α'} mit
dem K'-Eckwinkel α' = 0.787172, so gibt Satz P-I(i): sup_{D̄}|Π_N f_K − f_K| = O(N^{κ + 3/2 − 5.447α'}) = **O(N^{κ − 2.79})**. Mit dem gemessenen
κ ≈ ½ also O(N^{−2.29}) (grob; wahr ist N^{−4.8}). Der Exponent 3/2 ist Buchhaltung (ℓ²→sup und 1/dist) und lässt sich verbessern.

**Mit Gs Summenkriterium (G10, endliche Identitäten, gilt für f_K):** G_N(t) ≤ G_{n₀−1}(t) + Σ_{n₀≤n≤N}|𝒟_n(t)|. Also genügt
**sup_{t∈γ}|𝒟_n(t)| ≤ C_Δ für alle n ≥ n₀** (dann κ = 1 < 2.79), erst recht C_Δ n^{−1/2} (κ = ½).

### 1.4 Geometrie für f_K (Skizze mit Begründung; Zahlen aus vorhandenen Dateien)
- Lage der Schlitze: in L hat ∂D die äussere Normale n₁ (Richtung arg L = 76.6°), f(D̄) lokal ≈ L + μ(D̄ − L) die äussere Normale in Richtung
  76.6° + arg μ = 114.9°. Der Schlitz L(1+u) zeigt in Richtung n₁: Komponente längs der f(D̄)-Normalen cos 38.3° = 0.785 > 0 ⇒ ausserhalb f(D̄);
  ausserhalb D̄ trivial; die Linse ist in L tangential (Dicke O(τ^p)) ⇒ trifft den Schlitz nicht. dist(L(1+u), K_γ) ≥ 0.78·ρu (lokal).
- f(ℬ) ⊂ D nahe L: für x in der Linse ist x − L ≈ τ·(Tangente 166.6°), μ(x − L) zeigt nach 204.9°, Komponente längs der inneren Normalen (256.6°)
  cos 51.7° > 0 ⇒ nach innen; fern von L ist |f| < ρ auf Γin mit Abstand und die Linse dünn. Also K_γ = K' ∪ ℬ̄.
- m = f'∘f stetig bis L: f(D̄ ∪ ℬ̄) trifft den Schlitz von f' nicht (gleiche Rechnung), Randwert f'(L) = √L (K3).
- Eckwinkel von K_γ in L = Eckwinkel von K' (141.69°, P §9.3), denn die Linse ist tangential mit Dini-stetigem Tangentenwinkel (τ^{p−1}).
- Offen (numerisch zu bestätigen, Vs Liste B4.3): Lipschitz-Konstante von g_γ, Univalenz-/Abstandsfakten wie in P-II.

### 1.5 Numerik k01 (kubisches Modell, arb-Lösung, float-Normen; `k01_G_kurve.py`, Logs `k01_kubisch_*.log`)
Linse: s(τ) = Λ(τ)/x₀ mit Ws Kreisprofil Λ (W §18), τ = Winkel ab arg L Richtung Γin. ln G_N:

| Punkt | N = 100 | 200 | 300 | 400 | 600 | Rate 400→600 |
|---|---|---|---|---|---|---|
| Γout 0° / 30° / 60° / 70° | −0.597 / −0.466 / 0.074 / 0.574 | gleich | gleich | gleich | gleich | 0 |
| Γout 76° (0.6° vor L) | 1.392 | 1.610 | 1.707 | 1.755 | 1.790 | +0.0002 |
| L | 1.491 | 1.795 | 1.974 | 2.102 | 2.282 | +0.0009 (≈ ½ ln N) |
| Linse x₀ = 0.7, τ = 2° / 4° | 1.706 / 1.782 | 2.105 / 2.081 | 2.322 / 2.196 | 2.455 / 2.251 | 2.608 / 2.295 | +0.0008 / +0.0002 |
| Linse x₀ = 0.7, τ ≥ 15° | 1.84…1.99 | 1.93…2.02 | 1.95…2.02 | 1.95…2.02 | 1.95…2.02 | 0 |
| Linse x₀ = 0.95 (dünner), τ = 2° / ≥ 15° | 1.719 / 1.94…2.11 | | | 2.505 / 2.11…2.18 | 2.680 / 2.12…2.18 | +0.0009 / 0 |
| Kontrolle: Γin-Kreis 120° / 180° | 2.523 / 2.585 | 2.881 / 2.938 | 3.107 / 3.162 | 3.343 / 3.419 | 4.520 / 4.925 | +0.0059 / +0.0075 |

Fortsetzung N = 800 / 1000 / 1200 (Rechenzeit 55 / 131 / 251 s):

| Punkt | 800 | 1000 | 1200 | Rate 1000→1200 |
|---|---|---|---|---|
| Γout 0°…74° | wie oben, auf 3 Stellen unverändert | | | 0.00000 |
| Γout 76° | 1.795 | 1.794 | 1.792 | −0.00001 |
| L | 2.410 | 2.509 | 2.591 | +0.00041 (≈ N^{0.45}) |
| Linse x₀ = 0.7, τ = 2° / 4° / ≥ 8° | 2.694 / 2.308 / 1.95…2.02 | 2.749 / 2.311 / gleich | 2.787 / 2.313 / gleich | +0.00019 / 0.00001 / 0 |
| Linse x₀ = 0.95, τ = 2° / 4° | 2.785 / 2.472 | 2.859 / 2.488 | 2.916 / 2.500 | +0.00029 / +0.00006 |
| Kontrolle Γin-Kreis 120° / 180° | 6.419 / 7.148 | 8.409 / 9.432 | 10.419 / 11.735 | +0.0100 / +0.0115 |

⇒ **(H*_G) gilt numerisch auf γ bis N = 1200 mit κ ≈ 0.45** (Maximum in L, ≈ 0.39√N wie O §33.2); auf fast ganz γ ist G_N sogar **konstant in N**
(die Abschnittsantwort r_t konvergiert — vermutlich gegen s_t = J_∞⁻¹k_t, siehe §2). Auf dem Kreis über Γin wächst G mit der Eckrate (0.0115 → 0.012).
Beide Linsen (x₀ = 0.7 und die dünnere x₀ = 0.95) tragen.

**Kneser** (`kneser_z` = Ps zertifizierte Koeffizienten-Mittelpunkte, Radien nicht propagiert; `k01_kneser_200_600.log`): N = 200 / 400 / 600 —
jede Zelle gleich der kubischen auf ≤ 0.01 (Γout 0°…70°: −0.597 / −0.466 / 0.075 / 0.574 konstant; L: 1.795 / 2.102 / 2.282; Linse τ ≥ 15°:
konstant 1.95…2.17; Γin-Kreis 180°: 2.943 / 3.419 / 4.898, Rate 0.0074). **(H*_G) gilt für Kneser numerisch genauso.**

## 2. Nebenbefunde und verworfene Wege (11.10. 20:34)
- **Identität (bewiesen):** r_t = Π_N s_t mit s_t := J_∞⁻¹k_t (existiert für t ∉ D̄ nach O7: J_∞⁻¹: H²(D) → E²(int K')), denn J_N⁻¹T_NJ_∞s_t = J_N⁻¹T_Nk_t.
  Also ist (H*_G) bei t selbst eine (LP)-artige Aussage für die Eingabe s_t (die bei L die Kneser-artige Singularität (x−L)^{s₀}, Re s₀ = 5.447, trägt).
- **Bootstrap verworfen:** Der Versuch sup_{γ₂}G ≤ A + ε_N sup_{γ₁}G (Satz K1 auf die Eingaben s_t angewandt) scheitert: s_t hat einen Pol bei t,
  dessen Antwort Π_Nk_t/m(t) dieselbe Grösse wie G hat (für t nahe K_{γ₁}), und G auf einer inneren Kurve ist gegenüber einer äusseren nur mit
  e^{N·(Green-Abstand)} kontrolliert. Keine Kontraktion.
- **Multiplikative Linksgewichte ω = e^{Nφ} (BRETT 20:08) verworfen:** T_N[e^{Nφ}·w] verstärkt glatte Anteile der Eingabe bei arg ≈ π ebenso wie die
  Eckmoden; nur zusammen mit Phasenraumlokalisierung brauchbar — die liefert die Komposition C_k (IDEEN K3), nicht die Multiplikation.

## 3. Umformulierung von (H*_G): Stabilität der Abschnitte in der konformen Gleichungskoordinate (11.10. 20:54)

**Lemma K1′ (bewiesen).** Sei k: 𝔻 → Int γ konform mit k(0) = 0 (γ = Γout ∪ ∂ℬ ∪ konj. ist in L, L̄ C¹ und Dini-glatt, also k' stetig und ≠ 0
auf ∂𝔻), und C_k^{(N)}: P_N → P_N^{(w)}, b ↦ T_N^w[b∘k] (Taylor in w bis Grad N; untere Dreiecksmatrix, weil k(0) = 0 die Filtration erhält).
Setze Ĵ_N := C_k^{(N)} J_N : ℓ²_ρ(x) → ℓ²(w). Dann
  **‖Ĵ_N⁻¹‖ ≤ A_γ √(N+1) · max_γ G_N   und   max_γ G_N ≤ C_k √(N+1) · ‖Ĵ_N⁻¹‖.**
*Beweis.* (⇒) Für b ∈ P_N und jedes F ∈ A(Int γ) ∩ C mit F ≡ b mod x^{N+1}A(Int γ) ist J_N⁻¹b = J_N⁻¹T_NF, also nach Lemma K1.0 ‖J_N⁻¹b‖ ≤
A_γ max_γG·sup_γ|F|. Wegen x^{N+1}A(Int γ) = w^{N+1}A(𝔻) (k(w)/w ist eine Einheit) ist inf_F sup_γ|F| = dist_{L^∞(∂𝔻)}(b∘k, w^{N+1}H^∞)
= ‖L_N(C_k^{(N)}b)‖ (Carathéodory–Fejér/Schur: Norm der unteren Dreiecks-Toeplitz-Matrix der ersten N+1 Koeffizienten) ≤ √(N+1)‖C_k^{(N)}b‖_{ℓ²}.
(⇐) r_t = Ĵ_N⁻¹C_k^{(N)}T_Nk_t und C_k^{(N)}T_Nk_t = T_N^w[k_t∘k]; für t = k(w_t) ∈ γ hat 1/(t − k(w)) einen einfachen Pol auf ∂𝔻, also sind die
ersten N+1 Koeffizienten ≤ C(|k'(w_t)|⁻¹)√(N+1) in ℓ². ∎
**Bedeutung.** (H*_G) ⇔ (bis auf √N) die quadratischen Abschnitte werden **polynomial stabil, sobald man die Gleichungen in der konformen
Koordinate des um die Linse vergrösserten Kreises schreibt** (Unbekannte unverändert in x). Die Eckmoden mit σ ~ e^{cN} sind also ein Artefakt
der Gleichungskoordinate x: die Gleichungsseite ihrer Störung (Steg, Index 0.95N, Winkel π) wird in der w-Koordinate um e^{0.95N·s(π)} ≥ e^{cN}
verstärkt. Das ist die intrinsische Form der Kompositions-Vorkonditionierung (IDEEN K3, k = konforme Abbildung) und die natürliche Ausgangslage für
eine **semiklassische Parametrix** (h = 1/N) von Ĵ_N: fern von L elliptisch, an L ein Eckmodell mit √N-Verlust. Anders als bei einer Parametrix
für die Eckmode selbst (W §27.5: Amplitude global) muss eine Parametrix für Ĵ_N⁻¹ die Eck-Amplitude nicht kennen — sie ergibt sich aus der
Neumann-Reihe.
**Numerische Nebenbefunde (k03):** Fern von L konvergiert r_t^{(N)} koeffizientenweise gegen s_t = J_∞⁻¹k_t mit Fehler ≈ (5…60)·‖Schwanz von s_t‖,
und der Schwanz fällt wie k^{−7} (die Kneser-artige Singularität von s_t bei L) — die Konstanz von G fern von L ist also selbst eine (H₀)-artige
Aussage mit kleinem Verstärkungsfaktor, kein exponentieller Effekt. Nahe L (t = L, Linse 2°) konvergiert r_t nicht (s_L ∉ H²), G ~ √N.
Die Eck-Amplitude ist schwach global (k02: jet±0.01 ändert G(L)/√N um ±0.3 %).

## 4. Stand und nächste Schritte (11.10. 20:54)
- Satz K1 von V verifiziert (v_pruefstand/BERICHT §1.3); Zielaussage der Phase ist (H*_G) (Plan A1 erledigt). Nötige Stärke: κ < 2.79 für (LP),
  **κ < 1.79 mit Familie über X4a** (V-Korrektur), gemessen κ ≈ 0.45; die Buchhaltungs-3/2 in Satz K1 sind zu drücken (A2/A2').
- Nächster Block (A2/A2'): semiklassische Struktur von Ĵ_N = C_k^{(N)}J_N — (i) Hauptsymbol fern von L und Elliptizität in der w-Norm,
  (ii) Eckmodell an L (Os (P1)–(P3)) als lokale Inverse mit √N-Verlust, (iii) Klebefehler; V misst parallel das Möbius-Modell.
- Dateien: k_kern.py (Modelle, J̃, Cauchy-Seiten), k01_G_kurve.py (+ Logs), k02_ecke_lokal.py (k02.log), k03_konvergenz_rt.py (k03.log).

## 5. Lemma K1′ entschieden: die Abschnitte sind in der konformen Gleichungskoordinate GLEICHMÄSSIG stabil (numerisch) (11.10. 21:03)

### 5.1 Präzise Fassung (für V; alle Teile bewiesen, elementar)
**Daten.** γ = Γout ∪ ∂ℬ ∪ konj., Linse ℬ mit Polarradius ρe^{s(α)} über Γin, s ∈ C¹ ausser in L, L̄, s(α) ≍ |α − arg L|^p, p > 1 (Dini-glatt).
k: 𝔻 → Int γ konform, k(0) = 0, k'(0) > 0; k ∈ C¹(𝔻̄), k' ≠ 0 auf 𝔻̄ (Kellogg–Warschawski). Schreibe k(w) = ρ w E(w), E = e^F, F analytisch.
**Räume.** Unbekannte: x ∈ P_N mit ‖x‖ := ‖x‖_{H²(D_ρ)} = (Σ|x_j|²ρ^{2j})^{1/2}. Gleichungsseite: b ∈ P_N ↦ C_k^{(N)}b := ([w^i](b∘k))_{i≤N} ∈ ℂ^{N+1}
mit der ℓ²-Norm; in ρ-skalierten Koeffizienten ist C_k^{(N)}[i][j] = [w^{i−j}]E^j (untere Dreiecksmatrix, Diagonale E(0)^j).
**Operator.** Ĵ_N := C_k^{(N)}J_N : (P_N, ‖·‖_{H²(D_ρ)}) → (ℂ^{N+1}, ℓ²). Es gilt Ĵ_Nx = T_N^w[(J_∞x)∘k] (Filtration: T_N^w C_k = T_N^w C_k T_N).
**Lemma K1″.**
 (a) Für η ∉ K_γ = Int γ‾ ∪ f(Int γ‾): ‖Π_Nk_η‖_{H²(D_ρ)} ≤ ‖Ĵ_N⁻¹‖·sup_γ|J_∞k_η| ≤ ‖Ĵ_N⁻¹‖·(1 + M_γ)/dist(η, K_γ).
     *Beweis:* Π_Nk_η = J_N⁻¹T_NW = Ĵ_N⁻¹T_N^w[W∘k] (W = J_∞k_η analytisch auf Int γ, stetig bis γ), ‖T_N^w[W∘k]‖_{ℓ²} ≤ ‖W∘k‖_{H²(𝔻)} ≤ sup_{∂𝔻}|W∘k|. ∎
 (b) ‖Ĵ_N⁻¹‖ ≤ c_k · sup_{‖u‖=1} ‖Y_u‖_{L²(γ,|dt|)} ≤ c_k ‖G_N‖_{L²(γ,|dt|)},  c_k = (sup_{∂𝔻}|k'|/2π)^{1/2}/min_γ|t|.
     *Beweis:* ⟨c, Ĵ_N^{−T}u⟩ = ⟨(C_k^{(N)})⁻¹c, y_u⟩ = (1/2πi)∮_{|x|=ρ}(c∘k⁻¹)(x)Y_u(x)dx/x = (1/2πi)∮_{∂𝔻}c(w)·Y_u(k(w))k'(w)/k(w) dw
     (Kontur von |x| = ρ nach γ verschoben: c∘k⁻¹ analytisch in Int γ, Y_u analytisch in ℂ∖{0}); Cauchy–Schwarz auf ∂𝔻, |dt| = |k'||dw|. ∎
 (c) G_N(t) ≤ |t|·Q_k·√(N+1)·‖Ĵ_N⁻¹‖ für t ∈ γ, Q_k = sup_{w,w'∈𝔻̄}|w − w'|/|k(w) − k(w')|.
     *Beweis:* r_t = Ĵ_N⁻¹T_N^w[k_t∘k], k_t∘k(w) = q(w)/(w_t − w), |q| ≤ Q_k, also T_N^w[k_t∘k] = L_N(q)·T_N^w[1/(w_t − ·)], ‖L_N(q)‖ ≤ ‖q‖_∞. ∎
**Folgerung (verbessertes Satz K1).** Aus (a) wie in §1.3: |E_N(η, z)| ≤ C·‖Ĵ_N⁻¹‖·N^{3/2}·e^{−(N+2)g_γ(η)} (√N von ℓ²→sup auf D̄, N vom Abstand 1/N in
Bernstein–Walsh). **Ist ‖Ĵ_N⁻¹‖ gleichmässig beschränkt, so ist der Eck-Polynomfaktor 3/2 ≤ 2** (Ziel der Leitung für N₀ ≤ 400), ohne den √N-Verlust
von G: der stammt nur aus der Eingabe T_N^w[k_t∘k] (Pol auf ∂𝔻), nicht aus dem Operator.
**Zielaussage, neu und schärfer als (H*_G):** **(S_γ) sup_N ‖Ĵ_N⁻¹‖ < ∞** — gleichmässige Stabilität der Abschnittsfolge Ĵ_N = T_N^w C_k J_∞|_{P_N}.

### 5.2 Messung (k04_konform.py, kubisch; Linse s = 0.0166·(1 − e^{−(τ/12°)^{1.3}}), Theodorsen M = 2^15, F exakt als Polynom in arb, 1.2N + 300 Bit)
| N | 100 | 200 | 300 | 400 | 600 | 800 |
|---|---|---|---|---|---|---|
| ‖J̃_N⁻¹‖ (ℓ²_ρ, Eckmode) | 1.366 | 1.852 | 3.698 | 9.10 | 70.71 | 619.1 |
| **‖Ĵ_N⁻¹‖** (konform) | **1.251** | **1.251** | **1.251** | **1.251** | **1.251** | **1.251** |
| 2./3. Singulärwert von Ĵ_N⁻¹ | 1.202 / 1.155 | gleich | gleich | gleich | gleich | gleich |
Der Vektor grösster Verstärkung ist ein Bulk-Vektor (Index-Schwerpunkt 0.01–0.05 N, Maximum bei 180°), keine Eckmode. Randkurve k(∂𝔻): min|k|/ρ = 1.000000,
max 1.0167, passiert L im Gitterabstand. ⇒ **Lemma K1′ trägt, und zwar ohne jeden √N-Verlust: der Engpass ist eine gleichmässige Stabilitätsaussage
für eine neue Matrixfamilie.** (N = 800: 121 s Aufbau + 367 s Inversen, 1260 Bit. Kontrolle mit zu dünner Linse und Kneser folgen.)

### 5.3 Kontrollen (k04b_konform_S.py: Linsenhöhe S variabel; zusätzlich der Eckanteil ‖Ĵ_N^{−T}u_i‖ für die grössten Moden u_i von J̃_N⁻¹)
| Lauf | N | ‖Ĵ_N⁻¹‖ | σ₁(J̃_N⁻¹) | Eckanteil ‖Ĵ^{−T}u₁‖ |
|---|---|---|---|---|
| kubisch, S = 0.010 (unter c/0.95 ≈ 0.0126) | 300 / 500 / 700 | 1.288 / 1.288 / 1.288 | 3.7 / 24.8 / 207 | 0.461 / 0.546 / 0.637 (Rate +0.0008/Grad) |
| **Kneser**, S = 0.0166 | 300 / 450 / 600 | **1.250 / 1.250 / 1.250** | 3.65 / 14.5 / 68.4 | 0.321 / 0.369 / 0.403 (≈ N^{0.32}) |
Befunde: (1) **(S_γ) gilt für Kneser numerisch genauso** (‖Ĵ_N⁻¹‖ = 1.250 bis N = 600). (2) Bei zu dünner Linse wächst der Eckanteil exponentiell
mit der vorhergesagten kleinen Rate (lokale Eckrate ≈ 0.0105 minus Dämpfung 0.95·S ≈ 0.0095 ⇒ ≈ 0.001): die Linse muss über dem Eckprofil liegen;
das einfache „strikte Expansion genügt“ (IDEEN K14) ist **falsch** — die globale Eckrate steckt in (S_γ) als Ungleichung c < 0.95·s(π) (bzw. Λ ≤ 0.95 s).
(3) Auch bei passender Linse wächst der Eckanteil noch **polynomial** (≈ N^{0.3}, Kneser), liegt aber bis N = 600 bei 0.4 ≪ 1.25 (Bulk).
Asymptotisch ist also eher ‖Ĵ_N⁻¹‖ ≍ max(1.25, c·N^{κ_Ĵ}) mit κ_Ĵ ≈ 0.3 zu erwarten — immer noch κ_Ĵ + 3/2 ≈ 1.8 ≤ 2. Woher die N^{0.3} kommen
(Linse nahe L zu dünn oder echte Eckgrösse), testen die Varianten in k04b_varianten.log (S = 0.025; Linse nahe L dicker: τ₁ = 5°, p = 1.15).

### 5.4 Woher kommt das polynomiale Wachstum des Eckanteils? Linsenvarianten (kubisch, k04b_varianten.log) und Widom-Struktur (k05.log)
| Linse | ‖Ĵ_N⁻¹‖ (N = 300 / 500 / 700) | Eckanteil ‖Ĵ^{−T}u₁‖ | lokaler Exponent |
|---|---|---|---|
| S = 0.0166, τ₁ = 12°, p = 1.3 (Standard) | 1.251 (alle N) | 0.322 / 0.382 / 0.423 | ≈ 0.33 → 0.30 |
| S = 0.025 (fern von L dicker) | 1.210 (alle N) | 0.267 / 0.317 / 0.346 | ≈ 0.34 → 0.26 |
| S = 0.0166, τ₁ = 5°, p = 1.15 (nahe L dicker) | 1.251 (alle N) | **0.204 / 0.224 / 0.235** | **≈ 0.18 → 0.14** |
⇒ Das restliche Wachstum des Eckanteils kommt aus der **Linse nahe L**, nicht aus dem fernen Γin-Bogen; eine nahe L dickere (aber weiter tangentiale,
p > 1) Linse drückt es auf ≲ N^{0.15}, abnehmend. Ob ein Rest ~N^{0.05} bleibt (Innenschicht: |Y(L)| ~ √N auf Breite N^{−0.894}, also L²(γ)-Norm
~ N^{0.05}), ist bei N ≤ 700 nicht zu entscheiden. Für das Ziel „Eck-Polynomfaktor κ_Ĵ + 3/2 ≤ 2“ reicht das reichlich.
**Widom-Struktur (N' = 600 als Ersatz für 𝒥⁻¹):** ‖Ĵ_N⁻¹ − P_NĴ_{N'}⁻¹P_N‖ ≈ 0.51 / 0.54 / 0.56 / 0.47 (N = 200 / 300 / 400 / 500). Die Differenz sitzt in
einem Kantenband der Breite ∝ N (bei N = 200 tragen die letzten 100 Spalten 0.51 gegen 0.02 im Rest). Also: Ĵ_N⁻¹ = (Abschnitt der unendlichen
Inverse) + (Kantenkorrektur von Norm O(1), semiklassisch: Indexverhältnis nahe 1) — wie Widoms T_N(a)⁻¹ ≈ P_NT(a⁻¹)P_N + W_NP_NT(ã⁻¹)P_NW_N, nur dass
die Kante hier ein Band ξ ∈ (1−δ, 1] ist. **Folgerung für die Parametrix (Plan-Schritt 2):** sie braucht drei Bausteine — Bulk (= 𝒥⁻¹ = J_∞⁻¹C_k⁻¹,
nach O7 beschränkt), Kantenband (explizit zu bauen; im x-Bild ist das der „Steg“) und Ecke L; die Kantenkorrektur ist NICHT klein, aber beschränkt.

### 5.5 Radiales Quadrat (k04b_radial.log; Linse S = 0.0166, τ₁ = 5°, p = 1.15)
‖Ĵ_N⁻¹‖ = 1.310 für N = 300 / 500 / 700 (‖J̃_N⁻¹‖ = 3.2 / 18.9 / 142.5); Eckanteil 0.178 / 0.201 / 0.212 (≈ N^{0.2 → 0.14}). **(S_γ) gilt auch im
radialen Modell** — es ist damit das Testmodell für die Parametrix (Möbius scheidet aus, V 21:15).

## 6. Der „Reserve-Weg“ (Vorschlag V/Leitung 21:50): genügt eine n-te-Wurzel-Schranke plus Reserve? — Nein, und warum (11.10. 21:55)

**Frage.** In der konformen Koordinate ist die Rate der Eckmoden auf der Linse negativ (Dämpfung e^{−S·N} gegen Wachstum e^{cN}, V: Faktor κ = 1.00).
Genügt dann eine Ratenaussage mit o(1)-Fehler (Potentialtheorie W11/W12) plus gleichmässig kontrollierter Bulk, um W §26.3 zu umgehen?
**Antwort: nein.** Zwei Stellen haben Rate genau 0, nicht negativ, und dort ist polynomiale Genauigkeit erzwungen:
1. **Die niedrigen Koeffizienten der Eckmoden.** Y_u(t) = Σ_k y_{u,k}(ρ/t)^k → y_{u,0} für t → ∞; auf jeder Linse ist der Anteil Σ_{k<x₀N} y_k(ρ/t)^k nicht
   gedämpft. Es gilt exakt (bewiesen, eine Zeile): **|y_{u,0}| = |⟨(C_k^{(N)})⁻¹e₀, y_u⟩| ≤ ‖Ĵ_N⁻¹‖** für jedes Einheits-u, d. h. (S_γ) ⇒
   ‖J_N⁻¹T_N[1]‖ ≤ ‖Ĵ_N⁻¹‖; allgemein ‖J_N⁻¹T_Ng‖ ≤ ‖Ĵ_N⁻¹‖·‖g∘k‖_{H²(𝔻)} für jedes auf Int γ analytische g. Nun ist J_N⁻¹T_Ng = Π_N(J_∞⁻¹g), und
   J_∞⁻¹g hat bei L generisch eine Kneser-artige Singularität (x−L)^{s₀} (Re s₀ = 5.447, die Mellin-Nullstelle; Landkarte A2). **(S_γ) enthält also (H₀)
   für alle glatten Eingaben zugleich** — und für diese verlangt W §26.3 polynomiale, nicht e^{o(N)}-genaue Kontrolle der Eckanregung. Gemessen:
   max_{k<0.5N}|y_k| ≈ 0.03 (W), Rate 0.
2. **Die Ecke L selbst.** γ muss durch L gehen (sonst liegen Pole η nahe L in K_γ, und der Schlitzbeitrag ∫u^{5.447}|E_N|du bleibt unkontrolliert);
   in L ist die Rate 0 (|Y(L)| ~ √N). Mit Ratenfehler ε_N ist der „schlechte“ Bereich |τ| ≲ (ε_N/a)^{1/q}; er schrumpft nur dann auf die Innenschicht
   N^{−0.894}, wenn ε_N = O(polylog N/N) — das ist bereits polynomiale Genauigkeit.
**Was der Reserve-Gedanke trotzdem liefert (für die Planung):**
- Der **Steg-Anteil** (Indizes ≳ x₀N, fern von L) ist mit Reserve kontrollierbar: |Y^{high}(t)| ≤ σ₁√N e^{−x₀N s(t)}, also braucht man dort nur die
  Ungleichung c < x₀s (Wachstumsschranke (I) als n-te-Wurzel-Aussage) — mit dicker Linse (s(π) ≫ c) sogar nur eine grobe. Reserve: Standardlinse
  S − c_∞ ≈ 0.0048/Grad (Vs Gesetz κ = 1.00), mit dicker Linse fern von L ≈ 0.05; in L und für die niedrigen Koeffizienten 0.
- Damit zerfällt (S_γ) sauber in: (Bulk) σ_{r+1}(J̃_N⁻¹) ≤ B · (Steg) c < s(π) [n-te Wurzel] · **(Anregung) niedrige Koeffizienten y_{i,k}, k < x₀N, und
  Werte nahe L polynomial** — letzteres ist die echte, unvermeidbare Eckaussage (W §27.5: „(H-I) ist selbst eine (E_L)-Aussage“).
- W11 kann auch den Steg-Teil nicht allein liefern: die tropische Bilanz ist verschiebungsinvariant (W13) und fixiert die absolute Rate c nicht; man
  bräuchte (I) aus einer anderen Quelle.
**Folgerung:** Der Reserve-Weg ist als Abkürzung **gescheitert** (präziser Grund: Rate 0 bei den niedrigen Koeffizienten und in L, und (S_γ) ⇒ (H₀) für
glatte Eingaben). Er bestätigt aber die Arbeitsteilung: was zu beweisen bleibt, ist die polynomiale Eckanregung durch glatte Daten — in der
Operatorform (S_γ), für die die Parametrix (Bulk-Inverse 𝒥⁻¹ + Kantenband + Eckbaustein) das richtige Werkzeug ist.
**Familie (Hinweis V):** einverstanden — (S_γ) muss gleichmässig für T_N𝒥_P (P ∈ 𝒮_N ∩ B_N) gelten; die Vorkonditionierung C_k wirkt auf die
Gleichungsseite und ist von P unabhängig, K1″(a) gilt wörtlich mit W^P = 𝒥_Pk_η (analytisch auf Int γ, solange P(Int γ‾) ∌ η). Messbar wie A3.

## 7. Gegenlesen von Vs (K2)-Zertifikat (v06_k2_rand.py, v07_k2_zusammen.py, p32_rand_ableitungen.py) — 11.10. 21:57
**Ergebnis: Logik richtig; drei formale Punkte ohne numerische Wirkung; eine Abhängigkeit (p32a) von mir nicht geprüft.**
- Fourier-Schluss: b_k = (1/2π)∫f(ρe^{iθ})e^{−ikθ}dθ, sechsmal partiell integriert (∂_θ^j f stetig für j ≤ 5, ∂_θ⁵f ~ |θ−θ_L|^{0.447},
  ∂_θ⁶f ~ |θ−θ_L|^{−0.553} ∈ L¹ aus der Strukturform mit Re(1+γ) = 5.447) ⇒ |b_k| ≤ k^{−6}‖∂_θ⁶f‖_{L¹}/(2π) ✓; 2·(388.96 + 0.002 + 1060.36 + 0.057)
  = 2898.77, C₆ = 461.36 ✓.
- Zellschranke ∫_Zelle ≤ 2h_b|∂⁶f(θ_m)| + h_b²·sup|∂⁷f| (Mittelwertsatz) ✓; sup über die Zelle aus Jets auf der Krawczyk-Scheibe (z_c, z_r) ⊇ slog(Bogen) ✓;
  Annahmeregel beeinflusst nur den Aufwand ✓; nan/inf ⇒ Halbierung (die p32-Lücke) ✓.
- Zweigtreue innerhalb der Abschnitte (Punkt-Krawczyk am gemeinsamen Endpunkt, Scheibe in BEIDEN Zellscheiben) und an den Abschnittsgrenzen (v07,
  Scheibenradien beim Einlesen konservativ verkleinert) ✓; Startzweig reell (θ = 0: slog(ρ) ≈ 0.31, θ = π: slog(−ρ) ≈ −1.75) ✓.
- rest_schranke (|θ − θ_L| < δ_min) nachgerechnet: |w| ≤ |x−L|/(1−κ₂); |q| = e^{Re γ·log|w| − Im γ·Im u + 2π Im θ_K} ⇒ q_t = C_q w_t^{Re γ} ✓;
  Schwarz-Lemma für G_k (k ≥ 2, G_k(0) = 0) ✓; Koeffizient von H₆ in ∂_θ⁶f ist A⁶·6!·P'(H) ✓; |w| ≥ c_w|τ| und Exponent −0.553 < 0 ⇒ Integral
  exakt ✓.
- **Formale Punkte:** (F1) Die Mittelpunkt-Scheibe (z_m, z_mr) wird nicht explizit auf Enthaltensein in (z_c, z_r) geprüft — sonst könnte c6 in
  Prinzip von einem anderen Zweig stammen; eine Zeile `assert abs(zm − zc) + zmr < zr` schliesst das (sehr wahrscheinlich erfüllt, Newton vom
  Zentrum aus). (F2) v07 summiert float(integral.upper()) in Gleitkomma — Rundung ≤ 1e-15 relativ; formal mit Faktor (1 + 1e-12) oder arb-Summe.
  (F3) wie Vs eigene Liste: p32a (exakte Strukturform θ_K, Ps Beweis) ist Voraussetzung im Strukturgebiet und von mir **nicht** gegengelesen.
- Bewertung: (K2) mit C₆ = 461.36 gilt als **verifiziert bis auf p32a**; die Reserve gegenüber dem Bedarf (C₆ ≤ 5e5) ist ≈ 1000.

## 8. Grenzoperatoren der Abschnittsfolge Ĵ_N und ein Stabilitätssatz — was die klassische Theorie leistet und wo sie bricht (11.10. 22:10)

### 8.1 Die drei Grenzobjekte, numerisch bestimmt (kubisch; k05, k06)
- **Bulk:** Ĵ_N⁻¹ → 𝒥⁻¹ = J_∞⁻¹C_k⁻¹ (starker Limes), ‖𝒥⁻¹‖ ≈ ‖Ĵ_{N'}⁻¹‖ = 1.2512. Beschränktheit von 𝒥⁻¹: O7 (Fredholm, Index 0) + Injektivität +
  Subordination für C_k⁻¹ (bewiesen bis auf Os zitierte Standardsätze).
- **Kante fern von L:** keine sichtbare Korrektur. Erklärung (Skizze): am Index ξ ≈ 1 wird der Kompositionskanal abgeschnitten (er streckt die Grade),
  übrig bleibt der Identitätskanal 𝓛 = M̂·C_Ψ — untere Dreiecksgestalt, invertierbar mit beschränkter Inverse (Subordination), auf dem Linsenbogen
  zusätzlich strikt expandierend. Der Kanten-Grenzoperator fern von L ist damit „trivial invertierbar“.
- **Ecke L (und L̄):** die Widom-Differenz D_N = Ĵ_N⁻¹ − P_N𝒥⁻¹P_N hat Singulärwerte in Paaren (L, L̄): 0.492/0.484, 0.22/0.21 (N = 150), 0.512/0.506
  (N = 200; sauber, da N' = 600 ≥ 3N), führende Singulärvektoren auf BEIDEN Seiten beim Winkel arg L (77–79°) und im Indexband ξ ∈ (0.6, 1) mit
  Spitze bei 0.8–0.9. Deutung: am Fixpunkt wirkt die Komposition auf „lokale Frequenzen“ κ (Pakete e^{κ(x−L)}) als komplexe Streckung κ ↦ μκ;
  das Abschneiden ist |κ| ≲ N. Das Band (0.6, 1) sind die ersten Glieder der Kette ξ, ξ/|μ|, ξ/|μ|², … (1/|μ| = 0.853). In log-polaren Koordinaten
  (t = log|κ|, ϑ = arg κ) ist das ein **halbseitiges Wiener–Hopf-Problem für I + μ⁻¹D_μ** (Translation um log μ, abgeschnitten bei t = log N); wegen
  |μ⁻¹| < 1 ist sein Symbol 1 + μ^{s−1}·(…) auf der L²-Linie nullstellenfrei mit Windungszahl 0 (vgl. O7a, O10: λ_L = 0.98988) ⇒ invertierbar
  (Skizze). Die Periodizität in ϑ zusammen mit der Translation log μ ergibt genau Ws Kegel mit Deckverschiebung 2πi − log μ.

### 8.2 Warum die klassische Abschnittstheorie (Böttcher–Silbermann, Hagen–Roch–Silbermann, Roch–Santos–Silbermann, Lindner) nicht direkt greift
Die klassische Aussage „stabil ⇔ Operator und alle Grenzoperatoren (Kante, Unstetigkeitsstellen) invertierbar“ gilt für Folgen in der Algebra der
BESCHRÄNKTEN Folgen, A_N = P_NAP_N mit A beschränkt und bandbeherrscht bzw. in einer C*-Algebra lokalen Typs (Toeplitz/Hankel/Multiplikation/Mellin).
Hier bricht sie an zwei Stellen — die zweite ist die wesentliche:
1. **Exponentielle Fehlskalierung (technisch).** ‖Ĵ_N‖ ~ e^{0.47N} (Kompositionskanal auf Γout); Ĵ_∞ ist unbeschränkt und nicht bandbeherrscht
   (Komposition mit Φ = f∘k/ρ und auch mit der konformen Abbildung Ψ = k/ρ streckt die Grade: Spalte j hat Masse bis zu Grad ≈ j·Streckung, daher das
   Kantenband der Breite ∝ N). Gemessen: die natürliche Parametrix P_N𝒥⁻¹P_N^w hat RECHTS-Residuum E_N mit ‖E_N‖ ≈ e^{0.47N} (k07: 1.1e19 bei N = 100,
   6.7e39 bei N = 200), und die Schwanzabschnitte des beschränkten Operators 𝒥⁻¹ sind exponentiell singulär (k08: σ_min(Q_N𝒥⁻¹Q_N) ≈ 1e-17, also kein
   Ausweg über die Jacobi-Identität Ĵ_N⁻¹ = T11 − T12T22⁻¹T21, obwohl ‖T12T22⁻¹T21‖ ≈ 0.46 durch Auslöschung klein bleibt). Normen, in denen Ĵ beschränkt
   ist (Unbekannte in E²(int K_γ), wie in O7), zerstören die diskrete Stabilität (der gestreckte Kompositionskanal fällt durch das Abschneiden heraus).
   Gleichmässige Beschränktheit ist nur auf der Ebene der Inversen Ĵ_N⁻¹ sichtbar.
2. **Lokale Invertierbarkeit genügt NICHT (grundsätzlich).** In der x-Koordinate (J̃_N, ohne Linse) sind ALLE lokalen Objekte invertierbar — Bulk J_∞
   (O7), Kante auf Γin (dort nur der Multiplikationskanal L(m), untere Dreiecksgestalt), Kante auf Γout (Kompositionskanal abgeschnitten, Rest wie eben),
   Ecke (Os Mellin-Modell, ‖T_L⁻¹‖ = 1.28) — und doch wächst ‖J̃_N⁻¹‖ wie e^{cN}. Die Instabilität ist **global**: der Kompositionskanal transportiert
   Wellenpakete makroskopisch durch den Phasenraum (Fourier-Integraloperator), und die Eckmode ist eine geschlossene Phasenraum-Schleife (Beine bei arg L
   über alle ξ ⇄ Steg bei ξ ≈ 0.95 über ganz Γin) mit Netto-Verstärkung e^{cN}; c ist eine globale Grösse (W: freies Randwertproblem, θ-Abhängigkeit
   nach M). Bandbeherrschte Operatoren kennen keinen makroskopischen Transport, deshalb reicht dort die lokale Invertierbarkeit.
   In der w-Koordinate bricht die Linse die Schleife (Dämpfung e^{−S·N} auf dem Steg), und genau dann wird die Folge stabil: **V misst die scharfe Schwelle
   S* = c_∞** (v_pruefstand §3.1). Die nötige Zusatzbedingung ist also eine **dynamische (Antonevich-artige) Schleifenbedingung**: der Spektralradius der
   gewichteten Phasenraum-Zirkulation muss < 1 sein — wie bei Antonevichs Invertierbarkeitstheorie gewichteter Verschiebungsoperatoren aI + bW_α, wo die
   Bedingung über geometrische Mittel der Gewichte längs der Bahnen von α formuliert ist.

### 8.3 Stabilitätssatz in der Form, die die Daten verlangen (Vermutung; bewiesen ist nur die Rückführung K1″)
**Vermutung S.** Seien J_∞ invertierbar (Bulk), das Eckmodell an L invertierbar (halbseitiges Wiener–Hopf-Problem für I + μ⁻¹D_μ), und die Linse erfülle die
Schleifenbedingung s(τ) ≥ (1 + η)·Λ(τ) auf Γin (Λ = n-te-Wurzel-Profil der Zirkulation, Ws freies Randwertproblem; insbesondere s(π) > c_∞).
Dann ist sup_N‖Ĵ_N⁻¹‖ < ∞ (oder ≤ CN^{κ} mit kleinem κ, Eckanteil).
Struktur eines Beweises (Plan, nicht ausgeführt): Ĵ_N⁻¹ = (Bulk-Inverse) + (Eck-Korrektur bei L, L̄) + Rest, wobei der Rest durch die Schleifenbedingung
exponentiell klein ist. Die Schleifenbedingung ist eine globale Ungleichung (Wachstumsschranke der Zirkulation), alles andere ist lokal.
**Neue Mathematik, die gebraucht wird:** eine semiklassische (h = 1/N) Abschnittstheorie für Kompositions- bzw. gewichtete Verschiebungsoperatoren auf
Hardy-Räumen mit Abschneiden im Grad, deren Stabilitätskriterium „lokale Invertierbarkeit (Bulk, Kante, Ecke) + Spektralradius der gewichteten
Phasenraum-Zirkulation < 1“ lautet. Der globale Teil ist eine Max-Plus-/Tropen-Aufgabe (Ws tropische Bilanz W11, Unterlösungen 27.3) — aber als
UNGLEICHUNG (obere Schranke der Zirkulationsrate), nicht als Asymptotik. Ob die tropische Theorie diese Ungleichung liefert, ist offen: W13
(Verschiebungsinvarianz) zeigt, dass sie die absolute Rate in der x-Koordinate nicht fixiert; in der w-Koordinate ist die Frage neu (die Linse ist ein
absolutes Gewicht), aber unbeantwortet. **[23:39: beantwortet in §9.2 — die Linse ist ein Gewicht der Auswertung, nicht der Bilanz; Tropik allein
reicht nicht. Vermutung S ist durch S′ (§9.1) ersetzt: die Ecke ist keine eigene Hypothese.]** Die Eck-Korrektur selbst (lokal, polynomial) bleibt in jedem Fall ein Problem der starken Ebene (§6).

## 9. Lange Strecke A2a–A2c (11.10. 22:19–23:39): zwei Teilfragen entschieden, Kopplungssatz mit Familie formuliert

### 9.1 A2a entschieden — die Eck-Korrektur ist lokal und universell (numerisch); ein eigenständiges Eckmodell mit Gradabschnitt gibt es nicht (bewiesen + numerisch)
**(i) Universalität (k09_ecke_global.py, k09.log, k09b.log; N' = 600, N = 200, Linse S = 0.0166 wie k04; D_N = Ĵ_N⁻¹ − P_NĴ_{N'}⁻¹P_N).**

| Modell | mit kubisch gemeinsam | σ₁(J̃_600⁻¹) (globale Eckamplitude, x) | ‖Ĵ_200⁻¹‖ | Singulärwerte von D_200 |
|---|---|---|---|---|
| kubisch | — | 70.71 | 1.251155 | 0.51234 0.50566 0.23452 0.22319 0.09172 0.08747 |
| jet4 +1e-4 | L, μ, f'', f''' (global anders, Δf bis 5e-3 auf D̄) | 70.46 | 1.251179 | 0.51233 0.50567 0.23454 0.22320 0.09171 0.08749 |
| jet4 −1e-4 | dto. | 70.96 | 1.251168 | 0.51234 0.50566 0.23450 0.22319 0.09173 0.08745 |
| Kneser | nur L und μ | 68.41 | 1.250337 | 0.51215 0.50549 0.23437 0.22314 0.09175 0.08754 |
| radial | L, |μ|; arg μ = 31.4° statt 38.3° | 51.19 | 1.310262 | 0.45694 0.45345 0.17586 0.16493 0.06813 0.06226 |

(jet4_ε: f + (ε/ρ⁷)(x² − 2Re L·x + ρ²)⁴, verschwindet in L, L̄ von 4. Ordnung; `jet4_koeff` in k_kern.py.)
Die globale Amplitude reagiert auf die flache Störung in erster Ordnung (∓0.35 %), das führende D-Paar um ≤ 2e-5 relativ, die Paare 2–3 um ≤ 2.3e-4.
Kubisch und Kneser (nur L und μ gemeinsam) stimmen in D_200 überein (führendes Paar 4e-4, alle sechs ≤ 8e-4) (globale Amplituden 3.3 % verschieden); das radiale Modell (anderes
arg μ) weicht deutlich ab. V bestätigt unabhängig (BRETT 23:14: kubisch = Kneser auf 3–4 Stellen; D_N fällt fern von L wie e^{−SN}; wird mit
dickerer Linse an L kleiner). **Ergebnis: Die Eck-Korrektur der konformen Abschnitte ist eine universelle Funktion von μ und der Linsengeometrie
bei L; der globale Schleifengewinn geht nicht ein.**

**(ii) Das eigenständige Eckmodell ist nicht stellbar (k10_affin.py, k10.log).** Natürlichster Kandidat: f_aff(x) = L + μ(x − L), m ≡ μ —
in z = x − L genau I + μ⁻¹D_μ (bis auf den Faktor μ), also das Modell aus PLAN A2a.
- *Bewiesen (Zeilen):* In der x-Basis ist J_aff obere Dreiecksmatrix (affine Komposition erhöht den Grad nicht), C_k untere. Daher
  Ĵ_N⁻¹ = (J̃⁻¹)_N(C_k⁻¹)_N und **D_N = −P_N J̃⁻¹ Q_N C_k⁻¹ P_N** (Produkt zweier Eckblöcke; numerisch ‖D − D_Formel‖/‖D‖ = 4e-20).
- *Aber (numerisch):* ‖Ĵ_200⁻¹‖ = 3.4e60, ‖Ĵ_600⁻¹‖ = 4.2e182, σ₁(J̃_600⁻¹) = 6.6e182 — **Rate 0.700 pro Grad**, in beiden Koordinaten.
- *Erklärung der Rate (Skizze mit Zahl):* In der Basis (x−L)^k ist J diagonal, J⁻¹(x−L)^k = (x−L)^k/(μ^k + μ). Die Entwicklung
  1/(μ^k + μ) = Σ_j (−1)^j μ^j μ^{−k(j+1)} erzeugt Terme (φ^{−(j+1)}(x))^n; φ^{−j}(D̄) ist eine Kreisscheibe durch L, um −j·38.3° gedreht, mit
  max_{D̄}|φ^{−j}|/ρ = |1 − μ^{−j}| + |μ|^{−j} = 1.48, 1.82, 2.00, **2.021** (j = 4), 1.90, … ⇒ Rate ln 2.021 = 0.7036 (gemessen 0.700 bei N = 600).
- *Ursache auf Operatorebene (bewiesen):* J_aff ist auf E²(D ∪ φ(D)) **nicht injektiv**: (x−L)^{s_k} mit μ^{s_k−1} = −1, s_k = 1 + (2k+1)πi/log μ,
  Re s_k = 5.45 + 8.90k (k ≥ 0), liegen im Kern (der Zweig von log(x−L) ist auf dem Sektor S = H ∪ μH der Öffnung 218° konsistent, wie in O7a;
  J_aff(x−L)^s = (μ^s + μ)(x−L)^s = 0). Die lokalen Mellin-Lösungen setzen sich beim affinen Modell global fort; bei kubisch/Kneser verhindert das die
  globale Struktur (O7 + Injektivität D4).
- **Folgerung:** Eine Abschnittsaufgabe „Ecke allein“ gibt es nicht: der Gradabschnitt ist eine globale Bedingung, und die Dilatation um L mit
  komplexem μ dreht Teile von D̄ nach aussen. Lokal stellbar ist nur die Operator-Ecke (O7a: Mellin auf Re s = −½ im unendlichen Sektor,
  ‖T_L⁻¹‖ = 1.281) — und die ist bereits Teil der Bulk-Invertierbarkeit (O7). **Damit ist die Ecke keine eigene Hypothese mehr:**

**Vermutung S′ (ersetzt Vermutung S, §8.3).** (B) C_kJ_∞ ist auf Operatorebene invertierbar (O7 + Injektivität; enthält die Ecke über O7a) und
(SB) die Linse erfüllt die Schleifenbedingung s ≥ (1 + η)Λ auf Γin (9.2). Dann sup_N‖Ĵ_N⁻¹‖ < ∞, und D_N ist bis auf o(1) die universelle Eck-Korrektur
zu (μ, Linsengeometrie bei L).
*Konsistenz mit allen drei Beispielen:* affin — (B) verletzt ⇒ instabil (e^{0.70N}); x-Koordinate — (B) erfüllt, (SB) verletzt (keine Linse) ⇒ e^{cN};
w-Koordinate mit Linse — beides erfüllt ⇒ stabil (bis N = 2000, V). Vs Kerben-Test (BRETT 23:14) zeigt zusätzlich, dass (SB) punktweise scharf ist.

### 9.2 A2b entschieden — die Schleifenbedingung folgt NICHT aus Ws Tropik allein; sie folgt aus W12′ plus EINER Eck-Normierung
**Formulierung (n-te Wurzel).** u_γ(t) := limsup_N (1/N) ln G_N(t) (Ws Hüllkonfiguration 27.3(b), Quellen = Einheitsvektoren, Quellrate 0 ausserhalb D).
Nach K1″(b, c) ist ‖Ĵ_N⁻¹‖ bis auf den Faktor √N gleich max_γ G_N, also (Hartogs-Lemma) limsup (1/N) ln‖Ĵ_N⁻¹‖ = max_γ u_γ.
**(SB):** max_γ u_γ ≤ 0; Profilform mit Λ(τ) := inf{s : u_γ(ρe^{s'+iτ}) ≤ 0 für alle s' ≥ s}: **(SB) ⇔ s(τ) ≥ Λ(τ) auf Γin**. Gemessen (V): Λ(π) = c_∞
(scharfe Schwelle, κ = 1.00), Λ(60°) ≈ 0.0115 (Kerbe), Reserve der Standardlinse 0.0048/Grad.

**Satz K3 (bewiesen, elementar; verallgemeinert W13).** Gegeben eine exakte lineare Identität Σ_i A_i(Y) = Q (A_i linear in der Unbekannten, Q die Quelle)
und Grenzraten (a_i; q) mit Ws tropischer Bilanz (an jeder Stelle wird das Maximum mindestens zweimal angenommen). Sei Σ die Menge der Stellen, an
denen das Maximum GENAU von Q und einem einzigen A_i angenommen wird. Ist Σ = ∅, so erfüllt (a_i + δ; q) für jedes δ ≥ 0 wieder die Bilanz.
*Beweis.* Raten linearer Bilder verschieben sich gemeinsam um δ. Ist q < max a_i, so wird das Maximum von zwei A-Termen angenommen und bleibt es; ist
q = max a_i, so sind wegen Σ = ∅ mindestens zwei A-Terme maximal; q > max a_i ist durch die Bilanz ausgeschlossen. ∎
**Lesart:** Obere Schranken für den Abstand „Lösungsrate − Quellrate“ kann die Tropik nur an Σ gewinnen (dort, wo die Quelle allein ausgleicht).
- *x-Koordinate:* Σ ⊂ Umgebung von L (W 27.3(c), numerisch am Modell) — die Tropik bestimmt höchstens Gestalt und Anstieg (W12′), das Niveau sitzt in L.
- *w-Koordinate (die in §8.3 offene Frage):* Die Linse ändert die Identitäten nicht. ‖Ĵ_N⁻¹‖ ist nach K1″ dieselbe x-Hüllkonfiguration, nur auf γ
  ausgewertet; Σ bleibt bei L. **Die Linse ist ein absolutes Gewicht der Auswertung, nicht der Bilanz.**

**Ergebnis.** (SB) folgt aus der Tropik allein in keiner Koordinate. Genauer gilt die Zerlegung
  **(SB) ⇐ W12′ (Gestalt von u_γ auf dem Schleifengebiet bis auf eine additive Konstante eindeutig) + Eck-Normierung u_γ(L) = 0 + Profilvergleich s ≥ Λ_{W12′}**,
und der Schleifengewinn c = Λ(Steg) − Λ(L) ist dann eine reine Gestaltgrösse. Gemessen ist u_γ(L) = 0 (|Y(L)| ~ √N, W). Die Eck-Normierung ist (H₀) im
einzelnen Punkt L auf n-ter-Wurzel-Ebene, also eine Aussage über die Ecke — nach 9.1(i) universell. Offen bleibt W12′ (Ws Vermutung).
Das bestätigt Ws Befund 27.5 auch für die Linsenform und bindet A2b an A2a: **das Niveau kommt aus der Ecke, die Gestalt aus der Tropik.**

### 9.3 A2c — Kopplungssatz, mit Familie (formuliert; Familie numerisch geprüft, k11_newton.py)
**Klasse.** Zweikanal-Operatoren 𝒥 = C_P + M_m (Komposition mit P, Multiplikation mit m), P, m analytisch auf einer Umgebung von Int γ‾ ohne L, L̄,
stetig bis L, L̄, mit Fixpunkten L_P, L̄_P auf γ und |P'(L_P) − μ| ≤ ε_E. Abschnitte Ĵ_N^𝒥 := T_N^w C_k 𝒥|_{P_N} mit DEMSELBEN C_k für die ganze Klasse.
**Vermutung K-A2c (Kopplungssatz).** Aus (B_𝒥) ‖(C_k𝒥)⁻¹‖_{H²(𝔻)→E²} ≤ B und (SB_𝒥) s ≥ (1 + η)Λ_𝒥 auf Γin folgt sup_N‖(Ĵ_N^𝒥)⁻¹‖ ≤ C(B, η, ε_E),
mit einer Konstante, die nur von B, η, ε_E abhängt.
**Zur Familie (Auftrag der Leitung: „C_k hängt nicht von P ab, damit ist die Familie gleich mit erledigt“).** Das stimmt für die Vorkonditionierung,
aber nicht für die Hypothesen in einer Kugel. Zwei Befunde:
- **Satz K4 (Strukturidentität; bewiesen, eine Zeile).** Für jedes P mit P' ≠ 0 nahe 0 gilt P'∘P = (eˣ + R_P')/P', R_P := P∘P − eˣ. Also
  **T_N𝒥_P = T_N[C_P + M_{eˣ/P'}] + T_N[M_{R_P'/P'}]**. Der erste Operator hat keinen Überlauf (eˣ/P' wird in x ausgewertet, nicht in P(x) ∉ D̄);
  der zweite ist der Residuenterm. (Die Sekantenoperatoren 𝒥_P v = v∘P + M_P v mit M_P = ∫₀¹f'∘(f + t(P − f))dt haben diesen Überlauf
  nicht — dort wird das echte f' ausgewertet; er betrifft die Tangente J_N(P) aus (H-i).) Für Fast-Lösungen P ∈ 𝒮_N (ARCHITEKTUR §4.2: T_N R_P ∈ ℂ·x^N) ist R_P' ≡ cNx^{N−1} mod x^N, der Residuenterm
  wirkt also nur auf die Zeilen N−1, N (Rang ≤ 2, Kopf) — das ist genau die Stelle, die (H-iii) (Kopf-Stabilität) kontrolliert.
- **Numerik k11 (Kneser, Ks Linse; Spalte = ‖Ĵ_N⁻¹‖ der echten Newton-Matrix T_N𝒥_P, P'∘P exakt aus dem Polynom).**

| N | lin (J_N bei f_K) | P = p_N | P = Π_N f_K | P = T_N f_K | seg 1e-9 | seg 1e-6 | seg 1e-3 | seg 0.5 |
|---|---|---|---|---|---|---|---|---|
| 200 | 1.250337 | 1.250337 | 1.250337 | 13.08 | 1.250337 | 1.250337 | 1.250337 | 6.55 |
| 300 | 1.250337 | 1.250337 | 1.250337 | 5.87e6 | 1.250337 | 5.87 | 5.87e3 | 2.93e6 |
| 400 | 1.250337 | 1.250337 | 1.250337 | 4.43e12 | 4.43e3 | 4.43e6 | 4.43e9 | 2.22e12 |
| 500 | 1.250337 | 1.250337 | — | 2.95e19 | | | | |
| 600 | 1.250337 | (V: 1.250337 bis N = 1000) | — | 2.33e25 | | | | |

  (seg λ: P = p_N + λ(T_N f_K − p_N).) Überlauf im Multiplikationskanal bei T_N f_K: ρ^N|[x^N](P'∘P − f'∘f)| = 3.7e29, 1.3e49, 3.5e69, 1.5e90, 1.1e110
  (Rate 0.47 = die Rate von ‖Ĵ_N‖ auf Γout). Die Inverse wächst bei T_N f_K wie **e^{0.135N}** (13 → 2.3e25), längs der Strecke **linear in λ**:
  ‖Ĵ⁻¹‖ ≈ max(1.25, λ·K_N), K_N = ‖Ĵ⁻¹(T_N f_K)‖. Bei p_N und Π_N f_K ist die Inverse dagegen auf 7 Stellen die von J_N (V: bis N = 1000).
- **Folgerung (Familie).** (1) Eine Kugel in einer Grössen-Norm, die T_N f_K enthält (etwa sup_{D̄}|P − f| ≤ R_N, denn |T_N f_K − f_K| ~ N^{−5.4}), trägt
  KEINE gleichmässige Stabilität: dort ist die Newton-Matrix exponentiell schlecht konditioniert. (2) Die Architektur (H1/H2) braucht das auch nicht:
  sie verlangt (H) nur auf 𝒮_N ∩ B_N, und T_N f_K ∉ 𝒮_N (sein Residuum hat alle Grade). (3) Auf 𝒮_N liefert Satz K4 genau die richtige Zerlegung:
  **Familienfassung von K-A2c** = Kopplungssatz für G_P := T_N[C_P + M_{eˣ/P'}] gleichmässig in P ∈ 𝒮_N ∩ B_N (Hypothesen (B), (SB) für C_P + M_{eˣ/P'},
  offene Bedingungen in P) + Kopf-Stabilität (H-iii) für den Rang-2-Residuenterm. Das deckt sich mit Vs Befund (Versuch 7): keine Grössen-Norm, sondern
  Kohärenz; Satz K4 ist der Mechanismus dafür („der Schwanz von f_K hebt η'∘f auf“ ⇔ R_P ≡ 0 mod x^N).
- **Offene Stelle der Familienfassung (benannt):** (B) für C_P verlangt P auf Int γ, also auf der Linse. Ein Polynom P ∈ 𝒮_N ∩ B_N weicht dort um bis zu
  e^{Ns}·N^{−6.4} von f_K ab (Taylor-Schwanz jenseits des Konvergenzradius ρ); mit S = 0.0166 ist das ≤ 1e-7 bis N ≈ 2000, aber für N → ∞ nicht klein.
  Ob die Kohärenz das auch auf der Linse aufhebt (Vs Daten bis N = 1000/2000 sprechen dafür), ist zu zeigen — das ist der Rest des Familienproblems.
  **[23:59, §10: entschieden — für (S_γ) an der Familie irrelevant (Ĵ_N^P hängt nur vom Keim ab; S = 0.1, N = 600: Abweichung 8e9 auf der Linse,
  Norm unverändert). P auf γ geht nur in Daten-Schlüsse (H-ii)/H1 ein, Kosten e^{NF₀}·Kopf statt e^{Ns}; Rest = Kohärenz-Aussage (Fam-koh).]**

### 9.4 Literaturabgleich (gelesen = Zusammenfassung/Einleitung selbst gelesen; nur zitiert = nur Titel/Abstract über Suche)
| Quelle | Status | Was sie leistet | Greift hier? Lücke |
|---|---|---|---|
| A. Antonevich, *Linear Functional Equations. Operator Approach* (Birkhäuser OT 83, 1996); Antonevich–Lebedev 1994 | nur zitiert | Invertierbarkeit gewichteter Verschiebungsoperatoren aI − bT_α; Spektralradius über invariante Masse von α (bei Fixpunkten: Dirac-Masse ⇒ Bedingungen an den Fixpunkten + Verträglichkeit längs der Bahnen) | α ist dort Selbstabbildung, kein grosser Parameter; hier ist der Schleifengewinn e^{cN} ein WKB-Wirkungsintegral längs einer Phasenraum-Schleife und hängt von f global ab (k09: flache Störung ändert σ₁ in erster Ordnung) — keine Fixpunktgrösse |
| Yu. Karlovich–V. Kravchenko, Mat. Sb. 1981 (SIO mit Shift, endlich viele Fixpunkte auf Ljapunow-Kurve) | Abstract gelesen | Fredholm-Kriterien, auch Fixpunkte höherer Vielfachheit | Operatorebene, Shift auf der Kurve selbst; kein Abschnitt |
| A. Karlovich–Yu. Karlovich–A. Lebre, arXiv 1009.5656 / 1010.5336 | Einleitung und Hauptsätze (Zusammenfassung) gelesen | (aI − bW_α)P₊ + (cI − dW_α)P₋ auf L^p(ℝ₊), Fixpunkte 0, ∞: Fredholm ⇔ Binome invertierbar (Bedingungen an |a| − |b|(α')^{−1/p} in den Fixpunkten) + Mellin-Symbol; Methode Allan–Douglas + Mellin-ΨDO | reelles α' im Fixpunkt; hier komplexes μ (Spirale: Mellin-Symbol 1 + μ^{s−1} windet sich unendlich oft, O7a braucht Sektorgewichte); kein Abschnitt im Grad |
| Kravchenko–Litvinchuk (Monographie SIO mit Shift) | nur zitiert | Standardtheorie Carleman/Nicht-Carleman | wie oben |
| Kriete–MacCluer–Moorhouse, JOT 58 (2007); K. Quertermous, arXiv 1205.5786 | Zusammenfassung gelesen | C*-Algebren aus Kompositionsoperatoren mit Randfixpunkt (linear gebrochene Selbstabbildungen, φ'(ζ) > 0) modulo Kompakte ≅ C₀([0,1]) ⋊ ℤ mit Potenzabbildungen (in log-Koordinaten Translation = Mellin/Dilatation); Index 0 | Selbstabbildungen, reelles φ'(ζ), keine Abschnitte (endliche Schnitte nur in einem Determinantenkriterium) |
| B. Silbermann, ZAA 1982 (lokale Theorie des Reduktionsverfahrens für SIO); Roch–Santos–Silbermann (ZAA 1997; Buch 2011); Hagen–Roch–Silbermann; Lindner (Grenzoperatoren) | nur zitiert | Abschnitte stabil ⇔ Familie lokaler/Grenzoperatoren invertierbar, für Folgen in Algebren beschränkter Operatoren (PC-Multiplikation, Faltung, Mellin, Flip, Hankel) | hier unbeschränkt (‖Ĵ_N‖ ~ e^{0.47N}), nicht bandbeherrscht, nicht-Carleman; lokale Invertierbarkeit reicht nachweislich nicht (§8.2) |
| Didenko–Silbermann, *Approximation of Additive Convolution-Like Operators* (2008) | nur zitiert | Stabilität von Näherungsverfahren für SIO mit Konjugation (antilinear) über reelle Erweiterungen von C*-Algebren | Konjugation/Carleman, keine Kompositionsdynamik |
| Elschner–Graham (Übersicht Mellin-Gleichungen; Kollokation auf Polygonen), Elschner–Jeon–Sloan–Stephan (Numer. Math. 1997) | Abstracts gelesen | Ecken: Stabilität über Mellin nach Gitter-Gradierung (nichtlineare Parametrisierung); trigonometrische (gradbeschnittene) Verfahren an Ecken gelten als schwierig | Analogie: die konforme Koordinate w + Linse spielt die Rolle der Gradierung — aber zum Brechen der dynamischen Schleife, nicht zum Glätten; keine Komposition |

**Lücke = der neue Teil.** Keine gefundene Quelle behandelt Abschnitte im GRAD für Operatoren mit Kompositionskanal, dessen Abbildung keine
Selbstabbildung ist. Drei neue Elemente: (1) der Gradabschnitt erzeugt einen semiklassischen Parameter h = 1/N mit makroskopischem Transport im
Phasenraum (Fourier-Integraloperator), Stabilität = lokal (Bulk inkl. Operator-Ecke) + globale Schleifenbedingung (Vermutung S′); (2) der Randfixpunkt
mit komplexem Multiplikator ist im Gradabschnitt NICHT abtrennbar (9.1(ii)) — die lokale Theorie lebt nur auf Operatorebene, die Eck-Korrektur ist
universell, aber nur innerhalb eines global injektiven Operators definiert; (3) die Familie ist über die Lösungsmenge zu führen (Satz K4), nicht
über eine Kugel.

### 9.5 Stand nach der Strecke
- **Entschieden:** A2a (Eck-Korrektur universell — numerisch; eigenständiges Eckmodell widerlegt — bewiesen + numerisch); A2b (Tropik allein reicht
  nicht — Satz K3 bewiesen; Zerlegung (SB) ⇐ W12′ + Eck-Normierung).
- **Neu bewiesen:** Satz K3 (Homogenität), Satz K4 (Strukturidentität P'∘P = (eˣ + R_P')/P'), Nicht-Injektivität von J_aff, exakte D_N-Formel im affinen Modell.
- **Formuliert:** Vermutung S′ (zwei statt drei Hypothesen), Kopplungssatz K-A2c mit Familienfassung auf 𝒮_N ∩ B_N; offene Stelle: Linsen-Abweichung
  der Familie für N ≫ 2000.
- **Nächster Schritt (Vorschlag):** Beweis von K-A2c in der einfachsten nicht-trivialen Lage — zuerst eine Aussage „(B) + (SB) ⇒ n-te-Wurzel-Stabilität“
  (exponentielle Ebene; das ist eine Parametrix-Aussage mit Fehler e^{−ηN}, ohne Ecke), danach die Eck-Korrektur als universelle O(1)-Grösse
  (Widom-Zerlegung mit dem Bulk-Operator als äusserer Parametrix).

## 10. Die Linsen-Spannung der Familie (Auftrag der Leitung; 11.10. 23:44–23:59) — entschieden

**Spannung.** Für (SB) muss die Linse überall über Λ liegen (max Λ ≈ 0.0118, N-unabhängig); Polynome P der Familie weichen auf der Linse aber um
≈ e^{N s(τ)}·N^{−5.8…−6.4} von f_K ab (Taylor-Schwanz jenseits des Konvergenzradius ρ). Frage: echtes Hindernis für (S_γ) an der Familie, oder irrelevant?

### 10.1 Was in Ĵ_N^P eingeht (bewiesen, elementar)
**(a) Keimabhängigkeit.** T_N(δ∘P) hängt für δ ∈ P_N nur von P mod x^{N+1} ab (P^j mod x^{N+1} = (P mod x^{N+1})^j mod x^{N+1}). Nach Satz K4 ist
P'∘P ≡ (eˣ + R_P')/P' mod x^{N+1}; das hängt nur von P mod x^{N+2} und R_P mod x^{N+2} ab, und R_P = P∘P − eˣ mod x^{N+2} nur von P in D
(P∘P braucht P nahe P(0) = a₀ ∈ D). **Also ist Ĵ_N^P = C_k^{(N)}T_N𝒥_P|_{P_N} eine Funktion von P|_D allein (endlich viele Taylor-Koeffizienten bei 0).
Werte von P auf der Linse gehen in die Matrix nicht ein.**
**(b) K1″(b, c) gelten wörtlich für jede Matrix der Familie.** Mit G_N^P(t) := ‖(J̃_N^P)⁻¹ẽ_t‖:
  ‖(Ĵ_N^P)⁻¹‖ ≤ c_k‖G_N^P‖_{L²(γ)},   G_N^P(t) ≤ |t|·Q_k·√(N+1)·‖(Ĵ_N^P)⁻¹‖  (t ∈ γ).
Beweis wie §5.1: benutzt werden nur die Matrix J_N^P (adjungierte Erzeugende Y_u^P ist ein Polynom in 1/t, analytisch in ℂ∖{0}) und die Daten
k_t∘k bzw. c∘k⁻¹ (analytisch in Int γ) — nie P auf Int γ. ∎ **Also (S_γ)^P ⇔ G_N^P polynomial auf γ: die Linse tritt nur als Ort der Daten auf
(Cauchy-Kern bei t auf der Linse), nie als Argument von P.**
**(c) Numerischer Beleg (k11_newton.py mit Linsenhöhe als 3. Argument, k11d.log, Kneser, dicke Linse S = 0.1, F₀ = 0.048):** bei N = 400 ist
|p_N − T_N f_K| auf dem ∂D-Teil von γ 6.2e-13, **auf der Linse 1.6e2**, und doch ‖(Ĵ^{p_N})⁻¹‖ = 0.965153 = ‖(Ĵ^{f_K})⁻¹‖ auf 6 Stellen
bei N = 600 ist die Abweichung auf der Linse **8.1e9** (auf ∂D 8.8e-14), die Norm wieder 0.965153 = ‖(Ĵ^{f_K})⁻¹‖. Dicke Linse verkleinert die
Norm sogar (0.965 statt 1.250). (Die von k11 zusätzlich gedruckte Rang-1-„Empfindlichkeit“ der (N,0)-Stelle ist NICHT aussagekräftig — sie ist bei p_N
riesig, ohne Wirkung; nicht verwenden.)

### 10.2 Wo P auf γ doch eingeht — und wie teuer das wirklich ist
P auf γ erscheint erst in Schlüssen, die DATEN auf γ brauchen: K1″(a) für die Familie (Pol-Antwort (H-ii): W^P = 𝒥_Pk_η) und die A-priori-Schranke in
H1 (Π^P_N f mit 𝒥_P f). Auch dort zählt nur T_N W^P (Keim), also darf P durch JEDEN Repräsentanten P̃ ≡ P mod x^{N+2} ersetzt werden.
**Carathéodory–Fejér (bewiesen, Standard):** min_{P̃} sup_{Int γ}|P̃ − f_K| ist (bis auf Faktoren) die Grösse der ersten N+2 w-Koeffizienten von
(P − f_K)∘k; ein Abweichungs-Koeffizient d_N bei x^N liefert (generisch, ohne Auslöschung) den w-Koeffizienten ≈ d_Nρ^N·E(0)^N = d_Nρ^N·e^{NF₀}; untere Schranke
sup ≥ max_i|c_i| (Cauchy), obere ≤ Σ_i|c_i| (Schur). Die Kosten sind also
**e^{N·F₀}·(Kopfabweichung)** mit F₀ = harmonisches Mittel der Linsenhöhe — nicht e^{N·s(τ)} punktweise. Zahlen: Standardlinse F₀ = 0.0084;
minimale zulässige Linse (s = Λ ≤ 0.0118 auf Γin, Γin ≈ 57 % des Kreises) F₀ ≲ 0.0067 (Schätzung). Mit Vs Kopfprofil (abgeleitet aus Vs Angabe „7 % des ℓ¹-Gewichts in k ≥ 0.9N“: ≈ 0.7·N^{−5.8} je Koeffizient) wird e^{NF₀}N^{−5.8} erst bei **N ≈ 6300 (Standardlinse) bzw. ≳ 7800 (minimale Linse)** von der Grösse 1 — für die schlechteste,
inkohärente Abschätzung. Die Schwelle der Leitung (2000–3000) gilt nur für den Repräsentanten P selbst.
**Die echte Grösse hinter grossen N ist linsenunabhängig.** Schätzt man die Datenkorrektur T_N[Δ𝒥 h] statt auf γ auf einer Hilfskurve der Höhe s
zwischen ∂D und der Linse ab, so kostet der Repräsentant e^{sN}, die Antwort G_N dort e^{(Λ−s)N} (Vs Gesetz κ = 1.00): Produkt ≈ e^{ΛN}, unabhängig von s.
Das ist genau Vs linsenunabhängige Lipschitz-Konstante für inkohärente hohe Grade (e^{0.010…0.013N}, S = 0.0166…0.1). Für die ECHTEN Fast-Lösungen
misst V dagegen eine polynomiale Änderung (0.015·N·ε_N bis N = 1000) — die Kohärenz hebt e^{ΛN} auf.

### 10.3 Entscheidung
1. **Für (S_γ) an der Familie ((H-i) in w-Form) ist die Linsen-Abweichung KEIN Hindernis** — bewiesen (10.1(a, b)), numerisch belegt (10.1(c)).
   Die Linse muss nicht mit N dünner werden; die Familie muss nicht „über f_K auf der Linse“ ausgewertet werden.
2. **Folgen für die Formulierung von K-A2c (Familienfassung):** nicht mit 𝒥_P als Operator auf Int γ formulieren (dort würde P auf der Linse gebraucht),
   sondern in G-Form bzw. Keimform: Hypothesen nur über die Matrizen J_N^P; der Kopplungssatz wird am Referenzoperator (f_K) bewiesen und auf die
   Familie über eine **Kohärenz-Aussage** übertragen.
3. **Was als echte, benannte Lücke bleibt (nicht die Linse):** (Fam-koh) Für P ∈ 𝒮_N ∩ B_N und Δ = P − f_K gilt
   ‖(J_N^P)⁻¹T_N[(𝒥_P − 𝒥_f)h]‖ ≤ poly(N)·ε_N·sup_γ|h| — d. h. die Datenkorrektur in (H-ii)/H1 ist polynomial, nicht e^{ΛN}. Ohne Kohärenz trägt die
   Abschätzung über Carathéodory–Fejér bis N ≈ 6300–7800 (inkohärent schlechtester Fall); mit Kohärenz (V: bis N = 1000 polynomial) für alle N.
   Struktur: in erster Ordnung ist (𝒥_P − 𝒥_f)h = Φ''(f)[h, Δ] (symmetrisch in h und Δ), und für Fast-Lösungen ist Δ = λ(Π^P_N f − f) + Kopf
   (ARCHITEKTUR §5.1) — die Eckanregung durch Φ''(f)[h, Δ] ist dann selbst eine (H₀)-artige Eckaussage. Das ist dieselbe Ecke, keine neue Front.
4. **Für V (Vorrang laut Leitung):** bei N = 1500/2000 (S = 0.0166 und S = 0.1) neben ‖(Ĵ^{p_N})⁻¹‖ auch max_{Linse}|p_N − T_N f_K| ausgeben
   (Zeile in k11_newton.py) — erwartet: Norm unverändert trotz Linsen-Abweichung ≫ 1; zusätzlich die Datenkorrektur aus Punkt 3 für h = k_η (η auf dem
   Schlitz) messen: ist sie polynomial (Kohärenz) oder ~e^{ΛN}?

## 11. A2c Schritt 1 („(B) + (SB) ⇒ n-te-Wurzel-Stabilität“) — in der geplanten Form endgültig nicht trennbar von der Ecke; exakte Kopf-Zerlegung als Ersatz (12.10. 00:05–00:25)

### 11.1 Warum Schritt 1 „ohne Ecke“ nicht geht (bewiesen, zwei Zeilen)
(a) **L liegt auf γ.** Nach K1″(c) gilt G_N(t) ≤ |t|Q_k√(N+1)‖Ĵ_N⁻¹‖ für t ∈ γ, also folgt aus n-te-Wurzel-(S_γ) insbesondere **G_N(L) = e^{o(N)}** —
eine Aussage über die x-Abschnitte im EINEN Punkt L, unabhängig von der Linse (deren Höhe in L null ist). Kein linsenbasiertes Argument kann sie
liefern; sie ist die Eck-Normierung aus §9.2. (Numerisch, k12: G_N(L) = 7.20, 8.18, 9.03, 9.80 für N = 300…600 ~ N^{0.44}, also wahr.)
(b) **Auch fern von L steckt die Ecke im Ziel.** Lemma K6 unten führt die Linsenpunkte auf das Grad-Antwortprofil zurück; dessen Teil für glatte
Daten (kleine Grade) ist Ws (H-I) = Anregung der Eckmode durch glatte Daten (W §27.5) — eine Eckaussage auf exponentieller Ebene.
**Folgerung:** Die geplante Reihenfolge „exponentielle Ebene per Parametrix ohne Ecke, Ecke danach“ ist falsch. Die Ecke ist auf BEIDEN Ebenen der
Engpass; die Linse nützt nur fern von L. (Das hätte man aus §6 Punkt 2 und §9.2 schon ablesen können: das Niveau sitzt in L.)

### 11.2 Lemma K6 (bewiesen, eine Zeile) — Linsenpunkte über das Grad-Antwortprofil
Mit R_N(k) := ‖J̃_N⁻¹e_k‖ (Antwort auf den Einheitsvektor vom Grad k) gilt für |t| = ρe^{s}, s ≥ 0:
  **G_N(t) ≤ Σ_{k≤N} e^{−ks} R_N(k)**  (ẽ_t = Σ_k (ρ/t)^k e_k, Dreiecksungleichung).
Ist ψ_N(ξ) := (1/N)ln R_N(ξN) ≤ ε auf [0, ξ₀] und ≤ c' auf [ξ₀, 1], so ist (1/N)ln G_N(t) ≤ ε + O(ln N/N) an allen Linsenpunkten mit s ≥ c'/ξ₀.
**Gemessenes Profil (k12, kubisch, x-Koordinate):**

| N | ψ(0) | ψ(0.1 … 0.8) | ψ(0.9) | ψ(0.95) | ψ(1) |
|---|---|---|---|---|---|
| 300 | −0.0022 | −0.0002 | −0.0002 | +0.0017 | −0.0706 |
| 400 | −0.0016 | −0.0001 | +0.0001 | +0.0031 | −0.0696 |
| 500 | −0.0013 | −0.0001 | +0.0006 | +0.0040 | −0.0688 |
| 600 | −0.0011 | −0.0001 | +0.0014 | +0.0045 | −0.0683 |

Die exponentielle Grösse von J̃_N⁻¹ kommt AUSSCHLIESSLICH aus Daten im Kantenband ξ ≈ 0.85–0.97 (dem Steg); glatte und mittlere Daten haben Rate 0,
der Kopf (ξ = 1) ist exponentiell klein. Damit ist Lemma K6 fern von L fast scharf: Schwelle s ≥ ψ(Band)/ξ_Band ≈ c/0.95 — genau die von V
gemessene Linsenschwelle (k04b: „c/0.95 ≈ 0.0126“). **Fern von L reduziert sich Schritt 1 also exakt auf zwei x-Aussagen:** (I) Wachstumsschranke
im Kantenband (ψ ≤ c' mit c' < ξ₀·s; bei dicker Linse genügt eine GROBE Schranke) und (H-I)_exp: ψ(ξ) ≤ o(1) für ξ ≤ ξ₀ (glatte/mittlere Daten).
Nahe L reicht K6 nicht (s → 0); dort braucht man die Phasenstruktur — die Ecke.

### 11.3 Satz K5 (bewiesen, elementar) — exakte Kopf-Zerlegung (LU) der Abschnittsinversen
Sind alle Abschnitte J̃_n invertierbar, so hat die unendliche Matrix J̃_∞ eine eindeutige Zerlegung **J̃_∞ = LU** (L untere Dreiecksmatrix mit
Einsen auf der Diagonale, U obere Dreiecksmatrix; alle Produkte endliche Summen), und J̃_n = L_nU_n für jedes n (geschachtelte LU = Gauss ohne Pivot).
Daraus, jeweils in einer Zeile:
1. **J̃_N⁻¹ = (U⁻¹)_N(L⁻¹)_N = Σ_{n≤N} c_n ℓ_nᵀ**, mit c_n := J̃_n⁻¹e_n (Spalte n von U⁻¹, „Primal-Kopf“) und ℓ_nᵀ := e_nᵀJ̃_n⁻¹/(J̃_n⁻¹)_nn
   (Zeile n von L⁻¹, „adjungierter Kopf“, ℓ_nn = 1). Der Kopfterm c_Nℓ_Nᵀ = J̃_N⁻¹ − (J̃_{N−1}⁻¹ ⊕ 0) ist die Rang-1-Ränderung.
2. **Konforme Koordinate:** C_k ist untere Dreiecksmatrix ⇒ Ĵ = (C_kL)·U ist die LU-Zerlegung von Ĵ mit **DEMSELBEN U** (unabhängig von k).
   Also Ĵ_N⁻¹ = Σ_{n≤N} c_n·(ℓ_nᵀC_k⁻¹): die Linse ändert nur die adjungierten Köpfe (Erzeugende Λ_n auf γ statt auf ∂D, wie K1″(b)).
3. **Die Abschnittsinversen sind die Partialsummen der Kopf-Entwicklung der Bulk-Inversen**, und (bei starker Konvergenz Ĵ_N⁻¹ → Ĵ⁻¹)
   **D_N = Ĵ_N⁻¹ − P_NĴ⁻¹P_N = −P_N U⁻¹ Q_N L̂⁻¹ P_N = −Σ_{n>N} (P_Nc_n)(P_Nr̂_n)ᵀ** (r̂_n := C_k^{−T}ℓ_n) — die allgemeine Form der affinen
   Hankel-Produkt-Formel aus §9.1(ii): die Eck-Korrektur sind die auf Grade ≤ N projizierten KÜNFTIGEN Köpfe.
4. **n-te Wurzel:** max_{n≤N}‖c_n⊗r̂_n‖ ≤ 2 max_{n≤N}‖Ĵ_n⁻¹‖ und ‖Ĵ_N⁻¹‖ ≤ Σ_{n≤N}‖c_n⊗r̂_n‖, also
   **limsup (1/N)ln‖Ĵ_N⁻¹‖ ≤ 0 ⇔ limsup (1/n)ln(‖c_n‖_{H²(D_ρ)}·‖ℓ_nᵀC_k⁻¹‖) ≤ 0** — Schritt 1 ist eine Aussage über EINE Folge (die Köpfe).
5. **Gleichmässig, hinreichend:** ‖Ĵ_N⁻¹‖ ≤ ‖U⁻¹W⁻¹‖·‖WL̂⁻¹‖ für jede Diagonalmatrix W (Diagonalgewichte vertauschen mit P_N) — (S_γ) folgt aus der
   Beschränktheit zweier UNENDLICHER Dreiecksoperatoren in einem passenden Gewicht, ohne N.

**Kopf-Numerik (k12, kubisch; Kopfterm-Norm = ‖c_N‖·‖ℓ_N‖):**

| N | ‖c_N‖ | ‖ℓ_N‖ (x) | ‖ℓ_NᵀC⁻¹‖ (w) | Kopfterm x | Kopfterm w |
|---|---|---|---|---|---|
| 300 | 6.3e-10 | 6.6e9 | 7.1e7 | 4.14 | 0.0448 |
| 400 | 8.3e-13 | 1.3e13 | 3.5e10 | 11.1 | 0.0287 |
| 500 | 1.1e-15 | 2.8e16 | 1.9e13 | 31.6 | 0.0212 |
| 600 | 1.6e-18 | 5.8e19 | 1.1e16 | 92.3 | 0.0177 |
| 800 | 3.3e-24 | 2.5e26 | 4.6e21 | 831 | 0.0152 |
| 1000 | 7.1e-30 | 1.1e33 | 2.1e27 | 7796 | 0.0149 |
(Profil ψ bei N = 800/1000: ψ(0.9) = 0.0029/0.0037, ψ(0.95) = 0.0050/0.0067, ψ(1) = −0.0676/−0.0671, ψ(ξ ≤ 0.8) ≈ 0; G_N(L) = 11.1/12.3.)

In x wächst der Kopfterm wie e^{≈0.011n} (die Eckrate), in w fällt er zunächst und **sättigt bei ≈ 0.015** (lokale Exponenten 1.55, 1.36, 1.00,
0.53, 0.10). Beide Faktoren sind exponentiell (‖c_n‖ ~ e^{−0.067n}, ‖ℓ_nᵀC⁻¹‖ ~ e^{+0.062n} bei N = 1000); Schritt 1 ist also die EXAKTE Auslöschung
zweier Raten — gemessen ist sie da (Köpfe beschränkt). Weil die Köpfe nicht summierbar sind (Σ‖Kopf‖ ~ 0.015·N), ist **gleichmässige** Stabilität
eine Fast-Orthogonalität der Köpfe über n — wie schon bei der Einheitsmatrix (Köpfe e_n⊗e_n der Norm 1); an der Ecke heisst das: Mellin-Struktur
der künftigen Köpfe (D_N = projizierte künftige Köpfe ist universell). Das ist die polynomiale Eckaussage in Kopfform.

### 11.4 Entscheidung zu Schritt 1 und Wendung
- **Schritt 1 in der geplanten Form ist endgültig gescheitert** (11.1): die exponentielle Ebene enthält in L und für glatte Daten die Ecke.
- **Was bewiesen ist:** Satz K5 (Kopf-Zerlegung, k-unabhängiges U, D_N als projizierte künftige Köpfe, n-te Wurzel ⇔ Köpfe subexponentiell),
  Lemma K6 (Linsenpunkte über das Grad-Antwortprofil).
- **Wendung (Vorschlag):** A2c neu um die Ecke organisieren, auf beiden Ebenen in Kopfform:
  (E₀) exponentiell: G_N = e^{o(N)} auf γ ∩ U(L) und (H-I)_exp ψ(ξ) ≤ o(1) für ξ ≤ ξ₀; fern von L dann Lemma K6 + grobe Kantenband-Schranke;
  (E₁) polynomial: Köpfe ~ 1/n und fast orthogonal über die Skalen (Mellin) ⇒ gleichmässig (S_γ).
  Beide sind Aussagen über die universelle Eckstruktur (§9.1(i)) — das legt nahe, sie im Eck-Grenzobjekt zu beweisen (Köpfe in der Skala ξ = Index/n
  bei L), nicht global.

## 12. Anschluss an G und W, Eck-Kopfform-Fenster (12.10. 00:26–00:34)

### 12.1 Übersetzung Satz K5 ↔ G9/G15–G18/G22 und W §25.1(c) (damit nichts doppelt gebaut wird)
| K (§11.3) | G (RUNDE3/4, runde5) | W §25.1 |
|---|---|---|
| c_n = J̃_n⁻¹e_n (Spalte n von U⁻¹) | p_nQ_n, Q_n = J̃_n⁻¹e_n/p_n monisch (letzte Komponente 1), p_n = (J̃_n⁻¹)_nn = 1/h_n | Q_n (exakter Interpolant, (b)) |
| ℓ_n = e_nᵀJ̃_n⁻¹/(J̃_n⁻¹)_nn (Zeile n von L⁻¹, ℓ_nn = 1) | A_n (monisch), Erzeugende 𝒜_n(x) = Σ_k A_{n,k}(ρ/x)^k | A_n |
| Kopfterm c_nℓ_nᵀ | Q_nA_nᵀ/h_n (G15/G18) | Q_nA_nᵀ/h_n ((c)) |
| Kopfterm angewandt auf die Cauchy-Daten ẽ_t | 𝒟_n(t)·(Richtung Q_n/‖Q_n‖), 𝒟_n = ‖Q_n‖𝒜_n/h_n (G17) | — |
| ‖Kopfterm‖ in x (ℓ²_ρ → ℓ²_ρ) | ‖Q_n‖‖A_n‖/h_n (RUNDE3 Z. 230: 1.728 bei 200, 7795.68 bei 1000) | 1.2/1.7/4.1/11.1 (N = 150…400) |
| ‖Kopfterm‖ in w (ℓ²(w) → ℓ²_ρ) | ≈ ‖𝒟_n‖_{L²(γ)} (K1″(b)-Gewicht |k'|/|k|) | — |
| J̃_N⁻¹ = Σ_{n≤N} Köpfe; G_N ≤ Σ|𝒟_n| (Dreieck) | G10/G18 (Summenkriterium) | (c) |
Die Normierungen: K arbeitet mit ρ-skalierten Koeffizienten und unnormierten Köpfen; G mit monischen Q_n, A_n und Pivot h_n; dieselbe Matrix J̃.
Meine Zahlen stimmen mit G und W überein (x-Kopf 7796 vs. G 7795.68 bei N = 1000; 4.135/11.12 vs. W 4.1/11.1 bei 300/400).
**Neu in K5 gegenüber G/W** (nur das ist eigener Beitrag): (i) die LU-Zerlegung der konformen Matrix Ĵ = (C_kL)U mit DEMSELBEN U — die Linse ändert nur
die adjungierten Köpfe (𝒜_n auf γ statt ∂D); (ii) D_N = −Σ_{n>N} projizierte künftige Köpfe (allgemeine Widom-Formel); (iii) n-te-Wurzel-(S_γ) ⇔
w-Köpfe subexponentiell; (iv) Gewichtskriterium; (v) Messung der w-Köpfe (sättigen bei 0.015 bis N = 1000). Alles Übrige (Biorthogonalität,
Summenkriterium, x-Kopfnormen, Eckprofil) ist bei G/W vorhanden und wird nicht neu gebaut.

### 12.2 Was G zur Ecke schon hat (runde5_szego, gelesen: BERICHT, ETAPPE3, ETAPPE4 §1–3; spiral_analysis.json ausgewertet)
Komplexe Kegelskala d_n = n^{−a}, a = 1/κ₁ (G35), normiertes Profil S_n(ζ) = 𝒟_n(w_L + d_nζ)/𝒟_n(w_L), Wiener-Abstand 400→800 = 0.0061,
Vorhersage bei 1000 auf 6e-7; Amplitude B_n = √n|𝒟_n(L)| = 0.863 (200) … 0.607 (1000), B_∞ unsicher; G37: aus B_n ≤ C_A und ‖S_n‖_{W,1.05} ≤ C_S
folgt |𝒟_n| ≤ C_AC_S/√n auf den Eckkreisen; G40: Verschiebungsrang 2; G39: gemischtes Symbol 1/(1−sF) + M/(1−sz).
**Damit ist mein Abbruchkriterium „nicht skalenfrei beschreibbar“ numerisch beantwortet (skalenfrei: ja). Offen: Amplitude, Profilrest, Fast-Orthogonalität.**

### 12.3 Fenster für eine beweisbare Teilaussage in Kopfform (00:26–00:34) — Ergebnis: keine
Geprüft (Skripte in meinem Notizordner, Daten aus Gs spiral_analysis.json, nur gelesen):
- *Ist S_∞ eine bekannte Spezialfunktion?* S_1000 ist ganz vom Exponentialtyp ≈ 1.6 (|c_k|k! ~ 1.5^k). Kummer-Test: r_k = k c_k/c_{k−1} müsste eine
  Möbius-Funktion von k sein — Fehler 2.5e-3 (Daten genau auf ~1e-6 gegen S_1000); auch rationale Ansätze höheren Grades treffen nicht. **Kein ₁F₁.**
- *Erfüllt S eine lineare Differential- oder Dilatations-Differentialgleichung (q ∈ {μ, 1/μ, μ^{±a}, e^{2πia}, μ^κ})?* Kleinster Singulärwert fällt
  glatt mit der Zahl der Unbekannten, keine Lücke ⇒ **kein Befund** (S_1000 weicht von S_∞ um ~3e-3 ab; für eine Gleichungssuche bräuchte man S_∞).
- *Exakt lösbares Eckmodell mit komplexem μ?* Affin (auch mit M = 0) ist J obere Dreiecksmatrix ⇒ adjungierte Köpfe trivial (𝒜_n = w^n), kein
  Eckprofil; Gs gelöste quadratische Familie hat den Fixpunkt innen. **Kein Modell, das Ecke und adjungierte Struktur zugleich trägt.**
- Exakte Fixpunkt-Identität aus G39 bei s = 1/z_L: liefert (μ + 1/μ)𝒜_n(w_L) = (Schwanzsumme) − (ferne Urbilder) — das ist Ws Identität in L
  (27.3(c)), keine neue Schranke.
**Entscheidung:** Im Fenster keine beweisbare Teilaussage. Nach Auftrag der Leitung folgt die zweite Ideenrunde nur zur Ecke (IDEEN.md, Abschnitt
„Runde 2 — nur Ecke“), gemeinsam mit V.

## 13. Erste Tests zu R2, R5, R12 (Auftrag Leitung; 12.10. 00:36–00:47)

### 13.1 R2 Renormierung — Test bestanden (numerisch): Gs Korrekturordnungen sind Eigenwerte eines linearisierten Skalierungsoperators
Ohne Grenzprofil: Δ_n := S_{2n} − S_n (Gs normierte Eckprofile, 15 Taylor-Koeffizienten; spiral_analysis.json, nur gelesen, Auszug
k_kern/k14_profile_G_kubisch.npz). Liegt eine linearisierte Abbildung ℛ: S_n ↦ S_{2n} nahe S_∞ vor, so gilt Δ_{2n} = 𝒯Δ_n mit den Eigenwerten
2^{−α_i}. Kleinste-Quadrate-Bündel aus den Paaren (Δ_n, Δ_{2n}) (k14_r2_renorm.py, k14.log):

| Paare | Eigenwerte | erwartet |
|---|---|---|
| n = 100, 150, 200 (3 Paare; Rest 2.8e-5) | 0.2890−0.0038i, 0.4955+0.0009i, 0.5451−0.0114i | 2^{−2a} = 0.2896−0.0102i, 2^{−1} = 0.5, 2^{−a} = 0.5382−0.0094i |
| zusätzlich n = 50 (4 Paare; Rest 4.3e-7) | 0.1906, 0.2749, 0.4970, **0.5406−0.0094i** | dto., vierter ≈ nächste Ordnung |

Der relative Rest zeigt, dass die Abweichungen in einem 3–4-dimensionalen invarianten Unterraum liegen; die Eigenwerte treffen 2^{−a} (mit dem
richtigen Imaginärteil), 2^{−1} und 2^{−2a} auf 0.2–1.3 %. Das ist eine UNABHÄNGIGE Bestätigung von Gs Korrekturordnungen (G hat feste Potenzen
angepasst; hier kommen die Eigenwerte aus einer freien linearen Abbildung). Freie Exponenten-Anpassung (k14b) trennt a und 1 im Bereich n ≤ 1000 nicht
(effektiv 0.916 und 1.49); die Amplitude B_n bleibt unscharf (Differenzenquotienten 0.56/0.44/0.43). Deutung: 2^{−1} ist die triviale
Umparametrisierung n → n + const (ζS'), 2^{−a} die erste echte Eckkorrektur, 2^{−2a} ihre zweite Ordnung — wie „irrelevante Richtungen“ an einem
Renormierungs-Fixpunkt. **Entscheidung: R2 weiterführen.** Nächster Schritt: ℛ analytisch aus einer Gleichung für die Köpfe gewinnen (Kandidat:
Gs Volterra-Rekursion G22 auf der Kegelskala), dann Kontraktion in W_{1.05} zeigen; die Amplitude kommt aus der Richtung der Konstanten.

### 13.2 R5 IIKS — verworfen in der einfachen Form (bewiesen, elementar)
Gs Rang-2-Struktur G40 ist gleichwertig zu J_n = 𝒦_n + T_n(M) mit der Krylov-Matrix 𝒦_n = [T_nF^j]_j der Multiplikation mit F und der unteren
Toeplitz-Matrix T_n(M). Bezüglich (T_F, S): T_F𝒦_n − 𝒦_nS = T_n(F^{n+1})e_nᵀ (Rang 1), aber T_FT_n(M) − T_n(M)S = T_n(φ), φ := M(F − z) — voller Rang;
bezüglich (T_F, T_F) ist es umgekehrt. Es gibt also kein Sylvester-Paar mit kleinem Verschiebungsrang für J_n; G40 ist der Kommutator
[S, T_n(φ) + T_n(p_{n+1})e_nᵀ] (zweite Ordnung). IIKS braucht erste Ordnung ⇒ in dieser Form nicht anwendbar (nur über eine 2×2-Block-Linearisierung,
teuer ⇒ zurückgestellt). **Nebenbefund:** φ = M(F − z) verschwindet genau in den Fixpunkten — die Ecke ist die Stelle, an der die natürlichen
Verschiebungsstrukturen beider Kanäle zusammenfallen (Hinweis für R1).

### 13.3 R12 Quadranten-Wiener–Hopf — Prämisse halb bestätigt, Weg teurer als gedacht (k13_koepfe_2d.py, k13.log)
Alle Köpfe bis N' = 300 aus einer Inversen per Rang-1-Abbau (nur Mittelpunkte weitergeführt). Kontrollen: x-Kopf 1.728 bei n = 200 (= G),
D_N aus den Köpfen: 0.4930/0.4848 bei N = 150 (= k05: 0.492/0.484) — Satz K5(3) numerisch bestätigt.
- **Skalenfrei:** e_N(λ) := ‖P_Nc_{λN}‖·‖P_Nr̂_{λN}‖ ist für λ ≥ 1.2 nahezu N-unabhängig (λ = 1.5: 0.0044/0.0045/0.0046/0.0047; λ = 2: 0.00040–0.00042
  für N = 80…150); die Indexschwerpunkte liegen fest bei ξ ≈ 0.81–0.86 (c) und 0.91–0.93 (r̂).
- **Aber:** Weil e_N(λ) nicht wie 1/N fällt, wächst Σ_{n>N}‖Kopf‖ proportional zu N (≈ 3 bei N = 150), während ‖D_N‖ ≈ 0.49 bleibt — D_N entsteht
  durch eine mit N wachsende Auslöschung. Der Kern in (ξ, λ) ist also OSZILLIEREND (Gs schneller Phasenträger ω ≈ −1.485 je Grad); ein
  Quadranten-Wiener–Hopf bräuchte erst eine Phasenreduktion (Fourier-Integral-Form).
- **Entscheidung: R12 zurückstellen** (oder als Teil von R2 führen, wo die Phase über Gs χ_n schon abgespalten ist).

## 14. R2 als Hauptstrang: der Weg über G22 stösst auf eine harte Stelle; Amplitude und Fast-Orthogonalität (Cotlar–Stein); R1 aufgenommen (12.10. 00:48–00:55)

### 14.1 G22 ist Vorwärts-Substitution mit dem LU-Faktor L (bewiesen, eine Zeile)
Aus J̃ = LU (Satz K5) folgt J̃c_j = Le_j, also γ_nj/h_j = L_nj, und Gs G22b ist 𝒜 = L⁻¹(w^m)_m: die Zeilen von L⁻¹ durch Vorwärts-Substitution.
In den Köpfen: 𝒟_n(w) = (q_n/h_n)w^n − Σ_{j<n}K_nj𝒟_j(w) mit K_nj = (q_nh_j)/(h_nq_j)·L_nj — DERSELBE Kern für alle w; die Quelle ist für |w| ≤ 1
exponentiell klein (|Quelle| = 5.7e-4, 5.4e-7, 6.3e-10 bei n = 100, 200, 300), die Rekursion also praktisch homogen.

### 14.2 Harte Stelle: der Kern ist exponentiell gross, die Köpfe entstehen durch Auslöschung (k15_volterra_kern.py, k15b_wkern.py)
| n | Σ_j|K_nj| (x) | Maximum bei t = j/n | Σ_j|K̂_nj| (w-Koordinate, M = C_kL) | Maximum |
|---|---|---|---|---|
| 150 | 2.1e9 | 0.847 | 7.1e8 | 0.853 |
| 200 | 2.8e12 | 0.850 | 6.6e11 | 0.855 |
| 300 | 5.3e18 | 0.847 | 6.0e17 | 0.853 |
Wachstum ≈ e^{0.14n} in BEIDEN Koordinaten; schon |K_{n,n−1}| ≈ 0.11·n. Die Köpfe (√n|𝒟_n(L)| = 1.159, 0.863, 0.748 — gleich Gs B_n) sind also das
Ergebnis einer exponentiellen Auslöschung; **eine Kontraktion kann aus G22 nicht kommen**, auch nicht nach der Linse. Die Lage t* ≈ 0.85 ist
skalenfrei, aber NICHT 1/|μ|: mit |μ| = 1.30 bzw. 1.08 (k15, gleiches arg μ) liegt sie bei 0.860 bzw. 0.840 (1/|μ| wäre 0.769 bzw. 0.926) —
eine globale Grösse (Steg-Band wie in §11.2), keine lokale Gradstreckung.
**Folgerung für R2:** Die Renormierungs-Eigenwerte (§13.1) sind eine Eigenschaft der LÖSUNG, nicht einer kontrahierenden Rekursion. Ein analytisches ℛ
muss aus einer stabilen Darstellung kommen (Kandidaten: die Rückwärts-/Zukunftsform D_N = −Σ_{n>N}(projizierte Köpfe), oder die lokalen
Momentbedingungen W §25.1(a) direkt auf der Kegelskala). Das ist die harte Stelle; R2 bleibt als Struktur (Eigenwerte) bestehen.

### 14.3 Wie die Amplitude eingeht, und Fast-Orthogonalität als Cotlar–Stein-Bedingung (k16_cotlar.py, k16.log)
- **Amplitude:** Für das unnormierte Profil P_n(ζ) = √n𝒟_n(w_L + d_nζ) = B_n·S_n(ζ)·χ_n ist B_n die Richtung der Konstanten am Fixpunkt; B_n beschränkt
  heisst „Eigenwert des linearisierten ℛ in dieser Richtung = 1 (marginal)“, B_n konvergent heisst zusätzlich „übrige Eigenwerte < 1“ (gemessen 2^{−a}, 2^{−1},
  2^{−2a}). Die Differenzenquotienten von B_n (0.56/0.44/0.43) passen nicht sauber dazu — die Amplitude trägt eine zusätzliche, langsamere Korrektur.
- **Fast-Orthogonalität:** mit den w-Kopftermen T_n = c_n ⊗ r̂_n gilt nach Cotlar–Stein
  ‖Ĵ_N⁻¹‖ = ‖Σ_{n≤N}T_n‖ ≤ (A₁A₂)^{1/2}, A₁ = sup_n Σ_m (|⟨c_n,c_m⟩|‖r̂_n‖‖r̂_m‖)^{1/2}, A₂ = sup_n Σ_m (‖c_n‖‖c_m‖|⟨r̂_n,r̂_m⟩|)^{1/2}.
  Gemessen (Köpfe bis 300): **A₁ = 8.7 / 11.3 / 12.6, A₂ = 3.92 / 4.09 / 4.14** für N = 100 / 200 / 300 (wahre Norm 1.2512; Σ‖T_n‖ = 40 / 53 / 59).
  Die adjungierten w-Köpfe sind fast orthogonal (|cos(r̂_n, r̂_{n−k})| = 0.06 bei k = 3, 0.003 bei k = 40); die primalen Köpfe c_n korrelieren stark bei
  Nachbarn mit gerader Differenz (|cos| = 0.75–0.94 bei k = 2, 0.47–0.66 bei k = 10) — das ist Gs schneller Phasenträger (ω ≈ −1.485 je Grad, zwei
  Schritte ≈ −π). A₁ wächst langsam (etwa logarithmisch), A₂ sättigt.
  **Damit hat die polynomiale Eckaussage eine prüfbare Form:** (S_γ) folgt aus „A₁, A₂ beschränkt“, und schon „A₁ ≲ log N“ gäbe ‖Ĵ_N⁻¹‖ ≲ (log N)^{1/2} —
  für Satz K1 (κ < 1.79) mehr als genug. Das ist eine Aussage über Paar-Skalarprodukte der Köpfe, also über ihre Eckprofile auf zwei Skalen n, m und
  die Phasen — genau das, was R2 (Profile auf der Kegelskala) liefern soll.

### 14.4 R1 aufgenommen — Präzisierung
Im Log-Koordinaten z = z_Le^u, s = s_Le^v ist σ ≈ 1/(1 − e^{v+μu}) + μ/(1 − e^{v+u}) (exakt in der Koenigs-Koordinate). Auf dem Torus (u, v ∈ iℝ) ist die
Identitäts-Singularität die Diagonale v + u ∈ 2πiℤ, die Kompositions-Singularität v + μu ∈ 2πiℤ trifft den Torus wegen μ ∉ ℝ NUR im Eckpunkt.
Reelles Aufblasen dort: homogen vom Grad −1, Winkelfaktor g(φ) = i[1/(sin φ + μ cos φ) + μ/(sin φ + cos φ)] mit reellem Pol in Identitätsrichtung und
einem komplexen Pol tan φ = −μ. Die Ecke ist damit eine Quadranten-Ecke (n, n) des Indexquadrats mit homogenem Symbol — Malyuzhinets-Typ. Exakt gilt
a = 1/κ₁ = 1 − log μ/(2πi), d. h. n^{−a} = n^{−1}·exp(log μ·log n/(2πi)). **Nächster Schritt R1:** die Winkelgleichung aufstellen und prüfen, ob ihr erster
Exponent a ist; erst dann lohnt der Vergleich mit Gs Profil.
**Vorbehalt (00:55):** σ(z, s) ist die ERZEUGENDE Funktion der Matrixeinträge (J̃[k][j] = [z^k s^j]σ), kein Faltungssymbol — J̃ ist keine
zweidimensionale Faltung auf dem Quadranten. Die Singularität im Eckpunkt bestimmt die langsam abfallenden Anteile der Einträge in ALLEN (k, j), die
Ecke (n, n) entsteht erst durch den Abschnitt. Das Quadranten-Bild ist daher nur heuristisch; vermutlich reproduziert die Winkelgleichung lediglich das
bekannte κ₁ aus Os/Ws Kegelmodell. R1 bleibt zweiter Strang mit niedriger Priorität, bis V bewertet hat.

## 15. Hauptstrang Cotlar–Stein für die Köpfe (12.10. 00:56–01:07): Herkunft der Abnahme, schwächste Aussagen, zwei bewiesene Lemmata, Grenze

### 15.1 Woher die Abnahme der Paar-Skalarprodukte kommt (k16–k20; Köpfe kubisch bis 300)
- **Adjungierte w-Köpfe r̂_n — echte Abnahme in |n − m| aus einer FESTEN, nicht-stationären Phase.** Mit G_n(θ) = Σ_i r̂_{n,i}e^{−iiθ} gilt
  G_n ≈ a_n(θ)e^{inΦ(θ)}: das Phaseninkrement arg(G_n/G_{n−1}) ist n-unabhängig (Fenster n ∈ [200,230] und [270,300]: Abweichung ≤ 4e-4 rad),
  Φ' ∈ [−0.961, −0.557] auf dem ganzen Träger (keine stationären Punkte; mittlere Steigung −0.905), die Amplituden a_n = G_ne^{−inΦ} sind glatt
  (95 % der Fourier-Masse bei |k| ≤ 5, 99.9 % bei |k| ≤ 60) und ändern ihre Form langsam (‖â_250 − â_300‖ = 0.089). Gemessene Dekorrelation:
  max|cos(r̂_n, r̂_{n−k})| = 0.65, 0.10, 0.03 (k = 1, 2, 3), Hülle ≈ 2/k bis k ≈ 60, danach < 0.03.
- **Primale Köpfe c_n — KEINE Abnahme in |n − m|.** Sie sind kohärent mit einem schnellen Phasenschritt (|cos(c_n, c_{n+2})| = 0.75–0.94, also
  c_{n+2} ≈ −c_n), ihre Korrelation fällt nur mit dem Skalenverhältnis m/n (0.47 bei k = 40, 0.23 bei k = 80 für n = 250).
- **Eckanteil:** Die adjungierten Amplituden tragen eine Eckspitze bei θ_c = ±81.1° (Massenanteil in ±3n^{−Re a}: 3.3 % → 6.0 % für n = 150 → 300,
  wachsend); dort dekorrelieren die Köpfe im Absolutbetrag nur über das Skalenverhältnis. Die Summe der Eckanteile löscht stark aus:
  ‖Σ_{n≤N}T^E_n‖ = 0.52 → 0.63 bei Σ‖T^E_n‖ = 4.9 → 8.0 (N = 150 → 300) — eine VORZEICHEN-Auslöschung, keine absolute Fast-Orthogonalität.
  Ein Oszillator mit einer oder zwei Frequenzen erklärt sie nicht (Defekt ≈ 1.3–2.0, k17, k20): Primal- und Adjungiertenphase mischen sich.
**Antwort:** Die Abnahme kommt adjungiert aus der Phase (Dephasierung), primal aus der Skalenverschiebung; an der Ecke trägt nur die Vorzeichen-
Auslöschung, die Cotlar–Stein mit Absolutbeträgen nicht sieht.

### 15.2 Bewiesen
**Satz K7 (Cotlar–Stein-Kopfkriterium; Cotlar–Stein-Lemma, eine Zeile).** Für die w-Kopfterme T_n = c_n ⊗ r̂_n, τ_n = ‖T_n‖:
  ‖Ĵ_N⁻¹‖ ≤ (A₁A₂)^{1/2}, A₁ = sup_n Σ_{m≤N}(τ_nτ_m|cos(c_n,c_m)|)^{1/2}, A₂ = sup_n Σ_{m≤N}(τ_nτ_m|cos(r̂_n,r̂_m)|)^{1/2};
insbesondere (mit |cos(c_n,c_m)| ≤ 1) **‖Ĵ_N⁻¹‖ ≤ τ̄·((N+1)·B)^{1/2}, B := sup_n Σ_{m≤N}|cos(r̂_n, r̂_m)|^{1/2}, τ̄ = sup τ_n.**
**Lemma K8 (Dephasierung; partielle Integration).** Zerlege G_n = a_ne^{inΦ} + E_n mit fester Phase Φ ∈ C²(𝕋), |Φ'| ≥ c₀ > 0 auf 𝕋,
‖a_n‖_{C²} ≤ A‖r̂_n‖ und ‖E_n‖_{L²} ≤ √f_n·‖r̂_n‖. Dann
  |cos(r̂_n, r̂_m)| ≤ C(A, Φ)/(c₀²(n − m)²) + √f_n + √f_m + √(f_nf_m).
*Beweis.* Parseval (r̂_n sind endliche Koeffizientenvektoren): ⟨r̂_n, r̂_m⟩ = (1/2π)∫G_nḠ_m dθ. Glatter Teil: e^{i(n−m)Φ} = ∂_θe^{i(n−m)Φ}/(i(n−m)Φ'),
zweimal partiell integrieren auf dem Kreis (keine Randterme), Schranke durch ‖a_nā_m/Φ'‖_{C²}; Mischterme mit Cauchy–Schwarz. ∎
**Folgerung (schwächste Fassungen):**
- κ = 1 (reicht für Satz K1 und die Familie, κ < 1.79): **nur (H1) τ_n ≤ τ̄** — Dreiecksungleichung (Gs Summenkriterium in L²(γ)-Form).
- κ = 1/2 (reicht für das N₀ ≤ 400-Ziel κ + 3/2 ≤ 2): (H1) + Lemma K8 mit Σ_m f_m^{1/4} beschränkt (Eckanteil der adjungierten Köpfe klein) —
  **nur Phasen- und Lagekontrolle der ADJUNGIERTEN Köpfe, nichts über die primalen, keine Amplitudengrenzwerte.**
- Numerisch erfüllt: Φ fest und nicht-stationär (s. o.), τ̄ = 0.67 (bei n = 4), A₂ = 3.92 / 4.09 / 4.14 (N = 100 / 200 / 300).

### 15.3 Grenze (Entscheidung): absolute Fast-Orthogonalität reicht asymptotisch nicht unter κ = 1
Die Eck-Massenanteile f_n wachsen (3–6 %, nicht summierbar klein), und τ_n sättigt (0.015 bei N = 800–1000, k12); dann wachsen die Eck-Eck-Beiträge
zu A₂ und die kohärenten Primalbeiträge zu A₁ linear in N, und (A₁A₂)^{1/2} ist asymptotisch linear (mit kleiner Konstante ≈ τ·f^{1/2} ≈ 0.003).
**Cotlar–Stein mit Absolutbeträgen ist damit als Weg zu κ < 1 erschöpft**; die beobachtete Beschränktheit ist eine phasenaufgelöste Auslöschung an der
Ecke. Was dafür nötig wäre: die gemeinsame Phase von Primal- und Adjungiertenkopf an der Ecke als WKB-Kontinuum e^{inΨ}, mit langsam
veränderlichen Hüllen und ohne Resonanz Ψ ≡ 0 mod 2π auf dem Träger — schwächer als starke Asymptotik (keine Amplitudengrenze), aber eine volle
Phasenaussage; nicht bewiesen.
**Einordnung:** Für den Hauptsatz genügt κ = 1, also allein (H1). Die Fast-Orthogonalität ist eine Verfeinerung (für N₀ ≤ 400 bzw. (S_γ)); der Kern
bleibt (H1): ‖c_n‖·‖r̂_n‖ ≤ C, die exakte Auslöschung der beiden Kopfraten.

### 15.4 Nachtrag 01:16: Köpfe bis N = 600 (k13c, k16b, k21b) — die Grenze bestätigt sich, und warum
- Kontrolle: D_200 aus den Köpfen bis 600 = 0.5123 0.5057 0.2345 0.2232 0.0917 0.0875 = k09 (D_200 = 0.51234 0.50566 …) — Satz K5(3) auf 4 Stellen.
- Cotlar–Stein bis 600: A₁ = 8.7 / 11.3 / 12.6 / 13.4 / 14.0 / 14.4, A₂ = 3.92 / 4.09 / 4.14 / 4.17 / 4.19 / 4.20 (N = 100 … 600) — im gerechneten Bereich
  sehr langsam wachsend, (A₁A₂)^{1/2} ≈ 7.8 bei N = 600 (wahre Norm 1.2512). Die Zuwächse sind klein, weil τ_n hier noch fällt (0.045 → 0.018).
- **Warum es asymptotisch nicht hält — Laplace-Konzentration an der Ecke:** die Wachstumsrate je Grad der adjungierten Köpfe über dem Winkel,
  Λ(θ) = log|G_n/G_{n−1}|, hat ihr Maximum 0.0656–0.0659 bei θ ≈ ±79° (Ecke), gegen ≈ 0.060 längs Γin; der massengewichtete Mittelwert steigt
  0.0611 → 0.0617 → 0.0632 → 0.0642 (n ≈ 215 … 585) gegen das Maximum. Der Eckanteil wächst also wie e^{≈0.005n} (passt zu 3.3 % → 6.0 % zwischen 150 und 300),
  und mit ihm wird der adjungierte Kopf kohärent: |cos(r̂_n, r̂_{n−k})| bei k = 2 / 10 / 20 steigt von 0.05 / 0.15 / 0.09 (n = 250) auf 0.64 / 0.65 / 0.43
  (n = 550). Sobald τ_n sättigt (0.015 ab N ≈ 800), wachsen A₁ und A₂ linear — Cotlar–Stein mit Beträgen wird asymptotisch so gut wie die
  Dreiecksungleichung (κ = 1).
- **Kopf-Gleichgewicht an der Ecke:** die lokalen Normraten nähern sich der Eckrate: ‖r̂_n‖ +0.0638, ‖c_n‖ −0.0657 (n = 500 … 600), aus k12 zwischen
  800 und 1000: +0.0651 / −0.0653, Summe −0.0002. (H1) ist damit asymptotisch genau die Aussage „Eckrate des adjungierten Kopfes = Abklingrate des
  primalen Kopfes, mit beschränktem Vorfaktor“ — eine reine Eckaussage, gleichwertig mit Gs Amplitudenschranke B_n ≤ C_A samt Profilschranke.

## 16. H-A, erster Test (12.10. 01:23–01:27): die Kopf-Paarung ist KEINE invariante Dichte an der Ecke — die Köpfe sind in komplementären Sektoren homogen

**Abelkoordinate (bewiesen, KERNPROBLEM §3.1, gilt für jedes f mit lokaler Abelfunktion α an L):** φ := 1/(α'∘f) erfüllt φ∘f = m·φ, also
J_∞(φu) = (1/α'∘f∘f)·(u∘f + u). **Invariante Paarung (bewiesen, elementar):** Erfüllen δ und Λ nahe L die homogenen Eckgleichungen δ∘f = −mδ und
ℒΛ + mΛ = 0 (ℒΛ(y) = yΛ(g(y))g'(y)/g(y), g = f⁻¹ lokal), so ist die 1-Form ω = m·Λ·δ·dy/y unter f invariant (in s: ½-periodische Dichte).
Koordinatenfreier Test ohne Koenigs-Funktion (also ohne deren Genauigkeitsproblem): Q(y) = ω(f(y))/ω(y) ≈ 1? (k22_invariante.py, arb mit 1.3n + 400 Bit,
Köpfe direkt gelöst: J̃_nc = e_n, J̃_nᵀℓ = e_n; k22.log, k22b.log).

| Lage y = L + r·e^{iθ} | primaler Defekt |δ∘f + mδ|/|mδ| | adjungierter Defekt |ℒΛ + mΛ|/|mΛ| | Q |
|---|---|---|---|
| innen, θ = 240°, r = 0.12 (n = 200/300/400) | 1.6e-3 / 6.2e-6 / 1.7e-8 | ≈ 1.00 | 0.99 / 0.81 / 0.67 (Phase dreht) |
| innen, θ = 200° bzw. 280° | 1e-3 … 1 | 1.0 … 3.9 | 0.0000 … 6e6 |
| aussen, θ = 60°/100°, r = 0.12 (n = 300) | 81 / 1.0 | 0.046 / 0.008 | 80 / 0.009 |
| nahe L, r = 0.03 (alle θ) | 0.3 … 12 | 0.3 … 11 | 0.06 … 6 |

**Befund:** Der primale Kopf ist im INNEREN Ecksektor eine sehr genaue homogene lokale Lösung (Defekt 1.7e-8 bei n = 400 — die Quelle (x/ρ)^n ist dort
exponentiell klein und der Überlauf löscht sie aus); der adjungierte Kopf ist nur im ÄUSSEREN Sektor homogen (innen dominiert seine Quelle (ρ/t)^n);
im Eckbereich selbst (r = 0.03) ist keiner homogen. Es gibt also keinen Bereich, in dem beide die Eckgleichungen erfüllen — die Paarung ist dort keine
invariante Dichte, Q weicht um Faktoren 10^{±3} von 1 ab. **Nach dem Abbruchkriterium (Abweichung > 30 % bei n = 300–600) ist H-A in der Form
„Erhaltung von u·v auf dem Torus“ gescheitert.**
**Was bleibt (Strukturbefund):** Die Ecke ist eine GRENZFLÄCHE: innen lebt der primale Kopf als Summe lokaler Eckmoden (Lösungen von u(s+½) = −u(s), in
x: (x − L)^{s_k}-Typ), aussen der adjungierte; (H1) ist die Aussage, dass ihre Paarung über die Grenzfläche (den Kreis durch L, auf dem beide Quellen
Betrag 1 haben) gleichmässig nicht entartet. Das ist ein Transmissions-/Anpassungsproblem — dieselbe Stelle wie Os Eckformel und Ws starke Ebene;
die Abelkoordinate macht beide Seiten zu Problemen mit konstanten Koeffizienten, die Kopplung über die Grenzfläche aber nicht trivial.

## 17. Weg (a) „explizite lokale Parametrix“, Schritt a1 (12.10. 01:30–01:38): Stokes-Struktur der Köpfe an der Ecke und ihre Aussenexponenten

**Entscheidung (PLAN 01:30):** (a) zuerst. Das Vergleichsproblem ist Os (P1)–(P3) (o_operator §33.4), in Kopfform; es ist die lokale Hälfte von Ws (b).

### 17.1 Drei Stokes-Sektoren, skalenfrei auf der Kegelskala (k23_sektoren.py, k23.log; kubisch, n = 300 und 600, Köpfe in arb)
Defekte der homogenen Eckgleichungen auf Kreisen |x − L| = A·ρn^{−Re a} (A = 1, 4, 16), 24 Richtungen: **die Tabellen für n = 300 und n = 600 stimmen bei
gleichem A praktisch überein** — die Ecke der Köpfe ist auf der Kegelskala skalenfrei.
- A = 1, 4 (innere Zone): keiner der Köpfe ist homogen (Defekte 0.3 … > 10).
- A = 16: primal homogen in θ ∈ [195°, 300°] (Defekt bis 4e-7), adjungiert in θ ∈ [55°, 125°] (bis 9e-4); dazwischen „quellbestimmte“ Sektoren
  (Defekt ≈ 1.0: der Kopf ist dort Quelle/m) und ein „überlaufbestimmter“ Sektor um die Γout-Tangente (−45° … 45°, Defekt > 10: Komposition dominiert).
- Vorhergesagte Grenzkeile (A → ∞): primal S_in = (166.6°, 308.3°), adjungiert S_out = (24.9°, 166.6°), je Öffnung πα′ = 141.7° (α′ = 1 − arg μ/π),
  Lücke der Öffnung arg L um die Γout-Tangente. Gemessen bei A = 16 enger (110° bzw. 70°), zentriert passend.
**Deutung (Blickwinkel der Leitung bestätigt):** L ist ein Wendepunkt mit drei Stokes-Sektoren je Kopf — homogen / quellbestimmt / überlaufbestimmt; der
Kreis durch L (Tangente bei 166.6°/−13.4°) ist die Stokes-Linie, an der die Quelle (x/ρ)^n bzw. (ρ/t)^n zwischen exponentiell klein und gross wechselt.

### 17.2 Aussenexponenten in den homogenen Keilen — explizit durch L (k24_potenzen.py, k24.log)
Lokale Steigung von log|Kopf| längs des Mittelstrahls des homogenen Keils, A = 10 … 76:
| Kopf | gemessen (n = 600, A ≈ 10–13) | Vorhersage (Mellin-Wurzel der lokalen Gleichung) |
|---|---|---|
| primal, θ = 240° | Re s = −3.55 … −3.56, Im s ≈ −1.00 … −1.03 | s_{−1} = 1 − πi/log μ = 1 − 2πi/L = **2 − s₀** = −3.450 − 1.054i |
| adjungiert, θ = 90° | Re s = −6.40 … −6.42, Im s ≈ −1.14 … −1.54 | σ₀ = −2 − πi/log μ = **−1 − s₀** = −6.450 − 1.054i |
(s₀ = 1 + 2πi/L = Knesers Singulärexponent; lokale Gleichungen: δ(μz) = −μδ(z) ⇒ μ^{s−1} = −1; adjungiert μ^{−σ−2} = −1. Die Drift für grössere A
ist die Krümmung des Kreises und die nächste Mode.) **Die Anschlussdaten des Vergleichsproblems sind damit explizit**: der primale Kopf geht im
inneren Keil in z^{2−s₀} über, der adjungierte im äusseren Keil in z^{−1−s₀} — beide über log μ = L/2 durch L bestimmt (log L = L).
**Lemma K9 (Struktur homogener Keillösungen; Skizze eines elementaren Beweises):** Erfüllt U(z) auf einem Keil, dessen Öffnung grösser als arg μ ist,
U(μz) = −μU(z) und wächst polynomial, so ist U eine endliche Summe Σ c_k z^{s_k} (in w = log z ist U(w)e^{−s w} periodisch mit Periode log μ, setzt sich
zu einer ganzen periodischen Funktion fort, endlicher Typ ⇒ endliche Fourier-Summe). Entsprechend adjungiert mit σ_k.

### 17.3 Was das Vergleichsproblem jetzt ist (Formulierung für a2)
Innere Zone |z| ≲ 4 (Kegelskala): beide Köpfe inhomogen; die Quelle ist dort e^{±Kz}·(1 + langsam veränderlich) mit K ∝ n^{1−a}, 1 − a = log μ/(2πi)
(|K| wächst nur wie n^{0.106}), der Überlauf bzw. Rest T trägt denselben Faktor. Äussere Zone: homogene Keile mit den Exponenten 2 − s₀ (primal) und
−1 − s₀ (adjungiert), quellbestimmte und überlaufbestimmte Keile dazwischen. Gesucht: Lösung der Kegelgleichung (P1) mit Stokes-Quelle, die in den drei
Sektoren die richtige Form annimmt — Mellin in z macht die Gleichung diagonal (Faktor 1 + μ^{s−2}, Gammafunktionen), die Sektorbedingungen werden zu
Analytizitätsstreifen; das ist ein Malyuzhinets-Problem mit expliziten Daten.

### 17.4 a2, erste Schritte (bis 01:47)
- **Quelle und Rest heben sich auf (k25_rest_quelle.py, k25.log):** In der exakten adjungierten Kopf-Identität H = 𝒰 + T (H = adjungierter Operator auf Λ,
  𝒰 = (ρ/ζ)^n/ζ, T = Σ_{j>n}γ̃_j(ρ/ζ)^j/ζ) ist |H/𝒰| ≈ 1e-13 … 4e-9 auf der ganzen Kegelskala (A = 1, 4, 16; n = 300), gleich gross wie |Λ/(ζ𝒰)|: Quelle und
  Abschneiderest löschen sich auf die Grösse des Kopfes selbst aus. Die wirksame innere Gleichung ist daher keine Gleichung mit bekannter Quelle, sondern eine
  **Wiener–Hopf-Zerlegung nach Graden** (Λ: Grade ≤ n, H: Grade ≥ n, wegen der Biorthogonalität ⟨Λ, J_∞x^j⟩ = 0 für j < n).
- **Zwei Skalen, getrennt nur durch n^{0.106}:** auf der Kegelskala wird Λ von Graden k ~ n^{a} getragen (Os Befund: die grossen Koeffizienten bei k ≈ 0.95n
  löschen sich in L aus); die Abschneidekante liegt in der Frequenz λ = k·n^{−a} bei K = n^{1−a}, |K| ~ n^{0.106} → ∞ sehr langsam. Ein n-unabhängiges
  Vergleichsproblem ist der Grenzfall K → ∞ (scharfe Stokes-Schalter); die Messungen (k23: Tabellen n = 300 und 600 gleich; Gs schnelle Profilkonvergenz)
  sprechen dafür, dass die normierten Profile die K-Abhängigkeit nicht spüren.
- **Ordnung des Eckprofils (k26_ordnung_profil.py):** Gs S_1000 ist nicht vom Exponentialtyp, sondern von Ordnung > 1: der Ansatz c_k ∝ q^k/Γ(ak + b) mit
  a = 1/κ₁ macht die q-Schätzer über k = 4 … 38 auf 0.8 % konstant (a = 1: 1.2 %), frei angepasst a = 0.935 + 0.025i (Imaginärteil = Im(1/κ₁) = 0.0253),
  Streuung 0.13 %. Das passt zu Os WKB-Phase N·c₁ζ^{κ₁} (Ordnung κ₁ = 1/a ≈ 1.12). Eine einzelne Mittag-Leffler-Funktion E_{a,b} trifft die niedrigen
  Koeffizienten aber nicht (6–50 % Fehler) — das Vergleichsprofil ist eine Funktion dieser Ordnung, keine bekannte Einzelfunktion.
- Kontrolle: Kneser-Profil (n = 600) gegen kubisch (n = 1000) im Wiener-Abstand 3.7e-3 (Universalität).

### 17.5 Stand von (a) und nächster Schritt
Bewiesen: Abelkoordinate (konstante Koeffizienten), invariante Paarungsform, Lemma K9 (Skizze). Numerisch, belastbar: drei Stokes-Sektoren je Kopf,
skalenfrei auf der Kegelskala; Aussenexponenten **2 − s₀ (primal) und −1 − s₀ (adjungiert)**, s₀ = 1 + 2πi/L; Quelle/Rest-Auslöschung; Ordnung κ₁ des
Profils. **Nächster Schritt (Tage):** das Grenzproblem K → ∞ als Riemann–Hilbert-Aufgabe auf den drei Stokes-Strahlen durch L aufstellen (Sprünge aus
den Sektorrelationen: homogen ↔ quellbestimmt ↔ überlaufbestimmt, Dilatation μ), per Mellin lösen (Gammafunktionen, Malyuzhinets), Lösung gegen S_∞ und
die Aussenexponenten prüfen. Erst danach lohnt die Verbindungsformel mit Ws äusserer Lösung (b).

## 18. Zweiter K-Lauf: das Grenzprofil des adjungierten Kopfes ist explizit — eine Mellin–Barnes-Funktion mit dem Malyuzhinets-Faktor von 1 + μ^{s−2} (12.10. 01:50–02:14)

### 18.1 Herleitung des Kandidaten (Skizze, die Bausteine einzeln elementar)
Auf der Kegelskala gilt für den adjungierten Kopf y (Erzeugende Λ, ganz in der lokalen Variablen u) die lokale Identität
  **H(u) = μ y(u) + μ⁻¹ y(u/μ)**  (H = 𝒰 + T, ferne Urbilder und die Variation von m, ζ/x vernachlässigt; O §33.4, k25).
Mellin–Barnes-Ansatz y(u) = (1/2πi)∫_C Φ_y(s) u^{−s} ds ⇒ Φ_H(s) = μ(1 + μ^{s−2})Φ_y(s). Drei Forderungen:
(i) y ganz mit y(0) = 1 ⇒ Φ_y = Γ(s)·Ψ(s), Ψ ohne Pole links (die Pole von Γ in −ℕ₀ geben die Taylor-Reihe);
(ii) H ganz und im homogenen Keil ohne Potenzasymptotik (dort ist H nach k23/k25 klein gegen y) ⇒ Φ_H = Γ(s)·Ξ(s), Ξ ohne Pole rechts;
(iii) y fällt im homogenen Keil (keine wachsenden Exponenten).
Mit der Gamma-Zerlegung des Symbols, α := log μ/(2π) = L/(4π) (log μ = L/2):
  **μ(1 + μ^{s−2}) = K_L(s)·K_R(s),  K_L = 2πμ^{s/2}/Γ(½ − iα(s−2)),  K_R = 1/Γ(½ + iα(s−2)),**
K_R hat genau die Nullstellen s_m = 2 + (2m+1)·2πi/L (m ≥ 0, rechts), K_L genau die s_m mit m ≤ −1 (links). Aus Ξ = K_LK_RΨ folgt
Ξ/K_L = K_RΨ = E(s) ganz (linke Seite rechts analytisch, rechte Seite links analytisch) — Wiener–Hopf im Mellin-Raum. Die kleinste Wahl
E(s) = const·λ^{−s} ergibt den Kandidaten
  **Y(u) := (1/2πi)∫ Γ(s)Γ(½ + iα(s−2)) u^{−s} ds / Γ(½ − 2iα) = Σ_{k≥0} Γ(½ − iα(k+2))/Γ(½ − 2iα) · (−u)^k/k!,**
  **H_Y(u) = μY(u) + μ⁻¹Y(u/μ) = (2π/Γ(½ − 2iα)) Σ_{k≥0} (−u/√μ)^k / (k!·Γ(½ + iα(k+2))).**
Eigenschaften (aus der Formel): ganz von der (komplexen) Ordnung κ₁ = 1/a, a = 1 + iα = 1/κ₁ — denn Γ(s)Γ(iαs) hat dasselbe s log s-Wachstum wie
Γ(as) (Stirling) —, Taylor-Koeffizienten ~ q^k/Γ(ak + b) (das ist der k26-Befund), im homogenen Keil Y ~ c·u^{−s₀}, s₀ = 2 + 2πi/L (Residuum des
ersten Pols von Γ(½ + iα(s−2))), und H_Y ohne Potenzasymptotik (Φ_H hat rechts keine Pole). Die Skala λ ist lokal frei (die lokale Gleichung ist
dilatationsinvariant) — sie ist die Anschlusskonstante an die äussere Lösung.

### 18.2 Test gegen Gs Eckprofile (k27_gamma_test.py, k28_extrapol_vergleich.py; Daten spiral_analysis.json, nur gelesen)
Skaleninvarianter Test: q_k := r_k/r_{k−1}, r_k = k c_k/c_{k−1}; für Y gilt q_k = Γ(β − iαk)Γ(β − iα(k−2))/Γ(β − iα(k−1))² mit β = ½ − 2iα, OHNE freien Parameter.
| k | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| Modell | 1.0253−0.0087i | 1.0211−0.0067i | 1.0180−0.0054i | 1.0156−0.0045i | 1.0138−0.0039i |
| G, n = 400 | 1.0205−0.0086i | 1.0166−0.0067i | 1.0138−0.0054i | 1.0116−0.0045i | 1.0098−0.0038i |
| G, n = 1000 | 1.0232−0.0086i | 1.0191−0.0067i | 1.0161−0.0054i | 1.0138−0.0045i | 1.0120−0.0038i |
Imaginärteile stimmen schon bei endlichem n auf 4 Stellen, die Realteile laufen wie n^{−Re a} auf das Modell zu (Abstand 0.0048 → 0.0026 von 400 → 800,
Quotient 1.85 = 2^{0.894}; extrapoliert 1.0253).
**Grenzprofil:** Extrapolation je Koeffizient mit c + A n^{−a} + B n^{−1} + C n^{−2a} (+ D n^{−a−1}), n = 300…1000; λ aus c_1:
**λ = −0.203297 − 1.406007i** (|λ| = 1.420629, arg λ = −98.227°; in allen vier Extrapolationsvarianten auf 6 Stellen gleich).
| | k = 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| S_∞ (extrapoliert) | 0.16895+1.25576i | −0.79197+0.22429i | −0.15003−0.33165i | 0.10287−0.06725i | 0.02268+0.02495i |
| Y(λζ) | 0.16895+1.25576i | −0.79197+0.22429i | −0.15003−0.33165i | 0.10287−0.06725i | 0.02268+0.02495i |
Relative Abweichung k ≤ 10: ≤ 1.3e-4; **Wiener-Abstand (Radius 1.05) Modell − S_∞ = 2.5e-7** bei Norm 3.86 (beste Extrapolation; die hohen k ≥ 20 sind
in der Extrapolation unsicher, tragen zur Wiener-Norm aber < 1e-8 bei). Ein Parameter (λ) für 40 komplexe Koeffizienten.
**Konvergenz bei endlichem n:** ‖S_n − Y(λ·)‖_{W,1.05} = 3.14·n^{−0.894} … 3.40·n^{−0.894} (n = 50 … 1000, kubisch) bzw. 3.04 … 3.37·n^{−0.894}
(Kneser, n = 50 … 600) — reines n^{−Re a}-Gesetz mit langsam veränderlicher Konstante (nächste Ordnungen n^{−1}, n^{−2a}); Kneser und kubisch haben
dasselbe Grenzprofil (Universalität).

### 18.3 Test gegen die Sektorstruktur der Köpfe (k29_modell_sektoren.py, k29.log; Vergleich mit k23.log)
Abbildung der x-Richtung θ (x − L = re^{iθ}) auf u: arg u = θ + arg(−λ) − 2 arg L + Im(a)·ln n, |u| = |λ|·A (A wie k23). Modell-Defekt |H_Y|/|μY| gegen den
gemessenen adjungierten Defekt |ℒΛ + mΛ|/|mΛ| des echten Kopfes (k23, n = 600, exakte Operatoren des kubischen Modells):
| A | θ = −45° | 0° | 45° | 60° | 75° | 90° | 105° | 120° | 135° | 180° | 240° | 270° | 285° | 315° |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1, Kopf | 1.3 | 2.1 | 2.1 | 2.0 | 1.8 | 1.7 | 1.6 | 1.5 | 1.4 | 1.1 | 0.9 | 0.9 | 0.9 | 1.3 |
| 1, Modell | 1.3 | 2.1 | 2.1 | 2.0 | 1.8 | 1.7 | 1.6 | 1.5 | 1.4 | 1.1 | 0.9 | 0.9 | 0.9 | 1.3 |
| 4, Kopf | 7.4 | >10 | 5.6 | 3.2 | 2.0 | 1.4 | 1.2 | 1.1 | 1.1 | 1.0 | 0.9 | 1.2 | 1.8 | 7.4 |
| 4, Modell | 7.2 | >10 | 5.7 | 3.2 | 2.0 | 1.4 | 1.2 | 1.1 | 1.1 | 1.0 | 0.9 | 1.2 | 1.8 | 7.2 |
| 16, Kopf | >10 | >10 | 0.3 | 9e-3 | 1e-3 | 9e-4 | 3e-3 | 4e-2 | 1.4 | 1.0 | 1.0 | 1.0 | >10 | >10 |
| 16, Modell | >10 | >10 | 0.3 | 8e-3 | 9e-4 | 6e-4 | 2e-3 | 4e-2 | 1.5 | 1.0 | 1.0 | 1.0 | >10 | >10 |
Die ganze Sektorstruktur (homogen / quellbestimmt / überlaufbestimmt) und ihre Lage werden reproduziert; das Modell kennt nur μ, β und λ (aus Gs Taylor-
Koeffizienten bei 0, nicht an diese Tabelle angepasst). Lokale Exponenten des Modells im homogenen Keil: −6.447 − 1.058i = −s₀ exakt ab |u| ≈ 30 (gemessen
am Kopf −6.40 … −6.42 bei A ≈ 10–13, wegen der nächsten Mode und der Kreiskrümmung). Achtung Bezeichnung: −s₀(hier) = −1 − s₀(§17) mit s₀(§17) = 1 + 2πi/L.

### 18.4 Einordnung
- **Numerisch, belastbar:** Gs Grenzprofil S_∞ (adjungierter Kopf auf der komplexen Kegelskala) ist die Mellin–Barnes-Funktion Y(λζ), λ = −0.203297 − 1.406007i;
  die Konvergenz S_n → Y(λ·) folgt n^{−Re a}. Os Ψ_Y (Eckmode, nicht Kopf; mit eigener Quelle 𝔲 = O(y)) ist eine ANDERE Funktion (q_k = 1.075, 1.056, …) —
  die Eckmoden lösen die inhomogene Aufgabe (P1) mit Quelle, der Kopf die homogene.
- **Bedeutung für (a):** Die Grenzaufgabe K → ∞ des adjungierten Kopfes ist gelöst, und zwar explizit: Malyuzhinets-Faktorisierung des Symbols
  μ(1 + μ^{s−2}) in Gammafunktionen, die Lösung ist eine Fox–Wright-Funktion. Offen: die Formulierung, in der Y die EINDEUTIGE Lösung ist (§19), der
  primale Kopf, und die Anschlusskonstante λ (global).

### 18.5 Der primale Kopf ebenso explizit, und die Dualität λ_P = e^{iπa}λ_A (k30–k32; 02:15–02:32)
**Primale Kopfgleichung lokal:** P(μz) + μP(z) = Quelle + Überlauf; Mellin: μ(1 + μ^{−s−1})Φ_P; rechte Nullstellen t_m = −1 + (2m+1)·2πi/L (m ≥ 0) =
Pole von Γ(½ + iα(s+1)). Dieselben drei Forderungen (ganz; Quelle+Überlauf ohne Potenzasymptotik im INNEREN homogenen Keil; Abfall dort) geben
  **P(u) = Σ_{k≥0} Γ(½ + iα − iαk)/Γ(½ + iα) · (−u)^k/k!,  Φ_P = Γ(s)Γ(½ + iα(s+1))/Γ(½ + iα)  (β_P = ½ + iα).**
Messung (k30, primaler Kopf in arb, n = 200…800; Taylor bei L in x mit d_x = −(L²/ρ)n^{−a}, so dass ζ_x = ζ_G in linearer Näherung; Pipeline-Kontrolle:
der adjungierte Kopf derselben Rechnung trifft Gs q_k auf alle Stellen):
- Das primale Profil ist auf der Kegelskala glatt und konvergiert (c_1 = −0.41639−0.95096i bei n = 200 … −0.41473−0.94841i bei 800).
- Extrapolation (k31, n = 200…800, Terme n^{−a}, n^{−1}, n^{−2a}, n^{−a−1}): **freier β-Fit β = 0.393586 + 0.025316i, Modell ½ + iα = 0.393586 + 0.025316i**
  (6 Stellen; adjungiert ebenso frei 0.712826 − 0.050632i gegen ½ − 2iα = 0.712828 − 0.050632i). Wiener-Abstand Modell − Grenzprofil 1.7e-7 (Norm 3.05);
  ‖S_n^P − P(λ_P·)‖_W = 1.37 … 1.46·n^{−0.894}.
- Skala **λ_P = 0.603424 + 1.165020i**, und **λ_P/λ_A = e^{iπa} auf 5.7e-7** (= Genauigkeit der λ). Mit u := λ_Aζ heisst das: der primale Kopf ist P(e^{iπa}u),
  der adjungierte Y(u). Wegen (e^{iπa}u)^{κ₁} = −u^{κ₁} sind die WKB-Exponenten entgegengesetzt (Sattelpunkt: Y ~ u^{κ₁(β_A−½)}e^{−aC u^{κ₁}},
  P(e^{iπa}u) ~ u^{κ₁(β_P−½)}e^{+aC u^{κ₁}}, C = (iα)^{−iακ₁}): **das lokale, exakte Kopf-Gleichgewicht** — wo der adjungierte Kopf wächst, fällt der
  primale mit genau derselben Rate. Die Amplitudenexponenten addieren sich zu κ₁(β_A + β_P − 1) = −iακ₁ = κ₁ − 1, also ist das Produkt
  P·Y·du ~ const·u^{κ₁−1}du = const·dv in der uniformisierenden Koordinate v = u^{κ₁}: eine flache Paarungsdichte („Wronski“ von e^{v}, e^{−v}) —
  die in H-A gesuchte erhaltene Paarung, aber asymptotisch im WKB-Bereich und nicht punktweise an der Ecke (Vermutung, Sattelpunktrechnung; §19).
- Sektortest (k32): die primale Defekt-Tabelle des echten Kopfes (k23, n = 600) wird reproduziert — A = 1 identisch, A = 4 auf ±0.2, A = 16:
  Kopf 195°: 6e-2, 210°: 3e-4, 225°: 6e-6, 240°: 7e-7, 255°: 6e-7, 270°: 5e-6, 285°: 2e-4 gegen Modell 4e-2, 2e-4, 6e-6, 1e-6, 1e-6, 7e-6, 3e-4.
**Stand:** Beide Köpfe sind an der Ecke explizit (numerisch, belastbar); die einzigen freien Grössen sind EINE komplexe Skala λ_A (Anschluss an die äussere
Lösung) und die Amplituden. Die Grenzaufgabe K → ∞ ist damit gelöst; was fehlt, ist die Satzform (§19).

## 19. Satzform der Eck-Grenzaufgabe: Existenz (explizit), Eindeutigkeit (Mellin-Liouville), Dualität der Köpfe (12.10. 02:40–03:00)

Bezeichnungen: μ = √L (Hauptzweig), log μ = L/2, μ^w := e^{wL/2}; α := log μ/(2π) = L/(4π); a := 1 + iα = 1/κ₁ (κ₁ = 2πi/(2πi − log μ));
β_A := ½ − 2iα, β_P := ½ + iα. Mellin-Konvention Φ(s) = ∫_0^∞ y(u)u^{s−1}du, y(u) = (1/2πi)∫Φ(s)u^{−s}ds; dann y(u/μ) ↔ μ^sΦ(s), y(μu) ↔ μ^{−s}Φ(s).
Lokale Kopfoperatoren (Kegelskala, K → ∞): adjungiert **E_A y := μy + μ⁻¹y(·/μ)** ↔ μ(1 + μ^{s−2})Φ; primal **E_P p := p(μ·) + μp** ↔ μ(1 + μ^{−s−1})Φ.
Rechte Nullstellen: s_m^A = 2 + (2m+1)·2πi/L, s_m^P = −1 + (2m+1)·2πi/L (m ≥ 0); die übrigen (m < 0) liegen links (Re ≤ −2.45 bzw. −5.45).

### 19.1 Lemma K10.1 (Gamma-Zerlegung der Symbole) — bewiesen
  **μ(1 + μ^{s−2}) = 2πμ^{s/2} / [Γ(β_A + iαs)·Γ(1 − β_A − iαs)],   μ(1 + μ^{−s−1}) = 2πμ^{(1−s)/2} / [Γ(β_P + iαs)·Γ(1 − β_P − iαs)].**
*Beweis.* Ergänzungssatz Γ(z)Γ(1 − z) = π/sin πz mit z = β_A + iαs = ½ + iα(s − 2): sin πz = cos(πiα(s − 2)) = cosh(πα(s−2)) = cosh((s − 2)log μ/2)
= (μ^{(s−2)/2} + μ^{−(s−2)/2})/2, also Γ(z)Γ(1−z) = 2πμ^{(s−2)/2}/(1 + μ^{s−2}); mit μ multiplizieren. Primal ebenso mit z = ½ + iα(s + 1). ∎
Γ(β_A + iαs) hat genau die Pole s_m^A, m ≥ 0 (rechts), Γ(1 − β_A − iαs) genau die s_m^A, m < 0 (links); primal entsprechend. Das ist die
Malyuzhinets-/Wiener–Hopf-Zerlegung K = K_L·K_R des Symbols im Mellin-Raum.

### 19.2 Satz K10 (Existenz, explizite Lösung) — bewiesen bis auf die Standard-Asymptotik (iii)
Für β ∈ {β_A, β_P} sei **F_β(u) := Σ_{k≥0} Γ(β − iαk)·(−u)^k / (k!·Γ(β))**.
(i) F_β ist ganz von der Ordnung 1/Re a = 1.1191 (Stirling: log|Γ(β − iαk)/k!| = −Re(a)·k log k + O(k), weil −iα = 1 − a). F_β(0) = 1.
(ii) **E_A F_{β_A}(u) = (2π/Γ(β_A)) Σ_k (−μ^{−1/2}u)^k / (k!·Γ(½ + iα(k+2)))**, **E_P F_{β_P}(u) = (2π√μ/Γ(β_P)) Σ_k (−√μ u)^k / (k!·Γ(½ + iα(k−1)))**,
    beide ganz. *Beweis:* Taylor-Koeffizient von E_A F bei u^k ist μ(1 + μ^{−k−2})c_k; Lemma K10.1 bei s = −k. Primal ebenso. ∎ (k29: Identität auf 1e-55.)
(iii) Mellin–Barnes: F_β(u) = (1/2πi)∫_𝒞 Γ(s)Γ(β + iαs)u^{−s}ds/Γ(β), 𝒞 trennt −ℕ₀ von den Polen von Γ(β + iαs) und verläuft asymptotisch in Richtung
    e^{i(π/2 − arg a)}·ℝ (dort fällt der Integrand nur exponentiell; längs vertikaler Geraden gäbe es Glieder ∓Re(α)|τ|log|τ|). Folgerungen (Standard-
    Asymptotik von Barnes-Integralen, Braaksma 1964, hier mit komplexen Steigungen 1 und iα — **nicht im Einzelnen nachgeprüft, numerisch bestätigt**):
    in einem Spiralsektor 𝒮_hom (Ränder arg u = const + tan(arg a)·log|u|, tan(arg a) = 0.0283) gilt F_{β_A}(u) ~ Σ_{m≥0} r_m u^{−s_m^A} und
    E_A F_{β_A}(u) = O(|u|^{−M}) für alle M; ausserhalb wächst F wie (s*)^{β−½}e^{−as*}, s* = (u(iα)^{−iα}e^{2πik})^{κ₁} (Zweig k je Sektor).
    k34: Restglied der 3-Term-Entwicklung 4e-23 … 5e-26 bei |u| = 80 in arg u ∈ (−60°, 60°); |E_AF|/|F| ≤ 3e-27 in (−15°, 45°); Sattelformel trifft log|F| in
    den Wachstumssektoren auf ±0.1 (k = −1 für arg u ∈ [120°, 180°], k = 0 für arg u ∈ [−150°, −120°]).

### 19.3 Satz K11 (Eindeutigkeit, Mellin-Liouville) — bewiesen
Sei 𝒮(θ₁, θ₂) := {u = e^{a(t+iθ)} : t ∈ ℝ, θ₁ < θ < θ₂} (Spiralsektor; in v = u^{κ₁} ein gerader Sektor). Die Dilatation u ↦ u/μ bildet den Spiralstrahl
zu θ auf den zu θ − ϑ_μ ab, **ϑ_μ := 2π(Re κ₁ − 1) = 42.55°** (denn −log μ = a·2πi(1 − κ₁); bewiesen, eine Zeile). Sei y ganz, y(0) = 1, und
(E1) es gibt θ₂ − θ₁ > ϑ_μ, so dass gleichmässig auf abgeschlossenen Teilsektoren von 𝒮(θ₁, θ₂) gilt: y(u) = O(|u|^{−δ}) (ein δ > 0) und E_Ay(u) = O(|u|^{−M}) ∀M;
(E2) (Randbedingung im Mellin-Raum) mit einem λ ∈ ℂ*, q ≥ 0: |Φ_y(s)| ≤ C(1 + |s|)^q·|Γ(s)Γ(β_A + iαs)λ^{−s}| für s ausserhalb der Kreise vom Radius ¼
    um die Pole von Γ(s)Γ(β_A + iαs).
Dann ist **Φ_y = Π(s)·Γ(s)Γ(β_A + iαs)λ^{−s}/Γ(β_A)** mit einem Polynom Π vom Grad ≤ q, Π(0) = 1, also y = Π(−u d/du)F_{β_A}(λu); für q < 1 ist **y = F_{β_A}(λu)**.
*Beweis.* (1) Φ_y(s) := ∫_{Σ_θ} y u^{s−1}du längs eines Spiralstrahls Σ_θ ⊂ 𝒮 ist analytisch in 0 < Re(as) < δ·Re a und von θ unabhängig (Cauchy; Bögen
verschwinden wegen (E1) bei ∞ und wegen Re(as) > 0 bei 0). (2) Links: y − Σ_{k<K}y_ku^k = O(u^K) bei 0 ⇒ Φ_y − Σ_{k<K}y_k/(s+k) analytisch in Re(as) > −K·Re a;
also meromorph links mit einfachen Polen höchstens in −ℕ₀. (3) Rechts: für θ ∈ (θ₁ + ϑ_μ, θ₂) liegen Σ_θ und Σ_θ/μ = Σ_{θ−ϑ_μ} (Substitution v = u/μ) in 𝒮, also
Φ_{E_Ay}(s) = μ(1 + μ^{s−2})Φ_y(s) im Streifen; Φ_{E_Ay} ist wegen E_Ay = O(|u|^{−M}) und Analytizität in 0 für alle Re(as) > 0 analytisch. Damit ist
Φ_y = Φ_{E_Ay}/(μ(1 + μ^{s−2})) rechts meromorph mit einfachen Polen höchstens in {s_m^A : Re(as_m^A) > 0} = {s_m^A : m ≥ 0} (Re(a s_0^A) = 5.73,
Re(a s_{−1}^A) = −2.16). (4) E(s) := Φ_y(s)λ^s/(Γ(s)Γ(β_A + iαs)) ist ganz (Γ hat keine Nullstellen; seine Pole decken −ℕ₀ ∪ {s_m^A}_{m≥0} ab). (5) (E2) ⇒
|E| ≤ C(1+|s|)^q ausserhalb der Kreise, nach dem Maximumprinzip überall (mit grösserem C); Liouville ⇒ E Polynom vom Grad ≤ q. Residuum in s = 0: E(0) = 1/Γ(β_A). ∎
**Primal (Satz K11′)** wörtlich mit E_P, β_P, s_m^P (Pol-Kontrolle rechts: Re(a s_0^P) = 3.05 > 0, Re(a s_{−1}^P) = −4.84 < 0) — bewiesen wie oben.
*Bemerkungen.* (a) Ohne (E2) ist die Lösung NICHT eindeutig: jede Überlagerung ∫F(tλu)dν(t) von Dilaten mit gemeinsamem homogenem Sektor erfüllt (E1)
(z. B. E(s) = (μ^{s−s₀} − 1)/(s − s₀): nur die Potenz u^{−s₀} bleibt, y ist ganz) — die lokale Gleichung ist dilatationsinvariant. (E2) ist die
Anschlussbedingung an die äussere Lösung (Exponent → λ, algebraischer Vorfaktor → q < 1). (b) **Lemma K12 (Übersetzung, Skizze):** (E2) mit q = 0 folgt,
wenn y auf dem abgeschlossenen Spiralsektor, auf dem F(λ·) polynomial beschränkt ist, ebenfalls polynomial beschränkt ist und in jedem Wachstumssektor
y = c_j·(Sattelform von F(λ·))(1 + o(1)) mit DEMSELBEN Exponenten und Vorfaktor gilt. Lücke: die Phragmén–Lindelöf-Abschätzungen von Φ_y in der ganzen
s-Ebene (links über die Taylor-Koeffizienten, oben/unten über die Sektorränder) — Standard, nicht ausgeführt. (c) Damit ist die Grenzaufgabe K → ∞
**lösbar (explizit) und, gegeben die Anschlussdaten, eindeutig**. Dass die Grenzprofile der echten Köpfe die Hypothesen erfüllen, ist numerisch
belegt (§18), nicht bewiesen — das ist genau die Anpassung innen/aussen.

### 19.4 Dualität der beiden Köpfe — die Gammafaktoren sind Reflexionspartner (Identität bewiesen; λ-Relation numerisch)
Wegen β_P + iα = 1 − β_A gilt **Φ_{F_P}(1 − s) = Γ(1 − s)·Γ(1 − β_A − iαs)/Γ(β_P)**: der primale Kopf trägt nach der Parseval-Spiegelung s ↦ 1 − s genau die
KOMPLEMENTÄREN Faktoren (Γ(1−s) zu Γ(s), den linken WH-Faktor Γ(1 − β_A − iαs) zu Γ(β_A + iαs)). Mit der gemessenen Skalenrelation
  **λ_P = e^{iπa}·λ_A = −μ^{−1/2}·λ_A**  (k31: auf 5.7e-7; e^{iπa} = e^{iπ}e^{−πα} = −μ^{−1/2})
wird das Parseval-Produkt der beiden Mellin-Transformierten ELEMENTAR (Lemma K10.1 und Ergänzungssatz):
  **λ_A^{−s}Φ_{F_A}(s)·λ_P^{s−1}Φ_{F_P}(1 − s) = 4π²i / [μ·λ_P·Γ(β_A)Γ(β_P)·(1 − e^{−2πis})·(1 + μ^{s−2})].**
Deutung: die Ecke ist ein skalares Wiener–Hopf-Problem, dessen zwei Faktoren auf die beiden biorthogonalen Partner verteilt sind (wie bei Toeplitz:
Orthogonalpolynom ↔ äusserer Faktor, Cauchy-Transformierte ↔ der andere); die λ-Relation ist die lokale Form des Ratengleichgewichts g_P + g_A = 0 bei L
(entgegengesetzte WKB-Exponenten, (e^{iπa}u)^{κ₁} = −u^{κ₁}). Offen: eine Herleitung der λ-Relation aus der Biorthogonalität (Vermutung: sie folgt aus
⟨ℓ_n, J̃c_n⟩ = 1 im Anpassungsbereich).

## 20. Kopplung innen/aussen: was an τ_n lokal ist — und ein grundsätzliches Hindernis aus der Spirale (12.10. 02:35–03:10)

### 20.1 Zerlegung von τ_n (numerisch; k33, k35, k38)
τ_n = ‖c_n‖·‖r̂_n‖/|⟨ℓ_n, J̃c_n⟩| ist skalenfrei. Wo die beiden Normen sitzen (kumulierte Norm-Anteile in Kegeleinheiten A = |x − L|/(ρn^{−Re a}), beide Ecken):
| n | primal (Kreis) A ≤ 4 / 8 / 32 | adjungiert (Ks Linse) A ≤ 8 / 16 / 32 / 128 | max_γ|Λ|/|Λ(L)| (Lage) | τ_n |
|---|---|---|---|---|
| 300 | 0.746 / 0.860 / 0.968 | 0.251 / 0.392 / 0.502 / 0.641 | 2.54 (A = 9.6) | 0.0447 |
| 600 | 0.820 / 0.930 / 0.996 | 0.329 / 0.578 / 0.826 / 0.917 | 3.59 (A = 14.4) | 0.0180 |
| 1000 | 0.855 / 0.956 / 0.999 | 0.253 / 0.497 / 0.840 / 0.998 | **5.4 (A = 20.3)** | 0.0153 (k12/V: 0.0149) |
- **Primale Norm ist lokal und explizit:** ‖c_n‖ ≈ |δ_n(L)|·n^{−Re a/2}·Q_P(ψ_n), Q_P = (I_P(ψ_n)/(π|λ_P|))^{1/2}, I_P = ∫_ℝ|F_P(te^{iψ_n})|²dt längs der Tangente,
  ψ_n = −90° + arg λ_P − arg L + Im(a)·ln n. Modell 1.78 / 1.665 (n = 400 / 800, k33) gegen gemessen ‖c_n‖/(|δ(L)|n^{−0.447}) = 1.796 / 1.657 (k12/k30).
- **Adjungierte w-Norm ist NICHT kegel-lokal:** sie sitzt im Anpassungsbereich A ≈ 10–40, wo der Kopf längs der Linse (Γin-Seite) auf das 2.5–3.6-fache des
  Eckwertes anwächst. Das Eckprofil F_A wächst längs der Γin-Tangente (k33: |F_A| längs t < 0 bei t = 32: 10^{1.03} / 10^{1.39} / 10^{2.28} / 10^{4.05} für
  n = 400 / 1000 / 10⁴ / 10⁶) — die Tangente liegt am Rand des Wachstumssektors und wandert mit n hinein.

### 20.2 Der Mechanismus: komplexes κ₁ ⇒ log-Spirale an L (Skizze; jeder Baustein numerisch geprüft)
(a) Das Eckprofil ist in der KOMPLEXEN Skala n^{−a} stationär (G, §18); relativ zu festen Richtungen dreht es sich um Im(a)·ln n = 1.45°·ln n. Im
physikalischen Raum heisst das (eine Zeile: die ln n-Terme heben sich auf, weil tan(arg a)·Re a = Im a): die äussere Rate des adjungierten Kopfes bei L ist
die Kegelfunktion Re(c₁(x − L)^{κ₁}) mit log-spiralförmigen Null-Linien arg(x − L) = θ₀ + 0.0284·ln(1/|x − L|) — wie in Ws Kegel (§19 dort), dessen
Null-Linie bei ψ = 89.0° / 89.1° / 90.0° (r = 0.03 / 0.06 / 0.12) liegt, also genau an der Γin-Tangente (ψ = 90°), mit Drehung nach aussen für kleineres r.
(b) Folge: für r < r_c (≈ 0.1) liegt die Kreistangente auf der Wachstumsseite, mit einer Rate ≈ c'·r^{κ₁}·ln(r_c/r) (c' ≈ |c₁|·Re κ₁·0.0284 ≈ 0.02). Eine
tangentiale Linse (Höhe ~ r^p, p > Re κ₁ = 1.118) dämpft nur ~ r^p und verliert gegen r^{κ₁}ln(1/r) nahe L ⇒ die Netto-Rate längs γ ist in einem Bereich
r ∈ (0, r*) positiv (Schätzung für Ks Linse p = 1.3: r* ≈ 0.02, Maximum ≈ 1e-4) ⇒ **τ_n wächst asymptotisch wie e^{εn} mit kleinem ε > 0.**
(c) Eine Keil-Linse (Höhe c·r) gewinnt nahe L, bis die Spirale um den Keilwinkel weitergedreht ist (Δln(1/r) = Keilwinkel/0.0284) — das verschiebt das
Versagen zu sehr grossem N, verhindert es aber nicht: die Spirale windet sich unendlich oft um L, jede feste Kurve durch L liegt für geeignete kleine r auf
der Wachstumsseite. Die zulässige Richtung von γ in L ist aber durch die Geometrie beschränkt (Schlitz von f_K; f(ℬ) ⊂ D: Keil ≤ arg μ = 38.3°).

### 20.3 Belege
| Prüfung | Ergebnis |
|---|---|
| Modell (Eckprofil F_A längs γ, k36/k37; Annahme: Anpassung gilt bis A ≪ n^{Re a}) | max|F_A| längs Ks Linse (p = 1.3): 2.0 / 2.7 / 3.6 / 765 / 1.1e9 (n = 300 / 600 / 10³ / 10⁴ / 10⁵); Vs Linse (p = 1.15, τ₁ = 5°): 1.5 (10³), 3.4 (10⁴), 382 (10⁵); Keil c = 0.5: ≈ 1 bis 10⁸; Keil 0.25: 6.4e5 bei 10⁸ |
| **echte Köpfe, dünne Linse p = 2.0 (k38)** | **τ_n = 0.0494 / 0.0339 / 0.0488 (n = 300 / 600 / 1000) — Umkehr, τ wächst**; max_γ|Λ|/|Λ(L)| = 3.83 / 7.58 / 17.5 (A = 12 / 21 / 32); adjungierte Masse zu 99 % in A ≤ 64 |
| echte Köpfe, Ks Linse (k35, k12) | max_γ|Λ|/|Λ(L)| = 2.54 / 3.59 (n = 300 / 600; Modell 2.0 / 2.7); lokale Exponenten von τ_n 1.55 → 0.10 bei N ≈ 1000 (Minimum naht) |
| Vs VE1 (IDEEN) | tangentiale Linse: ‖D_N‖ = 0.477 … 0.584 wachsend (N = 100…500); Keil-Linse: sättigt bei 0.215 |
Achtung, Werkzeugfehler (behoben): die Modellspalte von k35/k38-Läufen vor 03:00 benutzte α in Doppelgenauigkeit; die Reihe von F_A löscht bei |u| > 25
um bis zu 10^{18} aus, die Modellwerte dort waren falsch (echte Kopfwerte unberührt). k36/k37 rechnen mit mpmath-α (richtig).

### 20.4 Folgerungen (Vermutung mit Belegen; Entscheidung bei V/Leitung)
1. **(S_γ), (H1) und (H*_G) mit FESTER Linse sind für N → ∞ vermutlich falsch** — nicht wegen des globalen Schleifengewinns, sondern wegen der Spirale
   an L: die w-Köpfe wachsen ab einem N_max(γ) wie e^{ε_γN}. N_max ≈ 10³–10⁴ für Ks Linse (Vorhersage, prüfbar), ≈ 10³ für p = 2 (gemessen), astronomisch
   für Keil-Linsen maximaler Öffnung. Die Messungen bis N = 2000 widersprechen dem nicht: ‖Ĵ_N⁻¹‖ = 1.25 wird vom Bulk bestimmt, der Eckanteil (τ ~ 0.02)
   ist darin unsichtbar.
2. Für den Hauptsatz heisst das: Satz K1 mit fester Kurve γ trägt nur bis N_max(γ); für „alle N" braucht man eine eck-angepasste Aussage — Kandidaten:
   (i) N-abhängige Kurven γ_N, die bis zur Kegelskala der Abklingseite der Spirale folgen (begrenzt durch den zulässigen Keil, reicht bis N ~ 10^{10…60});
   (ii) die Ecke nicht abschätzen, sondern mit dem expliziten Profil F_A AUSRECHNEN: in K1 trägt die Eckzone ∫_{γ∩U(L)} W(t) r_t dt bei, und r_t ist dort
   bis auf Amplitude bekannt; (iii) datenspezifisch statt Operatornorm (f_K regt die wachsende Eckrichtung nur über seine Singularität in L an).
3. Sofort umsetzbar: Standardlinse auf Keil-Linse (Steigung ~0.5) umstellen; τ_N für Ks Linse bei N = 1200/1500/2000 messen (Vorhersage: Minimum, dann
   Anstieg) — beides Aufgabe für V.

## 21. Welche Linse macht die Drehung mit? — Null-Spirale, Keil, Log-Spirale, mitgeführte Linse, und was K1″ dann hergibt (Auftrag Leitung 08:25; 08:25–08:37)

### 21.1 Die Null-Spirale des adjungierten Kopfes bei L (Gesetz bewiesen aus der Kegelform; Lage numerisch, k39)
Im Anpassungsbereich ist die n-te-Wurzel-Rate des adjungierten Kopfes relativ zu seinem Eckwert die Kegelfunktion g(x) = Re(−a·s*)/n = Re(C·(x − L)^{κ₁})
(Sattel k = −1 von F_A; C = −a[Λ₀(iα)^{−iα}e^{−2πi}]^{κ₁}, Λ₀ = −λ_Aρ/L²; Kontrolle Sattel gegen Reihe: konstanter Abstand 0.09 für arg u ∈ [100°, 180°]).
Ihre Null-Linien sind die Kurven **Im(κ₁·log(x − L)) = const**, also Strahlen in der uniformisierenden Koordinate v = (x − L)^{κ₁}, in Polarkoordinaten
  **θ(r) = θ₀ + t_s·ln(r/r₀),  t_s = −Im κ₁/Re κ₁ = 0.02833 rad = 1.623° je e-Faktor (3.74° je Dekade)** — sie drehen sich für r → 0 im Uhrzeigersinn,
von der Γin-Tangente (θ_c = 166.6°) zur äusseren Normalen (76.6°), also INTO die Linse hinein. Winkel φ_null(r) der Null-Spirale gegen die Tangente (nach aussen):
| r | 1e-1 | 3e-2 | 1e-2 | 1e-3 | 1e-5 | 1e-10 | 1e-20 |
|---|---|---|---|---|---|---|---|
| φ_null | 2.95° | 4.8° | 6.7° | 10.4° | 17.9° | 36.6° | 74.0° |
Formel: **φ_null(r) ≈ 6.69° + 3.737°·log10(0.01/r)**; auf der Kegelskala r_N = ρN^{−Re a}: **φ_null(r_N) ≈ 1.451°·ln N − 1.3°**. Für r ≳ 0.05 ist die
linearisierte Lage nur auf einige Grad genau (Nichtlinearität w = ρ/x: bei r = 0.2 Drehung um 9°; höhere Kegelglieder); W fand die Null-Linie seines
Eckmoden-Kegels bei 1° ausserhalb der Tangente (r = 0.03). Massgeblich ist das Drehgesetz, und das ist exakt (es folgt allein aus κ₁).
Prüfung am echten Kopf (k40, n = 600, log10(|Λ|/|Λ(L)|) längs θ_c − φ, Kopf/Modell): A ≤ 4 Übereinstimmung ±0.02; ab A ≈ 8 ist der echte Kopf systematisch
„wachsender“ (φ = 0: A = 16: +0.78/+0.45, A = 32: +0.38/−0.99; φ = 10°: A = 24: −2.07/−2.80) — endliche-n-Korrektur in Richtung Wachstum; n = 1000 siehe 21.5.

### 21.2 Die Linsenformen
Kriterium: γ ist auf der Skala r „sicher“, wenn sein Winkel φ_γ(r) gegen die Tangente grösser als φ_null(r) ist; unterhalb der Kegelskala r_N wächst nichts
(dort ist der Kopf das beschränkte Profil). Für N braucht man also φ_γ(r) > φ_null(r) auf [A₀r_N, r₁].
- **Tangentiale Linse** (Höhe ~ r^p, p > 1): φ_γ(r) → 0 ⇒ unsicher für alle r < r(φ_null = φ_γ); Ks Linse (p = 1.3, τ₁ = 12°): φ_γ(0.01) ≈ 1.7° < 6.7°.
  Netto-Rate längs γ bei r ~ 0.01 von der Ordnung 1e-4 ⇒ sichtbares Wachstum von τ_N ab N ~ 10³–10⁴; p = 2 (φ_γ(0.01) ≈ 0.2°): gemessen ab N ≈ 600.
- **Keil** (konstanter Winkel φ): sicher für r > r_×(φ), unsicher darunter; die Kegelskala erreicht r_× bei N_×(φ) := (ρ/r_×)^{1/Re a} (k39):
  | φ | 5° | 10° | 15° | 20° | 26° | 30° | 38.3° | 45° | 60° | 90° |
  |---|---|---|---|---|---|---|---|---|---|---|
  | r_× | 1e-2 | 4.6e-4 | 2.2e-5 | 1.0e-6 | 4.6e-8 | 2.2e-9 | 2.2e-11 | 2.2e-13 | 2.2e-17 | 2.2e-25 |
  | N_× | 2.5e2 | 7.7e3 | 2.4e5 | 7.4e6 | 2.3e8 | 7.1e9 | 1.2e12 | 2.1e14 | 6.4e18 | 5.7e27 |
  Jedes Grad Keilwinkel verschiebt N_× um den Faktor ≈ 1.95 (= e^{1°/1.451°}). Jenseits von N_× wächst τ_N wie exp(c·(N/N_×)^{≈1}) (die Netto-Rate bei
  r_×e^{−1/κ₁} ist ~ (r_×/r_N)^{κ₁}/N).
- **Log-Spirale** (Rand der Linse nahe L auf einer Kurve Im(κ₁log(x − L)) = const mit Abstand m zur Null-Spirale): macht die Drehung exakt mit, sicher auf
  allen Skalen — aber φ_γ(r) = φ_null(r) + m → ∞ für r → 0: die Kurve verlässt den zulässigen Bereich (21.3) bei r_min und muss dort abgeschnitten werden.
- **Mitgeführte (N-abhängige) Keil-Linse**: φ_N := φ_null(r_N) + m ≈ 1.451°·ln N + m′ — auf [r_N, ∞) sicher, weil φ_null(r) < φ_null(r_N) für r > r_N.
  Werte (m′ = 4°): N = 400: 12.7°; 2000: 15.0°; 10⁴: 17.4°; 10⁶: 24°; 10⁹: 34°; 10¹²: 44°; 10¹⁸: 64°. Gleichwertig mit der bei r_N abgeschnittenen
  Log-Spirale, aber mit geradem Eckwinkel (einfachere Geometrie, siehe 21.3); die Spirale hat auf den Skalen u ~ N^{−α} (wo der Kneser-Sprung für (LP) zählt)
  einen um ≈ 0.65°·ln N kleineren Winkel, also eine etwas bessere (LP)-Rate.

### 21.3 Trägt K1″ mit einer solchen Linse? (bewiesen, soweit angegeben; die Bedingungen sind elementar nachprüfbar)
Sei ℬ eine Linse über Γin, die bei L ein Keil mit Winkel φ ist (glatt an den Rest der Linse angeschlossen), konjugiert bei L̄.
1. **Schlitze** (Voraussetzung (A)): S_L zeigt in Richtung der äusseren Normalen (76.6°). ℬ liegt nahe L in den Richtungen [166.6° − φ, 166.6°], f(ℬ) in
   [204.9° − φ, 204.9°] (Drehung um arg μ = 38.3°), die Singularitätsstelle von m = f′∘f in Richtung 38.3°. Disjunkt genau für **φ < 90°** (f(ℬ): φ < 128.3°,
   m: φ < 128.3°). f(ℬ) ⊂ D nahe L gilt nur für φ ≤ arg μ = 38.3°, wird aber nicht gebraucht: für φ > 38.3° ändert f(ℬ) den Rand von K_γ bei L nicht (ℬ deckt
   die Richtungen schon ab).
2. **K1″(a)** (‖Π_Nk_η‖ ≤ ‖Ĵ_N⁻¹‖(1 + M_γ)/dist(η, K_γ)): braucht nur W = J_∞k_η analytisch auf Int γ und stetig bis γ, also (A) — **gilt**. M_γ = max|m| auf
   D̄ ∪ ℬ̄ ist gleichmässig in φ ≤ φ_max < 90° beschränkt.
3. **K1″(b)** (‖Ĵ_N⁻¹‖ ≤ c_k‖G_N‖_{L²(γ)}): an der einspringenden Ecke (Innenwinkel π + φ) ist k′ ~ (w − w_L)^{φ/π} beschränkt (→ 0), c_k < ∞ — **gilt**.
4. **K1″(c)** und Lemma K1′ (Q_k = sup|w − w′|/|k(w) − k(w′)|) — **gelten nicht** (k ist an der Ecke nicht bi-Lipschitz); sie werden für die Kette
   (S_γ) ⇒ (H₀) ⇒ (LP) aber nicht gebraucht, nur für die Rückrichtung G_N ≤ C√N‖Ĵ⁻¹‖ (Satz K1 in der G-Form). Satz V2 (Vs Konstanten-Kette) benutzt nur (a).
5. **Bernstein–Walsh / (LP)**: K_γ = (D̄ ∪ ℬ̄) ∪ f(D̄ ∪ ℬ̄) hat bei L den Aussenwinkel **π·α_γ = 141.69° − φ** (statt 141.69° für K′); g_γ ~ r^{1/α_γ} ist dort
   Lipschitz (Aussenwinkel < π), die 1/N-Umgebungsabschätzung bleibt. Satz P-I gibt mit beschränktem ‖Ĵ_N⁻¹‖:
   **sup_{D̄}|Π_Nf_K − f_K| = O(N^{3/2 − 5.447·α_γ})**. Mit Familie (X4a, ein Grad mehr) braucht man 5.447α_γ > 5/2 ⇒ **φ < 59.1°**; ohne Familie 5.447α_γ > 3/2
   ⇒ φ < 92° (dann bindet der Schlitz, φ < 90°).
6. **Mitgeführte Linse**: alle Konstanten in 1.–5. sind gleichmässig für φ_N ≤ φ_max < 90°; der einzige N-abhängige Teil ist der Eckexponent α_γ(N). Mit
   φ_N ≈ 1.451°·ln N + m′ wird die (LP)-Schranke **N^{3/2 − 4.288 + 0.0439·ln N + c} = N^{−2.79}·e^{0.044(ln N)²+...}**: fällt bis ln N ≈ 63, mit Familie
   (Bedarf < N^{−1}) bis **φ_N < 59° ⇔ N ≲ 10¹⁶ (Rand m′ = 4°) … 5·10¹⁷ (ohne Rand)**.

### 21.4 Ergebnis zur Frage der Leitung
- **Die Drehung macht genau die Log-Spirale Im(κ₁ log(x − L)) = const mit** (Steigung 1.623°/e-Faktor); praktisch gleichwertig und geometrisch sauberer ist
  der **mit N mitgeführte Keil φ_N ≈ 1.45°·ln N + 4°**. Ein fester Keil reicht bis N_×(φ) (Tabelle 21.2); die tangentiale Standardlinse nur bis ~10³–10⁴.
- **K1″(a)/(b) tragen** mit Keil- und mitgeführten Linsen (keine Glattheit von k nötig), K1′/K1″(c) nicht (unnötig). Der Preis ist der Eckexponent von g_γ:
  die Familienkette verträgt φ < 59°, also **N ≲ 10¹⁶–10¹⁷**, ohne Familie N ≲ 10²⁶–10²⁷ (Schlitz, φ < 90°).
- **Für „alle N“ reicht keine Linse**: der zulässige Winkelbereich an L ist endlich (Aussenwinkel 141.7° von K′, Schlitz), die Null-Spirale dreht unbeschränkt.
  Jenseits von ~10¹⁶ muss die Ecke in K1 anders behandelt werden — nicht über eine Kurve um K_γ abgeschätzt, sondern mit dem expliziten Eckprofil F_A (und der
  Kegelfunktion C(x − L)^{κ₁}) ausgerechnet: Vorschlag „K1 mit Eck-Einsatz“ — Cauchy-Darstellung auf γ ausserhalb U(L), in U(L) die Kopfentwicklung
  r_t = Σ_n c_nΛ_n(t) mit Λ_n = Λ_n(L)·F_A(u_n(t))(1 + Fehler) direkt gegen den Kneser-Sprung integriert. Das ist der nächste Baustein der Kopplung.
- Sofort für V: Standardlinse → mitgeführter Keil (für die Zertifikate N ≤ 2000 genügt ein fester Keil φ ≈ 15°, N_× ≈ 2·10⁵) und τ_N damit neu messen.

### 21.5 Prüfung am echten Kopf und eine zweite, linsenunabhängige Grenze (k40, k41; 08:30–08:37)
**k40 (adjungierter Kopf längs Strahlen θ_c − φ, n = 600 und 1000, log10(|Λ|/|Λ(L)|) Kopf/Modell):** die Abweichung Kopf − Modell bei festem A schrumpft von
600 auf 1000 um den Faktor ≈ 0.64 = (600/1000)^{Re a} (φ = 0, A = 24: 0.75 → 0.48; A = 32: 1.37 → 0.87; φ = 10°, A = 24: 0.73 → 0.47) — **die echten Köpfe
konvergieren auch im Anpassungsbereich wie n^{−a} gegen F_A**, und die Drift ist in beiden sichtbar (φ = 0, A = 24: Modell −0.07 → +0.61, Kopf +0.68 → +1.09;
φ = 5°, A = 16: Abklingen −0.44 → −0.12). Damit ist das Modell auch jenseits der Kegelskala der richtige Grenzwert.
**k41 (primaler Kopf längs des KREISES — die Unbekannten-Norm ℓ²_ρ hängt am Kreis, keine Linse wirkt dort):** log10(|δ|/|δ(L)|) im Modell F_P(λ_Pζ):
| n | Γin-Hälfte A = 16 / 32 / 100 | Γout-Hälfte A = 16 / 32 / 64 / 100 | arg u_P der Γout-Tangente |
|---|---|---|---|
| 1e3 | −0.9 / −2.0 / −7.0 | −2.4 / −3.9 / −4.9 / −5.6 | 86.0° |
| 1e6 | −2.6 / −5.0 / −6.9 | −0.7 / −2.1 / −4.8 / −5.5 | 96.0° |
| 1e8 | −3.8 / −5.2 / −6.9 | **+0.4 / +0.3 / −0.3 / −1.4** | 102.7° |
| 1e9 | −4.2 / −5.2 / −6.9 | **+0.9 / +1.5 / +2.3 / +2.9** | 106.1° |
| 1e10 | −4.2 / −5.1 / −6.8 | +1.5 / +2.7 / +5.0 / +7.2 | 109.4° |
| 1e12 | −4.0 / −5.1 / −6.8 | +2.6 / +5.1 / +10.2 / +15.8 | 116.1° |
Die Γout-Hälfte der primalen Tangente wandert (wie die Γin-Hälfte der adjungierten, nur ≈ 19° später) in den Wachstumssektor von F_P: **ab n ≈ 10⁸ ist die
primale Norm nicht mehr kegel-lokal und wächst super-polynomial — unabhängig von der Linse**, weil die Unbekannten auf dem Kreis gemessen werden. Das Gegenstück
der Dualität: die Ecke hat in jeder Richtung ein Wachstum von F_A oder F_P; Richtungen, in denen BEIDE nicht wachsen, sind zwei Fenster von ≈ 21° und ≈ 59°
(in u_A: (79°, 100°) und (−120°, −61°)), die mit 1.45°/e-Faktor durch die feste Geometrie wandern.
**Folgerung:** (S_γ) in der jetzigen Form (Unbekannte in ℓ²_ρ = Kreis, Gleichungsseite auf γ) ist für N ≳ 10⁸–10⁹ mit KEINER Linse zu halten (Modellaussage,
gestützt durch die n^{−a}-Konvergenz von k40). Der Zeiger für „alle N“ ist damit nicht eine bessere Linse, sondern eine Eckbehandlung, die beide Köpfe explizit
mitführt (Eck-Einsatz in K1, 21.4), oder eine Norm auch auf der Unbekannten-Seite, die der Drehung folgt — Letzteres kollidiert mit sup_{D̄} (die Punkte des
Kreises nahe L gehören zu D̄).

### 21.6 Antwort auf die Nachfrage der Leitung (08:40): Keil-Linse und Bessel-Formulierung (Vs Daten §3.13: Ks Linse τ_N-Minimum bei 1000, Anstieg auf 0.0204 bei 2000; Keil c = 0.5 fällt auf 4e-4)
**(a) Gelten K1″ und meine Sätze mit Keil?** Ja, soweit gebraucht:
- K10, K11/K11′, Lemma K10.1 sind lokal (Kegelskala) und linsenunabhängig. Satz K5 (LU, gleiches U), Satz K7 und Lemma K8 gelten für jede Linse (C_k bleibt
  untere Dreiecksmatrix, k(0) = 0). K1″(a) und (b) gelten mit Ecke (21.3); Satz K1 in der G-Form und K1″(c)/K1′ nicht (Q_k = ∞) — unnötig für (S_γ) ⇒ (H₀).
- Zahlen für Vs Keile (φ = arctan c): c = 0.25: φ = 14.0°, (LP)-Exponent 5.447·α_γ = 3.86 (statt 4.29), **N_× ≈ 1.2·10⁵**; c = 0.5: φ = 26.6°, Exponent 3.48,
  **N_× ≈ 4·10⁸**; c = 0.79: φ = 38.3°, Exponent 3.13, N_× ≈ 1.2·10¹². Für V2 (N₀ ≤ 400) ist c = 0.25 unbedenklich; für den asymptotischen Teil ist ein
  steiler Keil (c → 0.79) oder der mitgeführte Keil nötig. **„Feste Ecke trägt für alle N“ stimmt nicht**: jenseits N_× (adjungierte Seite) und unabhängig
  davon ab N ≈ 10⁸ (primale Kreisnorm, 21.5) wachsen die Köpfe wieder — Vs Messung bis 2000 sieht das nicht.
**(b) Macht das Grenzprofil die Bessel-Eigenschaft beweisbar?**
1. Exakt (bewiesen, Cauchy–Schwarz): Ĵ_N⁻¹ = Σ_{k≤N} τ_k ĉ_k r̃_k*, also ‖Ĵ_N⁻¹‖ ≤ (‖S_c‖·‖S_r‖)^{1/2} mit S_c = Σ τ_kĉ_kĉ_k*, S_r = Σ τ_kr̃_kr̃_k* —
   Vs Bessel-Bedingung ist sup_N ‖S_c‖, ‖S_r‖ < ∞.
2. **Adjungierte Seite:** fast orthogonal (feste nicht-stationäre WKB-Phase, Lemma K8; mit dem expliziten Profil ist die Phase bekannt) ⇒ ‖S_r‖ ≲ τ̄·log N.
3. **Primale Seite: KEINE Fast-Orthogonalität.** Die primalen Köpfe sind kegel-lokal mit Profil F_P auf der Skala k^{−a} und rotierender Phase χ_k; innerhalb eines
   Skalenbandes k ∈ [n, 2n] spannen sie im Wesentlichen den 2-dimensionalen Raum {Re p_n, Im p_n} auf (gemessen |cos(c_n, c_{n+2})| = 0.75–0.94). Daher
   (Lemma, Beweis durch Testen mit Re p_n, Im p_n; Voraussetzung die Profilkonvergenz): **‖S_c‖ ≳ Σ_{k∈[n,2n]} τ_k — gleichmässige Bessel-Schranke ⇔ τ_k ≲ C/k im Bandmittel.**
4. **Was das Grenzprofil für τ_k mit Keil liefert** (Skizze, Anpassung angenommen): primale Norm kegel-lokal ‖c_k‖ ≈ |δ_k(L)||d_k|^{1/2}Q_P; adjungierte w-Norm
   eckdominiert, wobei die einspringende Ecke (Innenwinkel π + φ, k′ ~ (w − w_L)^{φ/π}) die Eckzone im w-Mass streckt: ‖r̂_k‖ ≈ |Λ_k(L)||d_k|^{1/(2(1+φ/π))}N_Y(φ).
   Mit Gs B_k = √k‖c_k‖|Λ_k(L)|: **τ_k ≈ B_k·N_Y(φ)·k^{−1/2}·|d_k|^{1/(2(1+φ/π))} ~ B_k·k^{−0.89}** (φ = 26.6°; tangential ohne Wachstum wäre es k^{−0.95}).
   Grösse bei k = 2000: 0.6·1.2e-3·N_Y ≈ 7e-4·N_Y — **genau Vs „Sättigung“ 4e-4**: die Sättigung ist der Übergang in den eckdominierten Bereich, danach sollte τ_k
   wie k^{−0.89} weiter fallen (Vorhersage für N > 2000: τ(3000) ≈ 2.8e-4, τ(4000) ≈ 2.2e-4).
5. Folgerung: Bandsummen Σ_{[n,2n]}τ_k ~ B·n^{0.11} ⇒ **gleichmässige Bessel-Schranke asymptotisch knapp verfehlt** (sie bräuchte B_k → 0 wie k^{−0.11}), aber
   ‖Ĵ_N⁻¹‖ ≤ (‖S_c‖‖S_r‖)^{1/2} ≲ (N^{0.11}·log N)^{1/2}·sup B_k — **polynomial mit Exponent ≈ 0.06, für K1/K1″ (κ < 1.79) mehr als genug**, solange N < min(N_×, 10⁸).
   Der beweisbare Weg zu (H1) in Keil-Geometrie ist damit: (i) B_k beschränkt (Gs globale Eckamplitude — offen), (ii) Anpassungsfehler O(k^{−a}) für beide Köpfe
   im Bereich 1 ≪ A ≪ k^{Re a} (k40 zeigt die n^{−a}-Rate numerisch), (iii) ferne Teile der adjungierten w-Norm durch die Linse gedämpft (Vs Schleifenbedingung).
   Gleichmässiges (S_γ) für alle N ist mit keiner festen Norm zu haben (21.4–21.5).

### 21.7 Unabhängige Prüfung der Profilformel bei anderem μ (k42; kubisch, |μ| = 2, arg μ = arg √L ⇒ Im a = 0.110, Drehung 6.3°/e-Faktor)
Primaler Kopf, q_k bei n = 400 gegen das Modell mit α = log μ/(2π) und β_P = ½ + iα (KEIN freier Parameter):
0.9862−0.1120i / 1.0153−0.0785i / 1.0215−0.0539i / 1.0207−0.0386i gegen 0.9849−0.1140i / 1.0161−0.0812i / 1.0233−0.0559i / 1.0227−0.0401i (n = 200: 0.9867−0.1091i …,
Konvergenz aufs Modell). Profil stationär in der KOMPLEXEN Skala n^{−a}: Wiener-Abstand 100→200: 0.013 (komplex) gegen 0.223 (reell). Die Formel gilt also allgemein
mit α = log μ/2π, und die Drehung um Im(a)·ln n ist bei schnellerer Spirale direkt sichtbar. (Der adjungierte Kopf bei |μ| = 2 ist bei w_L von einer anderen
Struktur überdeckt — c₁ ≈ 0.005, q_k unregelmässig; nicht weiter verfolgt.)

## 22. L3 „K1 mit Eck-Einsatz“ — Reduktion, ein vierter Restpunkt, und eine Warnung für sehr grosse N (12.10. 08:47–08:55)

### 22.1 Lemma K13 (Kopfentwicklung mit kopfabhängigen Konturen) — bewiesen, elementar
Für W analytisch in einer einfach zusammenhängenden Umgebung von 0 (Keim bei 0 fest) und n ≤ N gilt ⟨ℓ_n, T_NW̃⟩ = Σ_{k≤n}ℓ_{n,k}ρ^kW_k (Λ_n ist ein Polynom in
1/t vom Grad n), also
  **J_N⁻¹T_NW = Σ_{n≤N} δ_n · I_n[W],   I_n[W] = (1/2πi)∮_{Γ_n} Λ_n(t)W(t)dt/t,**
mit JEDER Kontur Γ_n um 0, die in einem einfach zusammenhängenden Gebiet liegt, auf dem ein Zweig von W mit dem Keim bei 0 analytisch ist — **Γ_n und der
Schnitt dürfen von n abhängen** (I_n hängt nur von den ersten n+1 Taylor-Koeffizienten ab). Bei f_K darf Γ_n nahe L spiralförmig verlaufen; W wird dabei
auf weitere Blätter fortgesetzt, die singulären Glieder ändern sich je Umlauf um e^{±4π²Re(1/L)} = e^{±6.65}, bei der Drehung 1.45°·ln n ist das ein Faktor
≤ n^{0.027} (polynomial). [Beweis: Residuensatz auf kleinem Kreis, dann Deformation im Gebiet.]

### 22.2 Datenseite: der Eck-Einsatz trägt für ALLE n (Skizze mit benannten Annahmen)
Wähle Γ_n = (feste Kurve γ fern von L, mit Linse) ∪ (Eckstück, das auf allen Skalen r ∈ [r_n, ε] im Nicht-Wachstumsfenster von Λ_n liegt — die „mitgeführte“
Richtung, aber KOPFWEISE, also ohne gemeinsamen Kompromiss über n). Dann |I_n[W]| ≤ C·n^{p}·|Λ_n(L)|·sup|W|, sofern (ii) die Anpassung innen/aussen mit
Fehler O(n^{−a}) gilt (Λ_n folgt F_A im Fenster) und (iii) Λ_n auf dem fernen Teil von γ durch die Linse gedämpft ist (Schleifenbedingung). Das N_×-Problem
der festen Linse verschwindet: jeder Kopf hat seine eigene Ecke. Der Eckbeitrag ist sogar explizit (Mellin–Parseval): für Daten mit Singularität (x − L)^s
ist I_n^{Ecke} = Λ_n(L)·d_n^{s+1}·λ_A^{−s−1}·(Sprungfaktor)·Φ_{F_A}-Integral; für den Kneser-Exponenten ist er resonant (s₀ + 1 = s_0^A, d. h. die
Kneser-Singularität ist homogene Lösung des primalen Eckoperators, μ^{s₀−1} = −1), aber nicht null (Φ_H(s_0^A) = 2πμ^{s_0^A/2}Γ(s_0^A)/Γ(β_A) ≠ 0).

### 22.3 Unbekanntenseite: ein vierter Restpunkt, und er ist für grosse n verletzt
Aus K13: |J_N⁻¹T_NW(z)| ≤ Σ_n |δ_n(z)/δ_n(L)|·|δ_n(L)Λ_n(L)|·|I_n[W]|/|Λ_n(L)|. Für (H₀)/(LP) braucht man sup über z ∈ D̄ — also auch
**(iv) primale Lokalisierung: sup_{D̄}|δ_n| ≤ C n^{p}|δ_n(L)|.** z ∈ D̄ lässt sich nicht verschieben, keine Linse wirkt hier. Nach k41 verletzt das Grenzprofil
F_P (iv) längs der Γout-Hälfte des Kreises ab n ≈ 10⁸ (Wachstum log10|δ/δ(L)| = +2.9 bei A = 100 für n = 10⁹, +15.8 für 10¹²).
**Prüfung des Mechanismus bei schnellerer Spirale (k43, kubisch |μ| = 2, Drehung 6.3°/e-Faktor; λ_P aus n = 400):** echter primaler Kopf längs der Γout-Hälfte
gegen Modell, log10|δ/δ(L)| bei A = 2/4/8/16/32/64:
| n | echter Kopf | Modell F_P |
|---|---|---|
| 200 | +0.24 +0.36 +0.22 −0.38 −1.06 −1.85 | +0.25 +0.38 +0.23 −0.35 −0.87 −1.39 |
| 300 | +0.27 +0.42 +0.32 −0.24 −0.92 −1.56 | +0.29 +0.44 +0.34 −0.20 −0.78 −1.26 |
| 400 | +0.31 +0.49 +0.44 −0.05 −0.82 −1.39 | +0.31 +0.48 +0.43 −0.07 −0.73 −1.18 |
Übereinstimmung bis A = 16 auf ±0.04 und **dieselbe Drift** (A = 16: −0.38 → −0.24 → −0.05 im Kopf, −0.35 → −0.20 → −0.07 im Modell); das Modell sagt für
|μ| = 2 den Übergang ins Wachstum bei n ~ 10³–10⁴ voraus (n = 3000: +1.29 bei A = 16; 10⁴: +3.9 bei A = 64). (Der |μ| = 2-Kopf bei n = 600 ist numerisch
unbrauchbar — arb-Lösung mit 2300 Bit kollabiert —, verworfen; k44 bei |μ| = 2 zeigt globale x-Instabilität e^{≈0.18N}, uninformativ; bei √L konvergiert
die Abschnittslösung zu singulären Daten nahe L bis N = 400 ordentlich.)

### 22.4 Satz K14 (bedingte Reduktion; Beweis = Lemma K13 + Dreiecksungleichung + Satz P-I) 
Gelten für alle n ≥ n₀: (i) |δ_n(L)Λ_n(L)|/|⟨ℓ_n, J̃c_n⟩| ≤ C n^{β₁}; (ii)+(iii) |I_n[J_∞k_η]| ≤ C n^{β₂}|Λ_n(L)|/dist(η, K_γ) mit kopfweisen Konturen; (iv) sup_{D̄}|δ_n|
≤ C n^{β₃}|δ_n(L)|, so ist sup_{D̄}|Π_Nk_η| ≤ C N^{1+β₁+β₂+β₃}/dist(η, K_γ) — (H₀) mit polynomialer Konstante, für alle N, ohne N_×. Mit 1 + Σβ < 1.79 + …
(Familienbudget) folgt der Hauptsatz wie bisher. **Damit ist die Antwort auf die erste Frage der Leitung: (H1) bzw. (H₀) reduziert sich in der Eck-Einsatz-
Fassung auf vier, nicht drei Restpunkte; (i)–(iii) sind einzeln angreifbar, (iv) ist nach dem Grenzprofil für n ≳ 10⁸ FALSCH.**

### 22.5 Die Restpunkte einzeln
- **(i) B_n — legt K11 die Normierung fest? Nein.** Die Grenzaufgabe ist homogen; K11 bestimmt Gestalt und Skala λ, nicht die Amplituden. Was die Struktur aber
  festlegt, ist der n-Exponent: Anpassung der Sattelvorfaktoren (Y ~ u^{κ₁(β_A−½)}, P ~ u^{κ₁(β_P−½)}) an n-unabhängige äussere WKB-Amplituden gibt
  |δ_n(L)Λ_n(L)| ∝ |d_n|^{Re(κ₁−1)}·C_aussen, also **B_n ∝ n^{−(1−Re a)/2}·C_aussen = n^{−0.053}·C_aussen** (Vermutung). Gs Daten (n = 150…1000) unterscheiden das
  nicht von einem konstanten Grenzwert (beide Dreiparameter-Fits Rest 3e-3; mit n^{−0.053}: C = 0.86, ohne: B_∞ = 0.55–0.56). (i) reduziert sich damit auf die
  Endlichkeit der äusseren Paarungsnormierung (WKB-Normierungsintegral über den Paarungsbereich, in dem g_P + g_A = 0; nahe L ist die Dichte flach ~ dv,
  integrierbar) — eine Aussage über Ws äussere Lösung, angreifbar.
- **(ii) Anpassungsfehler** — lokal, mit expliziter innerer Lösung; numerisch O(n^{−a}) (k40); Standard-Programm der angepassten Asymptotik, braucht die äussere
  WKB-Lösung bis auf O(1) Amplitude.
- **(iii) Stegdämpfung** — Vs Schleifenbedingung s ≥ Λ auf Γin, numerisch punktweise scharf; braucht W12′-artige Gestaltaussage.
- **(iv) primale Lokalisierung** — für n ≲ 10⁷–10⁸ wahr (Fenster der Dualität), darüber nach dem Grenzprofil falsch.

### 22.6 Warnung: die Konvergenz nahe L selbst ist für sehr grosse N in Frage (Vermutung; Kette mit benannten Gliedern — KEIN Gegenbeweis; Entscheidung durch Test T1; Glied (c) präzisiert in §23.1: es ist genau „B_n polynomial beschränkt“)
Die Kopfentwicklung ist exakt, die Summanden δ_n(z)I_n[W] haben für n im primalen Wachstumsbereich bei z auf ∂D nahe L (Γout-Seite, |z − L| ~ 10⁻⁶)
die Grösse ~ e^{c(r_*/r_n)^{κ₁}}·r_n^{s+1}·B_n n^{−1/2} mit polynomial kleinem Eckfaktor r_n^{s+1} (für f_K s + 1 ≈ 6.45) und super-polynomial wachsendem
Primalfaktor; eine Auslöschung über n ist bei monoton wachsenden Beträgen und fester schneller Phase nicht möglich. Glieder: (a) Anpassung der echten Köpfe an
F_P bis A ~ 100 bei n ~ 10¹⁰ (gestützt durch k40/k43, nicht bewiesen), (b) I_n^{Ecke} ≠ 0 (Mellin-Rechnung 22.2), (c) keine Kompensation durch B_n.
**Folgerung, falls die Kette hält:** Π_Nf_K (und damit wohl p_N) konvergiert für N ≳ 10¹⁰–10¹¹ in Punkten von ∂D im Abstand ~10⁻⁶ von L nicht mehr gleichmässig —
die Vermutung in der Fassung „gleichmässig auf D̄“ wäre dann falsch, richtig bliebe sie z. B. auf D̄ ohne schrumpfende Umgebungen von L, L̄ bzw. für N unter
dieser Schwelle. Das ist der wichtigste offene Prüfpunkt nach diesem Block. Ein Test bei erreichbarem N braucht ein Modell mit schneller Spirale UND
gezähmter globaler Schleife (|μ| ≈ 1.6–2 mit passend dicker Keil-Linse) — Vorschlag für V.

## 23. Was bleibt wahr und beweisbar, falls „gleichmässig auf D̄“ nahe L kippt? (Auftrag Leitung 08:58; 09:00–09:10)
Kennzeichnung vorweg: §22.6 ist eine **Kette mit benannten Gliedern (a)–(c), kein Gegenbeweis**; entschieden wird sie durch Test T1 (V). Die folgenden Aussagen
sagen, was in jedem Fall gilt bzw. was man braucht, und was vom Ausgang von T1 abhängt.

### 23.1 Präzisierung der Warnung: Glied (c) ist genau „B_n polynomial beschränkt“
Mit Lemma K13 und Cauchy–Schwarz: sup_{D̄}|δ_n| ≤ √(n+1)·‖c_n‖_{ℓ²_ρ}, also
  **|δ_n(z)·I_n[W]| ≤ √(n+1)‖c_n‖·|Λ_n(L)|·|I_n[W]|/|Λ_n(L)| = B_n·(|I_n[W]|/|Λ_n(L)|)   für alle z ∈ D̄**   (bewiesen, eine Zeile).
Die primale Seite geht also NUR über Gs Eckamplitude B_n = √n‖c_n‖|Λ_n(L)| (mit ⟨ℓ_n, J̃c_n⟩ = 1) ein — Restpunkt (iv) aus §22 steckt vollständig in (i).
Satz K14 gilt daher in der schärferen Form: **(H₀) für alle N ⇐ (i′) B_n ≤ Cn^β, (ii) adjungierte Anpassung im Eck-Einsatz, (iii) Stegdämpfung.** Die Warnung
ist äquivalent zu: B_n wächst super-polynomial, sobald die primale Kreisnorm ‖c_n‖ (n ≳ 10⁸) nicht-lokal wird — ODER die globale Normierung gleicht das aus
(dann |δ_n(L)Λ_n(L)| kleiner). Welcher Fall eintritt, ist offen. Versuch, das bei |μ| = 2 zu messen (k45): unbrauchbar — dort sind die J̃_n für manche n fast
singulär (B_n springt 598 / 7778 / 200 / 1.2e12 / … für n = 100…350), das Modell ist für Normfragen pathologisch. Für T1 braucht V ein Modell mit gut
konditionierten J̃_n (|μ| etwa 1.3–1.6, B_n als glatte Folge prüfen); √L zum Vergleich: B_n = 1.159 / 0.863 / 0.748 / 0.697 (n = 100…400, = G).

### 23.2 Reduktionssatz in lokaler Fassung (bewiesen — Prüfung des vorhandenen Beweises von kn.satz.polynom_halbexp_reduktion)
Im Beweis wird „|p_N| ≤ B auf D“ nur zweimal benutzt: (1) Montel — dafür genügt LOKALE Beschränktheit auf D; (2) Cauchy-Schranke nahe 0 — genügt auf einer festen
Scheibe |x| ≤ r₀ < |L|. Also: **lokal beschränkt auf D + Korridor (Kompakta von V₀) + a₀^{(N)} ≤ ρ₀ ⇒ p_N → f_K lokal gleichmässig auf D.** Die Schranke auf D̄
ist überflüssig. ABER: die Eindeutigkeit (kn.satz.halbexp_beschraenkt, Schritt 4–5: einseitige Schranke für Re θ̃ bei ALLEN Höhen Im w, Bôcher) braucht den Korridor
bis an L heran (grosse Im w ↔ Nähe zu L). Die Ecke trägt also die Identifikation mit Kneser — periodische Störungen θ wachsen genau zur Ecke hin.

### 23.3 Die Kandidaten (a) und (b) gegen die Kette geprüft (Vermutungsebene; alle Aussagen „falls die Kette hält“)
Wo der primale Wachstumsbereich liegt (aus F_P und der Drehung, k41): im physikalischen Bild eine Zone um L vom Radius r_* ≈ 10⁻⁶ an der Γout-Seite des Randes,
die mit kleinerem r spiralförmig ins Innere von D dreht (bei r ~ r_*e^{−k} um 1.62°·k von der Γout-Tangente nach innen; die innere Normale erreicht sie erst bei
r ~ 10⁻³⁰, den Korridor V₀ bei r ~ 10⁻¹⁶). Folgen:
- **(a) lokal gleichmässig in D + beschränkt auf D̄:** die D̄-Schranke ist nicht robust (Zone am Rand). Auch „lokal gleichmässig in D“ ist nicht ganz robust: Kompakta
  von D, die L näher als r_* kommen, können Punkte der nach innen drehenden Zone enthalten.
- **(b) gleichmässig auf D̄ ohne Umgebung von L, L̄:** robust für eine FESTE Umgebung vom Radius r₀ ≳ r_* (mit Rand: r₀ = 10⁻³ ist sicher), nicht für schrumpfende.
- **Robuster Kern: gleichmässige Konvergenz auf K_{r₀} := D̄ ∖ (U_{r₀}(L) ∪ U_{r₀}(L̄)), insbesondere auf {|x| ≤ |L| − r₀}.** Daraus folgen alle Kernaussagen der Polynom-Wurzel-
  Vermutung: a_k^{(N)} → a_k (Cauchy auf einem kleinen Kreis um 0) und p_N(1) → T_e(½). Die Identifikation mit Kneser kommt dann aus der Referenz f_K in der
  Fehlerschätzung (B-Max auf K_{r₀}), NICHT aus dem Eindeutigkeitssatz — eine andere halbe e-Funktion f_θ unterscheidet sich von f_K schon auf {|x| ≤ 1.2} um O(|θ|).

### 23.4 Was die Beweisarchitektur dann braucht (Satz K15, bedingt; Beweis = K13 + Bernstein–Walsh punktweise in z + Cauchy über Schlitz und Γ_d)
Für z ∈ K_{r₀} gilt |δ_n(z)| ≤ C n^{p}|δ_n(L)| (v), solange K_{r₀} in der „ruhigen Phase“ der Spirale liegt (r₀ ≳ r_*; dort ist die äussere Primalrate ≤ 0, nahe L
der Potenzsektor mit Abfall |u|^{−3.45}, fern von L exponentiell klein). Dann
  **sup_{K_{r₀}}|Π_Nf_K − f_K| ≤ C·N^{κ−4.8}  ⇐  (i″) |δ_n(L)Λ_n(L)| ≤ Cn^{β} (Eckprodukt, NICHT ‖c_n‖), (ii) adjungierte Anpassung im kopfweisen Eck-Einsatz,
  (iii) Stegdämpfung, (v) ruhige Phase auf K_{r₀}.**
Gegenüber der Fassung „gleichmässig auf D̄“ entfallen: (S_γ), (H1), (H*_G), jede gleichmässige Aussage auf einer festen Kurve durch L, und die primale Lokalisierung
am Rand (iv). Die Ecke bleibt nur noch in zwei abgeschwächten, kopfweisen Formen: dem Eckprodukt (ein Skalar je n; vermutlich ~ n^{−0.105} nach der Anpassung,
mit Gs Daten verträglich) und der adjungierten Eckform für Daten mit Singularität in L (I_n^{Ecke} explizit über Mellin–Parseval, §22.2). (i″) ist robust gegen die
Zone: dort ist das Produkt der Köpfe nicht exponentiell klein, eine Zone kann das Eckprodukt über die Normierung höchstens verkleinern (Skizze, nicht bewiesen).
**H1/H2 mit Kugeln in sup_{K_{r₀}}:** formal tragen die Schritte (a-priori-Schranke → implizite Funktionen → Überlagerung → ein Blatt), weil P_N endlichdimensional
ist (jede Norm mit innerem Kompaktum ist äquivalent, 𝒵 bleibt kompakt) und T_N𝒥_P nur vom Keim abhängt (§10.1). **Neue Lücke:** die a-priori-Schranke braucht
Satz K15 für die Sekantenoperatoren 𝒥_P an Fast-Lösungen; der Eck-Einsatz verlangt dafür die Analytizität der Daten 𝒥_P f = f∘P + M_P f längs der kopfweisen Konturen
nahe L — also, dass P(x) für x in der Eckzone die Schlitze von f meidet. Mit der schwachen Kugel ist P nahe L nicht kontrolliert; nötig ist eine Eckbedingung an P
in der Kegelzone (P(L_P) = L_P, Multiplikator μ_P ≈ μ, P ≈ f auf der Skala r_N) — schwächer als sup_{D̄}, aber nicht trivial. H2 (σ-Homotopie) hat Daten x^M
(ganz): dort braucht man nur |c_M|·|δ^P_M(L_P)| (Eckwert des Schritts) und den Abfall ins Innere — die Zone stört nicht.
**Zur Frage der Leitung, ob das erklärt, warum alle Eckversuche scheiterten:** teilweise ja. Alle bisherigen Zielaussagen ((a) Operatornorm, (S), (H*_G), (S_γ), (H1),
Eckmoden auf dem Kreis) verlangten gleichmässige Kontrolle AN L bzw. längs ∂D nahe L in FESTER Geometrie — genau das, was die logarithmische Drehung des Eckmusters
(Periode e^{2π/Im a} in N) für alle N unmöglich macht. Was übrig bleibt, ist kopfweise und lokal auf der Kegelskala — dort ist die Ecke jetzt explizit (K10/K11).
Nicht erklärt ist damit, ob das Eckprodukt (i″) beweisbar ist; das ist die eine globale Zahl, die bleibt.

## 24. Das Eckprodukt ist eine Invariante: Rayleigh-Quotient des Eckkerns, P_n = aCL/(n d_n) (12.10. 09:10–09:14; k46)

### 24.1 Exakte Identität (bewiesen, eine Zeile)
K_n(x, t) := ε_xᵀ J̃_n⁻¹ ẽ_t = t·(J_n⁻¹T_nk_t)(x)  (ε_x = ((x/ρ)^k), ẽ_t = ((ρ/t)^k); Antwort des Abschnitts auf den Cauchy-Kern bei t, ausgewertet in x).
Aus der Rang-1-Ränderung (Satz K5) folgt **K_n(x, t) − K_{n−1}(x, t) = δ_n(x)Λ_n(t)**, insbesondere
  **P_n := δ_n(L)Λ_n(L) = K_n(L,L) − K_{n−1}(L,L),  K_n(L, L) = ẽ_L* J̃_n⁻¹ ẽ_L**   (weil |ρ/L| = 1, also ε_L = ẽ_L‾):
das Eckprodukt (mit der globalen Normierung ⟨ℓ_n, J̃c_n⟩ = 1) ist der Zuwachs EINES Rayleigh-Quotienten von J̃_n⁻¹ am Fourier-Vektor der Eckfrequenz. Die
freien Einzelamplituden kommen darin nicht vor — die Normierung durch den Cauchy-Kern ist fest.

### 24.2 Messung (k46, kubisch, arb; Kontrolle K_n − K_{n−1} = P_n auf alle gedruckten Stellen)
| n | P_n | P_n·n·d_n/(aL) mit komplexem d_n = ρn^{−a} | K_n(L,L) |
|---|---|---|---|
| 200 | 0.30337 − 0.02681i | −0.01021 − 0.59861i | 68.20 − 7.91i |
| 300 | 0.29075 − 0.02268i | −0.01018 − 0.59850i | 97.85 − 10.37i |
| 400 | 0.28212 − 0.01993i | −0.01016 − 0.59846i | 126.46 − 12.49i |
- **P_n ist glatt in n** (n = 200, 201, 202, 203: 0.30337, 0.30321, 0.30305, 0.30289 im Realteil) — die schnelle Phase ω ≈ −1.485 der Einzelköpfe hebt sich im
  Produkt auf.
- **P_n = a·C·L/(n·d_n)·(1 + O(n^{−2}))** mit **C ≈ −0.0102 − 0.5984i** (|C| = 0.5985, arg C = −90.97°; Konvergenz der Spalte ~ n^{−2}); die komplexe Skala
  d_n = ρn^{−a} ist nötig — mit |d_n| dreht die Phase genau um Im(a)·ln n. Also |P_n| = 0.535·n^{−0.106} (n = 400: 0.2834, gemessen 0.2828).
- Summiert: **K_n(L,L) = C·L/d_n + K₀ ≈ (0.580 − 0.149i)·n^{a}** — ein Fisher–Hartwig-artiges Gesetz mit KOMPLEXEM Exponenten a = 1/κ₁ für den Rayleigh-Quotienten
  an der Eckfrequenz (statt des Szegő-Wachstums ~ n/Symbol).

### 24.3 Herkunft (Skizze): Euler-Identität des Eckkerns
Auf der Kegelskala ist der Cauchy-Kern homogen vom Grad −1; ist K_n(L + d_nζ, L + d_nτ) ≈ (L/d_n)·𝒢(ζ, τ) mit festem Grenzkern 𝒢, so ist der Zuwachs in n bei festem x, t
die Euler-Ableitung: **P_n·F_P(λ_Pζ)F_A(λ_Aτ) ≈ (aL/(n d_n))·(1 + ζ∂_ζ + τ∂_τ)𝒢(ζ, τ)**. Weil (1 + ζ∂_ζ + τ∂_τ) den Cauchy-Kern annulliert, löst E𝒢 die homogenen
Eckgleichungen in ζ und τ und ist nach K11/K11′ ein Vielfaches von F_P(λ_Pζ)F_A(λ_Aτ): E𝒢 = C·F_P⊗F_A mit C = 𝒢(0,0). Das gemessene Gesetz ist genau dieser Fall bei
ζ = τ = 0. Integriert längs Strahlen: **𝒢(ζ, τ) = 𝒢_hom(ζ, τ) + C∫₀¹F_P(λ_Psζ)F_A(λ_Asτ)ds** (𝒢_hom homogen vom Grad −1, trägt den Pol) — die Kontinuumsform der
Christoffel–Darboux-Summe K_n = Σ_{m≤n}δ_mΛ_m. Den Wert von C legt diese Skizze nicht fest (er ist wie bei Fisher–Hartwig eine globale Konstante); für (i″) wird er
nicht gebraucht.

### 24.4 Folgerungen
1. **(i″) reduziert sich auf eine einzige skalare Aussage (R_L): |ẽ_L* J̃_n⁻¹ ẽ_L| ≤ C n^{β}** (dann |P_n| ≤ 2C n^β). Numerisch: = 0.5985·n^{0.894}(1 + O(n^{−1})).
   (R_L) ist die Pol-Antwort (H₀) im EINEN Punkt ζ = z = L — die Grenzlage, die die gleichmässige (H₀) nicht abdeckt —, also ein Fisher–Hartwig-Problem für einen
   Rayleigh-Quotienten. **Robust gegen die Spirale:** K_n(L,L) = Σ P_m enthält nur Eckwerte, nicht ‖c_m‖; die Zone bei r ~ 10⁻⁶ liegt für grosse n weit ausserhalb
   der Kegelskala. (G_N(L) dagegen enthält ‖c_m‖ und ist nicht robust.)
2. **Gs Eckamplitude erklärt:** B_n = √n‖c_n‖|Λ_n(L)| = √n·Q_n·|d_n|^{1/2}·|P_n| ∝ Q_n·n^{−0.053}; mit Q_n = 1.53 (n = 400, k45) ergibt das 0.695 (G: 0.697). Das
   in §22.5 vorhergesagte n^{−0.053}-Gesetz folgt also aus dem Eckproduktgesetz und der primalen Lokalität (Q_n → Q_∞); B_n → 0 langsam, falls Q_n konvergiert.
3. Für K15 genügt damit: (R_L), (ii) adjungierte Eckform, (iii) Stegdämpfung, (v) ruhige Phase. Alle Eckaussagen sind jetzt kopfweise-lokal bis auf (R_L), und
   (R_L) ist eine Zahl je n mit einem sauberen Gesetz.

### 24.5 Die Eckbedingung an die Sekanten-Daten (für H1 in sup_{K_{r₀}}) — Formulierung
Die Sekanten-Daten W^P_η = k_η∘P + M_P k_η, M_P = (f∘P − f∘f)/(P − f), f∘f = exp: weil P ein Polynom ist und f_K auf den benutzten Blättern nur die Verzweigungspunkte
L, L̄ hat, ist W^P_η nahe L analytisch bis auf den einen Punkt x_* mit P(x_*) = L (auf der Kegelskala, |x_* − L| ≲ |L_P − L|/|μ|) und die Polstellen P(x) = η.
Die kopfweise Kontur Γ_n kann im Nicht-Wachstumsfenster des adjungierten Kopfes laufen und dabei den primalen Wachstumssektor MEIDEN — die Dualität lässt dafür zwei
Fenster (u_A ∈ (19°, 100°) und (−120°, −61°)), der Übergang geschieht auf der Kegelskala, wo es keine Zone gibt. Gebraucht wird also nur:
**(EB_N) Für alle Fast-Lösungen P in der Kugel gilt |P − f| ≤ η_N auf 𝒲_N := {|x − L| ≤ A₀r_N} ∪ ⋃_{r_N ≤ r ≤ ε}(Ruhefenster bei Radius r) (und konjugiert),
mit η_N ≪ r_N und η_N·ln N ≪ 1** — dann liegen Fixpunkt L_P und Multiplikator μ_P von P so nahe an L, μ, dass die Eckprofile mit μ_P die von μ sind (Drehphase
Im(a_P)·ln N stabil), und die Daten sind längs Γ_n beschränkt. (EB_N) ist selbstkonsistent: P − f = λ(Π^P_Nf − f), und auf 𝒲_N gibt der Eck-Einsatz (gleiche
Abschätzung wie K15, aber an Punkten der Kegelskala und der Ruhefenster) |Π^Pf − f| ~ N^{−4.8+…} ≪ r_N.
**→ V (Messung):** für die echten Iterierten p_N (N ≤ 2000) sup|p_N − f_K| auf der Kegelscheibe |x − L| ≤ 16ρN^{−0.894} und längs der Ruhefenster-Strahlen (physikalischer
Winkel θ = u_A − arg Λ₀ − Im(a)ln N·(180/π), Λ₀ = −λ_Aρ/L², λ_A = −0.203297 − 1.406007i), verglichen mit N^{−4.8} und r_N.

## 25. E1: Die Eckkonstante C ist LOKAL bestimmt — Anpassung der Pol-Antwort an den äusseren Pol (12.10. 09:19–09:23; k46b, k47)

### 25.1 Herleitung (Skizze mit benannten Schritten)
K_n(L, L) = L·u_n(L) mit u_n = J_n⁻¹T_nk_L (Antwort auf den Eck-Pol). Aus der Christoffel–Darboux-Summe u_n(x) = Σ_{m≤n}δ_m(x)Λ_m(L), den Eckprofilen
δ_m(L + d_mζ) ≈ δ_m(L)F_P(λ̃_Pζ) (λ̃_P = −λ_Pρ²/L² in der Koordinate ζ = (x − L)/d_n, d_n = ρn^{−a}) und P_m = aCL/(m d_m) folgt durch Riemann-Summe (s = d_n/d_m,
ds = a(d_n/d_m)dm/m):
  **d_n·u_n(L + d_nζ) → C·V(ζ),   V(ζ) = ∫₀¹ F_P(λ̃_P sζ) ds = Σ_k Γ(β_P − iαk)(−λ̃_Pζ)^k/((k+1)!Γ(β_P))**.
Für |ζ| → ∞ im Potenzsektor von F_P ist V(ζ) ≈ Φ_{F_P}(1)/(λ̃_Pζ) mit Φ_{F_P}(1) = ∫₀^∞F_P = Γ(½ + 2iα)/Γ(½ + iα) (Mellin bei s = 1). Ausserhalb der Kegelskala muss u_n die
Bulk-Lösung s = J_∞⁻¹k_L sein, deren Eckteil der Pol c_s/(x − L) mit **c_s = −1/(μ + 1/μ)** ist (lokale Gleichung s(L + μz) + μs(L + z) = −1/z). Anpassen:
  **C = −λ̃_P·Γ(½ + iα) / ((μ + 1/μ)·Γ(½ + 2iα)) = −0.010118 − 0.598298i.**

### 25.2 Prüfung
- **Konstante:** aus k46/k46b (P_n·n·d_n/(aL), d_n komplex) C_n = −0.01021−0.59861i / −0.01018−0.59850i / −0.01016−0.59846i / −0.01014−0.59842i / −0.01015−0.59840i
  (n = 200/300/400/500/600); mit C_n = C + D/n aus 300/600 extrapoliert **C = −0.01012 − 0.59830i — Vorhersage −0.010118 − 0.598298i, Übereinstimmung auf 1e-5**.
- **Profil (k47):** d_n·u_n(L + d_nζ) gegen C·V(ζ) bei n = 400: ζ = 0: −0.00994−0.60097i / −0.01012−0.59830i; ζ = 1: −0.29785−0.50460i / −0.29565−0.50282i; ζ = 2i: 0.07858−2.11460i /
  0.07577−2.10980i; ζ = 2+2i: −1.7648−0.2247i / −1.7404−0.2123i; die Abweichungen (0.3–1.5 %) schrumpfen von n = 200 auf 400 um den Faktor ≈ 0.6 (Ordnung n^{−a}).

### 25.3 Folgerungen
1. **Das Eckprodukt (i″) ist keine globale Grösse mehr:** P_n = aCL/(n d_n)(1 + O(1/n)) mit dem lokal bestimmten C; eingehende Daten sind nur μ (über α, μ + 1/μ) und die Skala
   λ_P = −μ^{−1/2}λ_A der Köpfe — und λ_A ist die Anpassungskonstante an den äusseren Kegel (Ws c₁). Für eine SCHRANKE wird nicht einmal deren Wert gebraucht.
2. **(R_L) ist eine Aussage über EINE Antwort:** u_n = J_n⁻¹T_nk_L hat die Innenform (C/d_n)V(ζ) und die Aussenform J_∞⁻¹k_L; (R_L) heisst „angepasste Asymptotik dieser einen
   Antwort mit polynomialer Fehlerschranke“ — derselbe Typ wie (ii). Damit ruht K15 auf: (ii) kopfweise adjungierte Eckform, (R_L) Eck-Pol-Antwort, (iii) Stegdämpfung,
   (v) ruhige Phase — alles lokale bzw. äussere Aussagen ohne gleichmässige Randkontrolle.
3. **Beweisweg für (R_L) (Plan):** (a) Aussen: s = J_∞⁻¹k_L existiert (O7 + Injektivität) mit Eckentwicklung c_s/(x − L) + (homogene Eckterme) + regulär; (b) Innen: C·V/d_n
   explizit, C aus 25.1; (c) zusammengesetzte Näherung ũ_n und Residuum T_n[J_∞ũ_n − k_L]; (d) Fehler u_n(L) − ũ_n(L) = Σ_m δ_m(L)⟨ℓ_m, Residuum⟩ — kopfweise Paarungen
   (braucht (ii), (iii)) und Bootstrap mit |P_m| ≤ |K_m| + |K_{m−1}|. **Kritischer Punkt:** die Köpfe m ≈ n sehen das Residuum der Stufe n mit relativer Grösse ~n^{−a};
   die naive Summe über m ∈ [n/2, n] verliert einen Faktor n^{0.106}. Nötig ist dort die Phasen-/Skalenstruktur der Paarung (Mellin–Parseval zwischen F_A und dem
   Residuumsprofil) statt Beträgen — die Gammafaktoren sind bekannt, das ist eine explizite Rechnung.
4. Die „String-Gleichung“ n·P_n − a·K_n(L,L) = κ₀ + O(1/n) (gemessen κ₀ ≈ −0.48) ist die diskrete Form der Euler-Identität bei ζ = τ = 0; aus ihr allein folgt (R_L) durch
   Summation. Sie exakt herzuleiten wäre der kürzeste Weg (offen).

## 26. E1 weiter: Die String-Gleichung ist NICHT exakt — κ₀ ist global; Bootstrap für (R_L) mit glatten Zuwächsen (12.10. 09:24–09:27)

### 26.1 Befund zur String-Gleichung (numerisch, k46/k46b)
Mit dem lokal hergeleiteten C (§25) zerlegt: K_n(L,L) = C·L/d_n + K₀(n), gemessen K₀(n) = 0.5513 / 0.5557 / 0.5591 / 0.5619 / 0.5659 / 0.5691 / 0.5716 (n = 150 … 600), Drift wie
n^{a−1} = n^{−0.106}: K₀(n) ≈ K₀∞ − 0.253·n^{−0.106}, **K₀∞ ≈ 0.70** (Fit mit festem Exponenten; mit Exponent 0.08 … 0.2 passen alle Fits auf 1e-4 und
geben K₀∞ = 0.74 … 0.64 — der Grenzwert ist aus n ≤ 600 nur auf ±0.05 bestimmt). Entsprechend
  **n·P_n − a·K_n(L,L) = −a·K₀∞ + O(n^{a−1}),  κ₀ = −aK₀∞ ≈ −0.63**
(nicht −0.48 mit O(1/n), wie in §25.3 vermutet: die gemessene Folge −0.461 … −0.484 driftet weiter). K₀∞ = lim(Σ_{m≤n}P_m − CL/d_n) hängt an ALLEN frühen Köpfen —
eine globale Konstante; eine geschlossene Form ist nicht zu erwarten. Die String-Gleichung ist damit genau die diskrete Euler-Identität mit einem globalen
Achsenabschnitt; aus K5 + Euler-Identität folgt sie nur, wenn man die Konvergenz des Eckkerns auf der Kegelskala schon hat. **Der Weg „String-Gleichung exakt
⇒ (R_L) durch Summation“ trägt nicht.** (Für (R_L) ist κ₀ unwichtig — es geht nur um die Grössenordnung n^{Re a}.)

### 26.2 Bootstrap für (R_L) mit der richtigen Induktionsannahme (Plan, Rechnung explizit machbar)
Zusammengesetzte Näherung ũ_n = s + (C/d_n)V(ζ) − c_s/(x − L) (s = J_∞⁻¹k_L), Residuum r̃_n := T_nk_L − J̃_nT_nũ_n, Fehler e_n := u_n(L) − ũ_n(L) = Σ_{m≤n}P_m·ρ_{m,n}
mit den RELATIVEN Paarungen ρ_{m,n} := ℓ_mᵀr̃_n/Λ_m(L). Exakt P_n = L[ũ_n(L) − ũ_{n−1}(L)] + L[e_n − e_{n−1}], und L[ũ_n(L) − ũ_{n−1}(L)] = aCL/(n d_n)(1 + O(1/n)).
- Mit der groben Annahme |K_m| ≤ M/|d_m| verliert man bei m ≈ n einen Faktor n^{0.106} (|P_m| ≤ |K_m| + |K_{m−1}| verschenkt das 1/m).
- **Richtige Induktionsannahme (glatte Zuwächse):** (B_n) |P_m − aCL/(m d_m)| ≤ θ|aCL/(m d_m)| für m < n. Dann
  |e_n − e_{n−1}| ≤ |P_nρ_{n,n}| + Σ_{m<n}|P_m||ρ_{m,n} − ρ_{m,n−1}| ≲ (|C|/(n|d_n|))·ε(n)·(1 + ∫₀¹g(x)x^{−0.106}dx),
  sofern **Lemma ρ** gilt: |ρ_{m,n}| ≤ ε(n)g(m/n), |ρ_{m,n} − ρ_{m,n−1}| ≤ (ε(n)/n)g(m/n) mit ε(n) → 0 (erwartet n^{−a}, Profilfehler) und g integrierbar gegen
  x^{−0.106}. Für n ≥ n₀ (ε(n₀) klein) schliesst die Induktion, und (B_n) gilt für alle n ⇒ **(R_L)**, sogar mit dem Gesetz P_n ~ aCL/(n d_n).
- **Lemma ρ** ist eine explizite Mellin–Parseval-Rechnung: ρ_{m,n} ist die Paarung des adjungierten Eckprofils F_A auf der Skala d_m mit dem Residuumsprofil auf der Skala
  d_n (Gammafaktoren bekannt), plus der ferne Teil (Stegdämpfung). Es braucht also genau **(ii)** (Eckform der adjungierten Köpfe samt Anpassungsfehler) und **(iii)**.

### 26.3 Konsequenz für die Reihenfolge (ii), (iii), (v)
Alle vier Restpunkte — (R_L) über Lemma ρ, (ii), (iii), (v) — hängen an EINEM Kern: einem **Anpassungslemma für Köpfe** (der wahre Kopf ist auf der Kegelskala
Eckwert × F_{A/P} mit relativem Fehler O(n^{−a}·poly|ζ|) und folgt im Ruhefenster bzw. fern von L der äusseren WKB-Lösung mit Rate ≤ 0). Numerisch ist es belegt
(Defekt-Tabellen k23/k29/k32, k40: Kopf − Modell schrumpft wie n^{−a} auch bei A ≈ 24–32; Profile auf 1e-7). Zu beweisen ist es als Programm „innere Parametrix
explizit (K10/K11) + äussere WKB (Ws geschweisste Aufgabe) + Kleinheitsargument“ — das ist der verbleibende Forschungsschritt. Leichteste Reihenfolge: (ii) zuerst
auf der Kegelskala (rein lokal, explizite Funktionen), dann (v) (primal ins Innere, Potenzsektor — gleiche Technik), dann (iii) (global, Steg, Vs Schleifenbedingung).

## 27. E1 als Reduktion abgeschlossen: Satz K16 (Bootstrap) — (R_L) ist kein eigener Restpunkt; Integraldarstellung der Eckprofile (12.10. 09:27–09:38; k48)

### 27.1 Satz K16 (Bootstrap für das Eckprodukt) — bewiesen (bedingt auf (ρ1), (ρ2))
Seien ũ_n ∈ ℂ^{n+1} beliebige Vergleichsvektoren, r̃_n := ẽ_L − J̃_nũ_n, M_n := ε_Lᵀũ_n, p_n := M_n − M_{n−1}; ℓ_m = J̃_m^{−T}e_m, Λ_m(L) = ⟨ℓ_m, ẽ_L⟩_m ≠ 0
(⟨·,·⟩_m: Paarung der ersten m+1 Komponenten). **Exakt** (Rang-1-Ränderung K5: J̃_n⁻¹ = Σ_{m≤n} c_mℓ_mᵀ/ℓ_mm, aufgefüllt mit Nullen):
  **e_n := K_n(L,L) − M_n = ε_LᵀJ̃_n⁻¹r̃_n = Σ_{m≤n} P_m·ρ_{m,n},   ρ_{m,n} := ⟨ℓ_m, r̃_n⟩_m / Λ_m(L)**  (relative Paarung: Kopf m mit dem Residuum, geteilt durch Kopf m am Eck-Pol).
Annahmen: (ρ1) |ρ_{n,n}| ≤ ε_n ḡ; (ρ2) |ρ_{m,n} − ρ_{m,n−1}| ≤ (ε_n/n)·g(m/n) für m < n, mit ε_n ↓ 0, 0 ≤ g ≤ ḡ; (P) |p_n| ≥ c_p n^{−b} und
Σ_{m≤n}|p_m|g(m/n) ≤ G·n|p_n| (für p_n ~ n^{−b} heisst das G ≈ ∫₀¹x^{−b}g(x)dx < ∞).
**Aussage:** Mit B := max_{m<n₀}|P_m| und G' := 2G + Bn₀^{b}ḡ/c_p gelte ε_{n₀}ḡ ≤ ¼, ε_{n₀}G' ≤ ½. Dann für alle n ≥ n₀:
  **|P_n − p_n| ≤ ε_n(2ḡ + G')·|p_n|,  |K_n(L,L) − M_n| ≤ |e_{n₀}| + (2ḡ + G')Σ_{n₀<m≤n} ε_m|p_m|.**
*Beweis.* P_n = p_n + (e_n − e_{n−1}) und e_n − e_{n−1} = P_nρ_{n,n} + Σ_{m<n}P_m(ρ_{m,n} − ρ_{m,n−1}). Induktionsannahme |P_m| ≤ 2|p_m| für n₀ ≤ m < n:
|e_n − e_{n−1}| ≤ ε_nḡ|P_n| + (ε_n/n)[Bn₀ḡ + 2Gn|p_n|] ≤ ε_nḡ|P_n| + ε_nG'|p_n| (wegen Bn₀/n ≤ Bn₀^{b}|p_n|/c_p für n ≥ n₀). Also
|P_n|(1 − ¼) ≤ |p_n|(1 + ½), d.h. |P_n| ≤ 2|p_n|, und die erste Schranke folgt; die zweite durch Summation. ∎
**Lesart.** Für die zusammengesetzte Vergleichslösung (§25.3: aussen s = J_∞⁻¹k_L, innen C·V/d_n) ist M_n = L·s_reg(L) + CL/d_n explizit, p_n = aCL/(n d_n)(1 + O(1/n)),
b = 1 − Re a = 0.106. Gilt Lemma ρ mit ε_n = O(|d_n|) (Krümmungsterm f''(L), relativ O(d_nζ) auf der Kegelskala) und beschränktem g, so folgt
**P_n = aCL/(n d_n)(1 + O(n^{−0.894})) und K_n(L,L) = CL/d_n + O(log n)** — also (R_L) mit dem Gesetz, und die Konstanten (n₀, B, G', ε) sind für n₀ ≤ 400
explizit prüfbar (B aus arb-Werten P_m, k46). Die Daten zeigen sogar K₀(n) → const (§26.1), d.h. das log n tritt nicht auf.

### 27.2 Was Lemma ρ braucht (Skizze; Glieder benannt)
ρ_{m,n} ist ein Konturintegral: ⟨ℓ_m, r̃_n⟩_m = (1/2πi)∮ Λ_m(t) r̃_n(t) dt/t (Λ_m Laurent-Polynom in 1/t, r̃_n Polynom; jede Kontur um 0 erlaubt) — also die
KOPFWEISE Kontur Γ_m von Lemma K13 (Eckteil im Nicht-Wachstumsfenster des adjungierten Kopfes, ferner Teil über den Steg). Glieder:
- **(ii′)** auf dem Eckteil Λ_m(L + d_mτ)/Λ_m(L) = F_A(λ̃_Aτ)(1 + O(|d_m|(1+|τ|)^q)) — das ist (ii);
- **(iii′)** auf dem fernen Teil |Λ_m/Λ_m(L)| ≤ C m^{q′}·(exponentiell klein) — das ist (iii);
- **(K-Param)** eine POLYNOMIALE Vergleichslösung ũ_n, deren Residuum r̃_n längs aller Γ_m (m ≤ n) relativ O(|d_n|) ist. Die naive Wahl T_n[s + (C/d_n)V(ζ) − c_s/(x − L)] taugt
  nicht: V hat Ordnung 1/Re a = 1.12 und wächst in seinen Wachstumssektoren bis exp(c·n) am Rand der Kegelzone, T_n verschmiert das. (K-Param) ist die primale
  Parametrix der Eck-Pol-Antwort — dieselbe Art Aussage wie (v) (primaler Kopf im Inneren / Potenzsektor), nur für die eine Antwort u_n statt für die Köpfe.
- **Paarung (explizit):** mit (ii′) und r̃_n(L + d_nζ) = R₁(ζ)(1 + o(1)) ist ρ_{m,n} = d_n·𝓜((m/n)^a) + (fern), 𝓜(s) = (1/2πi)∫_Γ F_A(λ̃_Asζ)R₁(ζ)dζ, per Mellin–Parseval
  𝓜(s) = (1/2πi)∫Γ(w)Γ(β_A + iαw)(λ̃_As)^{−w}R̂₁(1 − w)dw/Γ(β_A); bei m ≈ n (s ≈ 1) ist das eine feste Zahl ⇒ (ρ1), (ρ2) mit ε_n = O(|d_n|) — der in §25.3
  befürchtete Verlust n^{0.106} tritt mit der Induktionsannahme „glatte Zuwächse“ nicht auf. Für m ≪ n entscheidet das Wachstum von R₁ längs Γ_m über g(x) bei
  x → 0; nötig ist nur ∫₀¹x^{−0.106}g(x)dx < ∞ (für R₁ beschränkt: g(x) ≲ x^{−0.894}, dann G logarithmisch — Satz K16 bleibt mit G' ~ log n und gibt (R_L) mit
  einem Faktor log n; das ist für K15 unerheblich).
**Folgerung für die Architektur:** (R_L) und damit (i″) folgen aus (ii), (iii) und (K-Param) ≈ (v). **K15 ruht auf drei Kopfaussagen: (ii) adjungierte Eckform,
(iii) Stegdämpfung, (v) primale Form im Inneren (einschliesslich der Eck-Pol-Antwort).** Alle drei sind „angepasste Asymptotik eines Kopfes mit polynomialer
Fehlerschranke“: innen explizit (K10/K11), aussen WKB (Ws Aussenlösung), dazwischen ein Kleinheitsargument auf dem Anpassungsring. Der globale Teil steckt in der
Aussenlösung — das ist die Stelle, an der W12 (11.10.) nicht geschlossen hat; ehrlich: der verbleibende Kern ist damit nicht kleiner geworden, aber er ist jetzt
nur noch EINE Art Aussage, und die Eck-Seite davon ist vollständig explizit.

### 27.3 Nebenergebnis: Integraldarstellung der Eckprofile — schliesst die Lücke in K10 (iii) für den Potenzsektor (bewiesen + numerisch, k48)
**F_β(u) = (1/Γ(β)) ∫₀^∞ t^{β−1} e^{−t − u·t^{−iα}} dt** für alle u ∈ ℂ (gliedweise Integration: ∫t^{β−iαk−1}e^{−t}dt = Γ(β − iαk); Majorante e^{|u|t^{Im α}}e^{−t} mit
Im α = 0.106 < 1). Also F_β(u) = „Erwartung“ von exp(−u·T^{−iα}) unter dem (komplexen) Gamma-Gewicht — jedes Eckprofil ist eine Überlagerung von Exponentialen
e^{−θu} mit θ = t^{−iα} auf einer Spirale, d.h. lokal von Monomen (x/L)^ν komplexen Grades. k48: Integral = Reihe auf 1e-45 (A) bzw. 1e-25 (P) an fünf Punkten.
**Potenzsektor (Watson-Lemma, Skizze vollständig):** v = u·t^{−iα}, t = (v/u)^γ, γ = i/α = 8.894 + 2.116i:
  **F_β(u) ~ (γ/Γ(β)) Σ_{j≥0} (−1)^j Γ(γ(β + j))/j! · u^{−γ(β+j)}**, Exponenten γ(β_A + j) = 2 + (2j+1)·2πi/L = s_j^A (genau die Pole aus §18.1), γβ_P = 3.447 + 1.058i.
  k48: drei Glieder gegen die Reihe bei |u| = 40, 80, arg u ∈ {−30°, 0°, 30°}: relative Abweichung 4e-12 … 7e-10 (A), 1e-16 … 2e-12 (P); bei |u| = 20 noch 4e-4 … 4e-3
  (A) — der Rest ist exponentiell klein (Sattelbeitrag), nicht algebraisch. Der Gültigkeitssektor ist der, in dem sich das Bild der t-Achse auf den Strahl arg v = 0
  drehen lässt, ohne dass e^{−t} verloren geht; der Wachstumssektor ist der Sattel t^a = iαu (Laplace-Methode, liefert die Formel von §19.2 (iii)). Damit ruht K10 (iii)
  im Potenzsektor nicht mehr auf Braaksma; für die Ränder der Spiralsektoren und den Wachstumssektor ist es eine Standard-Sattelpunktrechnung an EINEM Integral.

## 28. Eine Aussenhypothese (AH) „Eckdominanz des Kernzuwachses“ — (i″), (ii), (iii), (v) sind Faktoren EINER Ungleichung (Auftrag Leitung 09:40; 09:40–10:01; k49)

### 28.1 Exakte Umformung (bewiesen)
(a) **Kernzuwachs.** ΔK_n(z, t) := K_n(z, t) − K_{n−1}(z, t) = δ_n(z)Λ_n(t)/ℓ_nn (K5), K_n(z, t) = ε_zᵀJ̃_n⁻¹ẽ_t = t·(J_n⁻¹T_nk_t)(z) (Antwort in z auf den Cauchy-Kern
bei t). Lemma K13 lautet damit  **J_N⁻¹T_NW(z) = Σ_{n≤N} (1/2πi)∮_{Γ_n} ΔK_n(z, t)·W(t) dt/t** — jeder Summand ist ein Konturintegral des Kernzuwachses, und
K15 benutzt von den Köpfen nur |ΔK_n(z, t)| auf K_{r₀} × Γ_n.
(b) **Rang 1.** ΔK_n(z, t) = ΔK_n(z, L)·ΔK_n(L, t)/P_n mit P_n = ΔK_n(L, L), also |ΔK_n(z, t)| = |P_n|·|δ_n(z)/δ_n(L)|·|Λ_n(t)/Λ_n(L)|. Die vier Restpunkte von K15 sind
Einschränkungen EINER Funktion: (i″) auf (L, L), (v) auf K_{r₀} × {L}, (ii) auf {L} × (Eckteil der Kontur), (iii) auf {L} × (ferner Teil). Zudem
δ_n(z)/δ_n(L) = Δu_n(z)/Δu_n(L) mit der Eck-Pol-Antwort u_n = J_n⁻¹T_nk_L, und Λ_n(t)/Λ_n(L) = Δw_n(t)/Δw_n(L) mit w_n(t) = K_n(L, t).
(c) **Maximumprinzip = innere Hälfte in Schrankenform.** δ_n ist ein Polynom ⇒ max_{K_{r₀}} = max_{∂K_{r₀}}. Λ_n ist analytisch in t ≠ 0 ⇒ auf jeder Scheibe
|t − L| ≤ A₀r_n (r_n = ρn^{−Re a}) beherrscht der Kegelkreis C_n das Innere. Auf C_n sagt K10 den Wert voraus: max_{C_n}|Λ_n/Λ_n(L)| → M(A₀) = max_{|u|=|λ_A|A₀}|F_A(u)|
= 3.62 / 13.6 (A₀ = 1 / 2); gemessen (k49, n = 200 / 400 / 600) 3.64 / 3.63 / 3.63 bzw. 13.85 / 13.73 / 13.69 — auch in den Wachstumssektoren von F_A.
Die asymptotische Form der inneren Hälfte (Λ_n/Λ_n(L) → F_A) braucht man nur noch für das GESETZ von P_n (K16), nicht für die Schranken von K15.

### 28.2 Eine feste Kontur durch L für alle n (bewiesen aus der Kegelform; numerisch k49)
Die Ränder des Nicht-Wachstumsfensters im Profil sind arg u = c + T·ln|u|, T = tan(arg a) = 0.02833 (K10 (iii), §27.3). Physikalisch (u = λ(t − L)/d_n,
arg d_n = −Im(a)·ln n, |d_n| = ρn^{−Re a}): θ = c′ + T·ln r + (T·Re a − Im a)·ln n, und die letzte Klammer ist null. **Die Fensterränder sind feste Log-Spiralen
θ(r) = θ_j + 1.623°·ln(r/ρ), unabhängig von n** (dieselben wie Ws Null-Linien und §21.1; Vorzeichenfehler in §20.2 (a): richtig ist θ₀ − 0.0284·ln(1/r), Drehung im
Uhrzeigersinn für r → 0). Also kann der Eckteil ALLER kopfweisen Konturen von K13 eine feste Kurve sein:
  **Σ₁: θ(r) = u_A + 71.46° + 1.623°·ln(r/ρ), u_A ≈ 60–80° (ausserhalb D, zwischen S_L und Γin);  Σ₂: dasselbe mit u_A ≈ −90° (in D);  0 < r ≤ r₁**
(Vs Konvention arg Λ₀ = −71.46°, V §3.18). §20.2 (c) („jede feste Kurve durch L liegt irgendwann auf der Wachstumsseite“) gilt für Kurven mit fester Richtung in L —
Σ_j drehen mit. In §21 war die Log-Spirale als LINSE (K1-Rahmen, f(ℬ) ⊂ D) ausgeschieden; als K13-Kontur gibt es diese Bedingung nicht, nur die Analytizität der Daten
längs Σ_j (Fortsetzung über S_L wie in K13; für Σ₁ ausserhalb D die Überkonvergenz: V §3.18 misst p_N − f_K bei u_A = 85° bis r = 0.05 klein, nahe S_L nicht — daher
u_A ≈ 80°, r₁ = 0.05).
k49 (kubisch): längs Σ₁ (u_A = 40° / 60° / 80°) und Σ₂ (u_A = −70° / −90°), r ∈ [1e-7, r₁]: **max|Λ_n/Λ_n(L)| = 1.000 … 1.008 für n = 200, 400, 600** — längs beider
Spiralen ist der Kopf nirgends grösser als in der Ecke (das 1.008 liegt bei A ≈ 0.2). Am Fensterrand (u_A = 95° bzw. −110°): 1.35 bzw. 3.84 bei r = 0.05 (n = 200).

### 28.3 Die Hypothese (AH) „Eckdominanz“ — eine Grösse, ein Gebiet, eine Rate
Sei r₀ ∈ (r_*, 0.1] (r_* ≈ 1e-6 die Zone von §22.6), 𝒦 := K_{r₀} ∪ {L}, und Γ* := Σ₁ ∪ Σ₂ ∪ Γ*_fern ∪ (konjugiert bei L̄), Γ*_fern := Vs Kontur γ (Γout-Bogen + Keil-Linse
c = 0.25 über Γin) ausserhalb der r₁-Scheiben um L, L̄, mit kurzen Verbindungsbögen im Fenster.
  **(AH) Es gibt C, κ, n₀ mit  Q_n := max_{z ∈ 𝒦, t ∈ Γ*} |ΔK_n(z, t)| ≤ C·n^κ  für alle n ≥ n₀.**
Gleichwertig (Rang 1, 28.1 b/c): **Q_n = |P_n|·R_P(n)·R_A(n)**, R_P = max_{∂K_{r₀} ∪ {L}}|δ_n/δ_n(L)| ≥ 1, R_A = max_{Γ*}|Λ_n/Λ_n(L)| ≥ 1. In Worten: **die Köpfe sind (bis auf
n^κ) in der Ecke am grössten, und das Eckprodukt ist polynomial.**
- Grösse: ΔK_n = c_nℓ_nᵀ/ℓ_nn ist genau der Kopfterm von Vs τ_k, hier punktweise statt in Norm; Vs UL-Elimination liefert alle n ≤ N in einem Gang.
- Gebiet: zwei Kurven (∂K_{r₀} und Γ*), dank des Maximumprinzips.
- Rate (k49, kubisch, r₀ = 0.05; Γ*_fern nur grob als |t| = ρ(1 + 0.0166) über Γin und |t| = ρ über Γout, ohne 12° um L, L̄):

| n | |P_n| | R_P (Rand K_{0.05}; mit z = L) | R_A Spiralen | R_A fern (Lage) | Q_n |
|---|---|---|---|---|---|
| 200 | 0.3046 | 0.82 → 1 | 1.008 | 2.17 (180°, Steg) | 0.661 |
| 400 | 0.2828 | 0.51 → 1 | 1.008 | 1.16 (180°) | 0.328 |
| 600 | 0.2709 | 0.32 → 1 | 1.008 | 0.58 (180°) | **0.273 = 1.008·|P_n|** |

  (R_P für r₀ = 0.1 / 0.02: 0.58 / 1.00, 0.28 / 0.84, 0.14 / 0.66 für n = 200 / 400 / 600; Lage auf dem Kreis 1–4° neben L bzw. L̄ oder auf dem Bogen um L.)
  **Ab n ≈ 600 liegt das Maximum in der Ecke: Q_n = 1.008·|P_n| = 0.54·n^{−0.106} (strenge Eckdominanz, κ = −0.106).** Die ferne Steg-Spitze bei 180° fällt
  exponentiell (Faktor ≈ 0.5 je 200 Grade), R_P fällt beschleunigt (lokaler Exponent −1.0 → −1.7 für r₀ = 0.1, −0.7 → −1.1 für r₀ = 0.05).

### 28.4 Folgerung und Antwort auf die Fragen der Leitung
- **[ZURÜCKGEZOGEN 10:23, siehe §29.1 — die Datenseite nahe L trägt nur bis N ≈ 1e8…1e11]** ~~Satz K17 (bewiesen, bedingt auf (AH)):~~ Unter (AH) gilt die Schlussfolgerung von K15 auf K_{r₀}. *Beweis:* Lemma K13 mit Γ_n := Γ* für alle n (zulässig, sofern die
  Daten längs Γ* analytisch sind — dieselbe Datenbedingung wie in K15 und (EB_N)); 28.1 (a) gibt |Summand n| ≤ Q_n·(1/2π)∫_{Γ*}|W||dt|/|t|; Summe über n ≤ N gibt den Faktor
  Σ_{n≤N}Q_n ≤ C N^{1+κ}; die Datenseite wie in K15. ∎ Weil Γ* fest ist, folgt summiert sogar |K_N(z, t)| ≤ C N^{1+κ} auf 𝒦 × Γ* — **(AH) ist die kopfweise Fassung von
  Ws (H₀), auf einer mitdrehenden Spiralkontur statt auf einer Kurve mit fester Richtung in L.**
- **Gehen (iii) und (v) auf dieselbe Hypothese zurück? Ja** — und (ii), (i″) auch: alle vier sind Faktoren von |ΔK_n| auf 𝒦 × Γ*. Es bleibt genau EINE benannte
  Hypothese, (AH). K16 (Gesetz von P_n) und K10/K11 (Werte auf C_n, Fenster) erklären und sagen die Zahlen voraus; logisch nötig sind sie für K17 nicht.
- **Ehrlich zur Herkunft:** exponentiell ist (1/n)ln|ΔK_n(z, t)| → H_P(z) + H_A(t). H_A ist Ws H (normiert auf H(L) = 0); (AH) verlangt H_A ≤ 0 auf Γ* — nach Ws Daten ja
  (Fenster nahe L, Linse fern). Zusätzlich verlangt (AH) H_P ≤ 0 auf K_{r₀} — die primale Rate hat W NICHT berechnet (Ws Grenzproblem betrifft die adjungierten Moden).
  Polynomial sind es die Vorfaktoren der starken Asymptotik (Ws Frage 2). Beweisgegenstand für später also: **starke Asymptotik beider Kopffamilien im Aussenbereich mit
  polynomialem Fehler, und Eckdominanz der Raten (H_A ≤ 0 auf Γ*, H_P ≤ 0 auf K_{r₀})**; nahe L ist alles explizit (F_A, F_P fallen im Potenzsektor, Fenster fest).

### 28.5 Messprotokoll für V (Vorschlag E4)
- Objekte: c_n, ℓ_n, ℓ_nn für alle n ≤ 2000 aus der UL-Elimination, kubisch und Kneser. Auswertung am Mittelpunkt mit ≥ 1.3n + 400 Bit (acb-Horner nahe |x| = ρ bläht die
  Radien auf, V §3.18).
- Mengen: ∂K_{r₀} für r₀ = 0.1 / 0.05 / 0.02 (Kreis |z| = ρ ohne die Scheiben; Bögen |z − L| = r₀, |z − L̄| = r₀ in D̄); Σ₁ mit u_A = 60° und 80°, Σ₂ mit u_A = −90°,
  r ∈ [1e-7, 0.05], θ = u_A + 71.46° + 1.623°·ln(r/ρ), konjugiert bei L̄; Γ*_fern = Vs γ mit Keil-Linse c = 0.25, zum Vergleich Ks Linse.
- Ausgabe je n: |P_n|, R_P, R_A (getrennt Σ / fern, mit Lage), Q_n, lokale Exponenten.
- **Vorhersage:** Q_n = 1.008·|P_n| = 0.54·n^{−0.106} für alle n von ≈ 600 bis 2000; R_A(fern) fällt exponentiell; R_P fällt. **Entscheidung:** gestützt, wenn Q_n/|P_n|
  beschränkt bleibt; widerlegt für diese Kontur, wenn Q_n/|P_n| mit steigendem lokalem Exponenten wächst (Kandidat nur der Steg bei zu dünner Linse).
- Zusatz für den späteren Beweis: primale Rate h_P(z) = (1/n)·ln|δ_n(z)/δ_n(L)| auf ∂D für n = 500 / 1000 / 2000 — konvergiert sie gegen ein H_P ≤ 0?

## 29. Hauptsatz-Dokument, Rücknahme von K17, exponentielle Ebene mit abgeschnittener Ecke (Auftrag Leitung 10:05; 10:05–10:23; k50)

### 29.1 Warum K17 nicht trägt (Prüfung beim Schreiben von HAUPTSATZ_KR0.md)
K17 benutzte „Datenseite wie in K15“. Die Daten sind Pol-Kerne mit Quellen ζ = L(1+u) auf dem Schlitz (Cauchy-Darstellung von f_K): W_ζ = J[k_ζ]_{>N}, und
[k_ζ]_{>N}(f(t)) = (f(t)/ζ)^{N+1}/(ζ − f(t)) ist auf der Kontur nur klein, wo |f(t)| < |ζ|, nahe L also im FESTEN Sektor f⁻¹(D̄) = {θ ∈ (128.3°, 308.3°)}
(Richtung von t − L; f dreht um arg μ = 38.3°). Das adjungierte Fenster Σ₁ (u_A ≈ 60–95°) dreht mit 1.623°/e-Faktor aus diesem Sektor heraus: bei
u_A = 80° unterhalb r ≈ 8.5e-7, bei u_A = 95° unterhalb r ≈ 9e-11. Für Quellen mit u ≲ c·r ergibt das einen Faktor e^{N·c·r}, gegen den der Sprung u^{5.447} nur bis
N ≈ 1e8 (u_A = 80°) bzw. 1e11 (95°) gewinnt. **K17 ist damit für alle N nicht gezeigt** — dieselbe Drehung, jetzt auf der Datenseite. (§22.2 „Datenseite trägt für
ALLE n“ hatte die Pol-Quellen nahe L nicht geprüft; das gilt ebenso für K15.) (AH) bleibt als gemessene Eckaussage stehen.

### 29.2 Exponentielle Ebene: die Ecke lässt sich im festen Abstand umgehen (numerisch, k50, n = 300)
G_n(t) := (1/n)ln|Λ_n(t)/Λ_n(L)| + ln(|f(t)|/ρ) (Datengewicht der Taylor-Rest-Darstellung W = J(f − T_Nf), |W(t)| ≲ (|f(t)|/ρ)^N·Σ_{k>N}|b_k|):
- Bögen |t − L| = δ, δ ∈ {0.003, 0.01, 0.03, 0.1}: **G < 0 in den Richtungen 130° … 265°** (z. B. δ = 0.01: −0.0022 … −0.0045 … −0.0011), ≈ 0 bei 280°, > 0 ab 295°.
  Mechanismus (Skizze): die Daten fallen linear in δ, der Kopf wächst höchstens wie δ^{κ₁}, Re κ₁ = 1.118 > 1 — auf exponentieller Ebene gibt es KEINE Drehung,
  der gute Sektor ist asymptotisch genau f⁻¹(D) ∩ D nahe L.
- Γin-Seite: G ≈ −0.2 … −5 (Steg bei 180°: −4.99) — der Steg ist durch die Daten weit gedämpft, ohne Linse.
- **Γout-Seite in D ∩ f⁻¹(D̄): min_r G = +0.03 … +0.38 (θ = 288° … 72°)** — mit dem Taylor-Gewicht gibt es dort keine gute Kurve. Dort braucht es die Darstellung
  über den Rahmen K′ (Green-Funktion statt |f|/ρ), wie in (H-ii).
- Primal: max (1/n)ln|δ_n/δ_n(L)| auf ∂D ∖ U_{0.1} = −0.0031, ∖ U_{0.03} = −0.0007 (H_P ≤ 0 dort, n = 300).
Folgerung für Auftrag 2 der Leitung: Auf exponentieller Ebene ist die ECKE harmlos (fester guter Sektor); offen ist die Γout-Seite, also genau Ws Bucht-Problem —
mit dem richtigen Gewicht (Green-Funktion von K′ bzw. W̄).

### 29.3 HAUPTSATZ_KR0.md (phase3/)
Hauptsatz unter (H)^{r₀} = Ws/Ps (H) mit D̄ → K_{r₀} (Pol-Fehler nur für z ∈ K_{r₀}, Kugel in sup_{K_{r₀}}); (EB_N) wird in dieser Fassung nicht gebraucht; Kette
P-I + H1 + H2 (V) mit Umstellung auf K_{r₀} (B), Zertifikat bis 400 (S), Konstanten und fünf Prüfstellen für V.

## 30. Reicht in der K_{r₀}-Fassung die exponentielle Ebene? (Frage 1 der Leitung 10:25; 10:25–10:37; k51, k52) — global ja (mit strikter Reserve), an der Ecke nein

### 30.1 Die primale Reserve auf K_{r₀} ist POLYNOMIAL, nicht exponentiell (numerisch, k51, kubisch, n = 200 … 800)
| Punkt z | ln|δ_n(z)/δ_n(L)| (n = 200/300/400/600/800) | lokaler Exponent (Potenz) bzw. Rate |
|---|---|---|
| 0.9L, 0.7L, 0.5L | −7.32 … −11.57 / −11.68 … −15.92 / −14.05 … −18.30 | Potenz **−3.06** (n = 300 … 800 auf allen Stellen gleich) |
| r₀ = 0.1 bzw. 0.05, Richtungen 220° … 320° (innen) | z. B. r₀ = 0.1, 256°: −6.16 … −10.41 | Potenz −3.06 |
| −ρ (Steg-Seite) | −2.47 … −9.74 | Rate **−0.0120** (exponentiell, konstant) |
| r₀ = 0.05, 168° (Γin-Tangente) | −0.22 … −1.65 | Rate −0.0023 (exponentiell) |
Der Exponent ist explizit: **Re(a·s₀^P) = Re((1 + iα)(3.447 + 1.058i)) = 3.053** — das zweite Glied der angepassten Entwicklung der Eck-Pol-Antwort (Potenzglied
(z − L)^{−s₀^P}·d_n^{s₀^P} von F_P), es trägt die Ecke in ganz D hinein. Im Inneren von D ist δ_n(z)/δ_n(L) ≈ g(z)·φ(n) mit EINER Folge φ(n) ~ n^{−3.05}
(gleiche Raten an 0.5L, 0.7L, 0.9L auf 5 Stellen). Folge: die primale Seite schluckt e^{o(n)}-Fehler NICHT (Σ n^{−0.106−3.05}e^{o(n)} divergiert, sobald o(n) ≫ log n);
sie schluckt Potenzen bis zum Gesamtbudget unten.

### 30.2 Datenseite ohne abgeschnittene Pole — Satz K19 (bewiesen, elementar)
Für jedes q_N ∈ P_N gilt exakt Π_Nf − q_N = J_N⁻¹T_NJ(f − q_N). Sei 𝒞 eine feste geschlossene Kurve um 0 durch L, L̄ (Eckteil = Spiralen Σ₁, Σ₂ aus §28.2), E ein
Kompaktum ⊃ D̄ ∪ 𝒞 ∪ f(𝒞) ∪ f(int 𝒞), auf dem f_K stetig und im Innern holomorph ist (Spiral-Schnitt von L aus in der Lücke ≈ (25°, 124°), die E frei lässt), M = max_𝒞|m|.
Mit Lemma K13 (Kontur 𝒞 für alle n) und der Dreiecksungleichung:
  **sup_{K_{r₀}} |Π_N f_K − f_K| ≤ (1 + c_𝒞(1 + M)·Σ_{n≤N} Q_n)·dist_E(f_K, P_N),   Q_n := max_{z∈K_{r₀}, t∈𝒞} |ΔK_n(z, t)|.**
Die Daten sind gleichmässig klein auf E; es wird NIRGENDS |f(t)| < ρ verlangt ⇒ das Drehproblem von §29.1 (abgeschnittene Pole) entfällt, 𝒞 darf ausserhalb D liegen
(Linse über Γin, Γout-Bogen), ein Bucht-Weg in D ist nicht nötig. dist_E(f_K, P_N) ~ N^{−5.447·α_E}, α_Eπ = Aussenwinkel von E in L ≈ 99° ⇒ **s_E ≈ 3.0**
(Rate: Standard vom Dzjadyk-Typ, für die Spiralecke nicht nachgeprüft; für Konvergenz ohne Rate genügt Mergelyan). Budget: Σ_{n≤N}Q_n darf wie N^γ wachsen, γ < s_E.
**Einschränkung:** K19 ist die lineare Aussage bei P = f_K. Für die Familie (Sekanten) brauchen die Daten (f − q_N)∘P den Wert von P auf 𝒞 ∖ K_{r₀} (Linse, Ecke) —
den kontrolliert die K_{r₀}-Kugel nicht. Im Hauptsatz (HAUPTSATZ_KR0.md) bleibt deshalb (H)^{r₀} die Annahme.

### 30.3 Adjungierte auf der fernen Kontur: strikte exponentielle Reserve (numerisch, k52)
h_n(t) = (1/n)ln|Λ_n(t)/Λ_n(L)| auf |t| = ρ(1 + s): **Γout-Bogen (≥ 10° von L, L̄), s = 0: max h = −0.053 (n = 300), −0.037 (n = 600)**; Γin (Steg 180°): s = 0: +0.017/+0.015,
s = 0.0166: +0.0016/−0.0009, **s = 0.03: −0.011/−0.013**. Also: 𝒞 = Γout-Bogen ∪ Linse s ≈ 0.03 über Γin ∪ Spiralen in L, L̄ hat ausserhalb U_{r₁}(L, L̄) eine
strikte Reserve η ≈ 0.013; die Spiralen in U_{r₁}(L) liegen im Potenzsektor von F_A (max 1.008, k49).

### 30.4 Antwort auf Frage 1
In Q_n = |P_n|·R_P(n)·R_A(n) zerfallen R_P und R_A je in einen fernen Teil (ausserhalb einer festen Umgebung U_{r₁}(L, L̄)) und einen nahen Teil:
- **fern — exponentielle Ebene genügt:** adjungiert auf 𝒞 ∖ U_{r₁} (η ≈ 0.013, k52) ~~, primal auf ∂D ∖ U_{r₁} (Rate ≤ −0.0023 … −0.012, k51)~~ **[KORREKTUR §31.1: die primale Rate am Rand ist vorasymptotisch; asymptotisch H_P ≡ 0, die Primale braucht eine globale polynomiale Schranke]**; e^{o(n)} wird von e^{−ηn} geschluckt.
  Das ist genau Ws Ebene (n-te Wurzel), jetzt OHNE Bucht: Linse über Γin und der Γout-Bogen selbst.
- **nah — polynomial nötig, aber lokal und explizit:** (P) |P_n| ≤ Cn^{β₁} (gemessen n^{−0.106}, Gesetz §25); (A-lok) Adjungierte auf Σ_j ∩ U_{r₁}: ≤ Cn^{κ_A}
  (gemessen ≤ 1.008); (P-lok) Primale auf ∂K_{r₀} ∩ U_{r₁} (Bögen im Abstand r₀ von L): ≤ Cn^{κ_P} (gemessen fallend, Potenzsektor n^{−3.05}). Alles in Nicht-Wachstums-
  sektoren der expliziten Eckprofile (K10) — eine lokale Anpassungsaussage, keine globale.
- **Budget:** 1 + β₁ + κ_P + κ_A < s_E ≈ 3 (gemessen: 1 − 0.106 + 0 + 0 = 0.89).
**Also: Die exponentielle Ebene reicht für alles Globale; polynomiale Genauigkeit bleibt nur in einer festen Umgebung von L nötig — dort, wo die Profile explizit sind.**
Ws Sackgasse §26.3 ist damit auf die Eckumgebung eingeschränkt. Sie entstand, weil in Ws Rahmen die globale Aussage selbst polynomial sein musste (Ecke und Fernbereich
nicht getrennt; Daten als abgeschnittene Pole mit f(Kontur) ⊂ D̄).
Offene Glieder: (i) Satz für s_E an der Spiralecke; (ii) Familie (siehe 30.2); (iii) Beweis der lokalen polynomialen Aussagen (A-lok), (P-lok), (P) — das ist (ii)/(v)/(i″)
in lokaler Fassung, jetzt nur als SCHRANKE mit Budget ≈ 2, nicht mit relativem Fehler; (iv) globale exponentielle Aussagen aus Ws Potentialtheorie (Auftrag 2).

## 31. Fernraten H_A, H_P auf n-ter-Wurzel-Ebene (Auftrag 2 der Leitung 10:40; 10:40–10:48; k51b) — die Adjungierte hat eine Reserve, die Primale NICHT

### 31.1 Primale Fernrate: H_P ≡ 0 auf ganz D̄ (relativ zur Ecke) — Korrektur zu §30.4
k51b (kubisch, n = 200 … 1000), ln|δ_n(z)/δ_n(L)| auf Strahlen zum Rand:
| Strahl | 0.5ρ … 0.95ρ | 0.98ρ | ρ (Rand) |
|---|---|---|---|
| 120° | Potenz −3.06 (n ≥ 300) | Exponentialphase bis n ≈ 600, dann Potenz −3.1 | Rate −0.0105 (n ≤ 1000), Amplitude O(1) |
| 150° | Potenz ≈ −3 (verrauscht) | Übergang bei n ≈ 800 | Rate −0.0117 |
| 180° | Niveau e^{−20 … −24} (Interferenz der Beiträge von L und L̄ auf der reellen Achse) | Übergang bei n ≈ 900 | Rate −0.0120 |
| 30° (Γout-Seite) | Potenz −3.06 schon ab n = 200, Amplitude e^{−14 … −19} | wie innen | Potenz −3.06 |
**Befund:** Im Innern gilt δ_n(z)/δ_n(L) ≈ g(z)·n^{−a s₀^P} (|n^{−a s₀^P}| = n^{−3.053}) mit EINER Funktion g. Am Rand fern von L gibt es zusätzlich einen exponentiell fallenden
Anteil der Grösse O(1) mit Rate ≈ −0.011; der Potenzanteil hat dort nur die Amplitude e^{−14 … −20}, deshalb sieht man bis n = 1000 am Rand die Exponentialphase.
**Lemma K20 (bewiesen, eine Zeile):** Ist die n-te-Wurzel-Grenze u_P = lim sup (1/n)ln|δ_n/δ_n(L)| auf einer offenen Menge U ⊂ D gleich 0, so ist u_P ≥ 0 auf ∂U
(oberhalbstetige Regularisierung der subharmonischen Grenze). Mit dem Innengesetz (numerisch) folgt **H_P = 0 auf ganz D̄**: die Primale hat auf K_{r₀}
asymptotisch KEINE exponentielle Reserve; die Exponentialphase am Rand ist vorübergehend.
**Korrektur zu §30.4:** Dort stand „primal fern: Rate ≤ −0.0023 … −0.012 ⇒ e^{o(n)} geschluckt“ — das ist vorasymptotisch und für alle N falsch. Die Primale braucht auf
K_{r₀} eine POLYNOMIALE Schranke R_P(n) ≤ C(r₀)·n^{κ_P}, und zwar global (Budget κ_P < s_E − 0.894 − κ_A ≈ 2.1; wahr ist n^{−3.05}).

### 31.2 Was die primale Schranke ist (exakt) und warum keine Dualität sie auf die n-te Wurzel zurückführt
Aus der Rang-1-Ränderung: R_P(z) = |ΔK_n(z, L)|/|P_n| = |L|·|u_n(z) − u_{n−1}(z)|/|P_n| mit der Eck-Pol-Antwort u_n = J_n⁻¹T_nk_L = Π_n s (s = J_∞⁻¹k_L, Pol in L).
**(G-P) ⟺ die Galerkin-Lösung Π_n s des Bulk-Problems mit Eck-Pol konvergiert auf K_{r₀} mit polynomial beschränkten Zuwächsen.** Gemessen: Zuwächse ~ n^{−3.16}, also
Π_n s → s + O(n^{−2.16}) auf K_{r₀}. Die Herkunft ist explizit: das Potenzglied ζ^{−s₀^P} des inneren Profils V (§25) erzeugt aussen die Korrektur
c·d_n^{s₀^P−1}·ψ(z), ψ die Bulk-Lösung mit Ecksingularität (z − L)^{−s₀^P}. Das ist eine Aussage vom selben Typ wie der Hauptsatz (Konvergenz auf K_{r₀}), nur für
das eine Datum s; K19 hilft dafür nicht direkt (s hat in L ∈ E einen Pol, dist_E(s, P_N) = ∞).
Dualität auf n-ter-Wurzel-Ebene kann H_P nicht liefern, weil H_P auf D̄ keine Information trägt (≡ 0); die Information steckt in der Amplitude g und im Exponenten
Re(a s₀^P) — beides Eckgrössen.

### 31.3 Adjungierte Fernrate: strikte Reserve; sie ist eine GESTALT-Aussage — genau Ws W12′
Gemessen (k52): Γout-Bogen −0.037, Linse s = 0.03 −0.013 (n = 600, relativ zur Ecke). Theoretisch ist H_A = u_A − u_A(L) ein Unterschied von Raten. Satz W13 (V-geprüft)
sagt: das Grenzsystem aus W11 ist verschiebungsinvariant und legt das absolute Niveau nicht fest — für uns unschädlich, denn gebraucht wird nur die RELATIVE Rate,
also die Gestalt. Die gesuchte Aussage
  **(FA) lim sup_n max_{t ∈ 𝒞 ∖ U_{r₁}(L, L̄)} (1/n)·ln|Λ_n(t)/Λ_n(L)| ≤ −η(r₁) < 0**
folgt aus **W12′ (Gestalt des Grenzproblems bis auf Konstanten eindeutig; offen)** zusammen mit der Modellrechnung (Ws Schweissmodell, Reserve im Modell ≈ 0.013 auf der Linse,
0.037 auf dem Γout-Bogen). Nötig sind dazu nur W11 (Kompaktheit, tropische Bilanz; V) und die Eckstetigkeit u_A(t) → u_A(L) längs der Spiralen (lokal, Kegel).
Ohne W12′: ein Vergleichsprinzip (Collatz–Wielandt, Ws offener Ansatz (2)) mit einer Superlösung Φ, Φ(L) = 0, Φ < 0 auf 𝒞 ∖ U_{r₁}, würde (FA) direkt geben.
Neu gegenüber W: Die Kontur liegt jetzt auf dem Γout-Bogen und auf einer Linse über Γin — die Bucht in D wird nicht gebraucht (K19); die Reserve dort ist im Modell gross.

### 31.4 Gesamtbild nach Auftrag 2
Für die lineare Aussage bei f_K (K19) ist die Hypothese **(AH_{r₀}): Σ_{n≤N} Q_n ≤ C·N^γ mit γ < s_E, Q_n = max_{z ∈ K_{r₀}, t ∈ 𝒞}|ΔK_n(z, t)|** (§28 ohne z = L,
mit der Datenseite aus K19). Sie zerfällt in:
- Adjungierte fern: n-te-Wurzel-Ebene mit Reserve (FA) ⇐ W12′ (+ Modell) — **der einzige Teil, den Ws Potentialtheorie liefern kann**;
- Adjungierte nah (Spiralen in U_{r₁}): lokal, polynomial, explizites F_A (≤ 1.008);
- Eckprodukt P_n: lokal, polynomial (K16, Gesetz §25);
- **Primale auf K_{r₀}: global polynomial (G-P), = Konvergenz der Eck-Pol-Antwort Π_n s auf K_{r₀}** — keine n-te-Wurzel-Reduktion möglich (K20).
Die Familie (Sekanten) bleibt eigene offene Teilfrage (PLAN).

## 32. Zirkelprüfung (G-P) und die Hypothesenliste der linearen Aussage auf K_{r₀} (Auftrag Leitung 10:50; 10:50–10:52)
**Urteil (Begründung REFLEXION 10:52):** (G-P) ist der Eckteil der alten Pol-Antwort-Hypothese an der Quelle ζ = L, eingeschränkt auf z ∈ K_{r₀}, einen Operator und eine Quelle;
ein Beweis über die äussere Entwicklung reduziert sich (K16-artige Induktion, z-gleichmässig — neu und tragfähig) auf Lemma ρ und damit auf eine polynomiale Parametrix
der Eck-Pol-Antwort — den Schritt, an dem die Eckversuche gescheitert sind. Nicht direkt angreifen; als benannte Hypothese führen.
**Alle polynomialen Aussagen betreffen die beiden Eck-Pol-Antworten.** u_n = J_n⁻¹T_nk_L = Π_n s (primal, Quelle in L) und w_n(t) = K_n(L, t) (adjungiert, Auswertung in L):
- (P) Eckprodukt P_n = L·Δu_n(L) = Δw_n(L): |P_n| ≤ Cn^{β₁} (gemessen 0.535·n^{−0.106});
- (G-P) primal fern: sup_{K_{r₀}}|Δu_n| ≤ Cn^{κ_P − 0.106}·|L|^{−1} (gemessen ~ n^{−3.16});
- (A-lok) adjungiert nah: sup_{t ∈ Σ_j ∩ U_{r₁}}|Δw_n(t)| ≤ Cn^{κ_A}|P_n| (gemessen ≤ 1.008·|P_n|).
Zusammen: **(EP) Eck-Pol-Hypothese — die Zuwächse der beiden Eck-Pol-Antworten sind polynomial beschränkt: Δu_n auf K_{r₀} ∪ {L}, Δw_n auf den Spiralen in U_{r₁}(L)**,
mit Budget 1 + β₁ + κ_P + κ_A < s_E. Dazu die eine n-te-Wurzel-Aussage **(FA)** (adjungiert fern, ⇐ W12′ + Modell) und die klassische Approximationsrate s_E (V misst).
**Lineare Aussage (bedingt, Beweis = K19 + Rang-1-Faktorisierung):** (EP) + (FA) ⇒ sup_{K_{r₀}}|Π_Nf_K − f_K| ≤ C·N^{1+β₁+κ_P+κ_A}·dist_E(f_K, P_N) → 0.
Gegenüber (H-ii)^{r₀}: statt aller Quellen ζ auf Schlitz ∪ Γ_d nur die eine Quelle L (und ihr adjungiertes Gegenstück), statt der D̄- bzw. K_{r₀}-Pol-Antwort für alle ζ
nur die Zuwächse zweier Funktionenfolgen; der globale adjungierte Teil auf n-ter-Wurzel-Ebene. Familie: offen wie bisher (PLAN 10:48).
**Hebel für später:** duale, gewichtete Parametrix (Residuum klein nur nahe L und am Steg, exponentiell gross erlaubt, wo (FA) Reserve gibt).

## 33. (FA) über ein Vergleichsprinzip mit Superlösung — der Weg scheitert in dieser Form; genaue Stelle (Auftrag Leitung 10:55; 10:55–11:00)

### 33.1 Das Gerüst, das trägt (Skizze; Lemma K21 als Gerüst)
Für eine Grenzkonfiguration (u, w, p) (Raten von Y = Λ_n, T, P in Ws Identität, normiert u(L) = 0) und Vergleichsfunktionen Φ_y, Φ_T, Φ_P setzt man
M := max(sup(u − Φ_y), sup(w − Φ_T), sup(p − Φ_P)). Ist M > 0 und wird es an τ* von einer Komponente X angenommen, so gilt:
(a) Wo Φ_X superharmonisch ist, ist X − Φ_X subharmonisch ⇒ X − Φ_X ≡ M auf der ganzen Zusammenhangskomponente (starkes Maximumprinzip).
(b) Die tropische Bilanz (W11) am Punkt σ, an dem X als Glied auftritt: X ≤ max(übrige Glieder) ⇒ ein anderes Glied Z hat Z − Φ_Z ≥ M + (Φ_X − Φ_Z)(σ).
    Dominiert Φ_X im Vergleichssystem an σ streng, ist das ein Widerspruch; sind Φ_X und Φ_Z dort gleich, wandert der Überschuss M zu Z.
(c) An L: u(L) − Φ_y(L) = 0 < M und u − Φ_y ist oberhalbstetig ⇒ der Überschuss kann L nicht erreichen.
Läuft die Übertragung (a)–(b) durch die ganze Konfiguration bis an L, folgt M ≤ 0, also u − u(L) ≤ Φ_y und mit Modellreserve (FA). Die Normierung in L macht Ws
Verschiebungsinvarianz (W13) unschädlich — dieser Teil ist neu und richtig.

### 33.2 Hindernis O1 (bewiesen, eine Zeile): es gibt KEINE in der Ebene superharmonische Vergleichsfunktion
u = lim(1/n)ln|Λ_n| hat Δu = 2π(ν − δ₀) mit der Nullstellenverteilung ν (Wahrscheinlichkeitsmass): Gesamtmasse 0, wie für jede Funktion auf Ĉ. Eine auf Ĉ ∖ {0} superharmonische
Φ_y mit derselben Singularität −ln|τ| in 0 hätte ΔΦ_y = −ν′ − 2πδ₀ mit Gesamtmasse ≤ −2π < 0 — unmöglich. **Jede Vergleichsfunktion muss also positive
Laplace-Masse ≥ 2π tragen** (sie muss die Nullstellen der Köpfe „enthalten“), und dort versagt (a). Der Überschuss muss an diesen Stellen über die nichtlokale
Bilanz (b) weitergereicht werden.

### 33.3 Hindernis O2: die Weitergabe verlangt die Dominanzstruktur der unbekannten Konfiguration
(b) gibt nur „das Maximum wird zweimal angenommen“, nicht WELCHES Paar balanciert. Damit der Überschuss an einer Naht von Φ (positive Laplace-Masse) weitergeht, muss die
Konfiguration dort dasselbe Paar balancieren wie Φ. Auf Ws fester Vergleichsfläche ℳ₀ ist die Klebung über B (durch f) verträglich — V := max(u, w∘f) im y-Blatt ist
dieselbe Funktion wie max(w, u∘h) im T-Blatt und subharmonisch über B —, über die Linsen-Nulllinie λ (Klebung durch die Identität) aber nicht:
max(u, w∘f) ≠ max(w, u∘h). Subharmonie über λ bräuchte u ≥ w, w∘f in einem Übergangsstreifen, also Kenntnis der Dominanzgebiete — das ist genau Ws §22.3
(„Kleben auf einer festen Vergleichsfläche allein schliesst nicht“).
Wäre die Dominanzstruktur bekannt, gäbe 33.1 mit V auf ℳ₀ und der Liouville-Eigenschaft kompakter Flächen sogar V = H + const — also W12′. **Der Vergleichsweg
ist damit gleichwertig zu W12′ (Eindeutigkeit der Schweisslösung bis auf Konstanten) bzw. zur Kenntnis der Nullstellenverteilung der Köpfe**; er ist kein
leichterer Zugang.

### 33.4 Folgerung
(FA) bleibt eine benannte Hypothese auf n-ter-Wurzel-Ebene — „Gestalt“-Aussage, gemessen mit Reserve (η ≈ 0.013 auf der Linse s = 0.03, 0.037 auf dem Γout-Bogen),
beweisbar aus W12′ (offen). Ein numerischer Test eines konkreten Φ ist nach O1 sinnlos, solange Φ in der Ebene superharmonisch sein soll; auf ℳ₀ hängt er an O2.
Der eine neue Gedanke, der bleibt: die Normierung in L statt am Sockel (33.1(c)) — falls jemand W12′ angeht, ist die einseitige Fassung „jede in L normierte
Konfiguration liegt unter dem Modell“ das Ziel, nicht die volle Eindeutigkeit.

## 34. Die Familie in der K_{r₀}-Fassung: Keimabhängigkeit + K4 führen F2/F3 NICHT auf f_K zurück — die Lücke sitzt am Steg (Auftrag Leitung 11:06; 11:06–11:16; k51c)

### 34.1 Was Keimabhängigkeit und K4 leisten (bewiesen, §10.1, §9 Satz K4 — gilt wörtlich in K_{r₀})
- Die Matrizen T_N𝒥_P und J_N(P) hängen nur vom Keim P mod x^{N+2} ab; Werte von P ausserhalb von D gehen in sie nicht ein.
- K4: J(P) = T_N[C_P + M_{eˣ/P′}] + (Rang ≤ 2 im Kopf) für Fast-Lösungen — der Überlauf im Multiplikationskanal ist weg; der Rang-2-Rest ist F4.
Damit sind F1 und F4 Aussagen über Matrizen, die numerisch mit denen bei f_K übereinstimmen (‖Ĵ⁻¹‖ auf 7 Stellen gleich bis N = 1000, V). **Als Satz** brauchen
F1/F2 aber die Stabilität der KÖPFE unter der Keimstörung P − f, und der Kompositionskanal C_P − C_f verstärkt Störungen mit e^{0.47N} (Zeilen mit f(x) ∉ D̄;
inkohärent: T_N f_K ⇒ ‖Ĵ⁻¹‖ ~ e^{0.135N}, k11). Eine Störungsrechnung um f_K trägt also nicht; nötig ist die Kohärenz P − f = λ(Π^Pf − f) — das ist (Fam-koh), wie bisher.

### 34.2 F3 (P auf der Kontur 𝒞 ∖ K_{r₀}) zerfällt in zwei Teile
**(a) [KORREKTUR 12:33, §37.1: falsch — Fenster 1 ist primal ein WACHSTUMSfenster (k51d); F3 nahe L ist nicht lokal] ~~Nahe L (Spirale Σ₁ ausserhalb D, r ≤ r₁): lokal reduzierbar (Skizze).~~** Aus P − f = λ(Π^Pf − f) und der Kopfentwicklung von Satz A, ausgewertet in z ∈ Σ₁ ∩ U_{r₁}, folgt
|P − f| ≤ (Σ_n Q_n(z))·dist_E mit primalem Faktor R_P(z) auf Σ₁ — Fenster 1 ist auch primal kein Wachstumsfenster (§24.5) —, also F3 nahe L ⇐ **(EP⁺)**: (G-P) zusätzlich
für z auf Σ₁ ∩ U_{r₁} (lokal, Kegelprofil F_P). Geschlossen über ein Stetigkeitsargument längs λ ∈ [0, 1] wie bei (EB_N). Gemessen passt das: (EB_N) am echten p_N,
η_N/r_N = 5e-5 bei N = 2000 (V).
**(b) Fern auf der Linse über Γin (am Steg): echtes Hindernis (numerisch, k51c).** Die primalen Köpfe wachsen ausserhalb D exponentiell:
| Punkt | Rate (1/Δn)·ln|δ_n/δ_n(L)| (n = 300 → 600 → 1000) |
|---|---|
| ρ·e^{iθ}, θ = 120° / 150° / 180° | −0.0106 / −0.0117 / −0.0121 (vorübergehend, K20) |
| 1.01ρ·e^{iθ} | −0.0013 / −0.0023 / −0.0025 |
| **1.03ρ·e^{iθ}** (Linse s = 0.03) | **+0.0170 / +0.0164 / +0.0162** |
| 1.05ρ·e^{iθ} | +0.035 / +0.035 / +0.035 |
also ≈ ln(1 + s) − 0.012. Auf derselben Linse hat die Adjungierte −0.013 (k52). Für z auf der Linse ist Q_n(z) ~ e^{+0.016n}: **Π^P f − f, also P − f, ist auf der
Linse exponentiell gross — F3 für P selbst ist für grosse N falsch.** Nur der Keim zählt (10.2): gebraucht wird ein Repräsentant P̃ ≡ P mod x^{N+2} mit P̃(Linse) ⊂ D_{ρ(1−ε)}
(dafür genügt sogar eine O(1)-Abweichung, weil |f| am Steg ≈ 0.15ρ ist). Nach Carathéodory–Fejér kostet der beste Repräsentant e^{NF₀}·(Kopfabweichung ~N^{−5.8}),
F₀ ≈ 0.0067 … 0.0084 (harmonisches Mittel der Linsenhöhe) ⇒ trägt bis **N ≈ 6300 … 7800**; für alle N braucht es Kohärenz.

### 34.3 Warum man die Linse für die Familie nicht vermeiden kann (Zange; Skizze mit numerischen Gliedern)
Die Kopfpaarung ⟨ℓ_n, T_N W⟩ hat EINEN Wert; abschätzen kann man sie nur mit einer Datendarstellung und einer Kontur:
- **gleichmässige Daten** (K19, Approximation auf E): auf dem Γin-Kreis hat die Adjungierte keine Reserve (Steg +0.015 … +0.012, k52) ⇒ Kontur muss ausserhalb D (Linse) ⇒ P auf der Linse;
- **Taylor-Daten** (f − T_Nf, Gewicht (|f(t)|/ρ)^N): auf Γin riesige Reserve (G bis −5, k50), aber auf der Γout-Seite in D G ≥ +0.03 und ausserhalb D divergent ⇒ scheitert auf Γout;
- **Mischen** (Datenzerlegung, Knotenmass ν mit Interpolationsfehler e^{−N(U^ν(w) − U^ν(ζ))}): nahe L verlangt die Γin-Seite U^ν wachsend in Richtung 166.6°/204.9°, die Γout-Seite
  und der Schnitt U^ν nicht fallend in Richtung −14° … 25° bzw. ≤ in Richtung 76.6° — die Bedingungen an ∇U^ν(L) widersprechen sich in erster Ordnung (Richtungsintervalle
  (76.6°, 256.6°) und [−65°, 76°] disjunkt).
**Folge:** Eine Kontur in D̄ mit einer Datendarstellung gibt es nicht; die Familie braucht P (bzw. einen Keim-Repräsentanten) auf der Linse am Steg. Das ist (Fam-koh), jetzt
genau lokalisiert: **(Fam-koh)^{r₀}: für P ∈ 𝒮_N ∩ B_N^{r₀} gibt es P̃ ≡ P mod x^{N+2} mit P̃(Linse ∖ U_{r₁}(L, L̄)) ⊂ D_{ρ(1−ε)}.** (Mit dem Schweissbild: die Kohärenz müsste
genau die Hochfrequenzanteile ~0.95N am Steg auslöschen, die die primale Linsenrate +0.016 erzeugen.)

### 34.4 Ergebnis
- Familie in K_{r₀}: F1, F4 Matrix-Aussagen (num. wie bei f_K); F2 und der ferne Teil von F3 hängen an (Fam-koh)^{r₀}; der nahe Teil von F3 an (EP⁺) (lokal).
- Die K_{r₀}-Fassung verschiebt die Familienlücke NICHT, sie lokalisiert sie: am Steg (Γin, Winkel ≈ π), fern von L.
- Messvorschlag für V: (i) sup_{Linse}|p_N − f_K| bis 2000 (bestätigt das exponentielle Wachstum?); (ii) Carathéodory–Fejér-Repräsentant: min über g von sup_{Linse}|p_N + x^{N+2}g − f_K|
  (über die w-Koeffizienten von (p_N − f_K)∘k, k: 𝔻 → Int γ) — wächst er wie e^{NF₀}N^{−5.8} oder bleibt er polynomial (Kohärenz)?

## 35. (S_E) — die Approximationsrate auf E mit Spiralecke (Auftrag Leitung 11:17; 11:17–11:20; k53)

### 35.1 Aussage
Sei E wie in HAUPTSATZ_KR0 (Fassung 2) gebaut, mit den Spiralfenstern Σ₁ (u_A = u₁) und Σ₂ (u_A = u₂). Nahe L ist E der Spiralsektor vom Strahl Σ₁ gegen den Uhrzeigersinn
bis f(Σ₂). E muss enthalten: Σ₁, Σ₂, das Innere von 𝒞 und dessen f-Bild; f dreht um arg μ = 38.3°. Die Aussenöffnung ist
  **Θ = u₁ − u₂ − 38.3°,   α_E = Θ/π.**
Der Schnitt von f_K läuft als Log-Spirale in dieser Lücke. Bei L̄ gilt dasselbe konjugiert; sonst ist der Rand stückweise glatt, und f_K ist dort holomorph.
**Korrektur zur Definition von E in HAUPTSATZ_KR0 Fassung 2:** Statt E ⊃ D̄ gilt E ⊃ D̄ ∖ U_{r₀}(L, L̄). Der Spiralsektor dreht sich mit 1.623° je e-Faktor aus der festen Γout-Tangente
heraus; bei (80°, −90°) ist die Reserve von 27.8° bei r₁ unterhalb r ≈ 2e-9 aufgebraucht. Mit r₀ > r_* ~ 1e-6 schadet das nicht.
**Satz K22 (bedingt auf einen zitierten Satz, siehe 35.2):**
  **dist_E(f_K, P_N) ≤ C_E·N^{−s_E},   s_E = (Re s₀ − δ_T)·α_E/(1 + T²) ≈ 5.417·α_E,**
mit Re s₀ = 5.447 (führender singulärer Exponent von f_K in L, s₀ = 1 + 2πi/L), T = 0.02833 (Spiralsteigung) und δ_T = 1.058·T = 0.030 (Verlust durch die Blattdrehung).
| Fensterwahl (u₁, u₂) | Θ | s_E | Bemerkung |
|---|---|---|---|
| (60°, −90°) | 111.7° | 3.36 | Fenstermitten wie in k49 |
| (80°, −90°) | 131.7° | 3.96 | Standard in Fassung 2 |
| (95°, −110°) | 166.7° | 5.02 | nahe den Fensterrändern; die Adjungierte auf den Spiralen steigt dort auf 1.35 bzw. 3.84 (k49c) |
Budget von Satz A: γ ≈ 0.89 (gemessen) < s_E — mit Reserve ≥ 2.4 für jede Wahl oben.

### 35.2 Beweis (Glieder einzeln gekennzeichnet)
1. **Direkter Satz für Quasischeiben** (Dzjadyk-Typ; V. I. Belyi, Mat. Sb. 1977; Andrievskii–Belyi–Dzjadyk, *Conformal Invariants in Constructive Theory of Functions of Complex Variable*, 1995;
   Dzjadyk–Shevchuk, *Theory of Uniform Approximation of Functions by Polynomials*, 2008) — **nur zitiert, aus dem Gedächtnis, nicht nachgeschlagen.**
   - Annahmen: E ist eine Quasischeibe (Rand eine Quasikreislinie); f ∈ A(E), f^{(r)} hat auf E einen Stetigkeitsmodul ω (Dini-artig).
   - Aussage: Es gibt p_N ∈ P_N mit |f(z) − p_N(z)| ≤ C·ρ_{1/N}(z)^r·ω(ρ_{1/N}(z)) für z ∈ E. Dabei ist ρ_{1/N}(z) der Abstand von z zur Niveaulinie |Φ| = 1 + 1/N der äusseren Riemann-Abbildung Φ.
   - **[Erledigt 11:37, §36.1: gelesen in Andrievskii–Pritsker–Varga, Satz 2; keine Zusatzbedingung an ω.]**
2. **E ist eine Quasischeibe (B, Skizze).** Nahe L ist ∂E das Bild eines geradlinigen Winkels unter der quasikonformen Spiralabbildung z ↦ L + (z − L)|z − L|^{iT}. Ahlfors' Dreipunkt-Bedingung
   gilt, weil zwei Punkte im Abstand r von L durch einen Teilbogen vom Durchmesser ≲ r verbunden sind (Log-Spiralen haben endliche Länge ~ r·√(1 + 1/T²)). Sonst ist der Rand
   stückweise glatt mit Winkeln ∈ (0, 2π). Die globale Konstruktion von E (Jordan-Gebiet, keine Selbstberührung) ist als Geometrie-Prüfung offen (V-Prüfstelle 2 in HAUPTSATZ_KR0).
3. **Glattheitsklasse von f_K auf E (B, aus (K3)/p32a, dort V).** Nahe L gilt f_K = P(√L·W·R(q)), W = P⁻¹(x), q = W^{2πi/L} = W^{s₀−1}. Also
   f_K = L + μW + Σ_{j≥2}a_jW^j + c·W^{s₀}(1 + …).
   Auf einem Blatt des Spiralschnitts ist |W^{s₀}| = |W|^{5.447}e^{−1.058·arg W} mit arg W = θ₀ + T·ln|W| + O(1); das ergibt |W|^{5.447−δ_T}.
   Folglich f_K ∈ C^{5, ω}(E) mit ω(h) ≲ h^{0.417}; fern von L, L̄ ist f_K holomorph auf einer Umgebung von E (aus (F)).
4. **Abstand der Niveaulinie an der Spiralecke (B).** In der Koordinate ζ = ln(z − L) ist die äussere Lücke ein schräger Streifen der Breite Θ/√(1 + T²). Die Abbildung auf einen geraden
   Streifen ist linear, danach exp(π·/Breite). Also gilt |Φ(z)| − 1 ≍ |z − L|^{π√(1+T²)/Θ}, und ρ_{1/N}(L) ≍ N^{−α_E/√(1+T²)}.
   An glatten Randpunkten ist ρ_{1/N} ≍ 1/N.
5. **Zusammensetzen.** sup_E|f − p_N| ≤ C·max_z ρ_{1/N}(z)^5·ω(ρ_{1/N}(z)) ≍ C·ρ_{1/N}(L)^{5.417} = C·N^{−5.417·α_E/√(1+T²)}. Der Faktor √(1+T²) = 1.0004 ist vernachlässigbar. ∎
   (Für Konvergenz ohne Rate genügt Mergelyan; Satz A braucht die Rate.)

### 35.3 Numerische Prüfung der Rate (k53; unabhängig von f_K, gerade Ecke statt Spirale)
g(z) = (z − L)^{s₀} mit Schnitt in der Lücke; E mit einspringender Ecke der Aussenöffnung Θ; nahezu beste Approximation per kleinsten Quadraten auf dichten Randpunkten
(Vandermonde mit Arnoldi); Fehler auf feinerem Gitter:
| Θ | vorhergesagt 5.417·Θ/π | lokaler Exponent N = 40 → 60 → 80 |
|---|---|---|
| 111.7° | 3.36 | 3.38 → **3.37** |
| 131.7° | 3.96 | 3.98 → **3.97** |
| 166.7° | 5.02 | 5.03 → **5.02** |
Die Dzjadyk-Rate trifft auf 0.01–0.02 (für gerade Ecken ist der Exponent 5.447·α, also 3.38 / 3.99 / 5.04 — ebenfalls innerhalb der Messgenauigkeit).
**Für V:** Die Messung von s_E am echten f_K auf E sollte 3.96 (Fenster 80°/−90°) bzw. 3.36 (60°/−90°) ergeben. Eine Abweichung hiesse: E ist anders gebaut, als gedacht, oder
f_K hat auf E eine weitere Singularität.

### 35.4 Folgerungen
- (S_E) ist kein offener Punkt mehr im Sinn eines Forschungsproblems. Offen sind der zitierte Satz (nachzuschlagen; Bedingungen an ω), die globale Geometrie von E und die
  numerischen (F)-Eingaben (Univalenz). Für HAUPTSATZ_KR0 Fassung 2: (S_E) „bedingt auf zitierten Satz + Geometrie von E“, s_E = 5.417·α_E.
- Die Fensterwahl ist ein Abwägen: breitere Lücke (u₁ gross, u₂ klein) ⇒ grösseres s_E, aber grössere Konstanten in (A-lok) und (FA). Mit (80°, −90°) ist die Budgetreserve s_E − γ ≈ 3.1.

## 36. (S_E): der zitierte Satz im Wortlaut, und die globale Geometrie von E (Auftrag Leitung 11:22; 11:22–11:37; k54)

### 36.1 Der direkte Satz — jetzt GELESEN (Wortlaut)
Quelle: V. V. Andrievskii, I. E. Pritsker, R. S. Varga, *Simultaneous approximation and interpolation of functions on continua in the complex plane*, arXiv math/0005032v1 (2000),
J. Math. Pures Appl. (2001). **Gelesen:** S. 1–5 (Einleitung, Definitionen, Sätze 1 und 2).
- **Mengenklasse.** E ⊂ ℂ ist ein Kontinuum mit diam E > 0 und zusammenhängendem Komplement Ω.
  - **E ∈ H**, wenn je zwei Punkte z, ζ ∈ E durch einen Bogen in E der Länge ≤ c|z − ζ| verbunden werden können ((1.2)).
  - **E ∈ H***, wenn zusätzlich die Moduln der Vierecke Q_Z(r, R) ⊂ Ω für benachbarte Primenden um höchstens c₁ voneinander abweichen ((1.3)).
  - Wörtlich: „H* includes domains with quasiconformal boundaries“ (dort mit Verweis auf ihre [1], [20]) sowie Dzjadyks Klassen B_k*.
- **Glattheitsmass.** ω_{f,k,E}(δ) := sup_{z∈E} E_{k−1}(f, E ∩ D(z, δ)), der k-te Stetigkeitsmodul nach Dyn'kin (lokale Bestapproximation durch Polynome vom Grad k − 1).
  Für k = 1 ist das die gewöhnliche Schwankung. A^r(E) ist die Menge der f ∈ A(E), die auf E r-mal stetig differenzierbar sind.
- **Abstand zur Niveaulinie.** ρ_δ(z) := dist(z, L_δ) mit L_δ = {|Φ| = 1 + δ}; Φ: Ω → {|w| > 1} konform, Φ(∞) = ∞, Φ′(∞) > 0.
- **Satz 2 (APV).** Sei E ∈ H*, f ∈ A^r(E), r, k ∈ ℕ. Dann gibt es zu jedem n ≥ k ein p_n ∈ P_n mit
  **|f^{(l)}(z) − p_n^{(l)}(z)| ≤ c·ρ_{1/n}(z)^{r−l}·ω_{f^{(r)},k,E}(ρ_{1/n}(z))  (z ∈ ∂E, l = 0, …, r)**,
  mit c unabhängig von n. Interpolationsknoten auf ∂E sind dort zugelassen; wir brauchen keine.
- **Bedingung an ω: keine weitere.** Die Dini-artige Bedingung (1.8) wird nur für die exponentielle Rate im Innern (1.9) gebraucht. Meine Annahme „Dini-artig“ in §35.2 war
  unnötig und ist gestrichen.
- APV schreiben, Satz 1/2 verallgemeinere für k > 1 die direkten Sätze von Belyi und Tamrazov (Quasischeiben) und Shevchuk (Klasse B_k*). Belyi 1977 und Belyi–Tamrazov 1978:
  **sekundär** (über APV). Andrievskii–Belyi–Dzjadyk 1995: **nur zitiert**.
**Anwendung (ersetzt 35.2 Glied 1):** r = 5, k = 1, l = 0. f_K ∈ A⁵(E) mit ω_{f^{(5)},1,E}(δ) ≲ δ^{0.417}. Nach dem Maximumprinzip (f − p_n ∈ A(E)) gilt die Randschranke auf ganz E,
also sup_E|f_K − p_n| ≲ max_{z∈∂E} ρ_{1/n}(z)^{5.417} ≍ n^{−5.417·α_E}. **Satz K22 ruht jetzt auf einem gelesenen Satz**, bedingt nur noch auf E ∈ H* (Geometrie, 36.2)
und auf die Glattheitsklasse von f_K auf E (aus (K3)).

### 36.2 Globale Geometrie von E (kubisches Modell als Vorlage, k54; für f_K dieselbe Prüfung mit f_K auf 𝒞 → V)
**Konstruktion.**
- 𝒞 = Σ₁ (u_A = 80°, r ≤ r₁ = 0.05) → Verbindungsstrecke → Linse |t| = 1.03ρ über Γin → konjugiert → Σ₂ (u_A = −90°) → Verbindungsstrecke → Γout-Bogen |t| = ρ.
- E := Abschluss von int 𝒞 ∪ f(int 𝒞) ∪ (D ∖ U_{r₀}(L, L̄)), mit r₀ = 0.02.
- Der Rand von E besteht nahe L aus Σ₁ und f(Σ₂), weiter aussen aus der Linse und f(Γout-Bogen).

**Ergebnisse (kubisch):**
| Prüfung | Ergebnis |
|---|---|
| (1) Γin = {\|f\| < ρ} auf \|x\| = ρ | Wechsel bei 76.6° und 283.3° = arg L, arg L̄ ✓; max\|f\|/ρ auf Γout = 1.62 |
| (2) f univalent auf int 𝒞 | Windungszahl von f(𝒞) um f(z) = 1.000 für 600 Stichproben ✓; min\|f′\| auf int 𝒞 = 0.30 |
| (3) E zusammenhängend, ohne Löcher | Raster 441²: eine Komponente, Komplement zusammenhängend ✓ |
| (4) Spiralschnitt in der Lücke meidet E | ja, für r ∈ [1e-6, 1] ✓ |
| (6) Aussenöffnung bei L | 132.3° / 132.1° / 132.0° / 131.9° / 131.9° bei r = 3e-2 … 1e-6 (Vorhersage Θ = 131.7°); die Lücke dreht von (13°, 145°) nach (357°, 129°) — Spirale ✓ |
| (5) Dreipunkt-Bedingung auf dem Randbogen durch L (Σ₁ ∪ f(Σ₂ ∪ V₂ ∪ Γout)) | max diam(Teilbogen)/Sehne = 1.07 ✓ |

Damit ist E (kubisch) eine Quasischeibe mit Spiralecke der Öffnung Θ ≈ 132°. Damit gilt E ∈ H*; die Bogen-Sehnen-Bedingung (1.2) folgt für Quasischeiben aus der linearen Verbundenheit.
**Für f_K (Aufgabe V, Kneser-Motor):** dasselbe Skript k_kern/k54_geometrie_E.py mit f_K statt der Kubik, genauer f_K auf 𝒞. Zu prüfen:
- f_K auf der Linse s = 0.03 und auf Σ₁ ausserhalb D: f_K univalent auf int 𝒞, und f_K′ ≠ 0 (wegen m = eˣ/f_K′);
- f_K(int 𝒞) ∪ int 𝒞 ohne Löcher;
- f_Ks echter Schnitt (Zweig k₀ = 2 nach V §3.15) liegt in der Lücke und meidet E;
- die Aussenöffnung ist ≈ 131.7°;
- Dreipunkt-Verhältnis beschränkt.
Erwartung: wie kubisch, denn die Ecke hängt nur an μ. Unterschiede gibt es nur in der Grösse von f_K(Γout): max|f_K|/ρ auf Γout ist bei Kneser kleiner (A-Bogen ⊂ f_K(D)).

### 36.3 Status von (S_E) nach diesem Block
**Satz K22: dist_E(f_K, P_N) ≤ C_E·N^{−5.417·α_E}.**
- APV Satz 2: **gelesen**.
- E ∈ H*: kubisch numerisch geprüft; für f_K steht die Prüfung bei V aus.
- f_K ∈ A⁵(E) mit Hölder-Exponent 0.417: **B**, aus (K3), dort V.
- Rate an der synthetischen Ecke: num, k53.
Offen ist damit nur noch die f_K-Geometrie von E — ein endlicher numerischer Check, keine Analysis.

## 37. Vortest (b): Lemma ρ für explizite Parametrizes der Eck-Pol-Antwort — trägt nicht; dazu eine Korrektur zu Fenster 1 (Auftrag Leitung 11:40; 11:42–12:33; k55, k55b, k55c, k51d)

### 37.1 Korrektur: Fenster 1 (Σ₁) ist ein PRIMALES Wachstumsfenster (numerisch, k51d; Rechnung)
F_P wächst für arg u ∈ (100°, 260°) und folgt im Sektor (−90°, 95°) der Potenz (Tabelle: log10|F_P| bei R = 50 von −5.5 bis +23.5). Mit u = λ̃_Pζ,
arg λ̃_P = 89.38°, ζ = (x − L)/d_n und arg d_n = −1.451°·ln n wächst der primale Kopf in den physikalischen Richtungen **θ ∈ (10.6° − 1.451°ln n, 170.6° − 1.451°ln n)**.
Das ist praktisch die ganze Aussenseite von D nahe L. Fenster 1 der Adjungierten (u_A ∈ (19°, 100°), also θ ≈ (81°, 162°) bei n = 600) liegt darin.
k51d, ln|δ_n/δ_n(L)| längs der Spiralen (n = 200 / 400 / 600):
| Spirale | r = 0.003 | 0.01 | 0.03 | 0.05 |
|---|---|---|---|---|
| Σ₁, u_A = 80° | +0.12 / +0.21 / +0.30 | +0.37 / +0.66 / +0.92 | +1.04 / +1.86 / +2.68 | **+1.70 / +3.17 / +4.65 (Rate +0.0074)** |
| Σ₁, u_A = 60° | … | … | … | +3.18 / +6.15 / +9.12 (+0.0149) |
| Σ₂, u_A = −90° | −0.08 / … | −0.24 / −0.45 / −0.66 | −0.74 / −1.61 / −2.82 | −1.37 / −3.77 / −5.08 |
| Σ₂, u_A = −110° | … | … | … | −2.41 / −4.57 / −5.93 |
**Folge:** Die Behauptung in §24.5 und §34.2(a), Fenster 1 sei „auch primal kein Wachstumsfenster“, ist falsch. Damit trägt (EP⁺) nicht: Die Primale wächst längs Σ₁
exponentiell, und F3 nahe L ist KEINE lokale polynomiale Aussage. Das passt zu Vs E2: p_N − f_K ist bei u_A = 85°, r = 0.05 und N = 2000 noch 2e-10 klein,
denn e^{0.0074N}·N^{−5} ~ 6e-11; für N ≳ 5000 wächst es.
Für Satz A (bei f_K) ändert sich nichts. Dort wird die Primale nur auf K_{r₀} ⊂ D̄ gebraucht, und die Spiralen tragen nur die Adjungierte.
**Familie, neu zusammengefasst:** Auf der Γin-Seite ausserhalb D (Ecke UND Steg) wächst die Primale. Dort braucht die Adjungierte aber die Kontur, denn nur
ausserhalb D hat sie dort Reserve. Diese eine Spannung ist die Familienlücke. (Fam-koh)^{r₀} gilt also für Linse UND Σ₁.

### 37.2 Lemma ρ für explizite Parametrizes (k55, k55b; kubisch)
ρ_{m,n} = (L⁻¹r̃)_m/(L⁻¹b)_m mit r̃ = b − J̃_nũ; LU ohne Pivot in hoher Genauigkeit. Kontrolle: Für die exakte Lösung ist ρ ≈ 1e-157 … 1e-193.
| Parametrix | n | max_m\|ρ_{m,n}\| (Lage m/n) | Profil \|ρ\| bei m/n = 0.1 / 0.5 / 0.9 | Fehler in L |
|---|---|---|---|---|
| P1: T_n s (nur aussen, Taylor) | 60 | 2.8e6 (0.80) | 1e-17 / 189 / 6e5 | 60 % |
| P2: T_n s_reg + T_n[(C/d_n)V] (Taylor) | 60 | 8.5e4 (0.80) | 4e-3 / 5.3 / 2e4 | 3.1 % |
| P3: kleinste Quadrate auf ∂E, ungewichtet | 120 | 1.44 (0.92) | 0.63 / 0.57 / 1.2 | 35 % |
| P3: dual gewichtet (Linse e^{−0.013n}, Γout e^{−0.037n}) | 120 | 6.5 (0.02) | 3.9 / 1.8 / 0.6 | 6.3 % |
Lemma ρ verlangt ρ → 0. **Keine der expliziten Parametrizes erfüllt es.**
- P1/P2: Der Taylor-Abschnitt erzeugt den Kompositionsüberlauf. Der Schwanz [s]_{>n} wird an f(Γout) ausserhalb D ausgewertet, Wachstum ~ (|f|/ρ)^n, max|f|/ρ = 1.62 (kubisch).
  Das ist dieselbe Zange wie in §34.3, nur auf der Parametrix-Seite.
- P3: Der Approximant auf E kann der Zusammensetzung nicht folgen; der relative Fit-Fehler auf ∂E ist 0.9.
  Zwei technische Fehler sind behoben: Koeffizienten von V mit zu wenig Stellen, und der Zweig von s_reg auf dem Kreis.

### 37.3 Warum: im primalen Wachstumssektor folgt die echte Antwort NICHT dem inneren Profil (k55c)
Vergleich der exakten Eck-Pol-Antwort u_n mit C·V(ζ)/d_n an Punkten der Kontur bzw. ihres f-Bilds (|Verhältnis|, n = 120 / 300):
| Ort | \|ζ\| ≈ 1 | ≈ 2.5 | ≈ 5–6 | ≈ 12 | ≈ 29 |
|---|---|---|---|---|---|
| Σ₂ (Fenster 2, in D, primal Potenzsektor) | 1.007 | 0.993 / 0.998 | 0.958 / 0.978 | 0.962 (n = 300) | — |
| Σ₁ (Fenster 1, ausserhalb D, primal **Wachstum**) | 1.018 | 1.042 / 1.015 | 1.223 / 1.101 | 1.397 (n = 300) | — |
| f(Γout) nahe L (ausserhalb D, Wachstum) | 1.018 | 1.059 / 1.021 | 1.386 / 1.184 | 3.66 / 2.00 | **22.7** (n = 300) |
- **Im Potenzsektor (Fenster 2)** stimmt das innere Profil auf ≤ 4 % bis |ζ| ≈ 12 und wird mit n besser. Die Anpassung innen/aussen ist dort harmlos.
- **Im primalen Wachstumssektor** gilt das Profil nur bis |ζ| ≈ 2–3. Danach wächst die echte Antwort schneller (Faktor 22 bei |ζ| ≈ 29), getragen von der
  äusseren (WKB-)Lösung mit global bestimmter Amplitude.
- Die Kontur muss aber gerade dort Daten haben: auf der Γin-Seite ausserhalb D (Σ₁) und in f(𝒞), denn nur dort hat die Adjungierte Reserve. Eine explizite Parametrix
  „innen V, aussen s“ ist dort falsch.

### 37.4 Ergebnis und nächster Hebel
**Lemma ρ trägt mit expliziten Parametrizes nicht.** Der Grund ist strukturell und nicht numerisch: Ausserhalb D nahe L (primaler Wachstumssektor) wird die Eck-Pol-Antwort
von der äusseren WKB-Wachstumslösung bestimmt. Deren Amplitude ist global (Anpassung im Wachstumssektor). Genau dort braucht die Paarungskontur Daten.
Die Lücke ist dieselbe wie bei der Familie (§37.1) und am Steg (§34). **Alle drei offenen Stellen sind eine einzige: das Verhalten der primalen Grössen ausserhalb D auf
der Γin-Seite, wo die Adjungierte ihre Reserve hat.**
**Nächster Hebel (Vorschlag):**
- (i) Die primale Antwort ausserhalb D nicht approximieren, sondern über die exakte Identität J u_n = T_nk_L + R_n aus ihren Werten in D̄ bestimmen:
  u_n(f(x)) = k_L(x) − m(x)u_n(x) + R_n(x) für x ∈ D̄, mit dem Abschneide-Rest R_n = (I − T_n)Ju_n (Grade > n).
  Damit wird „primal ausserhalb D“ zu „R_n auf dem Γout-Bogen“. Das ist die primale Entsprechung von Ws T und genau die Grösse, die Ws Schweissbild beschreibt (T ≈ y∘h).
  Prüfbar: Ist R_n auf Γout polynomial (dann trägt eine Parametrix „in D explizit + aussen über die Identität“) oder exponentiell (dann ist es W12′-Gebiet)?
- (ii) Andernfalls bleibt (EP) benannte Hypothese. Die Leitung entscheidet, ob K (i) prüft oder der Strang auf Dokumentation geht.


## 38. Hebel (i): der Abschneide-Rest R_n = (I − T_n)Ju_n auf dem Γout-Bogen ist EXPONENTIELL — (i) trägt nicht (Auftrag Leitung 12:34; 12:35–12:45; k56, k56b)

### 38.1 Vorhersage (vor der Messung, REFLEXION 12:35)
Galerkin gibt T_nJu_n = T_nk_L, also exakt **R_n(x) = u_n(f(x)) + m(x)u_n(x) − T_nk_L(x)** [bewiesen, eine Zeile]. Auf dem Γout-Bogen fern von L sind m·u_n und T_nk_L
beschränkt, also R_n ≈ u_n(f(x)) mit f(x) ausserhalb D. Vorhersage: exponentiell, Rate ≈ ln(|f(x)|/ρ) (Bernstein-Walsh), Maximum ≈ 0.485 bei φ = 0, gegen 0 nahe L und L̄.
Damit wäre die Identität nur eine Umschreibung derselben primalen Grösse ausserhalb D.

### 38.2 Messung (k56: u_n = J̃_n⁻¹ b exakt wie k55c, kubisch, n = 100/200/300/400; k56b: Raten und Koeffizientenprofil)
| φ (y = ρe^{iφ}) | \|y − L\| | \|f(y)\|/ρ | \|R_n\| bei n = 100 / 200 / 300 / 400 | Rate je Grad (300→400) | BW = ln(\|f\|/ρ) | Rate/BW |
|---|---|---|---|---|---|---|
| 76.37° | 0.006 | 1.003 | 24 / 58 / 97 / 141 | +0.0037 (fällt mit n: polynomial) | 0.0032 | — |
| 75.62° | 0.024 | 1.013 | 28 / 77 / 145 / 238 | +0.0049 (fällt) | 0.0126 | — |
| 74.62° | 0.048 | 1.025 | 35 / 119 / 305 / 786 | +0.0095 | 0.0251 | 0.38 |
| 71.62° | 0.120 | 1.063 | 78 / 1.0e3 / 1.6e4 / 2.4e5 | +0.0274 | 0.0614 | 0.45 |
| 66.62° | 0.240 | 1.126 | 5.4e2 / 1.3e5 / 3.3e7 / 8.7e9 | +0.0557 | 0.1187 | 0.47 |
| 36.62° | 0.940 | 1.455 | 8.4e7 / 1.5e16 / 2.9e24 / 5.8e32 | +0.1911 | 0.3747 | 0.51 |
| 0° | 1.704 | 1.624 | 2.0e10 / 1.6e21 / 1.4e32 / 1.3e43 | **+0.2521** | 0.4850 | 0.52 |
| −40° | 2.339 | 1.424 | 1.0e7 / 5.7e14 / 3.4e22 / 2.1e30 | +0.1794 | 0.3533 | 0.51 |
| −70° | 2.633 | 1.084 | 6.2 / 2.3e2 / 8.6e3 / 3.4e5 | +0.0367 | 0.0804 | 0.46 |
| −76.12° (nahe L̄) | 2.672 | 1.006 | 0.39 / 0.68 / 0.66 / 0.30 | beschränkt | 0.0063 | — |
- Die Raten sind über n auf 3 Stellen konstant (φ = 0: 0.2513 / 0.2519 / 0.2521). Genauigkeit: Lösung in 920 Bit, Auslöschung bei u_n(f) ≤ 41 Stellen, Kondition grob ≲ 1e105 (aus den Raten 0.47 und 0.135) ⇒ ≥ 170 Stellen Reserve.
- In allen Zeilen ist |R_n| ≈ |u_n(f(y))|. Der Zielwert fern von L, k_L(y) − m(y)u_n(y) ≈ s(f(y)), ist O(1) (0.10 … 4.3). **u_n weicht also an f(Γout) exponentiell von der
  echten Fortsetzung s ab.**
- Nahe L: lokale Rate ≈ 0.23·|y − L| (0.0274/0.120, 0.0557/0.240). Polynomial ist R_n nur auf der Eckskala |y − L| ≲ C·(ln n)/n. Für festes r₀ ist der ganze Bogen
  ausserhalb U_{r₀}(L, L̄) exponentiell (Rate > 0, am Rand von U_{r₀}(L) ≈ 0.23·r₀).
- Koeffizientenprofil (k56b): (1/n)ln|ũ_k| ≈ −0.003 … −0.017 bis k = 0.9n, nur der letzte fällt (−0.10). Die Koeffizienten sind also flach. Die halbe Rate kommt aus
  Phasen-Auslöschung, nicht aus abfallenden Koeffizienten.

### 38.3 Vergleich mit der Vorhersage
- **Richtig:** exponentiell; Maximum bei φ = 0; gegen 0 nahe L und L̄; |R_n| ≈ |u_n∘f|; die Identität ist eine Umschreibung (tautologisch).
- **Falsch um den Faktor 2:** Die Rate ist nicht ln(|f|/ρ), sondern **≈ ½·ln(|f|/ρ)** (Verhältnis 0.45–0.52 im Inneren des Bogens). Bernstein-Walsh ist also nicht scharf.
  Der Grund ist offen (Vermutung: die Galerkin-Lösung ist auf einer Menge grösser als D̄ polynomial beschränkt, deren Green-Funktion auf f(Γout) etwa die Hälfte von ln(|w|/ρ) ist).

### 38.4 Ergebnis
**(i) trägt nicht.** R_n ist auf dem Γout-Bogen exponentiell (bis +0.252 je Grad); auf demselben Bogen hat die Adjungierte nur −0.053/−0.037 (§36, n = 300/600).
Die Parametrix „D̄ explizit + Identität“ wird nicht gebaut: Sie müsste genau das R_n tragen, das exponentiell ist. **(EP) bleibt benannte Hypothese** (HAUPTSATZ_KR0 §2.1).
Der gemeinsame Kern in seiner schärfsten Form steht als Nachtrag in phase3/KERNPROBLEM.md §4.


## 39. Produkt primal × adjungiert am Steg auf n-ter-Wurzel-Ebene — es gibt dort KEINE Lücke; das Produkt ist explizit (Christoffel-Darboux-artig). Korrektur zu (KK) (Auftrag Leitung 12:52; Zeitfenster 12:54–13:06; k57, k57b, k57c)

### 39.1 Unterschied zu H-A und Vorhersage (REFLEXION 12:54)
H-A (§16) prüfte an der ECKE, ob beide Köpfe die homogenen Eckgleichungen erfüllen. Das scheiterte an den komplementären Sektoren. Am Steg gibt es keine Ecke.
Gefragt war nur die n-te Wurzel des Produkts Π_n(t) := δ_n(t)Λ_n(t)/(δ_n(L)Λ_n(L)) = ΔK_n(t,t)/P_n, also die Summe h_P + h_A der Raten relativ zur Ecke.
Vorhersage: radial konstant, positiv ≈ +0.003 ± 0.002 bei 180°; die Primale am Steg quellbeherrscht (A_P := m·δ_n·(ρ/t)^n = O(1)).
Abbruchkriterium: trägt, wenn h_Π im Fenster 450 → 600 < +0.001 und fallend ist (oder der Überschuss nur in ℓ_nn sitzt); scheitert, wenn ≥ +0.002 und ortsabhängig.

### 39.2 Messung k57 (kubisch, Köpfe exakt, n = 150/300/450/600; Prüfung der Zerlegung Π_n = A_PA_A/(m²ℓ_nnP_n) auf ≤ 1e-188)
Rate h_Π je Grad (Fenster 150→300 / 300→450 / 450→600):
| Ort | s = −0.03 | s = −0.01 … +0.03 | s = +0.05 |
|---|---|---|---|
| Γin-Seite, θ = 100° / 120° / 150° / 180° | +0.0005 / −0.003 … +0.008 / **+0.028 … +0.033** (wächst: Primale innen, K20) | **+0.0005 / +0.0003 / +0.0002 (alle θ, alle s)** | +0.0005 / +0.0003 / +0.0002 |
| θ = 85° (8° neben L) | +0.0004 / +0.009 / +0.026 | +0.0005 / +0.0003 / +0.0002 | +0.0004 / +0.0018 / **+0.0147** (Eck-Wachstumszone) |
| Γout-Seite, θ = 0° / 40° (Fenster 450 → 600) | −0.015 … −0.017 | −0.016 … −0.017 | −0.017 |
Skalare Raten (450 → 600): ln|δ_n(L)| −0.0646, ln|Λ_n(L)| −0.0958, ln|ℓ_nn| **−0.16024** (auf 5 Stellen konstant), ln|P_n| −0.0002.
- **Auf der Γin-Seite ist das Produkt im ganzen Ring s ∈ [−0.01, 0.05] polynomial.** h_Π fällt über die Fenster (0.0005 → 0.0002) und ist gleich der Rate des
  Kreismittels (nächster Punkt). Auf der Γout-Seite fällt es (−0.017).
- Gewachsen wird nur (a) innen (s = −0.03), weil dort die Primale in ihr Innengesetz übergeht (K20) und die Adjungierte nach innen wächst; (b) in einer Zone nahe L
  ausserhalb D ab s ≈ 0.05 (θ = 85°) — das ist der Eck-Wachstumssektor (§37.3).
- **Meine Vorhersage war in Teil 1 falsch** (+0.003 statt 0) **und in Teil 2 falsch**: Die Köpfe sind NICHT quellbeherrscht. A_P ~ e^{−0.077n}, A_A ~ e^{−0.084n}, und
  ℓ_nn ~ e^{−0.1602n} gleicht beide genau aus.
- **Woher die +0.003 in KERNPROBLEM §4 kamen:** Ich habe eine STEIGUNG der Primalen (300 → 600, k51c) mit einem n-te-Wurzel-NIVEAU der Adjungierten
  ((1/n)ln|Λ_n/Λ_n(L)| bei n = 300/600, k52) addiert. Gleichartig gemessen (Steigung 450 → 600) ist am Steg bei s = 0.03 die Primale +0.0162 und die Adjungierte −0.0160 (Vs v58: −0.0163).
  Die Summe ist +0.0002. **Die Aussage „Summe ≈ +0.003, Lücke unabhängig vom Radius“ ist falsch und zurückgezogen.**

### 39.3 Das Produkt ist explizit (k57b, k57c; numerisch)
**(a) Diagonal-Zuwachs.** V_n(t) := ΔK_n(t,t)·m(t) hängt auf der Γin-Seite NICHT von n und kaum vom Radius ab:
| θ | 90° | 100° | 120° | 150° | 180° |
|---|---|---|---|---|---|
| V_n(t), s = 0, n = 150 … 600 | 0.8670 + 0.0144i | 0.9034 + 0.0083i | 0.9351 + 0.0036i | 0.9522 + 0.0012i | 0.9564 |
| Streuung über n und s ∈ [0, 0.03] | ≤ 0.001 | ≤ 0.0005 | ≤ 0.0002 | ≤ 0.0001 | ≤ 0.0001 |
| κ(t) = (1/n)·∂ ln\|δ_n\|/∂ ln\|t\| (n = 600) | — | 0.9034 | 0.9352 | 0.9523 | 0.9564 |
Auf der Γout-Seite ist V_n ≈ 0. **Re V_n(t) ist auf 4 Stellen der lokale radiale Exponent κ(t) des primalen Kopfes** (die Adjungierte hat −κ auf 1e-3).
**(b) Christoffel-Darboux-Form (k57c).** K_n(t,t)·m(t) − t·δ_n′(t)/δ_n(t) =: Δ(t) ist von n unabhängig (n = 150 / 300: −0.11174 / −0.11179 bei 180°; −1.815 + 1.441i /
−1.813 + 1.443i bei 100°). Also **K_n(t,t)·m(t) = t·(log δ_n)′(t) + Δ(t) + o(1)** am Steg [numerisch; Vermutung: exakte Identität vom Christoffel-Darboux-Typ].
Der Kern, also die Summe aller Produkte primal × adjungiert, ist dort die logarithmische Ableitung des primalen Kopfes allein. Sie zählt nach dem Argumentprinzip
dessen Nullstellen. Das ist die gesuchte „Wronski-Paarung“, aber am Steg, nicht an der Ecke.
**(c) Exakte Identität [bewiesen, eine Zeile].** Das Kreismittel von ΔK_n(t,t) ist auf JEDEM Kreis gleich Σ_k c_{n,k}ℓ_{n,k}/ℓ_nn = tr J̃_n⁻¹ − tr J̃_{n−1}⁻¹
(nullter Laurent-Koeffizient; Summe über n gibt tr J̃_N⁻¹). Gemessen: 0.62279 / 0.62194 / 0.62166 / 0.62152 (n = 150 … 600), also tr J̃_N⁻¹ ≈ 0.6215·N.
**(d) Spiralen (k57b).** Σ₁ (u_A = 80°, aussen): h_Π = +0.0005 / +0.0004 / +0.0000 / +0.0003 bei r = 0.01 / 0.03 / 0.05 / 0.1 (450 → 600), also polynomial; **bei r = 0.2 wächst
es (+0.0009 / +0.0051 / +0.0284)**. Σ₂ (innen): fällt (−0.0016 … −0.043). Also: Das Produkt wächst längs Σ₁ nicht, obwohl die Primale dort wächst (+0.0074, k51d);
die Adjungierte gleicht genau aus. Für die Familie heisst das: r₁ ≤ 0.1 wählen, nicht 0.2.

### 39.4 Abgleich mit Vs Carathéodory-Fejér-Test (V §3.24, BRETT 13:00)
- **Verschiedene Grössen, kein Widerspruch.** V misst: bester Keim-Repräsentant von p_N − f_K auf der Linse (Rate ≈ 0.90·(S − 0.0147)) plus adjungierte Steg-Rate
  (≈ 0.012 − 0.95·S) ≈ **−0.0024**, unabhängig von S. Ich messe: primaler KOPF plus Adjungierte ≈ **0** (fallend). Zusammen: Der Repräsentant ist um ≈ 0.0024 je Grad
  kohärenter als ein einzelner Kopf; die Köpfe selbst sind am Steg genau ausgeglichen. Mein +0.003 war das Artefakt aus §39.2.
- **Folgerung (numerisch, aus Vs zwei Ausgleichsgeraden):** Für Linsenhöhen S ∈ (0.0126, 0.0147) haben BEIDE Reserve: der Repräsentant ist polynomial (CF-Rate ≤ 0),
  und die Adjungierte fällt am Steg (≤ 0). Also könnten (F3′) punktweise und (FA) zugleich gelten, mit einer Linse S ≈ 0.0137. Die Reserve je Seite ist nur ≈ 0.001.
  Vorschlag an V: CF und adjungierte Rate bei S = 0.0137 auf der GANZEN Linse (nicht nur am Steg).
- **Freiheit je Kopf [bewiesen, eine Zeile aus K13]:** Der n-te Summand sieht nur die Koeffizienten ≤ n von W = 𝒥_P(f − q), also nur den Keim von P mod x^{n+1}.
  Je Kopf darf man also einen eigenen Repräsentanten P̃_n ≡ P mod x^{n+1} und eine eigene Kontur nehmen. Für n < 0.95N liegt das Stegband m ≈ 0.95N von V ausserhalb
  dieses Keims. Die Bedingung (F3′) trifft also nur die Köpfe n ≳ 0.95N.
- Vorbehalt: In W steht die Komposition (f − q)∘P̃. Das Gewicht |Λ_n| hilft dort nicht linear; gebraucht wird P̃ punktweise nahe f_K, wo der Integrand lebt.
  Darum ist das ungewichtete CF bei S im Fenster der passende Test, das gewichtete der zweite.

### 39.5 Ergebnis (Abbruchkriterium: TRÄGT, Fall (a))
1. **Am Steg ist das Produkt primal × adjungiert auf n-ter-Wurzel-Ebene polynomial, sogar konvergent:** ΔK_n(t,t)·m(t) → κ(t), K_n(t,t)·m(t) = t(log δ_n)′ + Δ(t) + o(1).
   Die „Zange“ am Steg (§34.3) ist für das Produkt keine exponentielle Lücke.
2. **(KK) in der Fassung von 12:44 ist falsch** (KERNPROBLEM §4, HAUPTSATZ §2.6 werden korrigiert). Die exponentiellen Hindernisse sitzen nicht am Steg, sondern
   (a) in der Komposition der Daten (P̃ punktweise; laut V gibt es ein Linsenfenster), (b) in der Eck-Wachstumszone ausserhalb D nahe L (θ = 85°, s ≥ 0.05; Σ₁ ab r ≈ 0.2),
   (c) für den Parametrix-Weg zu (EP) in der Primalen an f(Γout) (§38).
3. **Nächster Kandidat:** eine exakte Christoffel-Darboux-Identität für die Abschnitte J̃_n zu suchen, also K_n(t,t)·m(t) − t(log δ_n)′(t) als n-unabhängige Grösse zu beweisen.
   Sie würde das Produkt am Steg explizit machen. Dazu Vs Test des Linsenfensters S ≈ 0.0137.


## 40. Eine exakte Christoffel-Darboux-Identität (CD1) für J_n — und warum sie das Steg-Produkt erklärt (Auftrag Leitung 13:09; Zeitfenster 13:13–13:18; k58)

### 40.1 Herleitung [bewiesen, elementar]
Bezeichnungen (ungeskaliert): S = untere Verschiebung auf P_n, X = J_n⁻¹, e(z) = (z^j)_j, ê(t) = (t^{−k})_k, also K_n(z,t) = e(z)ᵀXê(t) (= ε_zᵀJ̃_n⁻¹ẽ_t).
G40 ist zweiter Ordnung und lässt sich nicht direkt invertieren (§13.2). Ich nehme die Verschiebung ERSTER Ordnung des Kompositionskanals:
- p_{j+1} − x·p_j = (f − x)f^j (der Term m·x^{j+1} hebt sich weg). Also **E := J S − S J = T_n M_{f−x} C_f − (T_n f^{n+1})·e_nᵀ** (exakt).
- Für die Inverse: S X − X S = X E X.
- Einschieben: e(z)ᵀS = z·e(z)ᵀ − z^{n+1}e_nᵀ, S·ê(t) = t·ê(t) − t·e_0.
Ergebnis:
> **(CD1)  (z − t)·K_n(z,t) = q_{n+1}(z)·λ_n(t) − t·ν_n(z) + Ξ_n(z,t)**
- **q_{n+1}(z) = z^{n+1} − (J_n⁻¹T_n f^{n+1})(z)**. Er hat Grad n+1 und erfüllt T_nJq = 0. Also ist er der monische primale Kopf der Ordnung n+1 (= δ_{n+1}/ℓ_{n+1,n+1}).
- **λ_n(t) = e_nᵀXê(t) = ρ^{−n}Λ_n(t)**, der adjungierte Kopf der Ordnung n.
- **ν_n = J_n⁻¹·1 = K_n(·, ∞)**, die Antwort auf konstante Daten.
- **Ξ_n(z,t) = J_n⁻¹T_n[(f − x)·K_n(f(x),t)](z) = (1/2πi)∮ K_n(z,s)·(f(s) − s)·K_n(f(s),t) ds/s** (Reproduktionseigenschaft; jeder Kreis um 0).
Das ist die Christoffel-Darboux-Form (Kopf n+1 mal Kopf n, geteilt durch z − t) mit zwei Korrekturen: eine vom Rang 1 und der Kanalterm. Der Kanalterm ist eine
Dyson-Gleichung, quadratisch in K, getragen von (f(s) − s), also vom Abstand zum Fixpunkt. Ohne Komposition (J = M_m) ist Ξ = 0, ν = 1/m, und es bleibt die
Toeplitz-CD-Formel (ν = 1/m im Konvergenzgebiet). Dieselbe Rechnung mit dem Multiplikationskanal (J S − T_F J) gibt (CD2) mit (z − f(t)) links. Die Herleitung ist analog,
nicht nachgerechnet; zusammen mit (CD1) ergibt sie nur eine Tautologie.
**Numerische Kontrolle:** (CD1) gilt bei n = 60 / 150 / 300 auf 1e-133 / 1e-148 / 1e-176, abseits der Diagonale und als Nebenbedingung auf der Diagonale (k58).

### 40.2 Zerlegung am Steg (k58; kubisch; θ = 100°…180°, s = 0 und 0.03)
Ableitung von (CD1) nach z in z = t: K(t,t) = q′_{n+1}(t)λ_n(t) − tν_n′(t) + ∂_zΞ_n(t,t). Nebenbedingung: q_{n+1}(t)λ_n(t) = tν_n(t) − Ξ_n(t,t).
| Grösse (mal m(t)), θ = 180°, s = 0 | n = 60 | n = 150 | n = 300 |
|---|---|---|---|
| m·K(t,t) | 57.337 | 143.411 | 286.865 |
| t·(log q_{n+1})′(t) | 58.405 | 144.479 | 287.933 |
| Differenz | −1.068 | −1.068 | −1.068 |
| **B(t) := m·q_{n+1}(t)λ_n(t)/t** | 1.03970 | 1.03968 | 1.03968 |
| m·ν_n(t) | 0.470077 | 0.470077 | 0.470077 |
| m·Ξ_n(t,t)/t | −0.56962 | −0.56960 | −0.56960 |
| −m·t·ν_n′(t) | −0.07642 | −0.07642 | −0.07642 |
| m·∂_zΞ_n(t,t) | −3.311 | −6.725 | −12.416 |
Bei 100° ebenso: B = 0.98164 + 0.0873i, Differenz −2.719 + 1.434i (n = 150) bzw. −2.715 + 1.434i (n = 300).
- ν_n(t), Ξ_n(t,t) und B(t) hängen am Steg nicht von n ab (auf 1e-5). **Das Produkt „nächster primaler Kopf × adjungierter Kopf“ ist dort also ein fester Wert B(t)·t/m(t).**
- Nur ∂_zΞ wächst, und zwar genau mit −(B − 1)·t(log q_{n+1})′. Rest: C(t) := m∂_zΞ + (B − 1)t(log q_{n+1})′ = −0.992 bei 180°, gleich für alle n.
- **Umnormierte Form (numerisch):** Mit Ξ^reg(z,t) := Ξ_n(z,t) + (1 − 1/B(t))·q_{n+1}(z)λ_n(t) gilt
  > **(z − t)·K_n(z,t) = q_{n+1}(z)λ_n(t)/B(t) − t·ν_n(z) + Ξ^reg(z,t),**
  und Ξ^reg ist am Steg auch abseits der Diagonale n-unabhängig. Bei z = t·e^{i5°} und n = 150 → 300 ändert es sich um 2e-4 (180°), 4e-4 (150°), 2e-3 (120°), 7e-3 (100°).
  Der rohe Ξ und q_{n+1}(z)λ_n(t) drehen dagegen ihre Phase mit n. Näher an L konvergiert es langsamer.
- **Damit ist die Messung aus §39.3 erklärt:** m·K_n(t,t) = t·(log q_{n+1})′(t) − m t ν_n′(t) + C(t). Weil t(log q_{n+1})′ = t(log δ_n)′ + κ(t) + o(1) ist,
  folgt Δ(t) = κ(t) − m t ν′(t) + C(t). Bei 180°: 0.956 − 0.076 − 0.992 = −0.112 (gemessen in §39.3: −0.1118).

### 40.3 Was bewiesen ist und was nicht
- **Bewiesen:** (CD1) und die Dyson-Form von Ξ (§40.1). Beides gilt für jedes f, m (auch Kneser, mit T_n statt Bandstruktur).
- **Numerisch (kubisch, n ≤ 300):** n-Unabhängigkeit von ν_n(t), B(t) und Ξ^reg am Steg; die Umnormierung mit 1/B(t).
- **Offen:** ein Beweis dieser n-Unabhängigkeit. ν_n(t) ist die Galerkin-Antwort auf glatte Daten, ausgewertet am Steg (für s > 0 ausserhalb D). Auf ∂D (s = 0)
  folgt die Konvergenz bedingt aus Satz A (unter (EP), (FA), (S_E) für s₁ = J_∞⁻¹·1). Ausserhalb D braucht sie wieder Kohärenz. Ξ ist quadratisch in K. Eine Betragsabschätzung seiner Kontur-Darstellung
  scheitert an der Primalen an f(Γout) (§38), dieselbe Stelle wie (KK′).

### 40.4 Ergebnis (Abbruchkriterium)
- Formal **Teilerfolg**: (CD1) ist exakt, aber ∂_zΞ wächst mit n. **In der Sache mehr:** Das Wachstum ist genau EIN expliziter Term (Projektion von Ξ auf den Generator).
  Nach Umnormierung sind alle Korrekturen n-unabhängig. Die gemessene Identität ist damit vollständig erklärt.
- **Bedeutung für das Programm:** Der ganze Kern K_N(z,t) (die Summe aller Kopfprodukte) ist durch EIN Kopfpaar (q_{N+1}, λ_N), die glatte Antwort ν_N und den
  Kanalterm Ξ_N gegeben. Satz A und K19 brauchen dann keine Betragssumme über Köpfe mehr. Die offene Grösse ist nur noch Ξ_N, also der Kompositionskanal.
  Er verschwindet in den Fixpunkten (Faktor f(s) − s).
- **Nächster Kandidat:**
  1. Ξ_N als Dyson-Term behandeln. Kontraktion in einer Norm, in der (f(s) − s)·K(f(s), t) klein ist. Das ist die Lage NAHE den Fixpunkten, also wieder die Ecke.
  2. Oder B(t) exakt bestimmen. B = m(ν − Ξ(t,t)/t) misst den Kompositionskanal am Steg; im reinen Multiplikationsmodell ist B = 1.


## 41. Ξ_N als Dyson-Term: Die Fixpunkt-Abbildung von (CD1) ist für grosse n KEINE Kontraktion, in keiner Norm (Auftrag Leitung 13:21; Zeitfenster 13:23–13:37; k59, k59b, k59c)

### 41.1 Formulierung [bewiesen, elementar]
In Matrixform ist (CD1) die Riccati-Gleichung **[S, X] = X E X** (X = J_n⁻¹, E = J S − S J; skaliert S̃ = ρS).
- Der Operator 𝒮⁻¹ löst [S̃, Y] = R längs der Diagonalen. Er ist die Division durch z − t; Rand: die erste Spalte (ν) oder die letzte Zeile (λ_n).
- Die Dyson-Abbildung ist 𝒯(X) = 𝒮⁻¹_Rand(X E X) mit festem Rand. Ihre Linearisierung am echten X lautet L[δX] = 𝒮₀⁻¹(δX·E·X + X·E·δX).
- Varianten: E voll, oder nur der Kompositionskanal E_c (der Rang-1-Teil −σλᵀ, also das Kopfpaar, gilt dann als Datum). Rand vorn (ν gegeben) oder hinten (λ_n gegeben).
- **Kriterium:** Der Spektralradius ρ(L) hängt nicht von der Norm ab. Ist ρ(L) ≥ 1, macht keine Kontur und kein Gewicht 𝒯 lokal zu einer Kontraktion.
- **Bemerkung [bewiesen, eine Zeile]:** E = J S − S J hängt nur von f ab, denn M_m vertauscht mit S. Also löst JEDES X_m = (C_f + M_m)⁻¹ dieselbe Riccati-Gleichung.
  (CD1) als Gleichung allein sieht den Multiplikationskanal nicht; m kommt nur über die Randdaten hinein.

### 41.2 Spektralradius (k59: L als Matrix in doppelter Genauigkeit aus X in hoher Genauigkeit; k59b: Potenzmethode in arb)
Kontrollen: Riccati-Residuum ≤ 1e-71; X wird aus seiner ersten Spalte bzw. letzten Zeile auf 4e-10 (n = 40) rekonstruiert.
| n | E voll, Rand vorn | E_c, Rand vorn | **E_c, Rand hinten (λ_n gegeben)** | E voll, Rand hinten |
|---|---|---|---|---|
| 10 | 7.67 | 0.941 | 0.612 | 0.711 |
| 20 | 262 | 1.51 | 0.705 | 1.61 |
| 30 | 9.4e3 | 2.76 | 0.798 | 4.11 |
| 40 | 3.8e5 | 6.24 | 0.866 (arb: ≈ 0.87) | — |
| 50 | — | — | 0.947 | — |
| 60 | — | — | 1.016 (arb: 1.020) | — |
| 80 | — | — | arb: 1.107 | — |
| 100 | — | — | arb: 1.196 | — |
- Die beste Variante (nur Kompositionskanal, Kopfpaar als Datum, Rand hinten) überschreitet 1 bei n ≈ 58 und wächst dann wie **e^{≈0.004·n}**.
  Die anderen wachsen viel schneller (e^{0.08n} vorn, e^{0.37n} mit vollem E).
- **Wo der wachsende Vektor sitzt** (k59b, n = 80): Er liegt vollständig im Block der niedrigen Indizes, Lösungsindex und Datenindex beide < n/4. Als Kern ist δK(z,t) glatt
  und von ähnlicher Grösse an der Ecke, am Steg und auf Γout. **Es ist kein Ecken-Effekt.** Darum hilft der Faktor f(s) − s, der in L verschwindet, nicht.
- Vermutung geprüft und verworfen (k59c): Die dominanten Eigenvektoren liegen NICHT in der Tangente der m-Familie (Anteil ≈ 0.3, zufällig wäre ≈ 0.16), und L wirkt
  auf diese Tangente nicht wie die Identität (Abweichung 0.5–0.6). Warum ρ > 1 ist, bleibt offen.

### 41.3 Welche Daten man bräuchte (Frage 3 der Leitung)
- In der einzigen Variante, die überhaupt in die Nähe von 1 kommt, sind die Daten das Kopfpaar (q_{n+1}, λ_n). Das ist dieselbe Information wie (EP) (Köpfe an der Ecke), also nicht leichter.
- Mit der glatten Antwort ν als Rand (Rand vorn) ist die Abbildung schon ab n ≈ 12 expansiv. ν allein reicht nicht.

### 41.4 Ergebnis (Abbruchkriterium: GESCHEITERT)
**(CD1) als Dyson-Fixpunktgleichung trägt nicht.** Der Spektralradius der Linearisierung ist ab n ≈ 58 grösser als 1 und wächst exponentiell (langsam, ≈ 0.004 je Grad).
Da er normunabhängig ist, hilft keine Kontur und kein Gewicht. Der wachsende Anteil ist glatt (niedrige Grade), nicht an der Ecke. Der Faktor f(s) − s greift dort nicht.
Was bleibt:
- (CD1) selbst, bewiesen, als exakte Darstellung des Kerns durch ein Kopfpaar, ν und Ξ (§40).
- Die Erkenntnis, dass (CD1) allein m nicht kennt (§41.1).
- Mit beiden Kanälen ((CD1) und (CD2) zusammen) wird das System linear: Die Differenz gibt X(T_F − S) = X(E − E′)X mit E − E′ = (T_F − S)J, also nur X = J⁻¹ [eine Zeile].
  Ein Dyson-Argument muss also einen Teil von J als bekannt behandeln. Das ist wieder Störungsrechnung, und die scheitert an der Fehlskalierung (KERNPROBLEM §2).


## 42. K3-2: Der Eckkern ist NICHT lokal integrabel — der Dyson-Term ist auf der Eckskala von derselben Ordnung wie der Generator (Auftrag Leitung 13:52; 13:53–14:01; k60)

### 42.1 Messung (k60, kubisch, n = 150/300/600; (CD1) jeweils auf ≤ 5e-150 erfüllt)
Punkte auf Σ₂ (innen, u_A = −90°) bzw. Σ₁ (aussen, u_A = 80°) im Abstand r = c·|d_n| von L. Gemessen ist |Ξ_n(z,t)| / |q_{n+1}(z)λ_n(t)|:
| z (Σ₂, c) | t (Σ₁, c = 0.5) | t (Σ₁, c = 2) | t (Σ₂, c = 2) | t (Σ₁, c = 5) |
|---|---|---|---|---|
| c = 0.5 | 0.80 / 0.91 / 1.02 | 0.74 / 0.86 / 0.99 | 0.81 / 0.90 / 1.01 | 0.53 / 0.67 / 0.82 |
| c = 1 | 0.97 / 1.08 / 1.20 | 0.85 / 0.98 / 1.12 | 1.02 / 1.10 / 1.20 | 0.46 / 0.63 / 0.80 |
| c = 2 | 1.29 / 1.42 / 1.55 | 1.01 / 1.18 / 1.35 | — | 0.26 / 0.46 / 0.69 |
| c = 5 | 2.09 / 2.30 / 2.50 | 1.12 / 1.39 / 1.67 | 3.69 / 3.68 / 3.68 | 1.00 / 0.65 / 0.26 |
|tν|/|q_{n+1}λ_n| liegt bei 1.7 … 6.9. **Ξ ist auf der Eckskala nicht klein und wächst langsam mit n** (etwa +0.11 je Verdopplung). Die eine fallende Zelle (5, 5) ist eine
Auslöschung im Generator.
- **Diagonale (t = z, exakt: q_{n+1}(t)λ_n(t) + Ξ_n(t,t) = t·ν_n(t)):** t·ν_n ist nahe L n-unabhängig, bis r = 0.003 auf 6 Stellen
  (z. B. 0.3744 + 0.5552i bei n = 150/300/600). Der Generator q_{n+1}λ_n hat auf der Diagonale ein Eckprofil in c und driftet wie n^{−0.11} (wie P_n).
  **Ξ_n(t,t) gleicht beides genau aus.** Der Kompositionskanal trägt also auf der Eckskala dieselbe Struktur wie der Generator, mit umgekehrtem Vorzeichen.
- An festen Punkten (r = 0.003 … 0.03) ist Ξ_n NICHT n-unabhängig (Änderungen 0.1 … 0.3). Auch die schwächere Form („Ξ glatt nahe L“) gilt nicht.

### 42.2 Warum (Skizze)
Auf der Eckskala ist alles skaleninvariant (Mellin). Der Faktor f(s) − s ≈ (μ − 1)(s − L) ist von der Grösse d_n. Der Kern K(f(s),t) im Eckbereich ist aber von der Grösse 1/d_n
(K wie 1/(t − f(s)) mal Profil). Das Produkt ist homogen vom Grad 0, also O(1), genauso gross wie der Generator.
Das Verschwinden von f − s in L liefert darum keinen kleinen Parameter; es wird von der Polstruktur des Kerns genau aufgefressen. Das passt zu Gs Symbol G39:
In L treffen sich beide Kanäle (1 − sF und 1 − sz) in erster Ordnung. Die Ecke ist echt zweikanalig, ein Rang-2-IIKS-Kern mit einem Generatorpaar reicht dort nicht.

### 42.3 Ergebnis (Abbruchkriterium: GESCHEITERT, strenge und schwache Form)
- K3-2 trägt nicht. Ein lokal integrabler Eckkern mit dem Kopfpaar als einzigem Generator existiert nicht. Der Kompositionskanal ist an der Ecke von führender Ordnung.
- Für einen RH-Ansatz heisst das: Die lokale Parametrix an L muss beide Kanäle enthalten, mindestens 3×3 oder mit Verschiebung (K3-4). Die Generatoren q_{n+1}, λ_n reichen nicht.
- Nebenbefund (numerisch): Die glatte Antwort ν_n = J_n⁻¹·1 ist bis r = 0.003 an L heran n-unabhängig. Also konvergiert die Galerkin-Lösung für glatte Daten bis dicht an
  die Ecke, ohne Eckstruktur.


## 43. K3-3: Gerüst für den (F3′)-Bootstrap — es schliesst termweise nur auf der kritischen Linse mit Marge 0; die Marge 0.0024 ist Kohärenz der Paarungen (Auftrag Leitung 14:02; 14:05–14:06)

### 43.1 Gerüst (Logik bewiesen bis auf die eine Ungleichung; Kennzeichnung je Schritt)
Daten: Linsenhöhe S, Spiralschnitt r₁, Kontur γ = γ_{S,r₁}, konforme Abbildung k: 𝔻 → Int γ. P liegt auf dem σ-Homotopiepfad vom Grad N, und P − f_K = Π^P f_K − f_K.
- **CF-Grösse:** CF(P) := min_g sup_γ |P − f_K + x^{N+2}g| = ‖T_{N+2}(c)‖, c_m = [w^m]((P − f_K)∘k) (Carathéodory-Fejér/Schur) [klassisch].
- **Hypothese (H_δ):** CF(P) ≤ δ₀. δ₀ ist so klein, dass der Extremal-Repräsentant P̃ die Kontur γ nach E abbildet und M_{P̃} auf γ beschränkt ist.
- **Schritt 1 (Daten) [bewiesen unter (H_δ)]:** W̃ := 𝒥_{P̃}(f_K − q_N) erfüllt sup_γ|W̃| ≤ (1 + M)·dist_E(f_K, P_N) ≤ C·N^{−s_E}. Das ist K19 mit (S_E).
- **Schritt 2 (Darstellung) [bewiesen]:** P − f_K = (q_N − f_K) + Σ_{n≤N} δ_n^P·π_n/ℓ_nn^P, π_n := ⟨ℓ_n^P, T_nW̃⟩ = (1/2πi)∮_γ Λ_n^P(s)W̃(s)ds/s.
  Der Wert π_n hängt vom Repräsentanten nicht ab; je Kopf genügt der Keim mod x^{n+1} (K13 kopfweise).
- **Schritt 3 (Paarungen) [bedingt auf F2]:** |π_n| ≤ |Λ_n(L)|·R_A^γ(n)·C·N^{−s_E}, R_A^γ(n) = max_γ|Λ_n/Λ_n(L)|.
- **Schritt 4 (Band):** c_m = [w^m]((q_N − f_K)∘k) + Σ_n [w^m](δ_n∘k)·π_n/ℓ_nn. Gebraucht wird
  > **(BU)  ‖T_{N+2}( Σ_{n≤N} [w^·](δ_n∘k)·π_n/ℓ_nn )‖ ≤ C·N^{κ_B}·N^{−s_E},  κ_B < s_E.**
- **Schritt 5 (Stetigkeit) [Logik wie (EB_N), bewiesen]:** CF(P) ist längs des Pfads stetig. Bei N = 400 gilt CF ≤ δ₀/2; das braucht eine endliche Prüfung von p₄₀₀ − f_K auf γ,
  noch nicht streng. Aus (BU) folgt CF ≤ δ₀/2 für N ≥ N₀, also wird δ₀ nie erreicht, und (F3′) gilt für alle N.
- **Kein Kompositionstransport:** Primale Köpfe werden nur auf γ gebraucht (über ihre w-Koeffizienten), die Adjungierte nur auf γ, P̃ nur auf γ.
  f − q wird auf P̃(γ) ⊂ E ausgewertet; das ist Datenseite, keine primale Grösse an f(Γout). Das Abbruchkriterium (Werte ausserhalb der Linse im Betrag) greift NICHT.

### 43.2 Termweise reicht es nur auf der kritischen Linse, mit Marge 0
Mit der Dreiecksungleichung wird (BU) zu Σ_n ‖T(δ_n∘k)‖·|Λ_n(L)|R_A^γ(n)/|ℓ_nn| = Σ_n |P_n|·R_P^γ(n)·R_A^γ(n) (bis auf Konstanten).
Dabei ist R_P^γ(n) = sup_γ|δ_n/δ_n(L)| (CF-Norm ≤ sup). Mit den gemessenen Raten am Steg (k57, V v58): h_P(S) ≈ 0.957·S − 0.012 und h_A(S) ≈ −h_P(S),
weil das Produkt flach ist (§39). Also:
- termweise Rate = max(h_P(S), 0) + max(h_A(S), 0) = |0.957·S − 0.012| ≥ 0;
- **Null nur auf der kritischen Linse S₀ ≈ 0.0125**, wo Kopf und Adjungierte beide subexponentiell sind (Vs Nullstelle der adjungierten Steg-Rate: 0.0126).
  Dort entscheiden die Polynomfaktoren. Eine strikte Reserve gibt es termweise NIRGENDS, eben weil das Produkt flach ist.
- Längs Σ₁ wachsen die Köpfe wie e^{≈0.15·r·n} (k51d), die Adjungierte nicht (A-lok). Termweise muss der Spiralschnitt also auf der Eckskala liegen, r₁(n) ~ C/n,
  mit kopfweisen Konturen (K13 erlaubt das). Ein fester r₁ = 0.05 kostet asymptotisch e^{0.0075n}, mehr als Vs Marge.

### 43.3 Was die Marge 0.0024 ist
Vs CF-Rate ≈ h_P(S) − 0.0024 (V §3.24) bedeutet: Der Repräsentant ist um 0.0024 je Grad besser als die termweise Schranke.
Die Dreiecksungleichung in (BU) wirft genau die Phasen der Paarungen π_n weg. Jede Schranke, die nur |π_n| ≤ ε_n benutzt, kann die Marge nicht sehen [elementar]. Denn das Mittel über die Phasen von ‖c‖₂² ist Σ|π_n|²·‖a_n‖₂²
(a_n = w-Koeffizienten von δ_n∘k), und ‖T(c)‖ ≥ ‖c‖₂. Also erreicht eine ungünstige Phasenwahl die termweise n-te-Wurzel-Rate.
**Also ist (BU) mit Marge eine Kohärenzaussage über die tatsächlichen Paarungen π_n**, dieselbe Art Aussage wie (Fam-koh), jetzt aber genau benannt und messbar.

### 43.4 Ergebnis (Abbruchkriterium: nicht ausgelöst; Teilerfolg)
- Das Gerüst steht. (F3′) für alle N folgt aus EINER Ungleichung (BU), plus der endlichen Startprüfung bei N = 400, plus F2. Kein Kompositionstransport.
- Termweise schliesst (BU) nur auf der kritischen Linse S₀ ≈ 0.0125, mit Spiralschnitt auf der Eckskala, und mit Marge 0. Das ist grenzwertig:
  Vs Primalrand-Übergang bei n ≈ 2200 (v56) kann die Raten dort verschieben.
- Robust (Fenster S ∈ (0.0126, 0.0147), Marge 0.0024) ist nur die kohärente Form. Sie braucht die Phasen der π_n, nicht nur ihre Grösse.
- **→ V, zwei Messungen:**
  (1) auf der kritischen Linse S₀ = 0.0125 bis N = 3000: R_P^γ(n) und R_A^γ(n), also polynomial oder nicht, über den Übergang bei n ≈ 2200 hinweg;
  (2) das gewichtete CF als direkter Test von (BU) mit den echten π_n.
