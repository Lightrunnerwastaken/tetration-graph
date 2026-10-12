# Beweis-Buchhaltung — Abhängigkeitsbild des Gesamtbeweises (Strang V, Phase 3)

## Stand des Gesamtbeweises (für Janis, 11.10. 22:54, ergänzt 12.10. 02:17)

**Ziel:** Die Polynome p_N (p_N∘p_N = eˣ bis zum Grad N) nähern sich auf der Scheibe D̄ gleichmässig Knesers f_K.
**Aufbau:** bis N₀ rechnet der Computer streng nach; ab N₀ zeigt ein Newton-Argument, dass p_N in einer kleinen Kugel um f_K bleibt.
Dafür muss die Inverse des linearisierten Problems für alle N höchstens polynomial wachsen (Bedingung (H)).

**Fertig und von einem zweiten Strang geprüft:**
- f_K selbst ist streng bekannt: Koeffizienten bis 640 und eine Schranke für alle weiteren (P-V, p32a, (K2): |b_k| ≤ 461·k⁻⁶).
- Der Grenzoperator ist nach unten beschränkt: σ_min ≥ 0.641 für alle N (Ms 𝒯-Weg).
- Der Rechnerteil reicht bis N₀ = 400: Zweig 1 → 193 (P), 193 → 400 (V, unabhängig nachgeprüft). Die Anschlussbedingung (A) gilt mit Reserve 1e8.
- Die logische Hauptkette und die Umformung des Eckproblems in die Form (S_γ) (Ks Sätze K1 und K1″).

**Genau offen:**
1. **(S_γ):** In der konformen Koordinate mit einer festen Linse gilt ‖Ĵ_N⁻¹‖ ≤ B für alle N ≥ N₀. Die ganze Kette von dieser Zahl bis zur
   Vermutung ist mit Zahlen durchgerechnet (BERICHT §3.9, §3.11, Satz V2). **Mit streng geprüfter Geometrie und streng eingeschlossener Green-Funktion
   genügt B ≤ 16.4 ab N₀ = 400** (das Zweig-Zertifikat reicht bis 400; seit 12.10. 09:00 mit STRENGER Sprungschranke C_J = 1.0378 für u ≤ 0.48,
   BERICHT §3.15 — vorher 17.9 mit dem numerischen 0.95). Gemessen ist B = 1.25 — bis N = 2000 an f_K, an p_N und an allen geprüften
   Zwischenpunkten der Homotopie. Noch nicht streng in dieser Kette: Xs Lemmata (ohne zweiten Leser) und der Datenteil der Familie (Fam).
   (Stand 12.10. 09:30: Ks Zielfassung K15 — gleichmässig auf D̄ ohne feste Umgebungen von L, L̄ — braucht (S_γ) nicht mehr, k_kern §23–§24.) **Das ist der Engpass:** ein Beweis von (S_γ) mit einer Konstante unter 17.9.
2. **Familie:** Dasselbe muss gleichmässig für alle Operatoren gelten, die im Newton-Schritt vorkommen (P nahe p_N, nicht nur P = f_K).
   Eine reine Störungsrechnung um f_K reicht nicht für alle N, und keine Norm, die nur Grössen misst, leistet es (Messung V, LANDKARTE C11, C12).
   Gemessen gilt (S_γ) aber am echten p_N bis N = 1000 mit derselben Konstante wie bei f_K. Der Weg ist daher eine Homotopie innerhalb der
   Lösungsmenge (wie in Satz H1/H2); Ks Kopplungssatz soll das mit abdecken.
3. Kleinere Punkte: Xs Lemmata X5, X7, X8′ und X4a ohne zweiten Leser; einige Randkonstanten nur numerisch; Kopf-Stabilität (H-iii) nur numerisch.

Stand: 11.10.2026 20:01 (aus `date`). Wird laufend gepflegt. Quellen: ARCHITEKTUR.md §4, §9 (massgeblich), p_pruefung/PRUEFLISTE.md,
x_alternativen/KONSTANTEN.md, m_modell_eckmode/BERICHT.md §15–§17, o_operator/BERICHT.md §33–§34, w_grenzschicht/BERICHT.md §22–§27,
g_gpt_hilfssaetze/ (Runden 3–5).

**Statuswörter (streng getrennt):**
- **V** = verifiziert: von einem zweiten Strang Schritt für Schritt nachvollzogen (bisher fast immer P).
- **S** = selbstgeprüft rigoros: vom Autor rigoros gerechnet/bewiesen (arb), kein zweiter Leser.
- **B** = behauptet bewiesen: Autor sagt „bewiesen“, niemand hat gegengelesen.
- **num** = numerisch gestützt (endlich viele N), kein Beweis. · **offen** · **Fehler** · **Lücke**.
- „|X“ = unter Annahme X.

## 1. Hauptkette (logisch geschlossen, bedingt)

```
(K) ──L8.2–8.4 [V|(K)]──▶ (F)  ◀── b35b′ [num], r > 0.02 [num], Univalenz [num]
(F) + (H₀) ──4.1 = P-I(i) [V]──▶ (LP): sup_D̄|Π_N f_K − f_K| ≤ C' N^{κ−4.81}
(F) + (H) ──H1 [V|(H-i),(H-ii)]──▶ ∃! p*_N ∈ B_N
(H) + (A) ──H2 [V|(H),(A)]──▶ p_N = p*_N ∀N ≥ N₀ ──▶ VERMUTUNG
(H) := (H₀) [an P = f_K] + Familie (H-i)–(H-iii) für alle P ∈ 𝒮_N ∩ B_N
(A) := Zweig-Zertifikat bis N₀ + Kopplung sup_D̄|p_{N₀} − f_K| ≤ R_{N₀}, N₀ ≥ N₀(H)
```
(K) für b = e: ersetzt durch Satz P-V + p32a [S, P] (zitiert: Ahlfors–Bers, Stoilow, Beurling-Norm, Paulsen–Cowgill Prop. 2).

## 2. Bausteintabelle

| # | Baustein | Aussage (kurz) | Status | Konstanten | braucht | Fundstelle |
|---|---|---|---|---|---|---|
| 1 | P-V | sexp_K auf Streifen ± 1.3e-168; a_0…a_640, Fehler·ρⁿ ≤ 2.9e-127 | **V** (V, 11.10. 20:20: Schritte 3–6 + kc_rig gelesen, Logik richtig; formale Lücken L1–L4 ohne Wirkung, BERICHT §1.1); M: Ausgabe num. gegengeprüft | ε = 1.3e-168 | p26–p29, Paulsen–Cowgill | ARCH §9.8 |
| 2 | p32a / (K3) | sexp_K = P(e^{L(z+θ_K)}), Im z ≥ 0.36; f_K(L) = L, f_K'(L) = √L | **V** (V 21:29, BERICHT §1.7) | |θ_K − θ̃| ≤ 4.4e-167 | P-V | ARCH §9.9 |
| 3 | (K2) | |b_k| ≤ C₆ k^{−6} ∀k ≥ 1 | **V** (Rechnung V 21:13, Neulauf mit Ks Korrekturen F1/F2 22:15 identisch; zweiter Leser K 21:57; P-Bausteine von V geprüft) | C₆ = 461.36, Σ_{k>640}|b_k| ≤ 8.6e-13; Datei v_pruefstand/kneser_halbexp_b_zert_V.json | P-V, p32a | BERICHT §2.2 |
| 4 | (F) | f_K holomorph um K bis auf Schlitze, Sprung ≤ C_J|ζ−L|^{5.447} | **V|(K)** + num.: b35b′ (0.7345), r > 0.02, Univalenz (0.3456) | C_J ≤ 1.0378 (u ≤ 0.48, **streng**, v44; num. 0.94) | (K) | PRUEF §2 |
| 5 | 4.1 / P-I | (H₀^K) ⇒ (LP) mit Rate N^{κ−5.447α_K} | **V** | Schlitz 12.7 [Formel V, Zahlen num.]; Kontur d = 0.1: 37.9, d = 0.3: 17.1 [num.] | (F) | ARCH §9.4, §9.7(1) |
| 6 | P-II | (a′)_r ⇒ Rest-Schranke in K' | **V** bis auf Univalenz [num.] | |w_η| ≤ 5.8 N^{3/2} | (a′)_r | ARCH §9.4 |
| 7 | P-III, P-IV | Eck-Kriterium (pα > 1), Lokalisierung | **V** (P eigen, kein zweiter Leser!) → eigentlich **S** | | O6, (E_L⁺), (E_sing) | ARCH §9.4 |
| 8 | X1, X3.1, X3.2 | (LP) ⇐ (a″) + (Ẽ); Rechteck-Splitting | **V** (P) | c₀, C_R | (a″), (Ẽ), P-II | ARCH §9.5 |
| 9 | (a″)_∞ Kneser | σ_min(J_∞|P_N) ≥ 0.641 ∀N (Fenster 7/23) bzw. 0.5696 (6/14) mit Vs (K2)-Datei | **V** (M; Schlusskette M5, IMS, M12, m26 und Lemmata M7–M10, M13 = m16, m22–m25, mk_kneser_arb von V gegengelesen 22:31; Prüfzeilen der Teilläufe einzeln nachgesehen) | c₀ = 0.6553 (wahr 0.7319) | (K1) ✓, (K2) ✓, (K3) ✓ | M §15.1, §17; BERICHT V §1.8 |
| 10 | (a″)_{4N/5N} | σ_min(T_{5N}J_∞|P_N) ≥ c₀ − τ_N | **B** (X5: bew.|(K), c₁ Skizze) | τ_N ≤ 3e-12 (N = 193) | (a″)_∞, X5 | KONST §2 |
| 11 | X7 | |w_η|_{ℓ²_ρ} ≤ 6.3 √M auf ∂K' | **B** (Beweis gelesen von niemandem; Gitter-Minima c_in = 0.6741, c_out = 0.3456 num.) | C_R = 12.7, κ_R = 1 | Univalenz, f_K-Werte auf ∂D | X §3.9.1 |
| 12 | X8 / X8′ | ‖(𝒥_P − J_∞)δ‖ ≤ (A_θ N + …) R ‖J_∞δ‖ | **B** (Randwerte num.) | A_θ = 19.7 / 13.4 (X8′) / 9.2 (3N-Zeilen, braucht Sprung bis u = 0.46) | X5, c₀ | X §3.9.2, §3.10.1 |
| 13 | X4a | eckfreier Teil der Familie gleichmässig | **B** (modulo (a″)_∞, X5, X8, Univ.) | θ_N = A_θ N R_N ≤ ½ | 9, 10, 12 | X §3.7 |
| 14 | (Fam-Eck) | Eckmoden der 𝒥_P relativ O(N^a R) stabil | **Vermutung** (num. N ≤ 600; Bogenmodell ε ≈ 0.06) | | | X §3.7, §3.9.6 |
| 15 | (H₀) Eckteil | **(S_γ) sup_N‖Ĵ_N⁻¹‖ < ∞** (Ks K1″, **V** 21:17: (S_γ) ⇒ (H₀^{K_γ}) mit Faktor N^{3/2}); schwächer (H*_G) (K1, V) | **offen** (Engpass); **numerisch bestätigt von V unabhängig** (kubisch N ≤ 1500: 1.251155; Kneser N ≤ 639: 1.250337) für Linsen über Ws Profil Λ; versagt darunter (Schwelle S* ≈ c_∞ = 0.0118) | BERICHT §3.1 |
| 16 | W1'', W3, W5 | Reduktionen für den Eckteil | **V** (P, 04:50) | | (H*) auf γ*, (G*) | PRUEF §5 |
| 17 | Kneser-Version von W1''/D6 | Eckformel für nicht-polynomiales f | **gegenstandslos** seit Satz K1 (braucht kein D6) | | (E_sing) | ARCH §9.2(2), O §33.1 |
| 18 | Zweig 1 → 193 | σ-Homotopie-Zweig rigoros | **S** (P: p18 + p24, auf aicgolling reproduziert p30) | Box ·ρ^k ≤ 1e-120 | — | ARCH §9.7(3b) |
| 19 | Zweig 193 → 400 | σ-Homotopie-Zweig rigoros | **S** (V: zert_kern = Ps p24 + Xs NC; Kette von V unabhängig mit eigenem Code geprüft, 21:59) | 207 Übergänge, Reserve ≥ 0.4999 | 18 | BERICHT §2.1 |
| 20 | Kopplung (A) | sup_{D̄}|p_N − f_K| ≤ 3.6e-12 (N = 318) … 1.75e-12 (N = 400); P: N = 140…193 | **S** (V v15, P p31) | Reserve ≈ 1e8 gegen R_{N₀} | 1, 3, 18, 19 | BERICHT §2.1 |
| 21 | N₀(H) | kleinstes N mit θ_N ≤ ½ und allen Eckbedingungen | **bedingt**: 318 mit X8′, Kontur d = 0.3 (Hybrid), **Eck-Konstanten C_E = 1, κ_E = 2 ANGENOMMEN** | 318 (bew. Teile), 262 (3N), 93 (A* = 0.63), 47 (gemessen) | 5, 11, 12, 15 | KONST §000 |
| 22 | (H-iii) | Kopf-Stabilität | **num.** (Quotient 1.0000, N = 60; J(f) ↔ J(pad p_N) 1e-10) | | Familie | ARCH §5.3 |
| 23 | G8–G11, G15–G23 | Diagonalisierung, Rang-1-Zuwachs, Summenkriterium G19/G20, Kreis-Majorante | **V** (V 20:25: G8–G11, G22, G23 richtig, elementar; die All-n-Hypothese G19 bleibt offen) | C_G = 6.26 falls C_Δ = 3.13 ∀n ≥ 200 | — | G R3, R4 |
| 24 | G25–G44 | Modelle, Zwei-Kanal-Symbol, Rang-2-Identität, G44 (F = 2z + ηz²) | G39, G40, G43, G44 **V** (V 20:25; G44 ohne Abschneide-Effekt, da F(0) = 0); G25–G38 **B** | |𝒟_n| ≤ 20(2/3)^n in G44 | — | G R5/ETAPPE4 |
| 25 | O1, O2, O7a, O7g–i, O9, O10 | Faktorisierung, Herglotz, Fredholm, c_∞ ≥ c(λ), λ_L | **V** (P) | λ_L = 0.98988321 | | ARCH §9.5 |
| 26 | O11 | c_∞ ≥ 0.54 durch Konturverschiebung | **num.** | | | O §30 |
| 27 | O12, Majorante G | Eckformel in L; |Y_u| ≤ G für alle u | **B** (O; elementar) | |𝒰| ≤ √(N+1)/ρ | | O §33.1 |
| 28 | M6/M11/M12 radial | 𝒯 akkretiv, (a″)_∞ ≥ 0.6395 | **B** (M) | | | M §11, §17 |

## 3. Was für den Gesamtbeweis noch fehlt (Kurzliste, nach Art)

**Analytisch (unendlich viele N, K-Seite):**
1. (H₀) Eckteil für alle N ≥ N₀ mit expliziten C_E, κ_E (κ_E ≤ 2 bei C_E ≤ 1, damit N₀ ≤ 400 trägt). Das ist der Engpass.
2. Kneser-Fassung der Eckreduktion (W1''/D6 für A_f, (E_sing)) — oder ein Eckweg, der nur endliche Matrixidentitäten benutzt (G10 tut das).
3. Familie: (Fam-Eck) für den Eckteil (der eckfreie Teil ist mit X4a B). Messung V 22:48 (BERICHT §3.3, LANDKARTE C11): im K1″-Rahmen
   ist Lip_N von P ↦ (Ĵ^P)⁻¹ linear in N für glatte Störungen, aber ≈ e^{0.011N} für Grade k ≥ 0.8N ⇒ Störung um f_K reicht nicht für alle N.
4. (H-iii) Kopf-Stabilität (bisher nur num.).

**Endlich (V-Seite, rechnerisch erreichbar):**
5. (K2) mit p = 6 und C₆ ≤ 5e5 (statt 300) — entsperrt M (Nr. 9), (A) (Nr. 20), X5.
6. Zweig 193 → 400 (läuft) + Kettenprüfung + Kopplung bei N₀ ∈ [318, 400].
7. Zertifizierte Randkonstanten aus P-V-Daten: c_in, c_out, min|f_K'|, max|m|, max|m'|, max|m/f'|, M2, C_J (erledigt, v44), b35b′, Univalenz, Sprungschranke bis u ≈ 0.46 (erledigt bis 0.48, v44),
   Holomorphie innerhalb Γ_{0.3}.
8. Korrektur der Koeffizientendatei (Dezimalrundung in die Radien) bzw. Nachweis, dass jeder strenge Verbraucher sie addiert.

**Prüfung (zweiter Leser fehlt):**
9. P-V Schritte 3–6 (p27, p28), p29-Fortpflanzung, kc_rig, p32a — Ps 7 Prüfstellen (ARCH §9.8).
10. M: m22–m26, mk_kneser_arb, mk_ecke_kneser, M5, M7–M10, M13.
11. X4a/b, X5, X6 (entbehrlich), X7, X8, X8′, Kontur d = 0.3.
12. G8–G44 (mindestens G9/G10/G11/G22/G23 für den Eckweg, G44 als Modell).
13. p24/zert_kern (Strenge der Röhre: Knotenrundung, Lebesgue-Konstante, Majoranten), p31 (Kopplung, Radien der Koeffizientendatei).
