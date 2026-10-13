# Plan (Phase 3)

Gewählte Angriffe mit Verantwortlichem, Zeitfenster und Abbruchkriterium. Bei jeder Änderung Datum und Grund dazuschreiben.

## Vorschlag K (11.10. 20:37) — von V zu bestätigen oder zu ändern
| # | Angriff | Verantwortlich | Zeitfenster | Abbruchkriterium |
|---|---|---|---|---|
| A1 | **Fundament K1**: Satz K1 gegenlesen (V), Kneser-Numerik (erledigt k01, N ≤ 600), Teil C auf (H*_G) umstellen | V prüft, K pflegt | bis 11.10. 23:00 | V findet eine Lücke, die sich nicht schliessen lässt ⇒ zurück zu (Ẽ)/X2 |
| A2 | **(H*_G) im Modell mit Randfixpunkt**: lox. Möbius zuerst (V misst Ecke/κ₁/G auf γ; K sucht exakte Struktur von G_N bzw. 𝒟_n auf γ), sonst radiales Quadrat | K (Analyse), V (Messung) | 2 Arbeitsblöcke (≈ 4 h) ab Bestätigung | Möbius-Ecke nicht universell ODER keine explizite Struktur für G_N auf γ in Sicht ⇒ Wechsel zum radialen Modell bzw. zu A2' |
| A2' | **Ecke direkt**: (KB) für G auf γ ∩ U(L) als lokale Parametrix mit Quelle (O §33.4, Abel-Zylinder K8) | K | nach A2, 1 Block | keine Formulierung, in der die Amplitude nicht global eingeht ⇒ dokumentieren, als benannte Hypothese führen |
| A3 | **Prüfsteine**: c₂-Formel (V2/K7), (Fam-Eck)-Test (V9), Falsifikation G_N auf γ bis N ≈ 3000 (V10) | V | parallel, 1 Tag | c₂ hängt nachweislich von globalen Daten ab ⇒ V2 streichen |
| A4 | **Rechteck-Satz** (K2) als gesicherter Zwischensatz | K schreibt, V prüft | nach A2, ≈ 3 h | LS-Zweig braucht neue Idee ⇒ nur linearer Teil |
Parallel (V, unverändert): (K2) mit C₆ ≤ 5e5, Zertifikat bis 400, Randkonstanten, Prüfschlange B5.

## Bestätigung V (11.10. 20:52) — einverstanden mit A1–A4, mit diesen Änderungen
- **A1 erledigt:** Satz K1 von V gegengelesen und verifiziert (v_pruefstand/BERICHT §1.3). Korrektur: mit Familie über X4a braucht (H*_G)
  κ < 1.79 (nicht 2.79); für N₀ ≤ 400 Eck-Polynomfaktor κ + 3/2 ≤ 2 bei Konstante ≲ 1 — die „3/2“ zu drücken gehört zu A2/A2'.
  Die d-Frage (Landkarte C5, Kante) entfällt mit K1 ganz.
- **A2:** V misst im loxodromischen Möbius-Modell (Ecke c, κ₁, Profil, G_N auf γ, N ≤ 800) — Zeitfenster 2 h ab 21:00; Abbruch wie von K.
- **A3:** Falsifikation nur bis N ≈ 2000 (kubisch, aicgolling, Hintergrund; 3000 kostet ≈ 18 h). Neu in A3: **Familie in G-Form** —
  G^P_N(t) = ‖(T_N𝒥_P)⁻¹T_Nk_t‖·|t| auf γ für P = p_N (echte Lösung) und für glatte/zufällige Störungen gleicher Grösse (Xs Tabelle in G-Form);
  Abbruch: Verstärkung > poly(N)·R_N bei kohärenten P.
- **A5 (neu, V): endlicher Unterbau des bedingten Hauptsatzes** — (K2) über v06 (C₆ ≤ 5e5 genügt, v04), Zertifikat 193 → 400 (läuft seit 19:55),
  Kopplung (A) bei N₀ ∈ [318, 400], zertifizierte Randkonstanten inkl. K1-Geometrie (Linse ∩ Schlitz = ∅, f(ℬ) ⊂ D, Lipschitz von g_γ),
  Prüfung von Ms 𝒯-Weg (m22–m26) und X4a/X5/X8′. Zeitfenster bis 12.10. abends. Abbruch (K2): wenn die Zellen nahe L nicht endlich werden ⇒
  Rückfall p = 5 oder grösseres δ_min mit eigener Restschranke.
- Ergebnis, das jederzeit vorzeigbar bleibt: **bedingter Hauptsatz mit genau einer linearen Hypothese (H*_G) für J_∞ und ihrer Familienfassung**,
  alles andere streng (Ideen V7 + K2).

## Anpassung K (11.10. 21:45), nach Auftrag der Leitung (21:00) und Versuch 3
- **A1 erweitert, erledigt:** Lemma K1′/K1″ verifiziert (V 21:17), numerisch entschieden: (S_γ) sup_N‖Ĵ_N⁻¹‖ < ∞ gilt bis N = 800 (kubisch) / 600 (Kneser);
  V rechnet unabhängig bis 1500 bzw. 639 und einen Linsenhöhen-Sweep. **Zielaussage des Engpasses ist jetzt (S_γ)** (stärker und sauberer als (H*_G));
  Eck-Polynomfaktor κ_Ĵ + 3/2 mit gemessen κ_Ĵ ≲ 0.15 (Linse nahe L dicker: τ₁ = 5°, p = 1.15).
- **A2 (Parametrix):** Möbius abgebrochen (V 21:15). Neu: radiales Quadrat. Schritte: (a) (S_γ) im radialen Modell prüfen (läuft);
  (b) Bausteine Bulk / Kantenband / Ecke einzeln numerisch isolieren (Widom-Differenz zerlegen, Kantenband als eigener Grenzoperator);
  (c) erst dann analytische Parametrix. Zeitfenster bis 12.10. mittags; Abbruch, wenn das Kantenband keine N-unabhängige Gestalt (in ξ = Index/N) zeigt.
- **A2' (Ecke):** zurückgestellt, wird Teil von A2 (Eckbaustein).

## Anpassung K (11.10. 22:12), nach Auftrag der Leitung (22:00) und Versuch 5
- **A2 neu gefasst („Stabilitätssatz“):** klassische Abschnittstheorie greift nicht direkt (BERICHT §8.2: exponentielle Fehlskalierung; lokale
  Invertierbarkeit genügt nicht, x-Koordinate als Gegenbeispiel). Ziel ist Vermutung S (§8.3): lokal (Bulk, Ecke) + Schleifenbedingung ⇒ stabil.
  - **A2a (K):** Eckmodell an L explizit (halbseitiges Wiener–Hopf für I + μ⁻¹D_μ in lokalen Frequenzen), Invertierbarkeit beweisen, gegen D_N testen.
    Fenster: bis 12.10. vormittags. Abbruch, wenn das Modell D_N nicht qualitativ trifft (Lage, Indexband, Paarstruktur).
  - **A2b (K):** Schleifenbedingung als globale Ungleichung formulieren und prüfen, ob Ws tropische Theorie (W11, Unterlösungen) sie in w-Koordinaten liefert.
  - **A2c (K, danach):** Kopplungssatz „lokal + global ⇒ stabil“ — die eigentliche neue Theorie.
- **Für V (A3, neu):** (i) Universalität der Eck-Korrektur D_N (kubisch / Kneser / radial: Singulärwerte, Paarstruktur, Profil in ξ); (ii) Abhängigkeit
  von der Linsenform nahe L (sättigt ‖D_N‖ bei dickerer Linse?); (iii) „Kante fern von L trivial“ direkt prüfen (D_N mit Winkel-Abschneidung komprimieren);
  (iv) Schärfe der Schleifenbedingung profilweise (Vs Linsenfenster) — die Daten dafür hat V schon.

## Anpassung K (11.10. 23:40), nach der langen Strecke A2a–A2b (BERICHT §9, REFLEXION [K] Versuch 6)
- **A2a erledigt:** Eck-Korrektur universell (numerisch); eigenständiges Eckmodell mit Gradabschnitt widerlegt (affines Modell, bewiesen + numerisch).
  Abbruchkriterium „Modell trifft D_N nicht“ ausgelöst — aber mit klarem Grund: die Ecke ist im Gradabschnitt nicht abtrennbar. Ersatz: Vermutung S′.
- **A2b erledigt:** Satz K3 — Tropik allein liefert (SB) in keiner Koordinate; (SB) ⇐ W12′ + Eck-Normierung u(L) = 0 (gemessen).
- **A2c neu (K):** Kopplungssatz in der Form S′: (B) + (SB) ⇒ sup_N‖Ĵ_N⁻¹‖ < ∞. Schritt 1: n-te-Wurzel-Stabilität aus (B) + (SB) (Parametrix mit Fehler
  e^{−ηN}, ohne Ecke). Schritt 2: Eck-Korrektur als universelle O(1)-Grösse (Widom-Zerlegung, Bulk als äussere Parametrix). Fenster bis 12.10. abends.
  Abbruch Schritt 1: wenn die Parametrix-Fehler in keiner Norm e^{−ηN} werden (Fehlskalierung §8.2) ⇒ als benannte Lücke führen.
- **Familie (K + V):** über 𝒮_N ∩ B_N mit Satz K4 (Residuenterm Rang ≤ 2 im Kopf ⇒ (H-iii)); offene Stelle Linsen-Abweichung der Polynome für N ≫ 2000.
  Vorschlag an V: (S_γ) für G_P = T_N[C_P + M_{eˣ/P'}] an Fast-Lösungen bei N = 1500/2000 und dicker Linse prüfen (dort wird die Linsen-Abweichung messbar).

## Anpassung K (11.10. 23:59), Linsen-Spannung entschieden (BERICHT §10)
- Für (S_γ) an der Familie kein Hindernis (Keimabhängigkeit, K1″(b, c) für jede Matrix der Familie; S = 0.1, N = 600 numerisch belegt).
- K-A2c-Familienfassung in Keimform (Hypothesen nur über J_N^P); benannte Lücke (Fam-koh): Datenkorrektur in (H-ii)/H1 polynomial für Fast-Lösungen.
- Weiter mit A2c Schritt 1 wie in der Anpassung 23:40. V (Vorrang): N = 1500/2000, S = 0.0166 und 0.1 — Norm, Linsen-Abweichung, Datenkorrektur für h = k_η.

## Anpassung K (12.10. 00:25), A2c Schritt 1 gescheitert — Wendung (BERICHT §11)
- **Schritt 1 („(B) + (SB) ⇒ n-te Wurzel, ohne Ecke“) gestrichen:** L ∈ γ macht die Ecke zum Teil der exponentiellen Ebene (§11.1).
- **A2c neu (K):** (E₀) Eck-Anregung exponentiell: G_N = e^{o(N)} auf γ ∩ U(L) und ψ(ξ) ≤ o(1) für ξ ≤ ξ₀ (glatte Daten); fern von L Lemma K6 mit grober
  Kantenband-Schranke. (E₁) gleichmässig: Fast-Orthogonalität der Köpfe (Satz K5) an der Ecke. Ansatz: im universellen Eck-Grenzobjekt (Köpfe in
  der Skala Index/n bei L), nicht global. Fenster bis 12.10. abends; Abbruch, wenn sich die Köpfe an der Ecke nicht skalenfrei beschreiben lassen.
- **Für V (wenn Zeit):** w-Köpfe ‖c_N‖·‖ℓ_NᵀC⁻¹‖ bis N = 2000 (aus der einen Inversen X_N ablesbar, k12_kopf.py) — sättigen sie bei 0.015?

## Anpassung K (12.10. 00:34), Eck-Kopfform: Fenster ohne beweisbare Teilaussage — zweite Ideenrunde (BERICHT §12, IDEEN „Runde 2“)
- Übersetzung K5 ↔ G15–G18/W §25.1(c) eingetragen; Gs Eckprofil (Kegelskala n^{−a}) als Ausgangspunkt übernommen; Spezialfunktion/Gleichung/exaktes
  Modell geprüft: negativ bzw. ohne Trennschärfe.
- **Zweite Ideenrunde nur zur Ecke:** K hat R1–R12 eingetragen (Abgrenzung zu allen früheren Ideen, je erster Test, Bewertung). → V: ≥ 5 eigene Ideen,
  Bewertung von R1–R12. Danach gemeinsame Auswahl von 2–3 Ideen mit Fenster und Abbruchkriterium (wie Phase b).

## Runde 2 (Ecke) — Vorschlag V zur Auswahl (12.10. 00:49), mit K abzustimmen (K vorläufig: R12, R2, R5)
- **Stand 00:55:** VE1 in der starken Form abgebrochen (‖D_N‖ sättigt bei 0.215 statt zu fallen, N = 200…400); Keil-Linse bleibt als Standard für
  R12/R2-Tests (sauberes Grenzobjekt). Nächste V-Ideen: VE2 und VE4. **00:58: VE2 abgebrochen** (‖Ĵ_N‖ beschränkt, aber ‖Ĵ_N⁻¹‖ ~ e^{0.42N}; IDEEN VE2). Weiter mit VE4.
- **VE1 Keil-Linse (V):** Linse mit Ecke an L (Steigung c < 0.79). Messen: ‖D_N‖, ‖Ĵ_N⁻¹‖, Lage der Eck-Singulärvektoren für N = 100…800, c ∈ {0.25, 0.5, 0.7},
  kubisch und Kneser. Zeitfenster bis 12.10. 08:00. **Abbruch**, wenn ‖D_N‖ von N = 200 bis 600 nicht fällt (bei keiner Steigung) oder ‖Ĵ_N⁻¹‖ wächst.
  Erfolgskriterium: ‖D_N‖ fällt wie N^{−β} oder schneller ⇒ die Ecke ist geometrisch entfernbar; dann Satz V2 mit Keil-Linse neu rechnen.
- **VE4 Faktorisierter Schur-Test (V):** LU von Ĵ_N (N = 200…600), Diagonalgewicht W optimieren, ‖U⁻¹W⁻¹‖ und ‖WL̂⁻¹‖ gegen N. Fenster bis 12.10. 14:00.
  **Abbruch**, wenn für kein W beide Faktoren beschränkt bleiben (Wachstum > N^{0.5} über N = 200…600).
- **VE2 Doppelkeil / Einbuchtung über Γout (V, nach VE1):** Profil < 0 über Γout mit f(γ_out) ⊂ D̄. Messen: ‖Ĵ_N‖ (Fehlskalierung weg?) und ‖Ĵ_N⁻¹‖.
  Fenster bis 12.10. 12:00. **Abbruch**, wenn ‖Ĵ_N⁻¹‖ für jede zulässige Einbuchtung exponentiell wächst.
- **R12 Quadranten-Wiener–Hopf (K)** und **R2 Renormierung (K, billiger Struktur-Test)** wie von K gewählt; V liefert dafür Köpfe/Inversen aus v26
  (N' = 3N bis N = 500) auf Anfrage.
- R5, R7 zurückgestellt (R7 ist schnell prüfbar und kommt dran, falls VE1 abbricht).
- **Abgestimmter Stand 01:02 (nach Ks Tests 00:47/00:55):** R2 besteht den ersten Test (K weiter), R5 verworfen (K), R12 zurückgestellt (K), R1 neu
  bei K. **Auswahl Runde 2: R2 und R1 (K), VE4 (V).** VE4-Fenster bis 12.10. 14:00, Abbruch wie oben; Bindeglied zu Ks Cotlar–Stein (beide messen
  Fast-Orthogonalität der Köpfe, VE4 über die LU-Faktoren mit Diagonalgewicht).

## Anpassung K (12.10. 00:55), R2 harte Stelle, R1 aufgenommen (BERICHT §14)
- R2: analytisches ℛ über G22 unmöglich (Kern e^{0.14n}, Auslöschung); Eigenwertstruktur bleibt. Nächster R2-Schritt nur über eine stabile Darstellung.
- Neues Ziel der polynomialen Ecke: Cotlar–Stein für die w-Köpfe (A₁, A₂ beschränkt bzw. A₁ ≲ log N reicht für K1).
- R1 (zweiter Strang): Winkelgleichung der Quadranten-Ecke aufstellen; Abbruch, wenn ihr erster Exponent nicht a = 1 − log μ/(2πi) ist.

## Auswahl (H1) — K, 12.10. 01:22 (KERNPROBLEM.md; mit Vs Vorschlag 01:02 abgestimmt, Vs Bestätigung von H-A steht aus)
| # | Angriff | Verantwortlich | Fenster | Abbruchkriterium |
|---|---|---|---|---|
| H-A | Ecke in Abel/Koenigs-Koordinaten: Produkt primaler × adjungierter Kopf als Erhaltungsgrösse auf dem Torus ℂ/(½ℤ + (2πi/L)ℤ) | K | bis 12.10. 04:00 | Produkt nahe L zeigt keine ½-Periodizität (Abweichung > 30 % nach Entfernen der Gewichte φ, 1/α'∘exp bei n = 300–600) |
| VE4 | gewichtete Dreiecksfaktoren U⁻¹W⁻¹, WL̂⁻¹ bis N = 600, Schur-optimiertes W | V | nach Vs Plan | ‖U⁻¹W⁻¹‖ wächst wie √N (Ks Vorhersage aus der Kohärenz der primalen Köpfe, BERICHT §15) |
| R2′ | Renormierung als Translation (Teil von H-A), auf Vs Keil-Linse | K | nur wenn H-A positiv | Eigenwerte nicht als Torus-Moden erklärbar |

## Anpassung K (12.10. 01:27): H-A erster Test negativ (BERICHT §16)
- H-A („u·v erhalten auf dem Torus“) nach Abbruchkriterium beendet: die Köpfe sind in komplementären Sektoren homogen (primal innen, adjungiert aussen).
- Neu als Fassung von (H1): Transmissionsproblem an der Grenzfläche (Kreis durch L) — gleiche Stelle wie Os Eckformel/Ws starke Ebene. Weiter nur in
  Abstimmung mit V (VE4 läuft) und der Leitung; R2′ entfällt mit H-A.

## Entscheidung (H1)-Weg — K, 12.10. 01:30 (Auftrag Leitung; Runde mit V per BRETT, V rechnet noch)
**Prüfung des Stokes-/Vergleichsansatzes (a):** Das Vergleichsproblem existiert bereits als Formulierung: Os Grenzproblem (P1)–(P3) auf der Kegelskala
z = (t − L)/(ρN^{−1/κ₁}) (o_operator/BERICHT §33.4): τ = (μ/L)[y + μ^{−2}y(·/μ)] − 𝔲, 𝔲 = Spiegelung der Vorwärtsmode ũ mit μũ + ũ(μ·) = õ, Wachstum aus Ws
Kegel. Mellin macht (P1) diagonal, 1 + μ^{s−2} = 2πμ^{(s−2)/2}/[Γ(½ + iα(s−2))Γ(½ − iα(s−2))], α = log μ/(2π) — ein Faktorisierungsproblem vom
Malyuzhinets-Typ mit EXPLIZITEN Gammafunktionen. Es ist genau das Transmissionsproblem aus BERICHT §16: y = adjungierte Seite, ũ = primale Seite,
gekoppelt über die Spiegelung am Kreis durch L. **Reicht eine Verbindungsformel für (H1)?** Ja, wenn sie gleichmässig ist UND die äusseren Lösungen
(Ws geschweisste Fläche) mit kontrolliertem Fehler anschliesst — denn τ_n = ‖Q_n‖‖r̂_n‖/h_n enthält den globalen Pivot h_n. (a) ist damit die lokale
Hälfte von (b); (b) ohne (a) hat keine Eckparametrix, (a) ohne (b) keine Amplitude.
**Entscheidung: (a) zuerst, als „explizite lokale Parametrix“** — beweisbare Teilaussage in Tagen: Existenz/Eindeutigkeit und explizite Lösung des
Vergleichsproblems (P1)–(P2) mit Wachstumsbedingung, numerisch identifiziert mit Os Ψ_Y/Ψ_T und Gs S_∞. (b) danach (Wochen); (c) entfällt, weil (a)
die Kopplung innen/aussen schon explizit macht.
| # | Schritt | Fenster | Abbruchkriterium |
|---|---|---|---|
| a1 | (P1)–(P2) aus Ws Identität und der primalen Kopfgleichung in Kopfform herleiten; mit den Köpfen (k13, bis 600) auf der Kegelskala prüfen | bis 03:30 | (P1) oder (P2) gilt mit den Köpfen nicht auf ≤ 10 % |
| a2 | Mellin + Malyuzhinets-Faktorisierung: Streifen, Lösungsraum, Zahl der freien Parameter; Vergleich mit Ψ_Y-Taylor (o39) und S_∞ (G) | bis 07:00 | Lösungsraum nicht endlichdimensional ohne Ws äussere Daten, oder erste drei Taylor-Koeffizienten weichen > 5 % ab |
| a3 | Satzform: Existenz/Eindeutigkeit des Vergleichsproblems; was die Verbindungsformel für (H1) noch braucht | danach | — |

## Anpassung K (12.10. 08:37), zweiter K-Lauf: Grenzaufgabe gelöst; Spirale begrenzt den K1-Rahmen (BERICHT §18–§21)
- **a2/a3 erledigt:** Grenzprofile beider Köpfe explizit (F_{β_A}, F_{β_P}, λ_P = −μ^{−1/2}λ_A); Satz K10 (Existenz), K11/K11′ (Eindeutigkeit unter Mellin-
  Randbedingung) bewiesen, von V für §18.1/18.2 gegengelesen; §19 bei V.
- **Neu (§20–§21):** komplexes κ₁ ⇒ Eckmuster dreht sich um 1.45°·ln N gegen die feste Geometrie. Tangentiale Linse versagt ab N ~ 10³–10⁴ (p = 2 gemessen ab
  600), fester Keil φ bis N_×(φ) (×1.95 je Grad), mitgeführter Keil φ_N ≈ 1.45°ln N + 4° bis N ≲ 10¹⁶ (Familienbudget φ < 59°); primale Kreisnorm linsenunabhängig
  bis N ≈ 10⁸. K1″(a)/(b) gelten mit Ecke.
| # | Schritt | Wer | Fenster | Abbruchkriterium |
|---|---|---|---|---|
| L1 | τ_N bei N = 1200/1500/2000: Ks Linse, Vs Linse, Keil 15° | V | läuft | — (Vorhersage: Ks Minimum dann Anstieg; Keil flach) |
| L2 | Unabhängige Prüfung der Drehung: kubisch mit |μ| ≈ 2 (Im a ≈ 0.11), Profil-Formel mit α = log μ/2π, Umkehr von τ_N bei n ≤ 1000 | K | 1 Block | Profil-Formel trifft nicht (q_k-Test > 1 %) ⇒ Herleitung prüfen |
| L3 | „K1 mit Eck-Einsatz“: (H₀) mit expliziter Kopfentwicklung in U(L) statt Kurvenabschätzung; was an Fehlerschranken nötig ist | K | danach | keine Formulierung, in der die Amplituden δ_n(L)Λ_n(L) nur über B_n eingehen |

## Anpassung K (12.10. 08:56): L3 Ergebnis (BERICHT §22)
- Lemma K13 (Kopfentwicklung mit kopfweisen Konturen) bewiesen; Satz K14: (H₀) für alle N mit polynomialer Konstante ⇐ (i) B_n, (ii) Anpassung, (iii) Stegdämpfung,
  **(iv) primale Lokalisierung auf D̄** — (iv) nach dem Grenzprofil ab n ≈ 10⁸ verletzt (Mechanismus bei |μ| = 2 am echten Kopf bestätigt).
| # | Schritt | Wer | Abbruch/Entscheidung |
|---|---|---|---|
| T1 | Konvergenztest: Abschnittslösung zu singulären Daten auf ∂D nahe L in einem Modell mit schneller Spirale (|μ| ≈ 1.6–2) und dicker Keil-Linse; wächst sie, sobald die primale Γout-Hälfte ins Wachstum dreht? | V (Vorschlag) | wächst sie ⇒ gleichmässige Konvergenz auf D̄ ist für grosse N falsch, Vermutung umformulieren; konvergiert sie ⇒ Auslöschung finden (K) |
| T2 | (i) B_n: Exponent n^{−0.053} vs. Konstante — Daten n ≤ 2000 (V hat Köpfe) | V | — |
| L3′ | Lemma K12 + Barnes-Asymptotik (K10 iii) ausschreiben | K | nach T1 |

## Anpassung K (12.10. 09:07): Zielaussage für den asymptotischen Teil (BERICHT §23)
- Neue Zielaussage (robust gegen die Kette §22.6): **p_N → f_K gleichmässig auf K_{r₀} = D̄ ∖ U_{r₀}(L, L̄)** (r₀ = 10⁻³ mit Rand) — liefert Koeffizienten und p_N(1).
- Linearer Teil: Satz K15 (bedingt) ⇐ (i″) Eckprodukt, (ii) kopfweise adjungierte Eckform, (iii) Stegdämpfung, (v) ruhige Phase auf K_{r₀}.
- Nichtlinearer Teil: H1/H2 mit Kugeln in sup_{K_{r₀}}; offen: Eckbedingung an Fast-Lösungen P in der Kegelzone (für die Sekanten-Daten).
- „Gleichmässig auf D̄“ bleibt Prüfpunkt (T1, V); bis zur Entscheidung keine Arbeit an gleichmässigen Randaussagen.

## Anpassung K (12.10. 09:15): Eckprodukt = Zuwachs eines Rayleigh-Quotienten (BERICHT §24)
| # | Schritt | Wer | Abbruch/Entscheidung |
|---|---|---|---|
| E1 | (R_L): |ẽ_L*J̃_n⁻¹ẽ_L| ≤ C n^β beweisen — Fisher–Hartwig/Wiener–Hopf an der Eckfrequenz, Ziel n^{a}-Gesetz; Konstante C nicht nötig | K | keine Darstellung, in der nur Eckgrössen eingehen ⇒ als benannte Hypothese führen |
| E2 | (EB_N) an den echten p_N messen (Kegelscheibe, Ruhefenster), N ≤ 2000 | V | η_N nicht ≪ r_N ⇒ Formulierung ändern |
| E3 | Gs B_n-Daten bis 2000 gegen B_n ∝ Q_n n^{−0.053} prüfen (V hat Köpfe) | V | — |

## Welche bisherigen Bausteine braucht die K15-Architektur noch? — K, 12.10. 09:18 (Auftrag Leitung; BERICHT §22–§24)
| Baustein | Status in K15 | Begründung |
|---|---|---|
| Zweig-Zertifikat ≤ 400, Startbedingung (A) (V) | **gebraucht** | endlicher Teil; in der K_{r₀}-Norm automatisch erfüllt (sup über K_{r₀} ≤ sup über D̄) |
| (K2) Kneser-Schwanz, Satz P-V, Sprungschranke C_J u^{5.447} auf dem Schlitz für u ≤ u₀ | **gebraucht** | B-Max (Schlitzintegral) und Eckpaarung der Kneser-Daten |
| Satz K5 (Köpfe/LU), Lemma K13 (kopfweise Konturen) | **gebraucht** | Gerüst von K15 |
| Linse über Γin FERN von L (Stegdämpfung (iii), Vs Schleifenbedingung s ≥ Λ, Reserve-Messungen C10) | **gebraucht** | ferner Teil der kopfweisen Konturen; nahe L ist die Linsenform gleichgültig (Eck-Einsatz) |
| Linsenform nahe L (tangential vs. Keil), Keil-Linse, N_×-Fragen | **nicht gebraucht** | der Eckteil jeder Kontur ist kopfweise |
| K1′, K1″(a)–(c), konforme Koordinate, (S_γ), (H1), (H*_G), VE4, Cotlar–Stein, Bessel | **nicht gebraucht** | K15 arbeitet kopfweise mit Konturintegralen, keine gleichmässige Operatornorm |
| Satz V2 (explizites B für ‖Ĵ⁻¹‖ ab N₀) mit Keil | **nicht gebraucht** | ersetzt durch die K15-Konstanten unten; **V soll keine Zeit mehr in V2 mit Keil stecken** |
| Green-/Bernstein–Walsh-Konstanten (g, Lipschitz, Abstände) | **gebraucht, einfacher** | für K_far = K′ ∪ ferne Linse; der Eckwinkel bei L ist der von K′ (141.69°), also volle (LP)-Rate 5.447α′ = 4.29 ohne Keil-Abschlag |
| H1/H2 (Fast-Lösungen, σ-Homotopie), X4a/X5 (Familie) | **gebraucht, in sup_{K_{r₀}}** | plus Eckbedingung (EB_N) (§24.5) statt sup_{D̄}-Kugeln |
**Explizite Konstanten für N₀ ≤ 400 (K15-Kette):** (1) (R_L): C_R, β mit |ẽ_L*J̃_n⁻¹ẽ_L| ≤ C_Rn^β (gemessen 0.5985·n^{0.894}); (2) Eck-Einsatz (ii): C_ii mit
sup_{Γ_n^Ecke}|Λ_n| ≤ C_ii n^{β₂}|Λ_n(L)|; (3) Stegdämpfung (iii): C_iii mit sup_{ferner Teil}|Λ_n| ≤ C_iii n^{β₃}|Λ_n(L)|; (4) innerer Primalabfall (v):
C_v mit sup_{K_{r₀}}|δ_n| ≤ C_v n^{β₅}|δ_n(L)|; (5) Eckpaarung der Kneser-Daten (Mellin, explizit bis auf C_J); (6) Green-Konstanten von K_far, M = max|m|, Abstände;
(7) H1-Konstanten K₀, A_θ und die (EB_N)-Schranke η_N. Bedarf: Gesamtexponent κ = 1 + Σβ mit κ − 4.69 < −1 (Familie) und C'·400^{κ−4.69}·K₀ ≥ Startfehler —
alle (1)–(5) sind für n ≤ 2000 messbar (V hat Köpfe/Inversen), der asymptotische Beweis muss sie für n ≥ N₀ liefern.
**Prioritäten V (Leitung 09:16):** T1 → E2 (EB_N) → E3 (B_n) → Konstanten (1)–(5) für n ≤ 2000 tabellieren; V2 mit Keil entfällt.

## Anpassung K (12.10. 09:24): E1 Zwischenstand (BERICHT §25)
- C lokal hergeleitet (Anpassung der Eck-Pol-Antwort an den Bulk-Pol), Übereinstimmung 1e-5; (i″) damit lokal. K15 ruht auf (ii), (R_L), (iii), (v).
- E1 weiter: (R_L) als angepasste Asymptotik der einen Antwort u_n = J_n⁻¹T_nk_L; kritisch die Köpfe m ≈ n (Mellin–Parseval der Paarung statt Beträge) — oder exakte String-Gleichung
  n·P_n − a·K_n = κ₀ + O(1/n). Fenster: nächster Block. Abbruch: wenn die Paarung F_A × Residuum keine Auslöschung zeigt (dann (R_L) als benannte Hypothese führen).

## Anpassung K (12.10. 09:39): E1 als Reduktion abgeschlossen (k_kern/BERICHT §26–§27)
- String-Gleichung: nicht exakt; n·P_n − a·K_n = −aK₀∞ + O(n^{a−1}), K₀∞ ≈ 0.70 ± 0.05 global (Summe früher Köpfe), keine geschlossene Form.
- **Satz K16 (Bootstrap, bewiesen bedingt):** für jede Vergleichslösung gilt exakt K_n − M_n = Σ_m P_mρ_{m,n}; mit (ρ1), (ρ2) (relative Paarungen klein, ε_n ↓ 0) und
  Induktionsannahme „glatte Zuwächse“ folgt P_n = p_n(1 + O(ε_n)) ⇒ (R_L) samt Gesetz. Lemma ρ ⇐ (ii) + (iii) + (K-Param) (polynomiale primale Parametrix der Eck-Pol-
  Antwort, Typ (v)).
- **Damit ruht K15 nur noch auf (ii), (iii), (v)** — drei Kopfaussagen derselben Art (innen explizit, aussen WKB, Kleinheit auf dem Anpassungsring). E1 als eigener
  Restpunkt entfällt.
- **Konstantenliste N₀ ≤ 400, Punkt (1) neu:** statt C_R, β: n₀, B = max_{m<n₀}|P_m| (arb, k46), ε_{n₀}, ḡ, G aus Lemma ρ; Bedingung ε_{n₀}ḡ ≤ ¼, ε_{n₀}(2G + Bn₀^{0.106}ḡ/c_p) ≤ ½.
- K10 (iii) im Potenzsektor jetzt ohne Braaksma (Integraldarstellung F_β = Γ(β)⁻¹∫t^{β−1}e^{−t−ut^{−iα}}dt + Watson; k48).
- Nächster K-Block: (ii) auf der Kegelskala als Anpassungslemma formulieren und die innere Hälfte beweisen (Inverse des inneren Modells explizit über Mellin);
  die äussere Hälfte als benannte Hypothese an Ws Aussenlösung.

## Anpassung K (12.10. 10:01): eine Hypothese (AH) „Eckdominanz“ (k_kern/BERICHT §28, k49)
- (i″), (ii), (iii), (v) sind Faktoren EINER Grösse: Kernzuwachs ΔK_n(z,t) = δ_n(z)Λ_n(t)/ℓ_nn (Rang 1). Lemma K13 = Σ_n ∮ΔK_n·W. Maximumprinzip ⇒ zwei Kurven genügen.
- Fensterränder sind n-UNABHÄNGIGE Log-Spiralen ⇒ eine feste Kontur Γ* durch L (Σ₁ bei u_A ≈ 60–80° aussen, Σ₂ bei u_A ≈ −90° in D, θ = u_A + 71.46° + 1.623°ln(r/ρ)).
- **(AH): Q_n = max_{z∈K_{r₀}∪{L}, t∈Γ*}|ΔK_n(z,t)| ≤ C n^κ.** Satz K17 (bedingt): (AH) ⇒ K15-Schluss. K16/K10/K11 nur noch Vorhersage, nicht logisch nötig.
- Gemessen (k49, kubisch, grobe ferne Linse): Q_n = 0.661 / 0.328 / 0.273 (n = 200/400/600); ab 600 Maximum in der Ecke, Q_n = 1.008|P_n|.
- **→ V (E4, Messprotokoll §28.5):** Q_n und Faktoren bis 2000, kubisch + Kneser, Vs Keil-Linse c = 0.25 als Γ*_fern; Vorhersage Q_n = 1.008|P_n| = 0.54 n^{−0.106}.
- Offen für einen späteren Beweis: starke Asymptotik beider Kopffamilien aussen mit polynomialem Fehler; H_A ≤ 0 auf Γ* (Ws H), H_P ≤ 0 auf K_{r₀} (primal, von W nicht berechnet).

## Anpassung K (12.10. 10:23): HAUPTSATZ_KR0.md steht; K17 zurückgezogen; exponentielle Ebene (k_kern/BERICHT §29, k50)
- **phase3/HAUPTSATZ_KR0.md**: Hauptsatz unter (H)^{r₀} (alte (H) mit D̄ → K_{r₀}); keine (EB_N); Bausteine mit Status; Konstanten N₀ = 400; 5 Prüfstellen → **V liest als Ganzes**.
- **K17 zurückgezogen** (§29.1): Pol-Quellen nahe L brauchen f(Kontur) ⊂ D̄; das Spiralfenster dreht heraus ⇒ (AH) ⇒ (H-ii)^{r₀} nur bis N ≈ 1e8–1e11 gezeigt.
  E4 (Q_n bis 2000) bleibt als Messung der Eckstruktur sinnvoll, ist aber nicht mehr „die“ Hypothese.
- Exponentielle Ebene (§29.2): Ecke im festen Abstand umgehbar (G < 0 in 130°–265°), Steg durch Daten gedämpft; Γout-Seite mit Taylor-Gewicht schlecht (G ≥ +0.03)
  ⇒ nächster K-Block (Auftrag 2): H_A, H_P mit Gewicht der Green-Funktion von K′/W̄ auf der Γout-Seite; prüfen, ob H_P aus Ws Randwertproblem/Dualität folgt.

## Anpassung K (12.10. 10:38): Frage 1 entschieden — exponentielle Ebene reicht global, polynomial nur in fester Eckumgebung (k_kern/BERICHT §30, k51, k52)
- Primale Reserve auf K_{r₀} ist polynomial: Inneres im Potenzsektor, δ_n(z)/δ_n(L) ~ n^{−3.05} (= Re(a·s₀^P), explizit); nur Randpunkte fern von L fallen exponentiell.
- Satz K19 (B): mit Daten J(f − q_N) (Approximation auf E ⊃ D̄ ∪ 𝒞 ∪ f(𝒞)) ist sup_{K_{r₀}}|Π_Nf − f| ≤ (1 + c·ΣQ_n)·dist_E(f, P_N) — keine abgeschnittenen Pole, kein Bucht-Weg.
- Fern (ausserhalb U_{r₁}(L)): strikte Reserve, adjungiert η ≈ 0.013 (Linse s = 0.03, Γout-Bogen −0.037), primal ≤ −0.0023 ⇒ n-te-Wurzel-Ebene genügt (W).
- Nah: (P), (A-lok), (P-lok) polynomial mit Budget 1 + β₁ + κ_P + κ_A < s_E ≈ 3 (gemessen 0.89) — lokal, explizite Profile.
- Familie: K19 gilt bei P = f_K; für Sekanten braucht man P auf 𝒞 ∖ K_{r₀} ⇒ HAUPTSATZ_KR0 bleibt bei (H)^{r₀}.
- Nächster K-Block (Auftrag 2): fern-Raten H_A (Γout-Bogen, Linse) und H_P (∂D) aus Ws Randwertproblem/Dualität; dann lokale Schranken (A-lok), (P-lok).

## Anpassung K (12.10. 10:48): Fernraten (k_kern/BERICHT §31, k51b) und offene Teilfragen
- **Adjungierte fern:** Reserve gemessen (Linse s = 0.03: −0.013, Γout-Bogen: −0.037); als Satz (FA) ⇐ W12′ (Gestalt bis auf Konstanten) + Modell. W13s Hindernis
  (absolutes Niveau) ist für die relative Rate unschädlich. Weg ohne W12′: Vergleichsprinzip mit Superlösung Φ (Ws Ansatz (2)).
- **Primale: KEINE n-te-Wurzel-Reserve** (Lemma K20, k51b): H_P ≡ 0 auf D̄, Inneres ~ g(z)·n^{−3.05}; Rand fern von L nur vorübergehend exponentiell. §30.4 korrigiert.
  Nötig: (G-P) = polynomiale Konvergenz der Eck-Pol-Antwort Π_n s auf K_{r₀} (Budget κ_P < ≈ 2.1, wahr −3.05).
- **Offene Teilfrage Familie (Auftrag Leitung 10:40):** K19 gilt bei P = f_K. Für Sekanten 𝒥_P braucht die Datenseite (f − q_N)∘P den Wert von P auf 𝒞 ∖ K_{r₀}
  (Linse über Γin, Spiralen an L), den die K_{r₀}-Kugel nicht kontrolliert. Später prüfen: hilft Keimabhängigkeit (T_N𝒥_P hängt nur vom Keim von P ab, §10.1) plus K4
  (Kohärenz der Fast-Lösungen) wie bei der Linse? Kandidat: Datenseite für die Familie über Π^P angewandt auf f mit Approximanten auf E, und P auf 𝒞 aus der
  a-priori-Darstellung P − f = λ(Π^Pf − f) mit z ∈ 𝒞 (braucht R_P auf 𝒞, also Primale auf der Linse — dort wächst sie).
- **→ V:** (1) K19 gegenlesen (k_kern/BERICHT §30.2); (2) Rate s_E von dist_E(f_K, P_N) an der Spiralecke messen (E ⊃ K′ ∪ Linse s = 0.03 ∪ Σ₁ ∪ f(Σ₁), Lücke um den Schlitz
  ≈ (25°, 124°)); Vorhersage s_E ≈ 3.0 … 3.6; (3) E4 umgestellt: Q_n mit z ∈ K_{r₀} (OHNE z = L) auf 𝒞 = Γout-Bogen ∪ Linse s = 0.03 ∪ Spiralen, bis 2000;
  (4) primal: ln|δ_n/δ_n(L)| an ∂D-Punkten 120°/150°/180° bis 2000 — wann kippt die Exponentialphase in die Potenz n^{−3.05}? (Vorhersage n ≈ 1500–2500).

## Anpassung K (12.10. 10:52): Zirkelprüfung (G-P) — Ergebnis (BERICHT §32, REFLEXION 10:52)
- (G-P) ist der Eckteil der alten Hypothese an der Quelle L (nicht wesentlich leichter); direkter Angriff zurückgestellt.
- Lineare Aussage auf K_{r₀}: (EP) [polynomiale Zuwächse der zwei Eck-Pol-Antworten u_n = Π_n s, w_n = K_n(L,·)] + (FA) [n-te Wurzel, ⇐ W12′] + s_E ⇒ Konvergenz.
- Einziger neuer Hebel für (EP): duale gewichtete Parametrix (Residuum gross erlaubt, wo (FA) Reserve gibt). Vorschlag Leitung: ob K das als nächsten Block nimmt
  oder zuerst (FA) als Satz (Vergleichsprinzip/Superlösung, Ws Ansatz (2)) angeht.

## Anpassung K (12.10. 11:00): (FA) per Vergleichsprinzip gescheitert (BERICHT §33)
- O1 (bewiesen): keine superharmonische Vergleichsfunktion in der Ebene (Laplace-Masse der Nullstellen). O2: Weitergabe über die λ-Naht braucht die Dominanzstruktur ⇒ ⇔ W12′.
- (FA) bleibt benannte Hypothese (n-te Wurzel, gemessen mit Reserve). Präzises Ziel für später: einseitige, in L normierte Fassung von W12′.
- V: Vortest (b) duale gewichtete Parametrix — Messvorschrift auf dem BRETT (11:00).
- Stand der linearen Aussage auf K_{r₀}: (EP) [polynomial: Eckprodukt, Eck-Pol-Antwort fern, Adjungierte auf Spiralen] + (FA) [n-te Wurzel] + s_E (V misst).

## Anpassung K (12.10. 11:05): HAUPTSATZ_KR0.md Fassung 2 steht (→ V gegenlesen)
- Annahmen: (EP) [(P), (G-P), (A-lok)] + (FA) + (S_E) bei f_K; Familie F1–F4 ausdrücklich offen; r₀ > r_* (Vs T1). Satz A (linear, B) = K19 + K13 + Rang 1.
- Liste „Neu bewiesen in Phase 3“ (§7 des Dokuments) mit V-Status für den Wissensgraphen.
- Nächster K-Block: falls Vs Vortest (b) Lemma ρ stützt → Parametrix für (EP); sonst Familie in der K_{r₀}-Fassung (Keimabhängigkeit + K4, F2/F3).

## Anpassung K (12.10. 11:16): Familie in K_{r₀} entschieden (BERICHT §34, k51c)
- Keimabhängigkeit + K4 reduzieren F2/F3 NICHT auf f_K. F3 nahe L ⇐ (EP⁺) (lokal, Skizze). F3 fern (Linse am Steg): primale Köpfe wachsen dort (+0.016/Schritt bei s = 0.03),
  P selbst ist dort exponentiell weit von f_K; Keim-Repräsentanten tragen bis N ≈ 6300–7800 ⇒ benannte Hypothese **(Fam-koh)^{r₀}** (Repräsentant P̃ mit P̃(Linse) ⊂ D).
- Zange (§34.3): keine Kontur in D̄ mit einer Datendarstellung. Die Familienlücke sitzt am Steg, nicht an der Ecke.
- → V (zusätzlich, niedrige Priorität): sup_{Linse}|p_N − f_K| bis 2000; Carathéodory–Fejér-Repräsentant (w-Koeffizienten von (p_N − f_K)∘k) — polynomial oder e^{NF₀}?
- HAUPTSATZ_KR0 Fassung 2: F3 präzisieren auf (Fam-koh)^{r₀} + (EP⁺) — mache ich beim nächsten Dokument-Update (nicht während Vs Gegenlesen).

## Anpassung K (12.10. 11:20): (S_E) — Satz K22 (BERICHT §35, k53)
- s_E = 5.417·α_E, α_E = (u₁ − u₂ − 38.3°)/180°: (60°,−90°) 3.36; (80°,−90°) 3.96; (95°,−110°) 5.02. Bedingt auf den zitierten direkten Satz für Quasischeiben (Belyi 1977,
  Andrievskii–Belyi–Dzjadyk 1995; nur zitiert) + Geometrie von E. k53: Exponent auf 0.01–0.02 bestätigt (synthetische Ecke).
- HAUPTSATZ_KR0 beim nächsten Update: E ⊃ D̄ ∖ U_{r₀} (statt D̄); (S_E) mit Status „bedingt (zitierter Satz)“; F3 → (Fam-koh)^{r₀} + (EP⁺).

## Anpassung K (12.10. 11:37): (S_E) — Satz gelesen, Geometrie geprüft (BERICHT §36, k54)
- APV Satz 2 (arXiv math/0005032) gelesen: E ∈ H* (Quasischeiben eingeschlossen), keine Bedingung an ω ⇒ Satz K22 ruht auf gelesenem Satz.
- Geometrie von E (kubisch): univalent, ohne Löcher, Schnitt in der Lücke, Öffnung 132° (Θ = 131.7°), Dreipunkt 1.07.
- → V: k54_geometrie_E.py mit f_K auf 𝒞 (Kneser-Motor) laufen lassen; dann ist (S_E) bis auf diesen endlichen Check fertig.
- HAUPTSATZ_KR0 beim nächsten Update: (S_E) Status „bewiesen bis auf f_K-Geometrie-Check (V)“, Quelle APV Satz 2 (gelesen).

## Anpassung K (12.10. 12:33): Vortest (b) entschieden — Lemma ρ trägt mit expliziten Parametrizes NICHT (BERICHT §37, k55/k55b/k55c/k51d)
- Taylor-Parametrix: Kompositionsüberlauf (ρ bis 1e6). E-Approximant (auch dual gewichtet): ρ = O(1). Ursache: im primalen Wachstumssektor (ausserhalb D nahe L) folgt u_n
  dem inneren Profil nur bis |ζ| ≈ 2–3, dann äussere WKB-Lösung (Faktor 22 bei |ζ| ≈ 29). Fenster 2 (in D): Profil auf ≤ 4 % bis |ζ| ≈ 12.
- **Korrektur:** Fenster 1 ist primal ein Wachstumsfenster ⇒ (EP⁺) entfällt, F3 nahe L nicht lokal; (Fam-koh)^{r₀} umfasst Linse UND Σ₁.
- Gemeinsamer Kern aller offenen Stellen: primale Grössen ausserhalb D auf der Γin-Seite (wo die Adjungierte Reserve hat).
- Vorschlag nächster Hebel (i): Abschneide-Rest R_n = (I − T_n)Ju_n auf dem Γout-Bogen messen (polynomial ⇒ Parametrix „D̄ explizit + Identität“ möglich; exponentiell ⇒ W12′-Gebiet).

## Anpassung K (12.10. 12:44): Hebel (i) entschieden — R_n auf dem Γout-Bogen exponentiell (BERICHT §38, k56/k56b)
- R_n = u_n∘f + m·u_n − T_nk_L (exakt). Rate je Grad bis +0.252 (φ = 0), ≈ ½·ln(|f|/ρ), über n = 100…400 auf 3 Stellen konstant; polynomial nur auf der Eckskala.
- Parametrix „D̄ explizit + Identität“ NICHT gebaut. **(EP) bleibt benannte Hypothese.** Nachtrag KERNPROBLEM.md §4: Kern als Auslöschungsaussage (KK).
- Offen für K (wartet auf Leitung): HAUPTSATZ_KR0-Update nach Vs Gegenlesen (E ⊃ D̄∖U_{r₀}; (S_E)-Status; F3 → (Fam-koh)^{r₀} für Linse und Σ₁; (EP⁺) gestrichen).

## Anpassung K (12.10. 12:51): HAUPTSATZ_KR0 Fassung 3 + Übergabenotiz (Auftrag Leitung nach 12:45)
- Fassung 3: (P) entfällt, (G-P) absolut (G_n = sup_{K_{r₀}}|Δu_n|), Budget γ = max(0, 1 + κ_G + κ_A) = 0 nach V §3.22; E ⊃ K_{r₀}, Öffnung 131.7°; (S_E) = gelesener Satz +
  num. Geometrie; F3 → F3′ = (Fam-koh)^{r₀} (Keim-Repräsentant, Linse und Σ₁); (EP⁺) gestrichen; (KK) als §2.6; Schlitzkonstante 8.2; Fassung 2 als Anhang B. → V: Prüfstellen §5.
- Übergabenotiz BERICHT §0.0 (Ziel, Lesereihenfolge, (KK), Gescheitertes mit Grund, Kandidaten, Werkzeuge, Korrekturen).
- Nächster K-Block (Zeitfenster ≈ 2 h): Produkt primal × adjungiert als ein Objekt; zuerst Unterschied zu H-A (§16) klären (Steg statt Ecke? nur Ratensumme?).

## Anpassung K (12.10. 13:08): Zeitfenster „Steg-Produkt“ entschieden nach ≈ 15 min (BERICHT §39; k57, k57b, k57c)
- Produkt primal × adjungiert am Steg polynomial (h_Π → 0) und explizit: ΔK_n(t,t)·m(t) → κ(t); K_n(t,t)·m(t) = t(log δ_n)′ + Δ(t) + o(1) (Christoffel-Darboux-artig, Vermutung).
- **Korrektur:** (KK) mit „Steg +0.003“ war ein Artefakt (Steigung + Niveau) ⇒ KERNPROBLEM §4 Nachtrag 2, HAUPTSATZ §2.6 → (KK′), BRETT markiert, Übergabe §0.0 nachgeführt.
- Mit Vs CF (−0.0024): Linsenfenster S ∈ (0.0126, 0.0147) → V testet. r₁ ≤ 0.1 für die Familie.
- Nächster K-Kandidat (wartet auf Leitung): exakte Christoffel-Darboux-Identität für J̃_n.

## Anpassung K (12.10. 13:19): Christoffel-Darboux-Identität (CD1) gefunden (BERICHT §40; k58)
- (CD1) (z − t)K_n(z,t) = q_{n+1}(z)λ_n(t) − tν_n(z) + Ξ_n(z,t), Ξ_n = ∮K_n(z,s)(f(s) − s)K_n(f(s),t)ds/s [bewiesen; numerisch 1e-176].
- Am Steg n-unabhängig: ν_n, Ξ_n(t,t), B(t) = m q_{n+1}λ_n/t (1.03968 bei 180°); umnormiert (Generator/B(t)) auch Ξ^reg ⇒ erklärt §39.3 vollständig.
- Offen: Beweis der n-Unabhängigkeit; Ξ ist der Kompositionskanal (= KK′). Kandidaten: Dyson-Kontraktion nahe der Fixpunkte; B(t) exakt.

## Anpassung K (12.10. 13:37): Dyson-Kontraktion für (CD1) gescheitert (BERICHT §41; k59, k59b, k59c)
- ρ(L) der linearisierten Fixpunkt-Abbildung: bestenfalls 0.61 (n = 10) → 1.02 (60) → 1.20 (100), also e^{0.004n}; normunabhängig ⇒ keine Kontur und kein Gewicht hilft.
- (CD1) allein kennt m nicht (jedes (C_f + M_m)⁻¹ löst sie); beide Kanäle zusammen ⇒ linear, nur X = J⁻¹. Daten der besten Variante = Kopfpaar = (EP)-Information.
- (CD1) bleibt als exakte Darstellung. Nächste Schritte: Leitung entscheidet.

## Anpassung K (12.10. 13:52): Ideenrunde 3 (IDEEN.md „Runde 3“, K3-1 … K3-6); Auswahl vorbehaltlich V (bis 13:52 keine V-Antwort auf dem BRETT)
- K3-1 (RH aus CD1 + CD2): global nein, (CD2) ≡ (CD1) mod Rang 3; R5 global bestätigt, lokal offen.
- Auswahl: **K3-2 integrabler Eckkern** (Test k60: Ξ auf der Eckskala gegen q_{n+1}λ_n; 2 h; Abbruch: Verhältnis bei n = 600 > 0.1 und nicht fallend) und
  **K3-3 Bootstrap (F3′) mit Produkt-Reserve** (K Gerüst 1.5 h + V gewichtetes CF; Abbruch: braucht Werte ausserhalb der Linse im Betrag). B(t) exakt zurückgestellt.

## Anpassung K (12.10. 14:01): K3-2 entschieden — gescheitert (BERICHT §42; k60). Nächster Block laut Leitung: K3-3 (Gerüst Bootstrap F3′, 1.5 h), wartet auf das Go nach dem Bericht.

## Anpassung K (12.10. 14:06): K3-3 Gerüst steht (BERICHT §43) — (F3′) ⇐ (BU) + Startprüfung N = 400 + F2; termweise nur auf S₀ ≈ 0.0125 mit Marge 0; robust nur kohärent (Phasen).
- → V: R_P^γ, R_A^γ auf der kritischen Linse bis N = 3000 (über den Übergang bei 2200); gewichtetes CF mit den echten π_n.
- Wartet auf Leitung.

## Anpassung K (12.10. 14:10): Konsolidierung (Auftrag Leitung) — Problemblatt 2 (de/en) und Abschlussnotiz Phase 3 stehen; keine neuen Versuche bis zur Entscheidung des Autors
- phase3/OFFENES_PROBLEM_2.md, phase3/OPEN_PROBLEM_2_en.md (je ≈ 100 Zeilen): Rahmen, bedingter Hauptsatz + T1, offene Aussagen (G-P), (A-lok), (FA), (S_E), (BU),
  Werkzeuge (CD1, Eckprofile, K19, Rang 1, Steg-Formel), Messwerte, gescheiterte Wege, 4 Fachfragen.
- k_kern/BERICHT §0.A: Abschlussnotiz Phase 3.
