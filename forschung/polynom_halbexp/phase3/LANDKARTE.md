# Landkarte des Beweises (Phase 3)

## A. Analytische Seite (Strang K)

Stand 11.10. 20:25 (aus `date`). Kennzeichnung: **bew.** / **Skizze** / **num.** / **Verm.**; Fundstellen wie in ARCHITEKTUR (W §x = w_grenzschicht/BERICHT.md usw.).
Notation: J_∞δ = δ∘f + mδ, m = f'∘f, J_N = T_NJ_∞|_{P_N}, Π_N = J_N⁻¹T_NJ_∞, Eckmoden J_N⁻¹ ≈ Σσ_iu_iv_iᵀ mit σ_i ~ e^{cN}, c ≈ 0.012,
adjungiert y_i = σ_iv_i, Y_i(t) = Σ_k y_{i,k}t^{−k}, T_i(ζ) = ⟨y_i, ρ_ζ⟩ (Eckanteil der Pol-Antwort).

### A1. Das Gerüst in einem Bild
```
(K) Kneser-Struktur ──► (F) f_K um K holomorph bis auf Schlitze, Sprung ~u^{5.447}           [bew.|(K); P-V zertifiziert (K) für b = e]
(F) + (H₀) ──► Π_N f_K → f_K  (lineare Konvergenz, Rate N^{κ−4.29…5.07} je nach Rahmen)     [bew., P-I / X1 / C13.1]
(F) + (H) ──► H1 (genau eine Lösung in B_N)  ──(A)──► H2 (Zweig = diese Lösung) ──► Vermutung [bew., V-verifiziert]
(H) = (H₀) an allen Fast-Lösungen P ∈ 𝒮_N ∩ B_N (Familie) + (H-i) + (H-iii)                   [offen]
(H₀) = eckfreier Teil [bew. aus (a″): X1 Schritt (4)–(5)] + Eckteil (Ẽ) bzw. (H*) [OFFEN — der Engpass]
```
Was an (H₀) offen ist, ist **genau eine Art von Aussage**: Die Eckmoden haben Norm e^{cN}, werden aber von Eingaben, die bei L sitzen, nur
**polynomial** angeregt. In y-Normierung: |T_i(ζ)| ≤ CN^κ für ζ nahe L (Schlitz) und e^{−εNg} weiter weg. Das ist eine **Auslöschungsaussage**
(‖y_i‖ = σ_i ~ e^{cN}, aber ⟨y_i, ρ_ζ⟩ = O(N^κ)). Jede der vorhandenen Formulierungen ist eine Fassung davon:

| Fassung | Aussage | ⇒ | Gültig für Kneser? |
|---|---|---|---|
| (H*) W1'' | \|Y_i\| ≤ CN^κ auf einer festen Kurve γ* (Linse ∪ h(A-Bogen)) | (E_K) | nur polynomiales f (braucht Eckformel D6; Kneser: I_F − B_y bzw. (E_sing)) |
| G19/G10 | \|𝒟_n\| ≤ C n^{−1/2} (oder nur ≤ C) auf γ* für alle n ≥ n₀ | G_N ≤ C'√N bzw. C'N ⇒ (H*) für alle u | Summenkriterium ja (endliche Identitäten); Schritt zu (E_K) wieder über W1'' |
| (Ẽ) X1/X2 | kanonischer Eckteil c_N(k_ζ) = Π_Nk_ζ − Π̂k_ζ (Rang ≤ d): \|c_N\|₁ ≤ CN^κe^{−εNg'} | (LP) | **ja** (nur Rationalität in ζ, P-II, (a″)) |
| (KB) + aussen | Kreisbedingung auf \|t−L\| = AN^{−0.894} + (H*)_loc + Wachstum (I) | (E_K) | wie (H*) |

**Mein Vorschlag für Teil C:** Zielaussage ist **(Ẽ) in der einteiligen Form X2** (Rc_N(L) ≤ CN^κ, Θ^c_N ≤ 1−ε) — die einzige Fassung, deren
Reduktion für Kneser ohne Zusatzlücke bewiesen ist; (H*)/G19 bleiben die Arbeitsform im kubischen und radialen Modell. Die nötige Stärke ist
gering: **jedes κ < 4.29** genügt (gemessen κ ≈ 0.5–1.1, Reserve gross); e^{o(N)} genügt nicht (W §26.3).

### A2. Bewiesen auf der analytischen Seite (Status nach Teil B)
- **Hauptkette** (K) ⇒ (F) ⇒ 4.1 ⇒ H1 ⇒ H2 — V-verifiziert. Sie benutzt **keine** Operatornorm von J_N⁻¹, nur kohärente Pol-Antworten.
- **Reduktionen von (H₀):** P-I/P-II/P-III (Rahmen, Rest in K', Eck-Kriterium p·α > 1), C13.1–13.4 (Randform, Schärfe), X1/X2/X3 (Rechteck-
  Splitting, eckfreier Teil aus (a″)), W1'' + W3 + W5 + W11 + W13, O1/O7/O9/O10/O12, G9/G10/G22/G23 (Grad-Zuwächse). Damit ist der **eckfreie
  Teil vollständig auf (a″) zurückgeführt**, und (a″)_∞ ist für Kneser rechnergestützt bewiesen bis auf (K2) (M; Prüfung durch V offen).
- **Modelle:** G27–G29 (Dreieck/affin, exakt), G41/G44 (Innenfixpunkt, gleichmässig — aber mit Kontraktionsreserve, die am Randfixpunkt fehlt),
  M1/M6 (radiales Quadrat: einzige reelle Familie mit radialem Randbild; dort ist J̃_N = P_N[Op(s(φ)^k) + M_m]P_N ein Pseudodifferential-
  operator ohne Winkeltransport, (a″) bewiesen), X Satz 1.5 (Möbius = Toeplitz × GL(2)-Darstellung, Krawtchouk).
- **Lokale Struktur bei L (Skizze/num., konsistent über drei Stränge):** Kegel mit Deckverschiebung 2πi − log μ, Exponent κ₁ = 1.118 − 0.032i
  (W §19), Innenschicht der Breite 1/N mit Vorfaktor √N (O §32), Anpassungsskala N^{−1/κ₁} inkl. logarithmischer Drehung (G Etappe 3), lokaler
  adjungierter Eckraum eindimensional (O §33.3). Neu notiert: Knesers Singulärexponenten bei L sind **genau die Mellin-Nullstellen** des
  lokalen Symbols (μ^{s−1} = −1 ⇔ s = 1 + (2k+1)πi/log μ; wegen log L = L ist log μ = L/2, also s₀ = 1 + 2πi/L = 1 + γ) — der singuläre Teil
  von f_K liegt lokal im Kern von J_∞; seine Anregung der Eckmoden kommt nur aus nichtlokalen Korrekturen (Einordnung, kein Satz).

### A3. Wo und warum es gescheitert ist — die Mechanismen
1. **Normschranken (Kantorovich, Kugeln, Stabilität, Gewichtsradien, Polynomgewichte, ℓ²-Extremalität).** Mechanismus: J_N⁻¹ hat echte
   Eckmoden mit σ ~ e^{cN}; ‖J_N‖ ~ e^{0.47N} (Komposition reicht bis |f| = 2.23 > ρ). Jede Aussage, die den Eingang nur über eine Norm sieht,
   verliert e^{cN} — und zwar **genau dort, wo Knesers Singularität sitzt** (M §5). Diagonale Gewichte (in k) sehen den Winkel nicht; die
   Information „Eingang sitzt bei arg L“ geht verloren. **Lehre:** Ein Beweis muss winkel- und indexaufgelöst sein (Phasenraum).
2. **Lokale Argumente (Koeffizienten-Induktion, Eckzone per Maximumprinzip, lokales Wiener–Hopf, quellenfreies Innenproblem).** Mechanismus:
   (a) X §2.4 (bew.): auf jedem Kompaktum in D gibt es exakte Konkurrenz-Lösungen f_c → f_K; die Auswahl geschieht **nur am Rand bei L**.
   (b) W3 (bew.): Bucht-Schranke ⇔ T-Schranke auf dem A-Bogen — der Kreisschluss ist intrinsisch. (c) W §24.2: die Eckzone erbt am Aussenrand
   winkelaufgelöste, polynomial genaue Daten. (d) O §33.5: die Amplitude der Eckmode (Eckwert) ist nicht lokal bestimmt. **Lehre:** Die
   Eckaussage ist eine Aussage über die globale Lösung, ausgewertet an der Ecke.
3. **Potentialtheorie / n-te Wurzel (W11, W5, W12, Nullstellenverteilung).** Mechanismus: liefert nur e^{o(N)}-Genauigkeit; Knesers Singularität
   sitzt **in** L mit nur polynomialem Abfall, also verlangt schon die Existenz (Stufe (i)) polynomiale Eckgenauigkeit (W §26.3). Zudem ist das
   Grenzsystem verschiebungsinvariant (W13): es bestimmt die Gestalt, nicht die Normierung; die Normierung entsteht in L (W §27.5).
   **Lehre:** Man braucht **starke Asymptotik** (polynomial genau) in einer festen Umgebung von L.
4. **Klassische starke Asymptotik (Gonchar–Rakhmanov, Szegő, Deift–Zhou).** Mechanismus: keine exakte Orthogonalität vom Typ {s^jW}: jede
   Bedingung koppelt ein durch f transportiertes Moment in a₀ mit einem Moment in 0 bei gleichem Index (W §25.2); die geklebte Fläche ℳ, auf
   der es GR-Form hätte, existiert nur im Grenzwert. Kein festes Toeplitz-Symbol (G39: zwei Kanäle 1 − sF und 1 − sz). **Lehre:** Es fehlt ein
   Werkzeug — eine RH-/Parametrix-Theorie für eine **Schweiss-Aufgabe mit nichtlokaler Verschiebung** h (Fixpunkte L, L̄).
5. **Exakt lösbare Modelle.** Mechanismus (M §2): in allen geprüften Familien schliessen sich exakte Lösbarkeit und Eckmode aus (affin: kein
   Abschneide-Effekt; Pol im Fixpunkt: keine Abschneidekante; G44: Innenfixpunkt mit Kontraktionsreserve). **Lehre:** Ein Modell, das die Ecke
   trägt, muss Gradanstieg (nicht-affine Komposition), Randfixpunkt und Multiplikationsterm zugleich haben — und ist dann nicht trivial lösbar.

### A4. Wo neue Mathematik nötig ist
- **(N1) Finite Abschnitte von Operatoren C_f + M_m mit Fixpunkt auf dem Abschneidekreis.** Es gibt keine Theorie, die exponentiell wachsende
  Randzustände (Eckmoden) solcher Abschnitte beschreibt. Bekannt sind nur: Stabilität ⇔ Invertierbarkeit lokaler Grenzoperatoren (Roch–Santos–
  Silbermann) — hier verletzt, aber mit **exponentieller**, nicht polynomialer Folge. Gebraucht wird nicht die Wachstumsrate selbst, sondern eine
  **Anregungsschranke für eingangsseitig lokalisierte Daten**.
- **(N2) Starke Asymptotik biorthogonaler Polynome einer Zwei-Kanal-Bilinearform mit Komposition** (Q_n∘f + mQ_n = O(xⁿ), A_n = h_nJ_n^{−T}e_n).
  Vorhandene Bausteine: g-Funktion aus dem freien Randwertproblem (W §20, c = 0.0118 ohne Fit), Rang-2-Verschiebung (G40), Vierkanal-Transport mit
  zyklischer Monodromie (G8), Schweiss-Formulierung (W §27.1), Innenschicht (O §32–33), komplexe Skala (G Etappe 3). Fehlt: der Rahmen, in dem
  „Parametrix + kleiner Rest ⇒ Schranke“ ein Satz ist (G43 gibt die abstrakte Form).
- **(N3) Lokale Parametrix in L mit Quelle**, Malyuzhinets-Faktorisierung von 1 + μ^{s−2} (O §33.4), angepasst über Skalen, die sich nur wie
  N^{0.106} trennen.
- **(N4) Störungstheorie dieser Asymptotik** in f (für (H₀) → (H) an Fast-Lösungen und für kubisch → Kneser): die Ecke reagiert polynomial auf
  glatte, kaum auf kohärente Störungen (X §3.8.2), aber ohne (N2) ist das nicht beweisbar.

### A5. Was auf der analytischen Seite sonst noch zu bauen ist (ohne neue Theorie)
- **Rechteck-Satz (Vorschlag K):** der eckfreie Teil ist bewiesen (X1 Schritt 4–5) und die Familie dafür über X8 (A_θ) + X5 (Schwanz) +
  (a″)_∞ gleichmässig kontrollierbar, weil σ_min eines **Rechteck**-Abschnitts unter kleinen Störungen stabil ist. Das sollte einen
  **vollständigen** Satz für das Kleinste-Quadrate-Verfahren (N+1 Unbekannte, N+d Gleichungen) samt Zweig ergeben, bis auf (K2) und
  Gitterkonstanten — und zeigt sauber: Quadratisches und rechteckiges Verfahren unterscheiden sich **nur** um den Rang-d-Eckteil c_N.
- **Kante (a″)_{αN}:** ist eine Eckaussage (M §15.2) — gehört zu (H₀), kein eigener Baustein.
- **Familie (H₀) → (H):** im eckfreien Teil erledigt (X4a, B); im Eckteil braucht es (N4). Mit X1 ist die Familie für c_N eine Aussage über
  denselben Rang-d-Operator an benachbarten Linearisierungen.
- **Kneser-Eckformel:** W1'' braucht D6; für Kneser über I_F − B_y (G Runde 2) oder (E_sing) (P-IV). Wird überflüssig, wenn man (Ẽ) als Ziel nimmt.

### A6. Der Engpass und meine Einschätzung
**Der Engpass ist (Ẽ) — polynomiale Anregung der Eckmoden durch Daten bei L — und er ist eine Frage fehlender Technik, nicht fehlender
Evidenz** (numerisch mit Reserve ε ≥ 0.6, κ ≈ 0.5–1 bei Schwelle 4.29, bis N = 1000, Kneser = kubisch auf 2–3 Stellen). Alles andere ist
entweder bewiesen, endlich (V) oder hängt an derselben Technik (Familie, Kante). Drei Wege sind denkbar:
1. **Die Technik (N2)+(N3) bauen** — direkt am echten Problem zu teuer; zuerst in einem **Modell mit Randfixpunkt, das die Ecke trägt**
   (radiales Quadrat: Pseudodifferentialstruktur ohne Winkeltransport; loxodromisches Möbius-Modell: Krawtchouk-Asymptotik bekannt).
2. **Eine Umformulierung finden, in der die Auslöschung strukturell wird** (z. B. einseitige analytische Vorkonditionierung — T_N vertauscht mit
   Linksmultiplikation durch analytische ω, anders als bei der verworfenen beidseitigen Agmon-Konjugation; Grad-Zuwächse 𝒟_n mit Rang-2-Struktur).
3. **Schwächere echte Sätze** sichern: Rechteck-Satz (A5), bedingter Hauptsatz mit (Ẽ) als einziger benannter linearer Hypothese.
Die Ideenliste (IDEEN.md) bewertet diese Wege einzeln.

## B. Geprüfte und endliche Seite (Strang V)

Stand 11.10. 20:01 (aus `date`). Vollständige Bausteintabelle mit Konstanten und Abhängigkeiten: `v_pruefstand/ABHAENGIGKEIT.md`.
Statuswörter: **V** verifiziert (ein zweiter Strang hat Schritt für Schritt nachvollzogen) · **S** selbstgeprüft (rigoros gerechnet, aber nur vom
Autor gelesen) · **B** behauptet bewiesen, niemand hat gegengelesen · **num** nur endlich viele N · **offen**.

### B1. Was bewiesen und verifiziert ist (Fundament)
- **Die Logik der Hauptkette** (K) ⇒ (F) ⇒ 4.1 ⇒ H1 ⇒ H2 ⇒ Vermutung ist **V** (P, PRUEFLISTE §1; Identitäten am echten Problem auf 1e-50).
  Bedingt auf (K), (H) = (H₀) + Familie, (A) und drei numerische Geometrieaussagen (b35b′, Radius r > 0.02 der Eckformel, Univalenz von f_K auf D̄).
- **Reduktionen V:** P-I (rahmenunabhängig, Rate N^{κ−5.447α_K}), P-II, X1 (Kleinste-Quadrate-Eckteil), X3.1/3.2 (Rechteck-Splitting), C13.x,
  D8, L7.1–7.3, W1'' (EINE Hypothese (H*) auf EINER Kurve γ* genügt für den Eckteil — **nur für polynomiales f**), W3, W5, O1, O2 (Struktur),
  O7a, O7g–i, O9, O10 (λ_L = 0.98988321).
- **Nur S (P hat es selbst gerechnet, niemand hat es gegengelesen):** Satz P-V (Kneser zertifiziert, a_0…a_640), p32a ((K3): f_K(L) = L,
  f_K'(L) = √L), Zweig-Zertifikat 1 → 193 (p18 + p24, auf aicgolling reproduziert), Kopplung (A) für N = 140…193 (p31, bis auf den Schwanz
  S₆₄₀), sowie P-III/P-IV (Ps eigene Sätze). M hat die *Ausgabe* von P-V numerisch gegengeprüft (1.1e-43 gegen A, Funktionalgleichung bis 2e-138)
  und einen Formfehler gefunden: die Radien der Koeffizientendatei enthalten die Rundung der 140-stelligen Dezimalmitten nicht (bis ~1e-139·|a_n|).
  Für alle Anwendungen bedeutungslos, für die Strenge zu beheben.

### B2. Was nur behauptet ist (B — gehört in die Prüfschlange)
- **M:** (a″)_∞ für Kneser, σ_min(J_∞|P_N) ≥ 0.6553 für alle N (rechnergestützt, arb) — bis auf (K2); Beweisskripte m22–m26 und Lemmata M5,
  M7–M10, M13 ungeprüft. Radial: M6/M11/M12.
- **X:** X4a (eckfreier Teil der Familie gleichmässig), X5 (Schwanz, Konstante c₁ nur Skizze), X7 (Eingangsnorm 6.3√M), X8/X8′ (A_θ = 13.4),
  Kontur d = 0.3 (Konstante 17.1; Holomorphie nur Skizze; Sprungschranke nur bis u = 0.2 statt 0.46).
- **O:** O12 (Eckformel in L), SVD-freie Majorante G(x) = ‖J̃⁻¹ẽ_x‖ ≥ |Y_u(x)| für alle Einheits-u (elementar, aber ungeprüft).
- **G (extern):** G8–G44, darunter G9/G10 (Rang-1-Zuwachs und Summenkriterium: |𝒟_n| ≤ C/√n auf γ* ∀n ⇒ G_N ≤ C′√N), G11 (Kreis-Majorante),
  G22/G23 (funktionale Schur-Rekursion), G44 (gelöstes Spielzeugmodell F = 2z + ηz²). Einzelgrad-Intervallzertifikate nur **kubisch**, nicht Kneser.

### B3. Wo es gescheitert ist und warum (Fehler und Sackgassen, mit Lehre)
| Was | Warum es scheitert | Lehre |
|---|---|---|
| Gleichmässige Stabilität sup‖J_N⁻¹‖ < ∞ | Eckmoden bei L, σ₁ ~ e^{0.012N} | Nur kohärente Eingaben tragen; jede Normschranke von J_N⁻¹ ist wertlos |
| Kantorovich, Kugelargumente, Newton-Kugeln | Zufallsrichtungen ≳ 1.17^N verstärkt; Kopfantworten auf P(D̄) ~1.3^N (p22) | Der Beweis darf nur entlang kohärenter Richtungen laufen |
| (S)/(H*) „polynomial im ganzen Band“ | numerisch falsch ab N ≈ 420–650 | Nur lokale Pol-Antwort auf Schlitz ∪ Γ_d ist wahr |
| (a′) ⇒ (E_R), (E_L) ⇒ (E_K) (Phase 1) | A-Bogen liegt in f_K(D): Eingaben ~(ρ/|h|)^N, Urbildterme unkontrolliert | Rahmen K' bzw. gemischter Weg / X1 |
| Lokale Abschätzungen, Koeffizienten-Induktion (W §22–24) | Bucht-Schranke ⇔ T-Schranke (W3); Auslöschung zwischen Bulk und Randterm | Die Eckaussage ist intrinsisch global |
| Gonchar–Rakhmanov in der Ebene (W §25) | Komposition koppelt Momente in 0 und f(0); keine {s^jW}-Struktur | GR-Form nur auf der Grenzfläche ℳ |
| n-te-Wurzel-Asymptotik allein (W §26/27) | Kneser-Singularität sitzt IN L, Daten fallen nur polynomial; e^{o(N)} lässt den Schlitzbeitrag wachsen | Es braucht starke (polynomial genaue) Asymptotik an der Ecke |
| Extremalität in ℓ²_ρ, andere Gewichtsradien | Rayleigh 10^90; jedes ρ' ≠ ρ wächst exponentiell | ρ ist erzwungen |
| Lokales Wiener–Hopf-Modell (O §32.4), quellenfreies Grenzproblem | Skalen trennen nur wie N^{0.106}; Quelle 𝔲 so gross wie y | Innen- und Aussenlösung müssen gemeinsam gebaut werden |
| Betragsmajoranten | Auslöschung ~10^180 bei n = 1000 (G: ‖Q_1000‖ ≈ 2.5e40, Produkt klein) | Nur ausgewertete Funktionen dürfen majorisiert werden |
| Störung von (a″)_∞ für kleines α (Kante) | Überlaufblock wächst 1.48^N (radial) / 1.62^N (Kneser) | (a″)_{αN} ist eine Eckaussage (K_d), keine Störung |
| Datenfehler | flint acb_series kappt still bei 10 Termen (D-Tabellen falsch); float-Horner in D ist Rauschen (C „Halbmond-Divergenz“, P „+0.031“); W1-Exponent 1.27 zurückgezogen | Präzisionsregeln; jede Zahl gegen zweite Methode |
| (K2) über Kugeln im Strukturgebiet (p32) | Kugelauswertung über q-Kugeln überschätzt ≈ 4e7 je Einheit Radius | zentrierte Punkt-Jets nötig |

### B4. Was noch gebaut werden muss (endlich, rechnerisch erreichbar) — Verantwortung V
1. **(K2) — erledigt (V, 21:13):** |b_k| ≤ 461.36·k^{−6} für alle k ≥ 1 (rechnergestützt, repariertes p32 = v06/v07; benutzt p32a und
   Ps rest_schranke, die noch gegenzulesen sind). Datei für M/O: v_pruefstand/kneser_halbexp_b_zert_V.json. Vorher gemessen: Ms 𝒯-Weg verträgt
   C₆ ≤ 5e5 (Fenster 6/14); mit C₆ = 461 läuft Ms Gesamtlauf jetzt auf aicgolling (6/14 und 7/23).
2. **Zweig-Zertifikat 193 → 400**: läuft seit 19:55 auf aicgolling (12 Arbeiter, Skripte = Xs Kopie von p24, eigener Ordner
   `v_pruefstand/zert400/`, Erwartung ≈ 2.5 h); danach Kettenprüfung und Kopplung p31 für N₀ ∈ [318, 400].
3. **Zertifizierte Randkonstanten** aus den P-V-Daten (heute alle Gitterwerte): c_in, c_out = min|f_K'| auf ∂D, max|m|, max|m'|, max|m/f'|,
   M2, C_J, b35b′, Univalenz, die Sprungschranke bis u ≈ 0.46 und die Holomorphie innerhalb Γ_{0.3}. Sie gehen alle in N₀ ein.
4. **Koeffizientendatei**: Dezimalrundung in die Radien (oder Nachweis, dass jeder strenge Verbraucher sie addiert; mk_kneser_arb tut es).

### B5. Was noch geprüft werden muss (Reihenfolge nach Gewicht für den Beweis)
1. P-V Schritte 3–6 (p27 Taylor-Modelle/Defekte, p28 Beltrami/Stoilow/Paulsen–Cowgill-Prüfungen, p29, kc_rig) und p32a — ohne zweiten Leser
   steht der ganze Kneser-Datensatz auf einer Person.
2. Ms 𝒯-Weg für Kneser (m22–m26, M5, M7–M10, M13): einzige Quelle für c₀ = (a″)_∞, geht in X4a und N₀ ein.
3. zert_kern/p24 (Röhren-Krawczyk: gerundete Knoten, Lebesgue-Konstante, Majoranten) und p31 (Kopplung).
4. X4a/X5/X7/X8′ und Kontur d = 0.3 (Grundlage von N₀ ≈ 318).
5. Gs Sätze für den Eckweg: G9/G10 (Summenkriterium), G11 (Kreis-Majorante), G22/G23; G44 als Modell.

### B6. Gegenstimme: Lücken in der Buchhaltung, die leicht übersehen werden
1. **(A) ist kein unabhängiger endlicher Test.** N₀ hängt an den Konstanten von (H); die Eck-Konstanten C_E, κ_E sind Teil von (H₀) und nur
   gemessen. „N₀ ≈ 318 beweisbar“ heisst: *falls* (H₀) mit κ_E ≤ 2, C_E ≤ 1 bewiesen wird. Das Zertifikat bis 400 sichert genau dieses Szenario
   ab, nicht mehr. Ein (H₀)-Beweis mit schlechterer Konstante (z. B. κ_E = 2.5, C_E = 1 ⇒ N₀ ≈ 2243) bräuchte ein viel längeres Zertifikat
   oder schärfere Konstanten (A_θ → A* ≈ 0.63).
2. **Kneser ≠ kubisch an der entscheidenden Stelle.** W1'' (die Reduktion des Eckteils auf (H*) auf γ*) ist nur für polynomiales f bewiesen;
   für f_K fehlt die Eckformel D6 (mit I_F − B_y bzw. (E_sing)). Gs Summenkriterium G10 benutzt dagegen nur endliche Matrixidentitäten und gilt
   für Kneser direkt — aber G liefert G_N auf γ*, und der Schritt G_N ⇒ (H₀) geht wieder über W1''.
3. **Familie (H₀) → (H):** nur der eckfreie Teil ist (B, X4a) gleichmässig; der Eckteil steht auf der Vermutung (Fam-Eck). Ausserdem ist
   (H-iii) nur numerisch. Ein Beweis von (H₀) allein beweist die Vermutung nicht.
4. **Zwei konkurrierende Zerlegungen** des Eckteils sind im Umlauf (SVD-Eckmoden Q_N bei W/O/G vs. kanonischer Eckteil c_N bei X1/X3);
   ihre Bedingungen sind ähnlich, aber nicht identisch formuliert. Für den Beweis muss EINE Zielaussage fixiert werden (Teil C).
5. **Viele „rigorose“ Teile hängen an Gitterwerten** (B4.3); sie sind billig zu zertifizieren, aber bis dahin ist N₀ nicht streng.

## C. Engpass und Folgerungen (gemeinsam)

*Entwurf K (11.10. 20:08), von V zu ergänzen oder zu korrigieren.*

**C1. Engpass (Konsens aus A6 und B6).** Offen ist im Kern nur der **Eckteil von (H₀) für alle N ≥ N₀** mit expliziten Konstanten — polynomiale
Anregung der Eckmoden durch Daten bei L. Familie (Fam-Eck), Kante (a″)_{αN} und Kopf-Stabilität (H-iii) hängen an derselben Ecke: der Kopf
J_N⁻¹x^N ist die letzte Spalte Q_N/h_N, und Q_N liegt zu 0.99 im Eckraum (W §25.1) — (H-iii) ist also auch eine Eckaussage, allerdings mit
exponentieller Reserve (‖Q_N‖/h_N ~ e^{−0.065N}, G Etappe 2).

**C2. Eine Zielaussage.** Für den Gesamtbeweis gilt **(Ẽ) in der einteiligen Form X2**: Rc_N(L) := |c_N(k_L)|₁ ≤ C_E N^{κ_E} und
Θ^c_N ≤ 1 − ε für den kanonischen Eckteil c_N = Π_N − Π̂ (d = d(N) bzw. 2N), Rahmen K'. Begründung: Reduktion X1/X2 ist V-verifiziert und braucht
für Kneser keine Eckformel D6. Arbeitsformen in Modellen dürfen (H*) (W1'') oder G19 (Grad-Zuwächse) sein; ein Modellbeweis zählt für Kneser erst,
wenn er in (Ẽ) übersetzt ist.

**C3. Konstanten-Kopplung.** Die Zielaussage muss so stark sein, dass N₀(H) ≤ 400 bleibt (B6.1): bei C_E ≤ 1 heisst das κ_E ≤ 2. Ein qualitativer
Beweis (unbekannte C, κ < 4.29) beweist die Vermutung erst zusammen mit einem Zertifikat bis zum dann entstehenden N₀ — oder mit einer
Konstantenverschärfung (A_θ → A* ≈ 0.63, X).

**C4. Arbeitsteilung, die daraus folgt.** K: (Ẽ) bzw. ein Modell mit Randfixpunkt, in dem die Eck-Anregung bewiesen wird; Rechteck-Satz als
gesicherter Zwischensatz. V: (K2) mit C₆ ≤ 5e5, Zertifikat bis 400, zertifizierte Randkonstanten, Prüfschlange B5 (zuerst P-V und Ms 𝒯-Weg, weil
alle eckfreien Teile daran hängen), schnelle Experimente zu Ks Vermutungen, Gegenlesen jedes K-Satzes.

*Ergänzungen und Korrekturen V (11.10. 20:14), nach kritischem Lesen von Teil A:*

**C5. Die Wahl von d entscheidet, wie viele Eckaussagen es sind (wichtig für C2).** Satz X2 gilt für jedes d; aber (Xs Bericht §3.5c, x13/x15b):
- **d = 2N (X3, keine Kante):** (a″)_{2N} bzw. (a″)_{4N} folgt aus (a″)_∞ + X5. Dafür ist Θ^c_N **asymptotisch grenzwertig**: gemessen −0.039/−0.012/
  +0.020/+0.072 (N = 200/300/400/600, Trend +0.028 je 100 Grade), Grenzwert aus Ws Profil Θ_∞ ≈ 0.82 (mit 20 % Aufschlag 0.99) — Reserve 0.18
  bzw. 0.01. Ursache: der grosse Eckteil c^{(2N)} klingt auf der f_K(Γout)-Seite nicht ab, und der Γin-Exponent p ≈ 1.4–1.5 liegt nur knapp über 1/α' = 1.27.
- **d = d(N) klein (X2):** Θ_∞ ≈ 0.07–0.20 (Reserve 5–14), aber dann braucht (a″)_{d(N)} die **Kante** — und die ist nach M §15.2 selbst eine
  Eckaussage (Überlaufblock wächst 1.62^N). Es sind dann **zwei** Eckaussagen: Kante + Θ^c_klein.
- Dazwischen: die Kante ist eine Störung von (a″)_∞ erst ab α ≥ 2.74 (Kneser; mit der Sprungschranke bis u ≈ 0.46 ab α = 2), Θ^c hat dort
  das grenzwertige d = 2N-Verhalten.
⇒ „(Ẽ) in X2-Form“ ist nur dann EINE Zielaussage, wenn man d = 2N nimmt und die knappe Reserve akzeptiert (oder einen besseren Rahmen als K'
für den grossen Eckteil findet). Vorschlag V: Zielaussage = **(Ẽ_{2N}): Rc_N(L) ≤ C_EN^{κ_E}, Θ^{c,(2N)}_N ≤ 1 − ε** in K'; als Rückfall
(Ẽ_{d(N)}) + Kante. Ein Beweis muss sagen, welche Variante er liefert.

**C6. Die nötige Stärke von κ ist kleiner als 4.29.** κ < 4.29 genügt für (LP) (P-I in K'). Für den Hauptsatz mit Familie braucht X4a
θ_N = A_θ·N·R_N → 0, also **κ < 4.29 − 1 = 3.29** (X §3.7: „strenger als in (H)“), und für N₀ ≤ 400 explizit κ_E ≤ 2 bei C_E ≤ 1 (C3).
Die gemessenen Werte (κ ≈ 0.5–1.1) erfüllen beides; die Formulierung in A1/A6 („jedes κ < 4.29“) sollte das sagen.

**C7. Statuskorrekturen zu A1/A2** (nach den Statuswörtern von Teil B): „eckfreier Teil bew. aus (a″)“ ist bew.|(a″), (Univ); (a″)_∞ selbst ist
**B** (M, ungeprüft) und hängt an (K2). G9/G10/G22/G23, O12, W11, W13 sind **B** (niemand hat gegengelesen; P hat W1'', W3, W5 geprüft, nicht
W11/W13). Bis zur Prüfung dürfen sie Arbeitsgrundlage sein, aber nicht Fundament.

**C8. Was V als Erstes absichert, damit C2 trägt:** (i) (K2) — mit einer Präzisierung: die nötige Schranke hängt am Fenster von Ms Lauf. Mit
T₁ = 6, T₂ = 14 (Ergebnis c₀ ≥ 0.5646) genügt nach Ms Abschätzung C₆ ≤ 5e5; mit T₁ = 7, T₂ = 23 (c₀ ≥ 0.6553) liegt das Fensterende bei
|z − L| ≈ e^{−23}·|L − L̄| ≈ 3e-10, und **schon C₆ = 5e4 lässt m24 scheitern** (v_pruefstand/m_kopie/v01, v02: „Streifenzweig“). Xs Bilanz
rechnet mit c₀ = 0.54, also reicht die schwächere Variante. (ii) Ms 𝒯-Weg prüfen (einzige Quelle für c₀). (iii) X4a/X5/X8′ prüfen.
(iv) Zertifikat bis 400 + Kopplung.

*Ergänzung K (11.10. 20:38), nach Satz K1 (k_kern/BERICHT §1) — C2 damit ersetzt (von V verifiziert und bestätigt 20:52, PLAN.md):*

**C9. Neue Zielaussage (H*_G) — eine Funktion auf einer Kurve, ohne d, ohne Eckformel.** Satz K1: ist G_N(t) = ‖J̃_N⁻¹ẽ_t‖ ≤ C_G N^κ auf
γ = Γout ∪ ∂ℬ ∪ konj. (ℬ eine in L, L̄ tangentiale Linse über Γin), so gilt (H₀^{K_γ}) mit K_γ = K' ∪ ℬ̄ und Exponent κ + 3/2 (grob), also (LP)
mit Rate N^{κ−2.79}. Der Beweis ist eine Cauchy-Darstellung; er gilt für f_K wörtlich (schliesst B6.2), braucht kein (a′)_r, keine Eckmoden, kein d
(erledigt C5: es gibt keine Poisson–Jensen-Bedingung mit exponentiellen Randdaten mehr, also keine grenzwertige Reserve — der Rahmen ist „vergrössertes
Kompaktum mit polynomialen Randwerten“ wie Ws W1, aber in K'). Konstanten-Bilanz (zu C6): Familie über X4a braucht κ + 3/2 − 4.29 < −1, also κ < 1.79;
N₀ ≤ 400 braucht entsprechend kleine C_G. Gemessen: κ ≈ 0.45, Maximum in L ≈ 0.39√N (O), Kneser = kubisch auf 0.01 (k01).
Mit Gs G10 genügt **sup_{t∈γ}|𝒟_n(t)| ≤ C_Δ für alle n ≥ n₀** (κ = 1). ⇒ Vorschlag: **C2 := (H*_G)** (bzw. die 𝒟_n-Form); (Ẽ_{2N}) als Rückfall.
Numerisch ist G_N auf γ ausserhalb einer Umgebung von L sogar konstant in N; der ganze offene Inhalt sitzt in γ ∩ U(L) (= Os (KB)).

**C9 (V, 11.10. 20:52) — nach Satz K1 (von V verifiziert):** Zielaussage für den Eckteil ist **(H*_G): sup_{t∈γ} G_N(t) ≤ C_G N^κ** auf der festen
Kurve γ = Γout ∪ Linse(Γin) ∪ konj., für Kneser ohne Eckformel gültig. Damit sind C5 (Wahl von d, Kante) und B6.2 (Kneser-Lücke) gegenstandslos.
Nötige Stärke: κ < 2.79 für (LP), κ < 1.79 mit Familie (X4a), κ + 3/2 ≤ 2 bei kleiner Konstante für N₀ ≤ 400. Offen bleiben: (H*_G) selbst
(Engpass), ihre Familienfassung für T_N𝒥_P, (H-iii), und der endliche Unterbau (A5 im Plan).

**C10 (V, 11.10. 22:02) — Reserve-Fenster der Linse für (S_γ)** (v_pruefstand/BERICHT §3.2). Untergrenze: das Wachstumsprofil Λ(τ) der Eckmoden auf
Γin (Λ(π) = c ≈ 0.012, Λ ~ τ^{1.5} bei L; eigene Messung bis N = 1500, bei 4°–15° Grenzwert auf Faktor 1.6 unsicher). Obergrenze (nur für Kneser,
aus der K1-Geometrie f(ℬ) ⊂ D): U(τ) ≈ 0.79τ bei L, ≥ 0.5 ab 30°. Die Operatorseite setzt keine Obergrenze (dicke Linsen bis S = 0.3 stabil).
**Reserve U/Λ ≈ 14 im ungünstigsten Winkel (τ ≈ 6°), sonst 20–50, bei L unbeschränkt.** Ein Beweis über Wachstumsraten (n-te Wurzel) muss das
Eckprofil also nur bis auf einen Faktor ≈ 14 kennen. Ks Linse S = 0.0166 liegt bei 4°–8° möglicherweise an der Untergrenze; als Bezugslinse
empfiehlt V S = 0.05 (Faktor 3–5 Abstand nach unten, 5–10 nach oben). Nachtrag 22:44: Ks Linse hält bei N = 2000 unverändert (1.251155) —
das konservative Λ war bei 4°–8° zu gross, Ks Linse ist zulässig.
Einschränkung (Ks Einwand 21:55, von V akzeptiert): die Reserve deckt nur den Steg fern von L; in der Eckzone um L (Skala N^{−0.894}) sind Linse und
Λ beide ≈ 0, die Raten verschwinden — dort bleibt polynomiale Genauigkeit nötig. Was (S_γ) gewinnt: keine exponentiell wachsenden Moden mehr; der
Rest ist ein reines Eckproblem mit Raten 0.

**C11 (V, 11.10. 22:48) — Familie in w-Form: Lipschitz-Konstante von P ↦ (Ĵ_N^P)⁻¹** (v_pruefstand/BERICHT §3.3, REFLEXION V Versuch 6).
Kubisches Modell, Sekantenoperatoren T_N𝒥_P, P = f + η, ε = sup_{D̄}|η|. Die naive Neumann-Reihe in ℓ² scheitert (‖Ĵ⁻¹ΔĴ‖, ‖ΔĴĴ⁻¹‖ ~ ε·e^{0.18N},
weil Ĵ Einträge bis 1.6^N hat). Die wahre Änderung der Inversen ist aber linear in ε (Glied erster Ordnung, Reihe konvergiert mit Quotient ≈ 10ε).
Lip_N = ‖(Ĵ^P)⁻¹ − Ĵ⁻¹‖/ε: glatte und auf L gerichtete Störungen ≈ 0.045N (Ableitungsverlust); Störungen hoher Grade ≈ e^{0.010…0.013·N},
**unabhängig von der Linsendicke** (S = 0.0166…0.1); Gradbänder: exponentiell nur für k ≥ 0.8N, für 0.3N ≤ k < 0.8N beschränkt. Mit dem echten Profil von p_N − f_K (Zertifikatsboxen; ≈ 7 % des ℓ¹-Gewichts im Band
k ≥ 0.9N, Abfall ≈ N^{−5.5}) ist ε·Lip_N ≈ 1e-14 bei N = 400 und ≈ 1 erst bei N ≈ 5000. ⇒ Für die Familie (Nr. 3 oben) genügt eine Störung um
f_K nicht für alle N; nötig ist eine gewichtete Norm (δ auf f(Int γ)) und (S_γ) direkt für Ĵ^{p_N} oder die Struktur von p_N − f_K (Lip_N sättigt bis N = 1200 nicht).

**C12 (V, 11.10. 23:07) — (S_γ) am echten p_N; keine Grössen-Norm für die Familie; Eck-Korrektur universell** (v_pruefstand/BERICHT §3.4–3.6).
(1) Ĵ_N^{p_N} (Tangente, p_N'∘p_N exakt) hat bis N = 1000 dieselbe Inversen-Norm wie Ĵ^{f_K} (1.250337 Ks Linse, 1.107257 bei S = 0.05);
Abstand 0.015·N·ε_N, ε_N ~ N^{−4.8}. (2) Einzelgrad-Antworten der Kneser-Tangente: r_k ≈ e^{N g(k/N)}, g(1) = 0.161, unabhängig von der Linse —
keine diagonale Gewichtsnorm und auch nicht die Graphennorm erfüllt zugleich „Lipschitz polynomial“ und „p_N − f_K klein“. Getragen wird die
Kleinheit von der Kohärenz (Schwanz von f_K hebt η'∘f auf; Richtungen J⁻¹(glatt) haben Lipschitz linear in N). ⇒ Familie über die Lösungsmenge
(H1) führen. (3) Ks Eck-Korrektur D_N ist für kubisch und Kneser auf 3–4 Stellen gleich (Paare L/L̄, ξ ≈ 0.85–0.92); fern von L fällt sie wie
e^{−SN}; an L wird sie mit dickerer Linse kleiner (0.49 → 0.33 → 0.24), ohne bis S = 0.1 zu sättigen. (4) Kerben-Test (blind vorhergesagt, §3.7): eine Kerbe
unter Λ(τ₀) löst Wachstum mit Rate ≈ Λ(τ₀) − s_min aus (15°: 0.0038; 60°: 0.0058 ⇒ Λ(60°) ≈ 0.0115), eine schmale Kerbe (±4°) erst ab N ≈ 1100 ⇒
die Schleifenbedingung s ≥ Λ ist asymptotisch punktweise scharf.

**C13 (K, 11.10. 23:40) — Ecke nicht abtrennbar; Schleifenbedingung = Gestalt + Eck-Normierung; Familie über Fast-Lösungen** (k_kern/BERICHT §9).
(1) D_N ist universell (μ und Linse bei L; kubisch = Kneser auf 4e-4 im führenden Paar), aber ein eigenständiges Eckmodell mit Gradabschnitt gibt es nicht: das affine
Modell ist exponentiell instabil (e^{0.70N}, J_aff nicht injektiv). ⇒ Vermutung S′: (B) Bulk inkl. Operator-Ecke + (SB) Schleife ⇒ (S_γ).
(2) Satz K3: Ws Tropik liefert (SB) in keiner Koordinate (die Linse wirkt nur in der Auswertung); (SB) ⇐ W12′ + u(L) = 0.
(3) Satz K4: P'∘P = (eˣ + R_P')/P' — auf 𝒮_N ist der Residuenterm ein Rang-2-Kopfterm; die Newton-Matrix an p_N, Π_N f_K ist so stabil wie J_N,
an T_N f_K exponentiell instabil (2.3e25 bei N = 600). Grössen-Kugeln um f_K tragen die Familie nicht; 𝒮_N ∩ B_N (Architektur H1) schon.

**C14 (K, 11.10. 23:59) — Linsen-Spannung der Familie: kein Hindernis für (S_γ)** (k_kern/BERICHT §10). Ĵ_N^P hängt nur vom Keim von P ab; K1″(b, c) gelten für
jede Familienmatrix. Kneser, S = 0.1: p_N weicht auf der Linse um 1.6e2 (N = 400) bzw. 8.1e9 (N = 600) ab, ‖(Ĵ^{p_N})⁻¹‖ = 0.965153 wie bei f_K.
P auf γ nur in Daten-Schlüssen (H-ii)/H1: Kosten e^{NF₀}·Kopf (Carathéodory–Fejér), inkohärent ab N ≈ 6300–7800; Rest (Fam-koh) = Kohärenz.

**C15 (K, 12.10. 00:25) — A2c Schritt 1 „ohne Ecke“ unmöglich; Kopf-Zerlegung** (k_kern/BERICHT §11). L ∈ γ ⇒ n-te-Wurzel-(S_γ) enthält G_N(L) = e^{o(N)}
(linsenfrei). Satz K5: J̃_∞ = LU, Ĵ = (C_kL)U mit k-unabhängigem U; Ĵ_N⁻¹ = Partialsummen der Kopf-Entwicklung, D_N = projizierte künftige Köpfe;
n-te Wurzel ⇔ Köpfe subexponentiell. Lemma K6: Linsenpunkte über das Grad-Antwortprofil ψ; ψ = 0 für ξ ≤ 0.8, exponentiell nur im Kantenband.

**C14 (V, 12.10. 01:40) — (S_γ) an Fast-Lösungen bis N = 2000; explizite Schranke B_max; Datenkorrektur klein** (v_pruefstand/BERICHT §3.8–3.10).
(1) Tangente an das echte p_N (hochgenau per Newton-Sprung, zwei unabhängige p_2000 auf 2e-206 gleich): ‖Ĵ⁻¹‖ = 1.250337 / 1.107257 / 0.965153 für
S = 0.0166 / 0.05 / 0.1 bei N = 1500 und 2000 — gleich f_K, auch bei der dicken Linse, wo das Polynom auf der Linse weit von f_K abweicht.
Sekanten und alle Zwischenpunkte der H1-λ-Homotopie und der H2-σ-Homotopie bis N = 1000 (1500 Sekante): derselbe Wert. ⇒ Ks Keimabhängigkeit bestätigt:
die Linsen-Abweichung ist für (S_γ) an der Familie harmlos. (2) Datenkorrektur Π^{p_N}k_η − Π_Nk_η ≤ 2ε_N, relativ fallend (0.0027ε_N bei N = 2000).
(3) Satz V2: mit Ks Linse genügt B ≤ 30.8 ab N₀ = 318 bzw. 76.7 ab 400 (Geometrie streng, Green-Konstanten numerisch); gemessen 1.25.

**C15 (V, 12.10. 02:17) — Satz V2 mit strengen Green-Konstanten: B ≤ 17.9 ab N₀ = 400** (v_pruefstand/BERICHT §3.11). Vergleichsgebiet K* = D₁ ∪ D₂
(Kreise durch L, L̄) ⊃ K_γ streng geprüft ⇒ g_γ ≥ g_{K*} explizit (Sektor-Öffnung α = 0.756); Lipschitz über innere Scheiben ≤ 0.861; |f_K| auf der
Kontur ≤ 5.03 über f = exp∘f∘log, log∘f∘exp und exp∘exp∘f∘log∘log. Ab N₀ = 318 bleibt nur B ≤ 1.36 — das Zertifikat bis 400 trägt.


**C16 (V, 12.10. 08:38) — Linse mit Ecke trägt, tangentiale nicht (numerisch)** (v_pruefstand/BERICHT §3.12–3.13). τ_N = ‖c_N‖‖ℓ_NᵀC⁻¹‖: Ks Linse Minimum 0.0149 bei
N ≈ 1000, dann Anstieg bis 0.0204 (2000); Keil c = 0.5 fällt auf ≈ 4e-4 und sättigt. VE4-Produkt: Ks Linse 1.33 → 2.98 (wachsend ab N ≈ 900), Keil 1.332 konstant
bis 2000 (kubisch = Kneser). ⇒ Für (H1)/(S_γ) für alle N eine feste Linse mit Ecke an L (Steigung 0 < c < 0.79); Satz V2 muss für den Keil neu gerechnet werden.

**C17 (V, 12.10. 09:00) — Sprungschranke streng: |Δf_K(L(1+u))| ≤ 1.0378·u^{5.4469} für u ≤ 0.48 (bewiesen bis auf Ps Zweigwahl)** (v_pruefstand/BERICHT §3.15,
`v44_sprung_zert.py`). Aus Ps oberer Struktur (p32a) und kc_rig; numerisch 0.94 → 0.62 fallend in u. Damit Satz V2 für Ks Linse: B ≤ 16.4 ab N₀ = 400 (streng bis auf X-Lemmata, (Fam)).

**C18 (V, 12.10. 10:25) — (R_L) und Eckamplitude bis N = 2000 (numerisch)** (BERICHT §3.17, `v48_RL.py`, `v50_E3_norm.py`). Eine LU liefert K_n(L,L) für alle n:
|K_n|/n^{0.894} = 0.5990 (kubisch) / 0.5991 (Kneser) bei 2000, C = P_n·n·d_n/(aL) → −0.0101 − 0.5983i (5 Stellen stabil); B_n = Q_n·0.6262·n^{−0.053} auf 4 Stellen,
Q_n → Q_∞ ≈ 1.307 wie n^{−0.89}; Kneser = kubisch auf 0.3 %.

**C19 (V, 12.10. 10:25) — (EB_N) an den echten p_N bis 2000 (numerisch)** (BERICHT §3.18, `v49_E2_EB.py`). Auf der Kegelscheibe |x − L| ≤ 16r_N: in D̄ ≈ 2N^{−4.8},
ausserhalb η_N = 7.4e-8 (N = 2000), η_N/r_N = 4.8e-5 fallend; Ruhefenster in D̄ ≈ N^{−4.8}; schlitznah ausserhalb divergiert p_N ab r ≈ 0.03 (Überkonvergenzrand).

**C20 (V, 12.10. 10:25) — T1: primales Γout-Wachstum ist in sauberen Modellen real, seine Wirkung auf Π_N U bisher klein** (BERICHT §3.14, `v45_t1_modell.py`,
`v47_t1_kopf.py`). Ks |μ| = 2-Kubik ist global fast-singulär (f(Γin) ⊄ D). Saubere Quintiken mit gleicher Ecke: M3 (|μ| = 2): Kopf auf Γout +0.5 Dekaden, Π_N U − U
fällt überall bis N = 2400; M4 (|μ| = 2.6): Kopf +1, Stillstand der Konvergenz in einer Γout-Zone |x − L| ≈ 0.15 seit N ≈ 1000, Anstieg ×1.9 bei 2400 (2800/3200 laufen).

**C21 (V, 12.10. 11:53) — K19, §32, HAUPTSATZ Fassung 2 gegengelesen: richtig mit Anmerkungen** (v_pruefstand/BERICHT §1.11): E ⊃ D̄ ∖ U_{r₀}; Aussenwinkel 131.7° (nicht 99°);
(P) neben (G-P) entbehrlich bei absoluter Formulierung; Faktor 2 in ARCH §9.7 (1) ist Schreibfehler (K′: 8.2).

**C22 (V, 12.10. 11:53) — s_E gemessen = Satz K22; f_K-Geometrie von E wie kubisch (numerisch)** (BERICHT §3.19, §3.23): 3.9–4.1 (80°/−90°), 3.3–3.4 (60°/−90°) bei r₁ = 0.2;
mit r₁ = 0.05 bei N ≤ 300 nur 2.9–3.3 (Verbindungsbogen). E für f_K: univalent (min|f_K′| = 0.336), ohne Löcher, Schnitt in der Lücke, Öffnung 132°, Dreipunkt 1.07.

**C23 (V, 12.10. 11:53) — F3 mit Linse 0.03 falsch; E4′ ohne z = L: Q_n ~ |P_n|·n^{−3}** (BERICHT §3.20, §3.22): |p_N − f_K| wächst auf der Linse wie e^{+0.015N} (1.75 bei 2000);
Q_n(0.05)/|P_n| = 0.78 → 0.017 (n = 200 → 1800, kubisch = Kneser), R_A = 1.008 nur von den Spiralen, fern exponentiell; primal am Rand bis 1800 noch exponentiell (−0.012/n).

**C24 (V, 12.10. 13:00) — Ks Sätze der Phase 3 gegengelesen** (v_pruefstand/BERICHT §1.12): K2–K6, K13, K16, K19, K20, K22 (Glattheit), O1, Watson/Integraldarstellung, Spiralfenster,
Reduktionssatz lokal — richtig; HAUPTSATZ_KR0 §7 trägt die V-Vermerke.

**C25 (V, 12.10. 13:00) — (Fam-koh) per Carathéodory–Fejér: Repräsentant wächst, Produkt mit dem adjungierten Steg-Gewicht nicht** (BERICHT §3.24): CF-Rate ≈ 0.90(S − 0.0147)
(S = 0.03: 1.5e-3 bei N = 2000), adjungierte Steg-Rate ≈ −0.94(S − 0.0126); Summe ≈ −0.0024 für alle S. Primale Randrate kippt bei n ≈ 2200 in die Potenz (§3.25).
