# Reflexion (Phase 3)

Eintrag nach jedem Versuch: `## [K|V] dd.mm. HH:MM — Versuch` mit Erwartung, Ergebnis, Warum, was es über das Problem lehrt,
Konsequenz für den Plan. Meta-Reflexion alle paar Versuche: Arbeiten wir am Engpass? Wiederholen wir Sackgassen? Was hat sich
an unserem Verständnis geändert?

## [K] 11.10. 20:37 — Versuch 1: Cauchy-Darstellung statt Eckformel (Satz K1) und Messung von G_N auf γ (k01)
- **Erwartung:** Beim Sortieren der Ideen fiel auf, dass J_N⁻¹T_NW für jede auf Int γ analytische Eingabe W nur die Antworten auf Cauchy-Kerne
  r_t (t ∈ γ) braucht. Erwartet war eine Reduktion „(H₀) ⇐ Schranke auf einer Kurve“, aber mit Zusatzlücken wie bei W1'' (D6, Kneser).
- **Ergebnis:** Die Reduktion ist elementar und vollständig (k_kern/BERICHT §1): G_N ≤ CN^κ auf γ = Γout ∪ Linse ⇒ (H₀^{K_γ}) ⇒ (LP), κ < 2.79,
  für jedes analytische f, auch f_K, ohne D6/(E_sing)/(a′)_r/Q_N. Numerisch (kubisch bis 1200, Kneser bis 600): G_N auf Γout und auf der Linse fern
  von L **konstant in N**, nahe L ~N^{0.45}; auf dem Γin-Kreis exponentiell mit Rate → 0.012.
- **Warum:** G_N(t) = sup_u|Y_u(t)| ist die einzige Grösse, die die Pol-Antwort braucht; die Eckformel D6 war ein Umweg über eine Zerlegung
  in Urbildterme, die man für die Schranke gar nicht braucht. Bulk und Ecke stecken beide in G.
- **Was es lehrt:** (1) Der Engpass ist jetzt EINE skalare Funktion auf EINER Kurve, und sie ist ausser in einer Umgebung von L numerisch völlig
  harmlos (konstant). Die gesamte Schwierigkeit sitzt in γ ∩ U(L) — Os Kreisbedingung (KB). (2) Zwei schnelle Versuche, den „harmlosen“ Teil weich
  zu beweisen, scheitern (BERICHT §2): Bootstrap (Maximumprinzip zeigt in die falsche Richtung) und Schwanzabschätzung (T_NJ_∞ auf Schwänzen ist
  wegen |f| > ρ auf Γout exponentiell gross) — auch die Konstanz von G auf Γout ist eine Auslöschung, kein weicher Effekt.
- **Konsequenz für den Plan:** Zielaussage in Teil C wird (H*_G) statt (Ẽ); V14 und V8 entfallen; die Familie wird „G^P ≤ CN^κ auf γ“.
  Die Modellarbeit (A2) soll gezielt (H*_G) im Modell angreifen, und die Ecke L ist der einzige Ort, an dem neue Analysis nötig ist.

## [K] 11.10. 20:46 — Versuch 2: Ist die Eck-Amplitude von G_N lokal? (k02, kubisch / jet±0.01 / radial, N = 300…700)
- **Erwartung:** Gibt die Jet-Störung (L, μ fest, global anders; λ ändert sich um 0.4 %) in L dieselbe Amplitude G_N(L)/√N, dann wäre die Ecke in
  der G-Form lokal bestimmt — dann reichte eine lokale Parametrix (Gegenprobe zu W §27.5).
- **Ergebnis:** ln(G_N(L)/√N): kubisch −0.878/−0.906/−0.925, jet+0.01 −0.880/−0.909/−0.928, jet−0.01 −0.875/−0.904/−0.922; radial (arg μ 31.4°)
  −0.860/−0.884/−0.899. Fern von L (Γout 0°) ändert jet ±0.02, nahe L ±0.002–0.003, ohne mit N abzunehmen.
- **Warum / was es lehrt:** Die Eck-Amplitude hängt schwach, aber messbar global ab (relativ ≈ 0.3 % je ε = 0.01, vergleichbar mit λ); μ (radial)
  wirkt stärker. Das bestätigt W §27.5 auch in der G-Form: eine rein lokale Parametrix kann die Amplitude nicht liefern, nur ihre Grössenordnung √N.
  Für (H*_G) braucht man aber nur eine Schranke, keine Asymptotik — die globale Abhängigkeit ist klein und glatt.
- **Konsequenz:** A2' (Ecke direkt) muss die Amplitude als Eingang behandeln, z. B. als Unbekannte in einer kleinen Gleichung (Lyapunov–Schmidt
  wie V11), nicht als lokales Resultat. Keine Planänderung.

## [V] 11.10. 20:56 — Versuche 1–3: Schärfe von (K2), Diagnose p32, Gegenlesen P-V und K1
- **Erwartung:** (1) Ps Ziel C₆ ≤ 300 für (K2) ist nötig (so stand es in ARCH §9.9). (2) p32 scheitert im Strukturgebiet an „Überschätzung“.
  (3) P-V und K1 halten, aber mit Zusatzannahmen.
- **Ergebnis:** (1) Ms ganzer 𝒯-Beweis läuft mit C₆ = 5e5, p = 6 durch (Fenster 6/14): σ_min(J_∞|P_N) ≥ 0.5695; mit dem schärferen Fenster 7/23
  scheitert schon C₆ = 5e4 (Fensterende 3e-10 neben L). (2) p32 lieferte inf (Taylor-Modell-Radius 1.5 um |H₀| ≈ 1.9 ⇒ elffache exp-Kette) und nahm
  nan an, weil `nan > x` falsch ist; mit Zellen ≈ 8e-4 ist alles endlich, und die Reserve aus (1) erlaubt Überschätzung ×5. Ausserdem rechnete p32
  θ als float — für d < 1e-15 unbrauchbar (jetzt d als Variable). (3) P-V: Logik richtig, nur formale Lücken (τ in Prüfung B, „Lipschitz 10“,
  Zweigtreue). K1: richtig, ohne Zusatzannahmen ausser Geometrie; Korrektur κ < 1.79 mit Familie.
- **Warum:** Die Konstanten in ARCH §9.9 waren nie gegen die tatsächliche Empfindlichkeit von Ms Skripten getestet; Schranken wurden schärfer
  verlangt als gebraucht. Bei p32 war der Fehler technisch, nicht grundsätzlich.
- **Was es lehrt:** Endliche Bausteine zuerst auf ihre wirkliche Empfindlichkeit prüfen, bevor man scharf zertifiziert — das spart Grössenordnungen.
  Für den Engpass ändert sich nichts: alle endlichen Teile sind erreichbar, offen bleibt allein (H*_G) (und ihre Familienfassung).
- **Konsequenz:** (K2) mit grober Zellteilung auf aicgolling (läuft); Plan A5. Meta: Wir arbeiten am Engpass (K: (H*_G) im Modell), V hält den
  endlichen Unterbau bereit; keine Sackgasse wiederholt (K1 benutzt keine Normschranke von J_N⁻¹, nur die ausgewertete Grösse G_N).

## [V] 11.10. 21:15 — Versuch 4 (Plan A2): loxodromisches Möbius-Modell f = M⁻¹(√L·M) — trägt es Kneser-Ecke und (H*_G)?
- **Erwartung:** Das Modell hat in L dieselben lokalen Daten (μ = √L) und exakte Toeplitz × GL(2)-Struktur; falls G_N auf γ wie bei
  kubisch/Kneser polynomial bleibt, wäre es das ideale Spielzeug für (H*_G).
- **Ergebnis (v08_moebius_G.py, aicgolling, N = 100/200/300):** in L wie Kneser (ln G − ½ln N = −0.82/−0.86/−0.89), Γout konstant (wie Kneser).
  Aber auf der ganzen Γin-Seite wächst G exponentiell, am stärksten nahe L̄: Kreis −80° ln G = 5.0/12.8/20.7 (Rate 0.078 = Ds d24-Rate),
  auch auf der Linse (−80°: 4.1/11.0/18.0; +180°: 2.2/4.0/6.9; +100°: 2.2/2.6/3.5).
- **Warum:** Eine Möbius-Abbildung hat nur EINEN abstossenden Fixpunkt; L̄ ist anziehend (f'(L̄) = 1/μ) und liegt ebenfalls auf dem
  Abschneidekreis. Die starke Instabilität (0.078) sitzt dort, nicht in der abstossenden Ecke L.
- **Was es lehrt:** Der anziehende Randfixpunkt erzeugt eine eigene, viel stärkere Instabilität, die kein Linsenrand umgeht — das Modell testet
  eine andere Frage. Die Kneser-Ecke braucht zwei abstossende Randfixpunkte (wie kubisch/radial). Ds Zuordnung „Eckmode bei L, Rate 0.077“
  (und Xs Vermutung „Eckmode bei M = 0“) ist damit korrigiert: die 0.078 gehört zu L̄.
- **Konsequenz:** Abbruchkriterium von A2 (Möbius) erreicht ⇒ laut Plan Wechsel zum radialen Quadrat (K5); Möbius nicht weiter verfolgen.

## [K] 11.10. 21:45 — Versuch 3: Lemma K1′ entschieden (k04/k04b/k05; V prüft unabhängig bis N = 1500)
- **Erwartung:** Nach K1′ (Äquivalenz bis √N mit max_γ G ~ N^{0.45}) höchstens ‖Ĵ_N⁻¹‖ ≲ N^{0.95}; ob die Abschnitte in der konformen
  Gleichungskoordinate wirklich stabil werden, war offen.
- **Ergebnis:** ‖Ĵ_N⁻¹‖ = 1.251 konstant für N = 100…800 (kubisch) und 1.250 für N = 300…600 (Kneser), während ‖J̃_N⁻¹‖ bis 619 wächst. Grösste
  Verstärkung = Bulk bei 180°. Kontrollen: (i) zu dünne Linse (S = 0.010 < c/0.95) ⇒ Eckanteil wächst exponentiell mit Rate ≈ 0.001 wie vorhergesagt;
  (ii) bei passender Linse wächst der Eckanteil nur polynomial (N^{0.3}, mit nahe L dickerer Linse N^{0.15}) und bleibt weit unter dem Bulk;
  (iii) Widom-Struktur: Inverse = Abschnitt der unendlichen Inverse + Kantenkorrektur O(1) in einem Band der Breite ∝ N. V: K1″ (a)–(c) logisch verifiziert,
  eigener Code reproduziert 1.251155 (N = 100, 200).
- **Warum:** In w-Koordinaten wird die Gleichungsseite der Eckmoden (Steg bei 0.95N, Winkel π) um e^{0.95N·s(π)} gewichtet; liegt die Linse über
  dem Eckprofil (s ≥ Λ/0.95), verschwindet die exponentielle Verstärkung; übrig bleibt die Innenschicht an L, und die hängt an der Linsenform nahe L.
- **Was es lehrt:** (1) Die Eckmoden sind ein Artefakt der Gleichungskoordinate — das ist die bisher klarste Beschreibung des Phänomens. (2) Die globale
  Eckrate verschwindet aber nicht aus dem Problem: sie steckt als UNGLEICHUNG (Linse über dem Profil, c < 0.95 s(π)) in der Voraussetzung (S_γ).
  Ein Beweis von (S_γ) beweist also nebenbei die Wachstumsschranke (I). (3) „Strikte Expansion des dominanten Kanals genügt“ (IDEEN K14) ist widerlegt.
  (4) Die Kante ist nicht klein, sondern ein O(1)-Band — die Parametrix braucht einen expliziten Kantenbaustein.
- **Konsequenz:** Zielaussage jetzt (S_γ) mit einer Linse, die nahe L dicker ist (τ₁ = 5°, p = 1.15). Plan-Schritt 2 (Parametrix) im radialen Modell
  (Möbius ist laut V ungeeignet: der anziehende Randfixpunkt L̄ erzeugt Wachstum 0.078/Grad auf der ganzen Γin-Seite).

## [V] 11.10. 21:46 — Versuch 5: (S_γ) unabhängig nachgerechnet, mit Fallen-Tests und Linsen-Sweep
- **Erwartung:** Ks 1.251 hält entweder (dann ist es der wichtigste Befund der Phase) oder es ist ein Artefakt der Abbildung/Präzision. Ein echter
  Effekt muss beim Ausdünnen der Linse (S → 0) in das bekannte Wachstum von ‖J̃_N⁻¹‖ übergehen.
- **Ergebnis:** eigener Code (direkte Spalten statt C_k·J̃_N, eigene Theodorsen-Iteration, eigenes Kubik) gibt 1.251155 für N = 100…1500 (kubisch) und
  1.250337 für N = 300…639 (Kneser). Unempfindlich gegen gröbere Abbildung, gestörte F_n und mehr Präzision. Der Sweep zeigt ein sauberes Gesetz
  ‖Ĵ⁻¹‖ ≈ max(B, A e^{(c(N) − S)N}) mit κ = 1.00; eine vorab notierte Vorhersage (S = 0.008: Einsatz ≈ 800) traf ein. Ein anderes Linsenprofil mit
  gleichem S wächst.
- **Warum:** In der w-Koordinate wird die Gleichungsseite der Eckmoden (oberste Indizes, Γin-Seite) um e^{−N s} gedämpft; stabil ist es genau dort,
  wo die Linse über dem n-te-Wurzel-Profil Λ(τ) der Eckmoden liegt — das ist Ks K1 in Normform.
- **Was es lehrt:** (1) (S_γ) ist echt und kein Rechenartefakt. (2) Die Linse ist NICHT frei wählbar: sie muss über Λ liegen — Λ ist eine
  n-te-Wurzel-Grösse (Ws freies Randwertproblem, c = 0.0118). (3) In dieser Koordinate sind die Eckraten negativ; die alte Sackgasse („e^{o(N)} reicht
  an der Ecke nicht“) gilt dort nicht mehr — ein Beweis könnte Ecke (n-te Wurzel mit Reserve) und Bulk (gleichmässig) trennen. Das ist eine
  Vermutung über die Beweisstrategie, kein Satz.
- **Konsequenz für den Plan:** Zielaussage = (S_γ) für Linsen mit Reserve über Λ. Für K: (i) Abschnittsmethode für Ĵ_∞ = C_kJ̃_∞ (Grenzoperatoren,
  Roch–Santos–Silbermann) mit Ecke als exponentiell kleinem Rest; (ii) obere Schranke für die Eckrate entlang Γin (n-te Wurzel, W11) — braucht Λ
  nur als Ungleichung. Für V: Familie in w-Form messen (‖(Ĵ^P_N)⁻¹‖ für P = p_N), Kopplung (A), dann Ms Lemmata.

## [K] 11.10. 21:55 — Versuch 4: Reserve-Weg geprüft (Vorschlag V/Leitung) — scheitert, mit präzisem Grund
- **Erwartung (V):** In der konformen Koordinate sind die Eckraten negativ (c − S < 0), also könnte n-te-Wurzel-Genauigkeit plus Reserve genügen.
- **Ergebnis:** Negativ ist nur die Rate des Steg-Anteils (hohe Indizes, fern von L). Rate genau 0 haben (a) die niedrigen Koeffizienten der Eckmoden
  (Y(∞) = y₀) und (b) die Werte in L. (S_γ) erzwingt für beide polynomiale Schranken: |y_{u,0}| ≤ ‖Ĵ_N⁻¹‖ exakt, und allgemein ist (S_γ) ⇒ (H₀) für
  alle glatten Eingaben (J_∞⁻¹g hat die Kneser-Singularität (x−L)^{s₀}). Für diese gilt W §26.3 unverändert.
- **Warum:** Die Linse dämpft nur Anteile, deren Gleichungsseite bei hohem Index sitzt; die Anregung der Eckmoden durch glatte Daten ist ein
  Effekt bei niedrigem Index und in L — dort gibt es keine Reserve.
- **Was es lehrt:** (S_γ) ist nicht schwächer als (H₀), sondern eine gleichmässige Operatorfassung davon. Der Reserve-Gedanke zerlegt aber sauber:
  Steg (n-te Wurzel + Wachstumsschranke (I)), Bulk (σ_{r+1} ≤ B), Anregung durch glatte Daten (polynomial, echt offen).
- **Konsequenz:** Kein Umweg um die starke Ebene. Plan-Schritt 2 (Parametrix) bleibt der Hauptweg; der Steg-Teil darin ist der leichte Teil
  (exponentielle Reserve), die Arbeit steckt im Eckbaustein und im Kantenband nahe L.

## [K] 11.10. 21:55 — Meta-Reflexion nach vier Versuchen
- **Arbeiten wir am Engpass?** Ja: alle vier Versuche zielten auf die Eckanregung. Der Engpass ist jetzt so scharf formuliert wie noch nie:
  (S_γ) sup_N‖Ĵ_N⁻¹‖ < ∞ — gleichmässige Stabilität einer festen Abschnittsfolge, numerisch mit 7 Stellen konstant bis N = 1500 (V).
- **Wiederholen wir Sackgassen?** Bootstrap (Versuch 1) und Reserve (Versuch 4) waren Varianten von „weich um die Ecke herum“ — beide scheitern aus
  demselben Grund wie W §22–§26 (Maximumprinzip in die falsche Richtung bzw. Rate 0 in L). Neue Regel für uns: **jede Idee zuerst daraufhin
  prüfen, ob sie die Eckanregung durch glatte Daten polynomial kontrolliert; wenn sie nur Raten liefert, ist sie für den Engpass wertlos.**
- **Was hat sich am Verständnis geändert?** (1) Die Eckmoden sind ein Artefakt der Gleichungskoordinate; in der konformen Koordinate bleibt eine
  gleichmässig beschränkte Inverse, deren Struktur „Bulk-Inverse + O(1)-Kantenband“ ist (Widom). (2) Die globale Eckrate tritt nur noch als
  Ungleichung (Linse über dem Profil) auf. (3) Das Problem ist ein Stabilitätsproblem einer Zwei-Kanal-Abschnittsfolge — Werkzeuge: Parametrix bzw.
  lokales Prinzip mit Kanten-Grenzoperator; der Kanten-Grenzoperator ist semiklassisch (Band ξ ∈ (1−δ, 1]), das ist das eigentlich Neue.

## [K] 11.10. 22:12 — Versuch 5: Grenzoperatoren und klassische Abschnittstheorie (Auftrag Leitung 22:00; k06, k07, k08)
- **Erwartung:** Stabilität ⇔ Invertierbarkeit von Bulk, Kanten- und Eck-Grenzoperatoren (Böttcher–Silbermann/RSS/Lindner); die Dreiteilung passt.
- **Ergebnis:** Die Grenzobjekte lassen sich numerisch bestimmen — Bulk 𝒥⁻¹ (1.2512), Kante fern von L trivial, Ecke: O(1)-Korrektur (0.49 → 0.56)
  genau bei arg L im Indexband ξ ∈ (0.6, 1), in Paaren L/L̄. Aber die klassische Theorie greift nicht: (1) exponentielle Fehlskalierung — Parametrix-
  Residuen e^{0.47N} (k07), Schwanzabschnitte von 𝒥⁻¹ exponentiell singulär (k08); (2) grundsätzlich: in der x-Koordinate sind alle lokalen Objekte
  invertierbar und trotzdem wächst die Inverse wie e^{cN}.
- **Warum:** Der Kompositionskanal transportiert Pakete makroskopisch durch den Phasenraum; die Eckmode ist eine globale Schleife. Bandbeherrschte
  Operatoren (Domäne der klassischen Theorie) haben keinen solchen Transport.
- **Was es lehrt:** Ein Stabilitätssatz braucht neben der lokalen Invertierbarkeit eine **dynamische Schleifenbedingung** (Spektralradius der gewichteten
  Zirkulation < 1, Antonevich-Typ); in w-Koordinaten ist das die Linsenbedingung s ≥ Λ, deren Schärfe V gemessen hat (S* = c_∞). Damit ist klar, welche
  neue Mathematik fehlt: eine semiklassische Abschnittstheorie für gewichtete Verschiebungs-/Kompositionsoperatoren mit Abschneiden im Grad.
- **Konsequenz:** Plan-Schritt 2 neu fassen: (a) Eckmodell (Wiener–Hopf für I + μ⁻¹D_μ, halbseitig) explizit aufstellen und gegen D_N testen;
  (b) die Schleifenbedingung als eigene, globale Ungleichung formulieren (Kandidat: tropische/Max-Plus-Schranke der Zirkulationsrate); (c) die Kopplung
  „lokal + global ⇒ stabil“ beweisen — das ist der Kern der neuen Theorie.

## [V] 11.10. 22:46 — Versuch 6 (Plan A3): Familie in w-Form — trägt eine Störungsrechnung um f_K die Familie T_N𝒥_P?
- **Erwartung:** Nach v17 (‖(Ĵ^P)⁻¹‖ unverändert bis ε = 1e-3) sollte eine Neumann-Reihe um Ĵ_N genügen: q = ‖Ĵ⁻¹ΔĴ‖ ≲ Cε gleichmässig in N.
- **Ergebnis:** Nein in ℓ²: q_L, q_R ~ ε·e^{0.18N} (Ĵ hat Einträge bis 1.6^N, weil |f| auf Γout bis 1.6ρ reicht). Trotzdem ist die wahre Änderung
  der Inversen linear in ε und gleich dem Glied erster Ordnung; die Reihe Σ(−X)^kĴ⁻¹ konvergiert rasch (Quotient ≈ 10ε). Die Lipschitz-Konstante
  Lip_N von P ↦ (Ĵ^P)⁻¹ wächst bei glatten Störungen linear in N (Ableitungsverlust), bei Störungen hoher Grade aber wie e^{≈0.011N} —
  unabhängig von der Linsendicke (0.0166 … 0.1). Realistisch (Profil von p_N − f_K aus den Zertifikatsboxen) ist ε·Lip_N ≈ 1e-14 bei N = 400,
  hochgerechnet ≈ 1 erst bei N ≈ 5000.
- **Warum:** Die Dämpfung der grossen Spalten leistet der rechte Faktor Ĵ⁻¹, nicht eine kleine Norm von ΔĴ. Die Eckrate c ≈ 0.012 ist in Ĵ⁻¹
  selbst durch die Linse beseitigt, taucht aber in der Ableitung nach P wieder auf (hohe Grade von η regen die Eck-Indizes an).
- **Was es lehrt:** (1) Eine Familien-Erweiterung braucht eine gewichtete Norm (δ auf f(Int γ)), nicht ℓ². (2) Reine Störungsrechnung um f_K
  reicht für alle N nicht, solange Lip_N exponentiell wächst; sie reicht weit über alle gerechneten N hinaus. (3) Ein Beweis muss entweder die
  Struktur von p_N − f_K nutzen (Gewicht fällt in den obersten Graden, auf L gerichtet: dort Lip kleiner) oder (S_γ) direkt für Ĵ^{p_N} zeigen
  (Fortsetzung in N statt Störung um f_K). Gradbänder (22:49): exponentiell nur für ξ = k/N ≥ 0.8, Mittelbänder beschränkt. Keine Sättigung: hoch bei N = 1200 Lip = 1107 (Rate 0.0096).

## [V] 11.10. 23:06 — Versuch 7 (Auftrag Leitung 22:52): gewichtete Norm für die Familie, und (S_γ) am echten p_N
- **Erwartung:** Ein Gewicht, das die obersten Grade wie e^{0.010k} misst, macht die Lipschitz-Konstante polynomial, und p_N − f_K bleibt darin klein.
- **Ergebnis:** Nein. Für Kneser mit der echten Tangente (Multiplikator P'∘P) ist die Antwort auf einen Einzelgrad r_k ≈ e^{N·g(k/N)}, g(1) = 0.161
  (nicht 0.010 wie im kubischen Sekantenmodell), unabhängig von der Linse; mit ω_k ≳ r_k ist das echte p_N − f_K schon bei N = 200 nicht mehr klein.
  Auch die Graphennorm scheitert (der Schwanz von f_K fehlt in T_N J_∞). Dagegen: (S_γ) gilt für Ĵ^{p_N} mit derselben Konstante wie für f_K bis
  N = 1000, Abstand nur 0.015·N·ε_N, und alle kohärenten Richtungen J⁻¹(glatt) haben Lipschitz-Konstante linear in N.
- **Warum:** Die grosse Antwort kommt aus η'∘f (Ableitung der Störung, ausgewertet auf f(Γout) bis 1.6ρ). Beim echten p_N hebt der Schwanz von f_K,
  mit f komponiert, genau diesen Teil auf; die Kleinheit ist eine Struktur-Eigenschaft (Fast-Lösung), keine Grössen-Eigenschaft.
- **Was es lehrt:** Die Familie gehört über die Lösungsmenge geführt (H1-Homotopie, ARCHITEKTUR §5.1(ii)), nicht über eine Kugel in einer Norm.
  Für Ks Kopplungssatz heisst das: (S_γ) an Fast-Lösungen + Lipschitz längs kohärenter Richtungen. Eine Norm, die nur Grössen misst, kann das nicht.

## [K] 11.10. 23:40 — Versuch 6 (Plan A2a/A2b, lange Strecke): Ist das Eckmodell eigenständig? Liefert Ws Tropik die Schleifenbedingung?
- **Erwartung:** (a) Die Eck-Korrektur D_N hängt nur am Jet in L; dann gibt es ein explizites Eckmodell (I + μ⁻¹D_μ mit Gradabschnitt), dessen
  Invertierbarkeit man beweist. (b) In w-Koordinaten bricht die Linse Ws Verschiebungsinvarianz (W13), und die Tropik liefert die Schleifenbedingung.
- **Ergebnis:** (a) Halb richtig. D_N ist universell (kubisch = Kneser auf 4e-4, flache Störung 4. Ordnung ≤ 2e-5, radial mit anderem arg μ deutlich
  anders; k09). Aber das affine Modell f = L + μ(x − L) ist exponentiell instabil (‖Ĵ_N⁻¹‖ ~ e^{0.70N}, Rate exakt = ln max_j max_{D̄}|φ^{−j}|/ρ),
  weil J_aff auf E²(D ∪ φD) einen Kern (x−L)^{s_k} hat (k10). (b) Falsch: die Linse ändert die Identitäten nicht, nur die Auswertung (K1″); Satz K3
  (Homogenität) zeigt, dass Schranken nur aus der Stelle kommen, wo die Quelle allein ausgleicht — das ist weiter nur L.
- **Warum:** Der Gradabschnitt ist eine globale Bedingung; die Dilatation um L mit komplexem μ dreht D̄ teilweise nach aussen. Ohne globale Struktur
  setzen sich die lokalen Mellin-Lösungen fort. Lokal ist nur die Operator-Ecke (O7a), und die steckt schon in der Bulk-Invertierbarkeit.
- **Was es lehrt:** Die Ecke ist keine eigene Hypothese (Vermutung S′: Bulk + Schleife ⇒ stabil, Eck-Korrektur universell). Die Schleifenbedingung
  zerfällt in Gestalt (tropisch, W12′) und Niveau (Eck-Normierung u(L) = 0). Nebenbei (k11, Familie): Satz K4 P'∘P = (eˣ + R_P')/P' erklärt,
  warum die Newton-Matrix an Fast-Lösungen (p_N, Π_N f_K) exakt so stabil ist wie J_N, an T_N f_K aber exponentiell instabil (2.3e25 bei N = 600).
- **Konsequenz:** A2a und A2b abgeschlossen (BERICHT §9). Weiter mit A2c in der Form S′: zuerst „(B) + (SB) ⇒ n-te-Wurzel-Stabilität“ (exponentielle
  Ebene, Parametrix mit Fehler e^{−ηN}), dann die Eck-Korrektur als universelle O(1)-Grösse. Familie über 𝒮_N mit Satz K4; offene Stelle: Linsen-
  Abweichung der Polynome für N ≫ 2000.

## [K] 11.10. 23:59 — Versuch 7 (Auftrag Leitung): Kippt die Linsen-Abweichung der Familienpolynome den Weg für grosse N?
- **Erwartung:** Echte Gefahr — P weicht auf der Linse um e^{Ns}N^{−6} von f_K ab, ab N ≈ 2000–3000 von der Grösse 1.
- **Ergebnis:** Für (S_γ) nein. Ĵ_N^P hängt nur vom Keim von P ab (bewiesen, mit K4), K1″(b, c) gelten wörtlich; S = 0.1, N = 600: p_N weicht auf der
  Linse um 8.1e9 ab, ‖(Ĵ^{p_N})⁻¹‖ = 0.965153 wie bei f_K (k11d). P auf γ geht nur in Daten-Schlüsse ein ((H-ii), H1-A-priori), dort kostet der
  beste Repräsentant (Carathéodory–Fejér) e^{NF₀}·Kopf, F₀ = harmonisches Mittel der Linsenhöhe ⇒ Schwelle N ≈ 6300–7800 im inkohärenten Fall.
- **Warum:** Die Linse ist ein Werkzeug auf der Gleichungsseite (Ort der Daten), kein Ort, an dem P ausgewertet wird. Die wirklich linsenunabhängige
  Grösse ist der Eckgewinn e^{ΛN} für inkohärente hohe Grade (= Vs Lip e^{0.011N}); Fast-Lösungen sind kohärent (V: polynomial bis N = 1000).
- **Was es lehrt:** Familien-Aussagen in Keim-/G-Form formulieren, nie mit 𝒥_P als Operator auf Int γ. Der Rest ist eine Kohärenz-Aussage (Fam-koh)
  über die Datenkorrektur Φ''(f)[h, Δ] — dieselbe Ecke wie (H₀), keine neue Front.
- **Konsequenz:** A2c Schritt 1 kann beginnen; K-A2c-Familienfassung in Keimform (BERICHT §10.3). V prüft bei N = 1500/2000 Norm + Linsen-Abweichung + Datenkorrektur.

## [V] 12.10. 00:14 — Versuch 8 (Auftrag Leitung 23:20): Konstanten von K1 explizit, Familie längs H1/H2
- **Erwartung:** Mit (S_γ) als einziger offener Zahl lässt sich die ganze Kette K1″ → P-I → X4a → H1/H2 ausrechnen; die Familie längs der
  Homotopie hat dieselbe Konstante wie f_K.
- **Ergebnis:** Geometrie der Linse streng (v32), Green-Konstanten numerisch (v30), Kette (v31): θ_N ≤ ½ für B ≤ 30.8 ab N₀ = 318 und B ≤ 76.7 ab 400,
  gemessen B = 1.25. Familie: Sekanten, Tangenten an P(λ) und σ-Homotopie geben bis N = 1000 überall 1.250337 (S = 0.05: 1.107257). Zwei technische
  Fallen: (1) Kompositionen f_K∘P über Σ a_kP^k sind ab N = 400 unbrauchbar (Datenfehler der a_k um bis 1e144 verstärkt) — Taylor um f(x) nehmen;
  (2) die σ-Homotopie braucht p_N auf ≈ 0.4·N Stellen, sonst ist der Kopfrest c_(N+1) Rundungsrauschen.
- **Warum es trägt:** Die Matrizen hängen nur vom Keim ab (Ks §10.1); die Fast-Lösungen sind kohärent, ihre Abweichung ändert die Inverse nur polynomial.
- **Was es lehrt:** Die einzige echte offene Zahl ist B; zwischen gemessenem B und benötigtem B liegt ein Faktor 25–60. Offen bleibt (Fam):
  der Datenteil der Familie (W_P bzw. Π^P k_η), den K über Keimform und Kohärenz führen will; gemessen ist er klein (v35).

## [K] 12.10. 00:25 — Versuch 8 (Plan A2c Schritt 1): n-te-Wurzel-Stabilität aus Bulk + Schleife, ohne Ecke?
- **Erwartung:** Eine Parametrix mit Fehler e^{−ηN} liefert die exponentielle Ebene; die Ecke kommt erst auf der polynomialen Ebene.
- **Ergebnis:** Nein, endgültig. L liegt auf γ; n-te-Wurzel-(S_γ) enthält G_N(L) = e^{o(N)} — linsenfrei, eine reine Eckaussage. Fern von L reduziert
  Lemma K6 auf das Grad-Antwortprofil ψ; dessen Teil für glatte Daten ist Ws (H-I), wieder Ecke. Neu bewiesen: Satz K5 (Kopf-/LU-Zerlegung,
  U unabhängig von k, D_N = projizierte künftige Köpfe, n-te Wurzel ⇔ Köpfe subexponentiell). Numerik k12 (N ≤ 1000): ψ = 0 für ξ ≤ 0.8,
  exponentiell nur im Kantenband; w-Köpfe sättigen bei 0.015, x-Köpfe wachsen wie e^{0.011n}; G_N(L) ~ N^{0.44}.
- **Warum:** Die Linse hat in L die Höhe null; das Niveau sitzt in L (§9.2) — die Trennung „exponentiell zuerst, Ecke danach“ war ein Planungsfehler.
- **Was es lehrt:** Vor jedem Teilziel prüfen, ob es an einem Punkt mit Linsenhöhe 0 ausgewertet wird. Die Ecke ist auf beiden Ebenen der Engpass.
- **Konsequenz:** A2c neu: (E₀) exponentiell an der Ecke (G_N auf γ ∩ U(L), glatte Daten ψ ≤ o(1)) + Lemma K6 fern von L; (E₁) polynomial als
  Fast-Orthogonalität der Köpfe (Mellin). Beides im universellen Eck-Grenzobjekt angreifen.

## [K] 12.10. 00:34 — Versuch 9 (Auftrag Leitung): Eck-Kopfform — Anschluss an G, kurzes Fenster, zweite Ideenrunde
- **Erwartung:** Mit Gs skalenfreiem Eckprofil lässt sich eine Teilaussage (Profil als Spezialfunktion oder als Lösung einer Gleichung) beweisen.
- **Ergebnis:** Nein. K5 ist in der Substanz Gs G15–G18 bzw. Ws §25.1(c) (Übersetzung BERICHT §12.1); neu sind nur die w-Fassung (gleiches U),
  D_N als projizierte künftige Köpfe und die n-te-Wurzel-Äquivalenz. S_1000 ist kein ₁F₁, eine Gleichungssuche hat keine Trennschärfe (S_1000 liegt
  ~3e-3 neben S_∞), ein exakt lösbares Eckmodell mit adjungierter Struktur gibt es nicht (affin: 𝒜_n = w^n).
- **Warum:** G hat die direkte Route (Profil, Amplitude, G37) bereits ausgereizt; was fehlt, ist eine Idee, keine weitere Rechnung derselben Art.
- **Was es lehrt:** Vor eigenen Konstruktionen Gs Ordner auf dieselben Objekte prüfen (Köpfe = Q_nA_nᵀ/h_n). Das Fenster war kurz, weil die naheliegenden
  Tests sofort negativ waren.
- **Konsequenz:** Zweite Ideenrunde nur zur Ecke eröffnet (IDEEN.md „Runde 2“, R1–R12); vorläufig R12 (Quadranten-Wiener–Hopf aus D_N), R2
  (Renormierung, Gs Korrekturordnungen passen), R5 (IIKS aus Gs Rang 2). V ergänzt und bewertet.

## [K] 12.10. 00:47 — Versuch 10: erste Tests R2, R5, R12
- **Erwartung:** R12 am stärksten (eigene Formel), R2 und R5 als Ergänzung.
- **Ergebnis:** Umgekehrt. R2 besteht klar (Eigenwerte 2^{−a}, 2^{−1}, 2^{−2a} aus einer freien linearen Abbildung), R5 scheitert an der Struktur
  (kein Sylvester-Paar), R12 hat einen oszillierenden Kern (Auslöschung ∝ N).
- **Warum:** Die Köpfe tragen eine schnelle Phase; Normen allein sind skalenfrei, die Summe nicht. R2 arbeitet mit Gs phasennormierten Profilen.
- **Lehre:** Bei Kopfsummen immer zuerst die Phase abspalten. **Konsequenz:** R2 als Hauptstrang zur Ecke, sobald V bewertet hat; R1 mit dem
  Nebenbefund φ(L) = 0 als möglicher zweiter.

## [K] 12.10. 00:55 — Versuch 11 (R2 Hauptstrang): Renormierungsoperator aus G22?
- **Erwartung:** G22 auf der Kegelskala ergibt einen Mellin-diagonalen Operator; seine charakteristische Gleichung liefert a, 1, 2a; Kontraktion in W_{1.05}.
- **Ergebnis:** Nein. G22 ist Vorwärts-Substitution mit dem LU-Faktor L; der Kern wächst wie e^{0.14n} (x und w), Masse bei j ≈ 0.85n (global, nicht 1/|μ|),
  die Köpfe entstehen durch Auslöschung. Neu: Fast-Orthogonalität als Cotlar–Stein-Bedingung, gemessen A₁ = 12.6 (langsam wachsend), A₂ = 4.1 bei N = 300.
- **Warum:** Vorwärts-Rekursionen mit exponentiell grossen Dreiecksfaktoren sind instabil; die Struktur sitzt in der Lösung.
- **Lehre:** Für Kontraktionsbeweise nur stabile Darstellungen verwenden (Zukunftsform, lokale Momentbedingungen). **Konsequenz:** R2 bleibt als
  Struktur; R1 aufgenommen (Quadranten-Ecke mit homogenem Symbol, Winkelgleichung); die Cotlar–Stein-Form ist das neue Ziel der polynomialen Ecke.

## [K] 12.10. 01:07 — Versuch 12 (Hauptstrang Cotlar–Stein): woher die Fast-Orthogonalität kommt
- **Erwartung:** Die schnelle Phase erzeugt die Abnahme; Cotlar–Stein mit Phasen- und Lagekontrolle gibt A₁, A₂ ≲ log N.
- **Ergebnis:** Teilweise. Adjungierte Köpfe: feste, nicht-stationäre WKB-Phase (n-unabhängig auf 4e-4 rad, Φ' ∈ [−0.96, −0.56]) ⇒ echte Abnahme
  (Lemma K8 bewiesen). Primale Köpfe: kohärent, nur Skalenabnahme. Ecke: Vorzeichen-Auslöschung (‖Σ T^E‖ ≈ 0.6 bei Σ‖T^E‖ ≈ 8), die Absolutbeträge
  nicht sehen. Satz K7 + K8 geben κ = 1/2 nur bei kleinem Eckanteil; gemessen wächst er (3–6 %) ⇒ asymptotisch κ = 1.
- **Warum:** Cotlar–Stein summiert Beträge; die Ecke lebt von Phasen.
- **Was es lehrt:** Für den Hauptsatz reicht κ = 1 = nur (H1); die Fast-Orthogonalität ist Verfeinerung. Der Kern ist (H1), die exakte Auslöschung der
  beiden Kopfraten. **Konsequenz:** Cotlar–Stein (absolut) abgeschlossen; nächster sinnvoller Schritt wäre (H1) selbst — mit V in der Ideenrunde.

## [K] 12.10. 01:22 — Meta-Reflexion (Auftrag Leitung/Autor): Kernproblem-Synthese
- **Was ich sah:** Alle Fassungen der offenen Aussage sind Aussagen über die Köpfe (K5); für den Hauptsatz genügt (H1), und (H1) ist asymptotisch eine reine
  Eckaussage. Alle gescheiterten Ansätze behandeln f als beliebige Abbildung mit drehendem Randfixpunkt; keiner benutzt die Herkunft aus f∘f = exp über
  die Abelkoordinate, in der J_∞ exakt die Verschiebung I + S_{½} ist (KERNPROBLEM.md).
- **Lehre:** Die bisherigen Werkzeuge (Normen, Tropik, Rekursionen, Fast-Orthogonalität) sind an derselben Stelle gescheitert, weil sie die Ecke ohne ihre
  exakte Symmetrie behandeln. **Konsequenz:** H-A — Ecke im Torus, Kopf-Gleichgewicht als Erhaltung des Produkts u·v (Liouville).

## [K] 12.10. 01:27 — Versuch 13 (H-A, erster Test): ist die Kopf-Paarung an der Ecke eine invariante Dichte?
- **Erwartung:** Beide Köpfe sind nahe L genähert homogene Eckmoden; ω = mΛδ dy/y ist f-invariant, auf dem Torus konstant (Liouville) ⇒ (H1).
- **Ergebnis:** Nein. Primal homogen nur innen (Defekt bis 1.7e-8), adjungiert nur aussen; im Eckbereich keiner. Q weicht um Faktoren 10^{±3} von 1 ab.
- **Warum:** Jeder Kopf ist die Antwort auf eine oberste Quelle, die auf der einen Seite des Kreises durch L exponentiell klein, auf der anderen
  exponentiell gross ist; die Ecke liegt genau auf diesem Kreis.
- **Lehre:** Die Ecke ist eine Grenzfläche zwischen zwei Sektoren, nicht eine volle Umgebung; Liouville auf dem Torus greift nur für Objekte, die rund um L
  homogen sind. **Konsequenz:** H-A in dieser Form abgebrochen; übrig bleibt die Transmissionsform von (H1) (innen/aussen gekoppelt über den Kreis durch L).

## [K] 12.10. 01:47 — Versuch 14 (Weg (a), a1 + Anfang a2): Vergleichsproblem an der Ecke
- **Erwartung:** Os (P1)–(P3) ist mit der Kopfquelle direkt das Vergleichsproblem; Mellin-Lösung mit Gammafunktionen.
- **Ergebnis:** Struktur klar, Lösung nicht: drei Stokes-Sektoren je Kopf, skalenfrei; Aussenexponenten exakt 2 − s₀ und −1 − s₀ (numerisch auf ~3 %);
  Quelle und Rest löschen sich auf Kopfgrösse aus (Wiener–Hopf nach Graden); zwei Skalen nur durch n^{0.106} getrennt; Profil von Ordnung κ₁.
- **Warum:** Die Ecke koppelt eine Ortsskala (Dilatation μ, Kegel n^{−a}) mit einer Frequenzkante (Grad n); beide treffen sich nur bei K = n^{1−a} → ∞.
- **Lehre:** Das Vergleichsproblem ist der Grenzfall K → ∞ — eine Riemann–Hilbert-Aufgabe auf drei Strahlen, nicht eine Einzel-Spezialfunktion.
- **Konsequenz:** a2 weiter als RH-Aufgabe auf drei Strahlen (Tage); bewiesene Teilaussage in Reichweite: Lösbarkeit/Eindeutigkeit dieser Grenzaufgabe.

## [K] 12.10. 02:14 — Versuch 15 (zweiter K-Lauf, Weg (a) Schritt a2): Mellin–Barnes-Kandidat für das Eckprofil des adjungierten Kopfes
- **Erwartung:** Die Grenzaufgabe K → ∞ ist ein RH-Problem auf drei Strahlen; ihre Lösung liefert Gs S_∞ — wahrscheinlich nur nach langer Rechnung, und
  erst nachdem die Sprungbedingungen hergeleitet sind.
- **Ergebnis:** Abkürzung über den Mellin-Raum: die drei qualitativen Forderungen (y ganz; H ohne Potenzasymptotik im homogenen Keil; y fällt dort) legen
  Φ_y = Γ(s)Γ(½ + iα(s−2))·E(s) fest (Wiener–Hopf-Zerlegung des Symbols μ(1+μ^{s−2}) in Gammafunktionen, E ganz). Die kleinste Wahl E = λ^{−s} trifft Gs
  extrapoliertes S_∞ mit EINEM komplexen Parameter auf 2.5e-7 (Wiener-Norm, 40 Koeffizienten), die Konvergenz S_n → Y(λ·) ist reines n^{−0.894}, Kneser = kubisch,
  und die ganze Defekt-Tabelle von k23 (A = 1, 4, 16; alle Richtungen) wird Ziffer für Ziffer reproduziert.
- **Warum:** Auf der Kegelskala sind die Quellen ±∞; sichtbar bleibt nur die lokale Dilatationsgleichung plus Analytizität (ganz) plus „H rezessiv im homogenen
  Keil“. Im Mellin-Raum ist das ein skalares Wiener–Hopf-Problem mit explizit faktorisierbarem Symbol — der Malyuzhinets-Typ, den O vorgeschlagen hatte, aber
  für den KOPF (homogen), nicht für die Eckmode (mit Quelle).
- **Was es lehrt:** (1) Der Kopf ist das richtige Objekt: seine Ecke ist quellfrei und daher explizit; Os Eckmoden tragen eine eigene Quelle und haben ein
  anderes Profil. (2) „Ordnung κ₁“ (k26) und „Korrekturordnungen n^{−a}, n^{−1}, n^{−2a}“ (G, R2) sind Folgen derselben Gammastruktur. (3) Vor langen
  RH-Herleitungen lohnt es, die Mellin-Analytizitätsbedingungen direkt zu formulieren — sie sind kürzer als die Sprungmatrizen.
- **Konsequenz:** §19: Satzform (Existenz/Eindeutigkeit) der Grenzaufgabe in Mellin-Form mit präzisen Wachstumsbedingungen; primaler Kopf analog
  (Faktor Γ(½ + iα(s+1))?); dann die Paarung an der Ecke (Pivot) aus beiden expliziten Profilen — das ist der Teil, der für (H1) zählt.

## [K] 12.10. 03:00 — Versuch 16 (primaler Kopf, Dualität) und Versuch 17 (Satzform der Grenzaufgabe) [Zeit korrigiert: geschrieben vor 03:03]
- **Erwartung:** Der primale Kopf hat (wie Os Eckmode) eine 1/N-Innenschicht und ist auf der Kegelskala nicht glatt; Satzform nur als lange Skizze.
- **Ergebnis:** Der primale Kopf ist auf der Kegelskala glatt, sein Grenzprofil ist F_{β_P}(λ_Pζ) mit β_P = ½ + iα (freier Fit auf 6 Stellen), und
  λ_P = −μ^{−1/2}λ_A (auf 5.7e-7). Die Gammafaktoren beider Köpfe sind Reflexionspartner (s ↦ 1 − s); mit der λ-Relation ist das Parseval-Produkt
  elementar. Satzform: Lemma K10.1 (Gamma-Zerlegung), Satz K10 (Existenz, explizit), Satz K11/K11′ (Eindeutigkeit unter einer Mellin-Randbedingung,
  Liouville) — bewiesen; Übersetzung der Randbedingung in Sektorbedingungen (Lemma K12) nur Skizze.
- **Warum:** Die 1/N-Innenschicht gehört zu Os Eckmode (mit Quelle), nicht zum Kopf. Ohne Randbedingung ist die Grenzaufgabe nicht eindeutig
  (Dilatationsinvarianz ⇒ Überlagerungen von Dilaten) — die Eindeutigkeit kommt aus dem Anschluss nach aussen, wie bei jeder lokalen Parametrix.
- **Was es lehrt:** Die Ecke ist ein skalares Wiener–Hopf-Problem, die zwei Faktoren sitzen in den zwei Köpfen; das Kopf-Gleichgewicht ist lokal exakt.
- **Konsequenz:** Kopplung: welcher Teil von τ_n ist lokal? (k33/k35: primale Norm lokal, adjungierte nicht) → §20.

## [K] 12.10. 03:00 — Meta-Reflexion (nach gut 1 h im neuen Lauf) [Zeit korrigiert]
- **Am Engpass?** Ja: die Grenzaufgabe K → ∞ war der benannte nächste Schritt; sie ist (numerisch belastbar) explizit gelöst, ihre Eindeutigkeit unter
  Anschlussdaten bewiesen. Für (H1) zählt aber die Kopplung, und dort hat die explizite Lösung sofort etwas Unerwartetes gezeigt (§20, Spirale).
- **Sackgassen-Abgleich (KERNPROBLEM §2):** kein Rückfall in Normen/Tropik/Rekursionen; die Mellin-Analytizität ist neu (sie benutzt die exakte Gamma-
  Zerlegung, die alle früheren Ansätze nicht hatten). Gefahr: zu viel Gewicht auf Modell-Extrapolationen in den Anpassungsbereich legen — deshalb
  jeder Schritt dort mit echten Köpfen prüfen (k35, k38) und V um die entscheidenden Messungen bei N = 1500/2000 bitten.

## [K] 12.10. 08:30 — Versuch 18 (Kopplung, 02:35–03:05; nachgetragen nach technischem Abbruch): Wo sitzen die Kopfnormen, und trägt eine feste Linse?
- **Erwartung:** Mit beiden expliziten Eckprofilen lässt sich τ_n lokal ausrechnen (Normen an der Ecke, Pivot aus einer lokalen Paarung) — (H1) würde
  zu einer Aussage über das Grenzprofil.
- **Ergebnis:** Nur halb. Die primale Norm ist kegel-lokal und explizit (Q_P trifft ‖c_n‖ auf 1 %). Die adjungierte w-Norm sitzt im Anpassungsbereich
  längs der Linse, und dort wächst der Kopf: das Eckmuster dreht sich wegen komplexem κ₁ gegenüber der festen Geometrie um 1.45°·ln n, die Γin-Tangente
  wandert in den Wachstumssektor von F_A. Mit einer dünnen Linse (p = 2) ist das schon bei n ≤ 1000 sichtbar: τ_n = 0.0494 / 0.0339 / 0.0488 (Umkehr).
- **Warum:** Die äussere Rate des adjungierten Kopfes bei L ist eine log-spiralförmige Kegelfunktion Re(C(x−L)^{κ₁}); jede feste Richtung durch L liegt
  bei genügend kleinem r auf ihrer Wachstumsseite. Eine tangentiale Linse (Höhe ~ r^p, p > Re κ₁) dämpft nahe L schwächer als die Rate wächst.
- **Was es lehrt:** (S_γ), (H1) und (H*_G) mit FESTER Kurve sind asymptotisch wohl falsch; die bisherigen Messungen bis N = 2000 sehen das nicht, weil
  der Bulk ‖Ĵ⁻¹‖ = 1.25 bestimmt. Frühere Warnzeichen waren da (Eckanteil ~ N^{0.3} bei der Standardlinse, §5.3; Vs ‖D_N‖ tangential wachsend, VE4 1.33 → 1.49).
- **Konsequenz:** Analytisch klären, welche Linsenform die Drehung mitmacht und was K1″ dann noch hergibt (Auftrag Leitung 08:25) → §21.

## [K] 12.10. 08:37 — Versuch 19 (Auftrag Leitung 08:25): Welche Linse macht die Drehung mit, und trägt K1″?
- **Erwartung:** Eine log-spiralförmige oder mitgeführte Linse folgt der Drehung; K1″ bleibt mit Ecke gültig; damit wäre (S_γ) für alle N rettbar.
- **Ergebnis:** Die Drehung macht genau die Log-Spirale Im(κ₁log(x − L)) = const mit (1.623°/e-Faktor, φ_null(0.01) ≈ 6.7°); geometrisch sauber ist der
  mitgeführte Keil φ_N ≈ 1.45°·ln N + 4°. K1″(a)/(b) tragen mit Ecke, K1′/K1″(c) nicht (unnötig). ABER: (i) der (LP)-Eckexponent kostet 0.030 je Grad,
  mit Familie nur φ < 59° ⇒ N ≲ 10¹⁶–10¹⁷; (ii) unabhängig von jeder Linse wird die primale Norm (am Kreis gemessen) ab n ≈ 10⁸ nicht-lokal (k41) — die
  Dualität lässt nur zwei schmale Richtungsfenster, in denen beide Köpfe nicht wachsen, und die wandern durch die feste Geometrie. Die echten Köpfe
  konvergieren im Anpassungsbereich wie n^{−a} gegen das Modell (k40), die Extrapolation ist also gut gestützt.
- **Warum:** Komplexes κ₁ macht das Eckproblem log-periodisch in N (Periode e^{2π/Im a} ≈ 10¹⁰⁸); jede N-unabhängige Norm/Kurve sieht im Lauf einer Periode
  alle Phasen des Eckmusters, auch die ungünstigen.
- **Was es lehrt:** Der K1-Rahmen (gleichmässige Abschnittsstabilität in fester Norm) ist für N → ∞ zu starr; er reicht weit (10⁸ und mehr), aber nicht
  für alle N. Eine Beweisstrategie für alle N muss die Ecke ausrechnen statt abschätzen — die Mittel dafür (explizite Profile, Dualität, Satz K10/K11) liegen
  jetzt vor.
- **Konsequenz:** Meilenstein an die Leitung; nächster Baustein „K1 mit Eck-Einsatz“ (Cauchy-Darstellung ausserhalb U(L), Kopfentwicklung mit F_A/F_P in U(L)),
  und eine unabhängige Prüfung der Drehung mit einem Modell schnellerer Spirale (kubisch mit |μ| ≈ 2, Im a ≈ 0.11) — dort liegt die Umkehr bei n ≤ 1000.

## [K] 12.10. 08:44 — Versuch 20 (Nachfrage Leitung 08:40): Keil und Bessel; Meta-Reflexion
- **Erwartung:** Mit Keil trägt alles; die Grenzprofile liefern Fast-Orthogonalität beider Kopffamilien ⇒ Bessel ⇒ (S_γ), (H1).
- **Ergebnis:** K1″(a)/(b) und alle K-Sätze gelten mit Keil; fester Keil trägt bis N_× (c = 0.25: 1.2·10⁵, c = 0.5: 4·10⁸, c = 0.79: 1.2·10¹²), die primale
  Kreisnorm linsenunabhängig bis ≈ 10⁸. Bessel: adjungiert fast orthogonal, primal NICHT (zwei Dimensionen je Skalenband) ⇒ gleichmässig nur bei τ_k ≲ 1/k;
  das Grenzprofil gibt mit Keil τ_k ≈ B_k·N_Y·k^{−1/2}|d_k|^{1/(2(1+φ/π))} ~ k^{−0.89} (trifft Vs 4e-4 bei 2000) ⇒ polynomiale Bessel-Schranke N^{0.06}, für K1 reichlich.
  Unabhängig bestätigt: Profilformel mit α = log μ/2π bei |μ| = 2 (primal, ohne freien Parameter) und Stationarität in komplexer Skala.
- **Meta:** Am Engpass gearbeitet (Kopplung → (H1)); die Grenzaufgabe hat sich als Werkzeug bewährt (zwei Vorhersagen von V bestätigt: Umkehr bei Ks Linse,
  Keil fällt). Neue Gefahr: „für alle N“ ist mit festen Normen nicht zu haben — das ist kein Sackgassen-Rückfall, sondern eine Strukturgrenze (komplexes κ₁ ⇒
  log-periodisches Eckproblem). Der Hauptsatz braucht für sehr grosse N eine Eckbehandlung, die beide Köpfe explizit mitführt.

## [K] 12.10. 08:56 — Versuch 21 (L3, Eck-Einsatz): Reduktion auf Restpunkte — es sind vier, und der vierte bricht
- **Erwartung:** Mit kopfweisen Konturen (Lemma K13) verschwindet das N_×-Problem, und (H₀) reduziert sich auf B_n, Anpassung, Stegdämpfung.
- **Ergebnis:** Datenseite ja (jeder Kopf seine eigene Ecke, W darf auf weitere Blätter, polynomialer Preis). Aber sup_{D̄} verlangt auch (iv) primale
  Lokalisierung auf dem Kreis — und die verletzt das Grenzprofil ab n ≈ 10⁸ (k41); bei |μ| = 2 folgt der echte primale Kopf dem Modell samt Drift (k43).
  B_n: K11 legt die Amplitude nicht fest; Anpassung an n-unabhängige äussere Amplituden sagt B_n ∝ n^{−0.053} voraus (mit Gs Daten verträglich, nicht entscheidbar).
- **Warum:** Die Unbekannten werden auf D̄ ausgewertet; dort ist keine Kontur frei wählbar. Die Dualität lässt nur schmale Fenster, in denen beide Köpfe ruhen.
- **Was es lehrt:** Die Spirale trifft nicht nur die Beweistechnik, sondern womöglich die Aussage selbst (gleichmässige Konvergenz bis an L heran für N ≳ 10¹⁰).
  Das muss vor jeder weiteren Beweisarbeit geklärt werden — sonst beweisen wir eine falsche Aussage.
- **Konsequenz:** Meilenstein mit Warnung an die Leitung; entscheidender Test (V): schnelle Spirale + gezähmte Schleife, Konvergenz der Abschnittslösung zu
  singulären Daten am Kreis nahe L. Parallel (K): Lemma K12/Barnes-Asymptotik ausschreiben, sobald die Richtung feststeht.

## [K] 12.10. 09:07 — Versuch 22 (Auftrag Leitung 08:58): Was bleibt wahr, wenn „gleichmässig auf D̄“ kippt?
- **Erwartung:** (a) lokal gleichmässig + beschränkt auf D̄ ist robust und lässt die Eckaussage fallen.
- **Ergebnis:** Glied (c) der Warnung ist genau „B_n polynomial beschränkt“ (sup_{D̄}|δ_n| ≤ √n‖c_n‖). Die D̄-Schranke ist nicht robust; die Reduktion braucht sie
  aber gar nicht (lokale Beschränktheit genügt, Beweis geprüft). Robust ist gleichmässige Konvergenz auf D̄ ohne FESTE Umgebung von L, L̄ (Radius ≳ 10⁻⁶) —
  das liefert Koeffizienten und p_N(1). Dafür braucht (Satz K15, bedingt) nur das Eckprodukt (ein Skalar je n), die kopfweise adjungierte Eckform und die
  Stegdämpfung; (S_γ), (H1), (H*_G) und die Randlokalisierung entfallen. Neue Lücke: H1 für die Sekanten braucht Kontrolle von P in der Kegelzone bei L.
  |μ| = 2 taugt nicht als Testmodell für Normen (J̃_n fast singulär).
- **Warum:** Die Unbekannten auf D̄ nahe L sehen die drehende Zone; Punkte in fester Entfernung sehen nur die n-unabhängige Phase der Spirale.
- **Was es lehrt:** Die Ecke trägt zwei Rollen — Identifikation mit Kneser (über die Referenz f_K, nicht über Eindeutigkeit nötig) und den polynomialen
  Fehler; nur die zweite braucht eine gleichmässige Eckaussage, und die war die falsche (zu starke) Zielgrösse.
- **Konsequenz:** Ziel umstellen auf K_{r₀}; als nächstes (K): Eckprodukt (i″) und die Sekanten-Eckbedingung; T1 (V) entscheidet die Kette.

## [K] 12.10. 09:15 — Versuch 23 (Auftrag Leitung 09:08): Ist das Eckprodukt eine Invariante?
- **Erwartung:** Die Paarung aus K11 plus Anpassung legt das Produkt der Amplituden fest (Wronski-artig), vermutlich nur bis auf eine globale Konstante.
- **Ergebnis:** Ja, und schärfer: exakt ist P_n = δ_n(L)Λ_n(L) = K_n(L,L) − K_{n−1}(L,L) mit dem Rayleigh-Quotienten K_n(L,L) = ẽ_L*J̃_n⁻¹ẽ_L (Antwort auf den
  Cauchy-Kern bei L). Gemessen P_n = aCL/(n d_n)(1 + O(n^{−2})), C ≈ −0.0102 − 0.5984i, glatt in n (die schnelle Phase hebt sich auf); K_n(L,L) ≈ C·L/d_n ~ n^a.
  Erklärt durch eine Euler-Identität des Eckkerns (Zuwachs in n = Skalenableitung, E annulliert den Cauchy-Kern, K11 ⇒ E𝒢 = C·F_P⊗F_A). Gs B_n-Gesetz folgt daraus.
- **Warum:** Die Cauchy-Quelle hat eine feste Normierung; die Einzelamplituden der Köpfe sind frei, ihr Produkt ist ein Zuwachs einer festen Grösse.
- **Was es lehrt:** Die einzige globale Zahl der K_{r₀}-Fassung ist ein Fisher–Hartwig-Problem für EINEN Rayleigh-Quotienten mit komplexem Exponenten a — robust gegen
  die Spirale. Die Eckbedingung (EB_N) an die Sekanten lässt sich auf Ruhefenster beschränken, die die kopfweise Kontur ohnehin benutzt.
- **Konsequenz:** Ziel (R_L) |ẽ_L*J̃_n⁻¹ẽ_L| ≤ poly beweisen (Fisher–Hartwig/Wiener–Hopf an der Eckfrequenz); V misst (EB_N) an den echten p_N.

## [K] 12.10. 09:23 — Versuch 24 (E1): Ist die Eckkonstante global wie bei Fisher–Hartwig?
- **Erwartung:** C ist eine globale Konstante; für (R_L) braucht man eine schwere Fisher–Hartwig-Analyse.
- **Ergebnis:** Nein, C ist lokal: die Antwort auf den Eck-Pol hat das Innenprofil C·∫₀¹F_P(λ̃_Psζ)ds (Kontinuumsform der Christoffel–Darboux-Summe), und die Anpassung an
  den äusseren Pol −1/((μ+1/μ)(x−L)) gibt C = −λ̃_PΓ(½+iα)/((μ+1/μ)Γ(½+2iα)) = −0.010118 − 0.598298i; gemessen (1/n-extrapoliert) −0.01012 − 0.59830i. Profil auf < 1 %.
- **Warum:** Mellin bei s = 1 verbindet Innen (Gammafaktoren von F_P) und Aussen (Residuum des Bulk-Pols); die Skalierung der Köpfe macht die Summe über m zum Skalenintegral.
- **Was es lehrt:** Die letzte „globale Zahl“ von K15 ist lokal. Die ganze K_{r₀}-Fassung ruht auf lokalen/äusseren Aussagen mit expliziten Eckfunktionen.
- **Konsequenz:** (R_L) als angepasste Asymptotik einer Antwort beweisen; kritisch ist die Summe über m ≈ n (Phasenstruktur statt Beträge) oder eine exakte String-Gleichung.

## [K] 12.10. 09:27 — Versuch 25 (E1, Vorschlag Leitung): String-Gleichung exakt herleiten?
- **Erwartung:** n·P_n − a·K_n = κ₀ + O(1/n) folgt exakt aus K5 + Euler-Identität, κ₀ hat eine geschlossene Form.
- **Ergebnis:** Nein. K_n(L,L) = CL/d_n + K₀(n) mit K₀(n) → K₀∞ ≈ 0.70 wie n^{−0.106}; also κ₀ = −aK₀∞ ≈ −0.63 mit Korrektur O(n^{a−1}), und K₀∞ ist die Summe über alle
  frühen Köpfe — global. Die Gleichung ist die diskrete Euler-Identität mit globalem Achsenabschnitt, nicht exakt herleitbar ohne die Eckkern-Konvergenz.
- **Warum:** Euler-/Skalenidentitäten gelten nur asymptotisch auf der Kegelskala; der konstante Term kommt von Graden, die nicht auf der Kegelskala leben.
- **Was es lehrt:** Für (R_L) braucht man die Grössenordnung, nicht κ₀. Der richtige Bootstrap verlangt glatte Zuwächse als Induktionsannahme (sonst verliert man n^{0.106}).
- **Konsequenz:** (R_L) ⇐ Lemma ρ (Mellin–Parseval) ⇐ (ii) + (iii). Alle Restpunkte hängen an einem Anpassungslemma für Köpfe — das ist der nächste Block.

## [K] 12.10. 09:39 — Versuch 26: Mellin–Parseval-Bootstrap für (R_L)
- **Erwartung:** Mit einer zusammengesetzten Vergleichslösung und kopfweisen Paarungen lässt sich (R_L) direkt beweisen; das kritische Stück (m ≈ n) per Mellin–Parseval.
- **Ergebnis:** Der Bootstrap selbst ist sauber und bewiesen (Satz K16), sobald die relativen Paarungen klein sind; die Induktionsannahme muss „glatte Zuwächse“ sein,
  dann geht bei m ≈ n nichts verloren. ABER: Lemma ρ braucht die adjungierte Eckform (ii), die Stegdämpfung (iii) und eine polynomiale primale Parametrix (Typ (v)) —
  die naive Komposition scheitert, weil V in seinen Wachstumssektoren exponentiell wächst und T_n das verschmiert.
- **Warum:** (R_L) ist ein Eckwert der Inversen; jede Darstellung der Inversen, die nicht über die Norm geht, läuft über die Köpfe (CD-Summe). Ohne Kopfkontrolle keine Paarungen.
- **Was es lehrt:** (R_L) war nie unabhängig von (ii)/(v); die Architektur hat jetzt drei Restpunkte gleicher Art statt vier verschiedener.
  Nebenbei: die Integraldarstellung F_β = E_Γ[exp(−uT^{−iα})] macht die Sektorasymptotik elementar (Watson statt Braaksma).
- **Konsequenz:** Weiter mit (ii) als Anpassungslemma: innere Hälfte beweisen, äussere Hälfte sauber als Hypothese an W formulieren.

## [K] 12.10. 09:39 — Meta-Reflexion (nach Versuch 25–26)
- Zweimal hintereinander hat ein „kürzester Weg“ (exakte String-Gleichung, direkter Bootstrap) auf denselben Kern zurückgeführt: Kopfasymptotik mit Fehlerschranke.
  Das ist ein Muster: jede Eckgrösse, die eine Inverse enthält, führt über die CD-Summe auf Köpfe. Künftig zuerst fragen, ob ein Vorschlag die Inverse umgeht — wenn nicht,
  gleich auf das Anpassungslemma zielen.
- Gut lief: exakte Identitäten (K5, CD-Summe) vor der Analysis prüfen; schnelle numerische Gegenprobe (k48) bevor etwas in den Bericht geht.
- Zu vermeiden: Zeitstempel geschätzt eingetragen (09:30 statt 09:27) — korrigiert; immer erst `date`.

## [K] 12.10. 10:01 — Versuch 27: eine messbare Aussenhypothese für (ii), (iii), (v)
- **Erwartung:** (ii) als Anpassungslemma: innere Hälfte über K10/K11 beweisen, äussere Hälfte als Hypothese; (iii) und (v) vielleicht auf dieselbe Hypothese.
- **Ergebnis:** Besser als erwartet in der Form, ehrlich gleich schwer im Kern. (1) Rang-1-Struktur: alle vier Restpunkte sind Faktoren des Kernzuwachses ΔK_n auf
  K_{r₀} × Γ — EINE Ungleichung (AH). (2) Für Schranken ist die innere Hälfte nur das Maximumprinzip; K10/K11 sagen den Wert auf dem Kegelkreis voraus (3.63 / 13.7
  gemessen gegen 3.62 / 13.6). (3) Die Fensterränder sind n-unabhängige Log-Spiralen ⇒ EINE feste Kontur durch L für alle Köpfe; längs der Spiralen ist der Kopf
  nirgends grösser als in der Ecke. (4) Ab n ≈ 600 liegt das Maximum von |ΔK_n| genau in der Ecke.
- **Warum:** K15 braucht nur Beträge, und Beträge eines Rang-1-Kerns faktorisieren; das Eckmuster dreht sich, aber seine Fenster sind im physikalischen Raum fest.
- **Was es lehrt:** Vor einer Anpassungsrechnung prüfen, ob die Anwendung die asymptotische Form braucht oder nur eine Schranke — hier nur die Schranke.
  Mein Ansatz von 09:39 (K16, Lemma ρ, K-Param) war für K15 ein Umweg; er bleibt für das Gesetz von P_n.
- **Konsequenz:** (AH) an V (E4); Rest des Beweisprogramms = starke Asymptotik aussen + Eckdominanz der Raten; die primale Rate H_P ist noch nicht formuliert.

## [K] 12.10. 10:23 — Versuch 28: Hauptsatz-Dokument — und dabei K17 zurückgenommen
- **Erwartung:** Das Dokument fasst nur zusammen: (AH) + K17 + Kette ⇒ Hauptsatz; (EB_N) folgt aus (AH) oder ist zweite Annahme.
- **Ergebnis:** Beim Ausschreiben der Datenseite fiel K17 durch: Pol-Quellen auf dem Schlitz nahe L verlangen f(Kontur) ⊂ D̄, das feste Spiralfenster dreht
  heraus (Grenze N ≈ 1e8–1e11). Das Dokument stützt sich darum auf die alte, von V geprüfte Kette mit (H) auf K_{r₀} umgestellt — dort braucht es keine (EB_N).
  Exponentiell (k50) ist die Ecke harmlos, die Γout-Seite braucht das Green-Gewicht.
- **Warum:** In §22.2/§28.4 hatte ich die Datenseite „wie in K15“ übernommen, ohne die Pol-Quellen nahe L zu prüfen. Die Drehung, die ich auf der Kopfseite
  beseitigt hatte, sitzt auf der Datenseite wieder.
- **Was es lehrt:** Ein Satz ist erst fertig, wenn BEIDE Seiten (Köpfe und Daten) im selben Gebiet ausgeschrieben sind. Ein Dokument für den Autor zwingt dazu — gut so.
  Meldung 10:01 war in diesem Punkt zu stark; das muss die Leitung wissen.
- **Konsequenz:** Hauptsatz unter (H)^{r₀}; (AH) als Eckstruktur; nächster Block exponentielle Ebene mit Green-Gewicht (Γout-Seite).

## [K] 12.10. 10:38 — Versuch 29: Reicht in K_{r₀} die exponentielle Ebene?
- **Erwartung:** Fester Abstand r₀ ⇒ primale Köpfe auf K_{r₀} exponentiell klein ⇒ e^{o(n)} überall geschluckt ⇒ (H)^{r₀} ganz über Ws Potentialtheorie.
- **Ergebnis:** Teilweise. Die primale Reserve ist im Inneren nur polynomial (n^{−3.05}, Potenzsektor, Exponent explizit) — die Ecke reicht algebraisch in ganz D hinein.
  Aber: (1) mit Daten J(f − q_N) statt abgeschnittener Pole (K19) wird die Datenseite klassisch und die Kontur darf aussen liegen; (2) auf der fernen Kontur (Linse
  s = 0.03, Γout-Bogen) hat die Adjungierte eine strikte exponentielle Reserve (η ≈ 0.013). Damit genügt global die n-te-Wurzel-Ebene, polynomial nötig nur nahe L.
- **Warum:** Die Ecke ist ein algebraischer Punkt (Potenzsektoren beider Profile); alles fern von ihr ist exponentiell getrennt, sobald man eine Kontur mit Reserve findet,
  und die abgeschnittenen Pole hatten diese Kontur künstlich in D hineingezwungen.
- **Was es lehrt:** Die Datendarstellung bestimmt, wo die Kontur liegen darf. Taylor-Rest zwingt nach f⁻¹(D); Approximation auf einem grösseren Kompaktum gibt die
  Freiheit zurück. Erst messen (k51 in 1 Min.), dann die Frage beantworten — die naheliegende Erwartung „exponentielle Reserve“ war falsch.
- **Konsequenz:** Programm „global exponentiell + lokal polynomial“; Familie bleibt bei (H)^{r₀}; nächster Block: Fern-Raten aus Ws Problem.

## [K] 12.10. 10:48 — Versuch 30: Fernraten aus Ws Problem — und eine zweite Korrektur
- **Erwartung:** Beide Fernraten auf n-ter-Wurzel-Ebene mit Reserve, die Primale vielleicht über Dualität aus der Adjungierten.
- **Ergebnis:** Adjungiert ja (gemessen; als Satz ⇐ W12′, und W13 stört nicht, weil nur die relative Rate zählt). Primal nein: im Innern ~ n^{−3.05}, am Rand nur eine
  vorübergehende Exponentialphase; mit Oberhalbstetigkeit ist H_P ≡ 0 auf D̄. Meine Aussage von 10:38 („primal fern exponentiell“) war vorasymptotisch.
- **Warum:** Die Ecke strahlt über das Potenzglied des inneren Profils in ganz D hinein — algebraisch, nicht exponentiell. Für die Primale ist die Ecke ein Multipol.
- **Was es lehrt:** Endliche-n-Raten am Rand können vorübergehend sein; immer auch Punkte knapp innen messen (0.95ρ, 0.98ρ) — dort sieht man das Kippen.
  Zweimal heute hat eine „Reserve“ nicht gehalten (K17-Datenseite, primale Randrate); beide Male hat eine schnelle Messung die Frage entschieden.
- **Konsequenz:** (AH_{r₀}) zerfällt in einen n-te-Wurzel-Teil (adjungiert fern, W12′) und polynomiale Teile (lokal nahe L; global nur die Primale = Eck-Pol-Antwort Π_n s).

## [K] 12.10. 10:52 — Meta-Reflexion: Zirkelprüfung zu (G-P) (Auftrag Leitung 10:50)
**Ist (G-P) leichter als der Hauptsatz, oder dieselbe Aussage in anderer Kleidung?** Überwiegend dieselbe — genauer: (G-P) ist der ECKTEIL der alten Hypothese, an der
schlimmsten Quelle. Begründung in drei exakten Zeilen:
- Π_N k_L = J_N⁻¹T_NJk_L und Jk_L = (μ + 1/μ)k_L + (Pole an den fernen Urbildern von L) + Polynom, also Π_N k_L ≈ (μ + 1/μ)·Π_N s + (ferne Antworten). (G-P) ist damit
  (H-ii)^{r₀} im Grenzfall ζ → L (Quelle IN der Ecke), eingeschränkt auf eine Quelle und einen Operator.
- (G-P) ist der Schnitt t = L von (AH_{r₀}) (R_P = |ΔK_n(z, L)|/|P_n|) und, in ℓ¹ statt sup_{K_{r₀}}, Ps alte Bedingung R_N(L) ≤ CN^κ (P-I(ii)).
- Ein Beweis über die äussere Entwicklung läuft wie K16: Induktion in n mit z-gleichmässiger Annahme R_P(m) ≤ Cm^κ — DIESER Schritt ist neu und trägt
  (die Zuwächse des expliziten Teils s + c·d_n^{s₀−1}ψ sind ~ n^{−3.16}). Er braucht aber Lemma ρ, also eine POLYNOMIALE Parametrix ũ_n der Eck-Pol-Antwort mit
  kleinem Residuum (K-Param, §27.2). Genau dort sind die Eckversuche gescheitert: Parametrix aus der Bulk-Inversen (k07/k08, Fehlskalierung e^{0.47N}),
  T_n verschmiert die Wachstumssektoren des inneren Profils (§27.2).
**Was s leichter macht:** ein Operator statt der Familie; Quelle fest in L; inneres Profil C·V/d_n explizit (C auf 1e-5, §25); äussere Korrektur explizit bis auf ψ
(Exponent Re(a s₀^P) = 3.053 gemessen); Budget grosszügig (κ_P < 2.1, wahr −3.05). Das ist echte Vereinfachung — aber keine, die den gescheiterten Schritt umgeht.
**Was neu wäre:** nur eine Parametrix in DUALER, gewichteter Form — das Residuum darf dort exponentiell gross sein, wo die adjungierten Köpfe eine Reserve haben (FA,
η ≈ 0.013 auf Linse/Γout-Bogen), und muss nur nahe L und am Steg klein sein, wo die Profile explizit sind. Das ist der Deift–Zhou-Gedanke „Kleinheit in der
Norm, die die Köpfe sehen“. Ob das trägt, ist offen; es ist der einzige Hebel, den ich sehe, der nicht schon gescheitert ist.
**Urteil:** nicht überzeugend leichter. (G-P) wird nicht direkt angegriffen, sondern als die EINE benannte primale Hypothese geführt (sie ist schwächer als (H-ii)^{r₀}:
eine Quelle, ein Operator, nur fern von L). Erreichbarer Satz: „Adjungierte fern über Ws Potentialtheorie (FA ⇐ W12′) + lokale Eckformeln + K19 ⇒ lineare
Konvergenz auf K_{r₀}, unter (G-P)“ — mit den lokalen polynomialen Aussagen ehrlich als eigener Posten (siehe BERICHT §32).

## [K] 12.10. 11:00 — Versuch 31: (FA) über Vergleichsprinzip mit Superlösung
- **Erwartung:** Mit Normierung in L (statt am Sockel) umgeht ein Collatz–Wielandt-artiges Vergleichsprinzip Ws W13; Ws Modell liefert die Superlösung.
- **Ergebnis:** Das Gerüst (Überschuss-Weitergabe über starkes Maximumprinzip + tropische Bilanz, Abbruch in L) stimmt. Aber: (O1, bewiesen) in der Ebene gibt es
  keine superharmonische Vergleichsfunktion mit der Singularität der Köpfe (Laplace-Masse), also muss die Weitergabe an Nähten über die nichtlokale Bilanz
  laufen; (O2) dazu braucht man die Dominanzstruktur der Konfiguration — über λ (Klebung durch die Identität) ist die naheliegende Glättung max(u, w∘f) nicht
  verträglich. Das ist Ws §22.3; der Weg ist gleichwertig zu W12′.
- **Warum:** Die Köpfe haben Nullstellen (Masse 2π), und deren Lage IST die unbekannte Gestalt; jede Vergleichsfunktion muss sie vorwegnehmen.
- **Was es lehrt:** Massen-/Ladungsbilanz vor jedem Superlösungs-Ansatz prüfen (eine Zeile hätte die Ebenen-Fassung sofort ausgeschlossen).
- **Konsequenz:** (FA) bleibt benannte n-te-Wurzel-Hypothese (⇐ W12′); die einseitige, in L normierte Fassung von W12′ ist das präzise Ziel für später.

## [K] 12.10. 11:16 — Versuch 32: Familie in K_{r₀} über Keimabhängigkeit + K4
- **Erwartung:** Keimabhängigkeit (Matrizen nur vom Keim) + K4 (kein Multiplikationsüberlauf) führen F2/F3 auf f_K zurück, wie bei der Linse für (S_γ).
- **Ergebnis:** Nein. Für die Matrizen ja (wie §10), aber F3 braucht P auf der Linse am Steg, und dort wachsen die primalen Köpfe exponentiell (+0.016 bei s = 0.03, k51c) —
  P − f ist dort exponentiell gross. Keim-Repräsentanten (Carathéodory–Fejér) retten bis N ≈ 6300–7800. Eine lineare Datendarstellung ohne Linse gibt es nicht
  (Zange: gleichmässige Daten brauchen die Linse, Taylor-Daten scheitern auf Γout, Mischen widerspricht sich in L).
- **Warum:** Die Adjungierte braucht am Steg eine Linse (s > 0.0118), die Primale wächst auf derselben Linse um ln(1+s) − 0.012 > 0 — die Bilanz am Steg ist knapp
  positiv (+0.016 − 0.013), und nur Kohärenz kann sie kippen.
- **Was es lehrt:** Die K_{r₀}-Fassung entschärft die Ecke, nicht den Steg. Die Familienlücke ist eine Steg-Aussage fern von L.
- **Konsequenz:** (Fam-koh)^{r₀} als benannte Familienhypothese (Repräsentant auf der Linse); nahe L reicht (EP⁺). Weiter mit der Parametrix, sobald Vs (b) vorliegt.

## [K] 12.10. 11:20 — Versuch 33: (S_E) über den direkten Satz von Dzjadyk-Typ
- **Erwartung:** Klassische Theorie gibt s_E = Re s₀·α_E an der Ecke; zu prüfen sind die Voraussetzungen (Quasischeibe, Glattheitsklasse) für die Spiralecke.
- **Ergebnis:** Ja: E ist nahe L das quasikonforme Spiralbild eines Winkels (Quasischeibe), f_K ∈ C^{5,0.417}(E) aus (K3); s_E = 5.417·α_E, α_E = (u₁ − u₂ − 38.3°)/180°.
  k53 (synthetische Ecke) trifft den Exponenten auf 0.01–0.02. Zitierter Satz nicht nachgeschlagen (Belyi 1977, Andrievskii–Belyi–Dzjadyk 1995).
  Nebenbefund: E ⊃ D̄ stimmt nahe L nicht für alle Skalen (Spiraldrehung); richtig ist E ⊃ D̄ ∖ U_{r₀}.
- **Warum:** Die Rate ist eine reine Eckgrösse (Aussenwinkel × Singularitätsexponent); die Spirale ändert nur Faktoren 1 + O(T²).
- **Was es lehrt:** Klassische Bausteine erst mit einem kleinen synthetischen Test absichern, bevor man sie zitiert — kostet zwei Minuten.
- **Konsequenz:** (S_E) ist bedingt bewiesen (zitierter Satz + Geometrie von E); Budgetreserve s_E − γ ≈ 3.1 bei Fenster (80°, −90°).

## [K] 12.10. 11:37 — Versuch 34: zitierten Satz nachgeschlagen, Geometrie von E geprüft
- **Erwartung:** Belyi/ABD im Wortlaut finden; E als Quasischeibe numerisch bestätigen.
- **Ergebnis:** Wortlaut gefunden und gelesen in Andrievskii–Pritsker–Varga (arXiv math/0005032, Satz 2): für E ∈ H* (enthält Quasischeiben), f ∈ A^r(E):
  |f − p_n| ≤ c ρ_{1/n}^r ω_{f^{(r)},k}(ρ_{1/n}) auf ∂E — ohne Zusatzbedingung an ω (meine „Dini“-Annahme war unnötig). Geometrie (kubisch, k54): univalent, ohne Löcher,
  Schnitt meidet E, Öffnung 132° (Vorhersage 131.7°), Dreipunkt-Verhältnis 1.07. Ein Skriptfehler (Bogenreihenfolge) gab zuerst 3e7 — erkannt und behoben.
- **Warum:** Klassische konstruktive Funktionentheorie ist genau für solche Ecken gebaut; die Spirale ist quasikonform harmlos.
- **Was es lehrt:** Wenn ein Prüfwert absurd ist (3e7), zuerst die Prüfung prüfen (Reihenfolge/Duplikate), dann die Mathematik.
- **Konsequenz:** (S_E) hängt nur noch an der f_K-Version des Geometrie-Checks (V, endlich).

## [K] 12.10. 12:33 — Versuch 35: Vortest (b) selbst — Lemma ρ mit expliziter Parametrix
- **Erwartung:** Die Komposition „aussen s, innen C·V/d_n“, als Polynom vom Grad n realisiert (Taylor oder Approximant auf E, notfalls dual gewichtet), hat kleine relative
  Paarungen ρ_{m,n}; dann (EP) über K16.
- **Ergebnis:** Nein. Taylor: Kompositionsüberlauf (ρ bis 1e6). Approximant auf E: Fit unmöglich, ρ = O(1), auch gewichtet. Grund (k55c): im primalen Wachstumssektor
  (ausserhalb D nahe L, also auf Σ₁ und f(Γout)) folgt die echte Antwort dem inneren Profil nur bis |ζ| ≈ 2–3. Danach dominiert die äussere WKB-Wachstumslösung
  (Faktor 22 bei |ζ| ≈ 29). Nebenbefund mit Folgen: Fenster 1 ist primal ein Wachstumsfenster (k51d) — §24.5/§34.2(a) korrigiert, F3 nahe L ist nicht lokal.
- **Warum:** Die Adjungierte hat ihre Reserve nur ausserhalb D auf der Γin-Seite, die Primale wächst genau dort. Dieselbe Spannung steckt in der Familie (§34),
  am Steg und jetzt in der Parametrix.
- **Was es lehrt:** Ich habe zweimal ungeprüft „Fenster 1 ist primal ruhig“ übernommen (§24.5). Eine Zeile Sektorrechnung (arg λ̃_P + arg ζ) hätte es gezeigt.
  Behauptungen über Sektoren immer mit einer Punktmessung belegen.
- **Konsequenz:** Alle offenen Stellen sind eine: primale Grössen ausserhalb D auf der Γin-Seite. Vorschlag (i): über die exakte Identität J u_n = T_nk_L + R_n die Werte
  ausserhalb D aus D̄ und dem Abschneide-Rest R_n bestimmen und prüfen, ob R_n auf Γout polynomial ist.

## [K] 12.10. 12:35 — Vorhersage vor Messung (i): Abschneide-Rest R_n = (I − T_n)Ju_n auf dem Γout-Bogen
Exakt gilt R_n(x) = u_n(f(x)) + m(x)u_n(x) − T_nk_L(x) (Galerkin: T_nJu_n = T_nk_L). Auf dem Γout-Bogen fern von L sind die Glieder m·u_n ≈ m·s und T_nk_L beschränkt.
Also ist R_n(x) ≈ u_n(f(x)), und f(x) liegt ausserhalb D. Ein Polynom vom Grad n, das auf D̄ eine Funktion mit Pol in L nähert, hat Koeffizienten ~ ρ^{−k} und wächst
ausserhalb D wie (|w|/ρ)^n. k55c passt dazu: n = 300, |d_n·u_n| ≈ 2.1e4 (also |u_n| ≈ 4.7e6) an f(y), y auf Γout im Abstand 0.2 von L, während s dort nur ≈ |c_s|/0.244 ist.
**Vorhersage:** R_n wächst auf dem Γout-Bogen exponentiell, mit der lokalen Rate (1/n)ln|R_n(x)| ≈ ln(|f(x)|/ρ). Das Maximum über den Bogen ist ≈ ln 1.624 = 0.485
(kubisch, max|f|/ρ = 1.624, k54). Nahe L geht die Rate gegen 0 (|f(x)|/ρ → 1). Auf Γin ist R_n klein (f(x) ∈ D).
**Folge, falls bestätigt:** (i) trägt nicht. „Primal ausserhalb D“ ist über die Identität genau der Kompositionsüberlauf (‖J_N‖ ~ e^{0.47N}).
(EP) bleibt dann benannte Hypothese, und der Nachtrag in KERNPROBLEM.md wird geschrieben.

## [K] 12.10. 12:44 — Versuch 36: Hebel (i), Abschneide-Rest auf Γout — Vorhersage halb richtig
- Qualitativ getroffen (exponentiell, Maximum bei 0°, gegen 0 an L und L̄; Identität tautologisch). Die Rate war um den Faktor 2 zu hoch vorhergesagt:
  gemessen ≈ ½·ln(|f|/ρ), nicht Bernstein-Walsh. Lehre: BW ist für Galerkin-Lösungen nicht scharf; die Koeffizienten sind flach, die Hälfte kommt aus Phasen.
- Die Messung dauerte Sekunden (n ≤ 400). Die Vorhersage aus der Identität hätte die Entscheidung schon allein getragen; die Messung macht sie quantitativ.
- Meta: Viermal hintereinander (Steg §34, Fenster 1 §37.1, Parametrix §37, R_n §38) endet ein Hebel an derselben Stelle. Weitere Hebel derselben Art
  (eine Seite getrennt abschätzen) sind nicht sinnvoll. Ein neuer Ansatz müsste das Produkt primal × adjungiert als ein Objekt behandeln (KERNPROBLEM §4).

## [K] 12.10. 12:54 — Zeitfenster „Produkt primal × adjungiert am Steg“ (bis ≈ 14:55): Unterschied zu H-A, Vorhersage, Abbruchkriterium
**Was anders ist als H-A (§16).** H-A prüfte an der ECKE, ob beide Köpfe die HOMOGENEN Eckgleichungen erfüllen (dann wäre m·Λ·δ·dy/y invariant). Das scheiterte,
weil der primale Kopf nur innen und der adjungierte nur aussen homogen ist. Am STEG ist es umgekehrt: Dort sollten beide Köpfe von ihren QUELLEN beherrscht sein,
δ_n ≈ (t/ρ)^n/m(t) und Λ_n ≈ (ρ/t)^n/m(t). Das Produkt der Quellen ist exakt 1, ohne radiales Wachstum. Das erklärt (Skizze) die Aufhebung der Steigungen (+0.94/−0.93).
Gefragt ist also nicht eine Erhaltungsgrösse homogener Lösungen, sondern die Grösse der Abweichungen von den Quellen.
**Exakte Zerlegung (eine Zeile):** Π_n(t) := ΔK_n(t,t)/P_n = δ_n(t)Λ_n(t)/(δ_n(L)Λ_n(L)) = A_P(t)·A_A(t)/(m(t)²·ℓ_nn·P_n), mit A_P := m·δ_n·(ρ/t)^n, A_A := m·Λ_n·(t/ρ)^n.
Dazu gilt exakt: Das Kreismittel von δ_nΛ_n ist auf JEDEM Kreis Σ_k c_{n,k}ℓ_{n,k} (nullter Laurent-Koeffizient), also unabhängig vom Radius.
**Vorhersage (k57, kubisch, n = 150…600, θ = 120°/150°/180°, s = −0.03 … 0.05):**
1. Die n-te-Wurzel-Rate h_Π von Π_n am Steg ist radial konstant (|∂_s h_Π| ≤ 0.01) für s ∈ [0, 0.05], und positiv ≈ +0.003 ± 0.002 bei 180°.
2. A_P ist am Steg für s ≥ 0 O(1) (Rate 0 ± 0.002): Die Primale ist quellbeherrscht. Der Überschuss sitzt in A_A (≈ +0.004) oder in 1/ℓ_nn. Welches von beiden, weiss ich nicht.
**Abbruchkriterium (vorher festgelegt):**
- Der Weg TRÄGT (weiterverfolgen), wenn (a) h_Π → 0, d. h. die Rate im Fenster 450 → 600 ist < +0.001 und fällt über die Fenster, oder (b) der ganze Überschuss in ℓ_nn sitzt
  (ein Skalar, ortsunabhängig) und A_P, A_A am Steg polynomial sind.
- Der Weg ist GESCHEITERT, wenn h_Π(180°) im Fenster 450 → 600 ≥ +0.002 bleibt UND der Überschuss in einem ortsabhängigen Faktor sitzt (A_A oder A_P, Variation über
  θ ∈ [120°, 180°] > 0.002). Dann bringt „Produkt als ein Objekt“ am Steg nichts, was die getrennten Raten nicht schon sagen.
- Spätestens 14:55 wird abgebrochen, mit dem Stand.

## [K] 12.10. 13:08 — Versuch 37: Steg-Produkt — Vorhersage in beiden Teilen falsch, Ergebnis besser als erwartet
- Vorhersage 1 (+0.003) falsch: Das war ein Artefakt aus meinem eigenen Vergleich von Steigung und Niveau (12:44). Lehre: Raten nur gleichartig vergleichen
  (beide als Steigung im selben n-Fenster). Das hätte ich vor dem Nachtrag in KERNPROBLEM prüfen müssen. Jetzt korrigiert, an allen Stellen markiert.
- Vorhersage 2 (Primale am Steg quellbeherrscht) falsch: Beide Köpfe sind exponentiell kleiner als ihre Quellen, und ℓ_nn gleicht genau aus.
- Positiv: Das Produkt ist am Steg nicht nur polynomial, sondern n-unabhängig und explizit (κ(t), Christoffel-Darboux-Form). Damit hat „Produkt als ein Objekt“ getragen,
  nur anders als gedacht: nicht als Erhaltungsgrösse homogener Lösungen (H-A), sondern als Kern-Diagonale = logarithmische Ableitung des Kopfes.
- Das Zeitfenster war nach ≈ 15 min entschieden. Ich melde mich früher, wie von der Leitung erlaubt.

## [K] 12.10. 13:13 — Zeitfenster „exakte Christoffel-Darboux-Identität für J̃_n“ (≈ 3 h, bis ≈ 16:10): Ansatz, Vorhersage, Abbruchkriterium
**Ansatz (vor der Rechnung, Herleitung auf Papier).** G40 ist zweiter Ordnung; für die Inverse taugt sie nicht direkt (K §13.2). Ich nehme stattdessen die zwei
EINZELNEN Verschiebungen erster Ordnung, deren rechte Seiten nicht klein im Rang, aber strukturiert sind:
- (V1) J S − S J = T_n M_{f−x} C_f − (T_n f^{n+1}) e_nᵀ (Kompositionskanal),
- (V2) J S − T_F J = T_n M_{(x−f)m} − (T_n f^{n+1}) e_nᵀ (Multiplikationskanal).
Für X = J_n⁻¹ gilt S X − X S = X(V1)X. Eingeschoben zwischen e(z)ᵀ und ê(t) gibt das (Herleitung in §40.1):
  **(CD1) (z − t)·K_n(z,t) = q_{n+1}(z)·λ_n(t) − t·ν_n(z) + Ξ_n(z,t),**
mit q_{n+1} = monischer primaler Kopf der Ordnung n+1, λ_n = adjungierter Kopf der Ordnung n, ν_n = J_n⁻¹·1 und dem Kanalterm Ξ_n = J_n⁻¹T_n[(f − x)·K_n(f(x),t)].
Das ist eine Christoffel-Darboux-Formel (Kopf n+1 mal Kopf n durch Differenz) plus zwei Korrekturen.
**Vorhersage.** (1) (CD1) gilt auf Arbeitsgenauigkeit (es ist eine Rechnung). (2) Am Steg sind q_{n+1}λ_n, tν und Ξ alle O(1) und von n unabhängig bis auf o(1);
damit erklärt (CD1) die gemessene Form K_n(t,t)m(t) = t(log δ_n)′ + Δ(t). (3) Ξ ist NICHT klein; er trägt einen n-unabhängigen Teil von Δ.
**Abbruchkriterium.** Erfolg: (CD1) exakt, UND am Steg sind ν_n, Ξ_n(t,t), ∂_zΞ_n(t,t) zwischen n = 150 und 300 auf ≤ 1e-3 gleich. Dann ist das Produkt am Steg
durch eine exakte Formel mit n-unabhängigen Korrekturen erklärt. Teilerfolg: (CD1) exakt, aber Ξ wächst mit n. Dann ist der Kompositionskanal genau die offene Grösse.
Misserfolg: (CD1) gilt nicht (Herleitungsfehler, der sich nicht beheben lässt). Spätestens 16:10 abbrechen.

## [K] 12.10. 13:19 — Versuch 38: Christoffel-Darboux — Vorhersage 1 richtig, 2 halb, 3 richtig
- (1) (CD1) exakt: ja (1e-176). Die Herleitung über die Verschiebung ERSTER Ordnung des Kompositionskanals war der Schlüssel. G40 (zweite Ordnung) direkt zu invertieren
  hätte nicht funktioniert (§13.2). Der Hinweis der Leitung auf Gohberg-Semencul hat die Richtung gegeben; die Struktur ist aber „Toeplitz-CD + Dyson-Kanalterm“, nicht Rang 2.
- (2) Halb: ν, Ξ(t,t), B sind n-unabhängig, ∂_zΞ aber nicht. Das Wachstum ist genau eine Projektion auf den Generator; nach Umnormierung mit 1/B(t) ist alles n-unabhängig.
- (3) Ξ ist nicht klein (≈ −0.57·t/m am Steg): richtig.
- Nach dem Abbruchkriterium formal ein Teilerfolg, in der Sache eine vollständige Erklärung der Messung aus §39. Offen bleibt der Beweis der n-Unabhängigkeit.
  Ξ ist der Kompositionskanal, also wieder (KK′).
- Zeit: 13:13–13:18 für Herleitung und Messung (die Herleitung stand vor der ersten Rechnung). Früher Bericht, wie von der Leitung erlaubt.

## [K] 12.10. 13:23 — Zeitfenster „Ξ_N als Dyson-Term: Kontraktion?“ (≈ 3 h, bis ≈ 16:20): Ansatz, Vorhersage, Abbruchkriterium
**Ansatz.** In Matrixform ist (CD1) die Riccati-Gleichung [S, X] = X E X (X = J_n⁻¹), mit der Eichung „erste Spalte ν = Xe_0 gegeben“. Die Dyson-Abbildung ist
𝒯(X) = 𝒮_ν⁻¹(X E X); 𝒮⁻¹ löst [S, Y] = R spaltenweise (Summation längs Diagonalen = Division durch z − t). Ihre Linearisierung am echten X ist
L[δX] = 𝒮₀⁻¹(δX·E·X + X·E·δX) auf {δX : δX e_0 = 0}. Zwei Varianten: mit dem vollen E, und mit nur dem Kompositionskanal E_c (Rang-1-Teil als Daten).
**Erster, normunabhängiger Test:** der Spektralradius ρ(L). Ist ρ(L) ≥ 1, gibt es KEINE Norm, in der 𝒯 lokal kontrahiert; eine Konturwahl kann daran nichts ändern.
**Vorhersage:** ρ(L) > 1, und zwar wachsend mit n. Begründung (Skizze): [S, X] = XEX ist in J-Koordinaten LINEAR, [J, S] = E. Dessen Lösungen sind die ganze
Familie C_f + M_m (m beliebig). Die Dyson-Iteration in X-Koordinaten hat Richtungen, in denen sie fast die Identität ist (δX = X M X mit Toeplitz-M). Ausserdem
steckt in E der Kompositionskanal mit dem Überlauf e^{0.47n}.
**Abbruchkriterium:**
- Der Weg ist GESCHEITERT, wenn ρ(L) ≥ 1 für n = 10 … 40 (beide Varianten) und nicht fällt. Dann wird nur noch der Grund festgehalten.
- Er TRÄGT (weiter zur Norm-Konstruktion auf Konturen), wenn ρ(L) ≤ 1 − c mit c n-unabhängig.
- Ist ρ(L) < 1, fällt aber gegen 1: Teilerfolg; dann wird die Rate bestimmt. Spätestens 16:20 abbrechen.

## [K] 12.10. 13:37 — Versuch 39: Dyson-Kontraktion — Vorhersage in der Richtung richtig, im Detail falsch
- Vorhergesagt: ρ(L) > 1 und wachsend. Richtig für drei der vier Varianten. Die beste Variante (nur Kompositionskanal, Kopfpaar als Datum, Rand hinten) liegt aber
  bis n ≈ 58 UNTER 1 und wächst nur langsam (e^{0.004n}). Das war knapper als erwartet.
- Meine Begründung (m-Familie, Richtungen mit Eigenwert ≈ 1) erklärt nur die Variante mit vollem E. Für die beste Variante ist sie geprüft und verworfen (k59c).
  Der wachsende Anteil ist glatt (niedrige Grade) und kein Ecken-Effekt. Sein Grund bleibt offen.
- Abbruchkriterium: gescheitert, nach ≈ 15 min Rechnung. Lehre: Der Spektralradius als normunabhängiger erster Test hat lange Normsuchen erspart.

## [K] 12.10. 13:53 — K3-2 (integrabler Eckkern), Vorhersage vor k60
- Grössenordnung: Auf der Eckskala ist (z − t)K_n ~ |d_n|·n^{0.9} = O(n^{−0.1}). Die Generatorterme q_{n+1}λ_n, tν und Ξ müssen sich dazu aufaddieren.
  Am Steg waren alle drei O(1) und vergleichbar (B = 1.04, mν = 0.47, mΞ/t = −0.57).
- **Vorhersage:** Ξ_n ist auf der Eckskala NICHT klein gegen q_{n+1}λ_n (Verhältnis 0.3 … 1). Die strenge Form von K3-2 scheitert dann nach dem Abbruchkriterium.
- **Schwächere Form, die ich für möglich halte:** Ξ_n und ν_n sind nahe L glatt und n-unabhängig, also ohne Struktur auf der Eckskala. Die ganze Eckskala steckt dann in
  q_{n+1}(z)λ_n(t), und auf der Diagonale ist q_{n+1}(t)λ_n(t) eine feste glatte Funktion (wie B am Steg).
- Abbruch wie in IDEEN festgelegt (Verhältnis > 0.1 bei n = 600 und nicht fallend). Gilt dann die schwächere Form (n-Unabhängigkeit auf 1e-2), zählt das als Teilerfolg.

## [K] 12.10. 14:01 — Versuch 40: K3-2 integrabler Eckkern — wie vorhergesagt gescheitert, auch in der schwächeren Form
- Die Vorhersage (Ξ nicht klein, Verhältnis 0.3 … 1) traf zu; gemessen 0.3 … 3.7, langsam steigend. Die schwächere Hoffnung (Ξ glatt und n-unabhängig nahe L) traf nicht zu.
  Ξ trägt das Eckprofil, gespiegelt zum Generator.
- Der Grund ist einfach und hätte vorher auffallen können: Skaleninvarianz der Ecke. f − s ~ d_n wird von K ~ 1/d_n aufgefressen. Lehre: Ein „verschwindender Faktor“ an
  einem Fixpunkt ist nur dann ein kleiner Parameter, wenn der Rest dort nicht mit derselben Potenz singulär ist. Erst die Homogenitätsgrade vergleichen.
- Nebenbefund: ν_n bis r = 0.003 an L n-unabhängig.

## [K] 12.10. 14:06 — Versuch 41: K3-3 Gerüst — Teilerfolg; die Marge ist eine Phasen-Aussage
- Vorher erwartet: Das Gerüst schliesst mit Vs Marge 0.0024. Tatsächlich gibt die Dreiecksungleichung genau die Marge 0, weil das Produkt Kopf × Adjungierte am Steg flach ist (§39).
  Die 0.0024 stecken in den Phasen der Paarungen. Das habe ich erst beim Aufschreiben von Schritt 4 gesehen.
- Positiv: Kein Kompositionstransport im Gerüst (Abbruchkriterium nicht ausgelöst). Die offene Ungleichung ist genau benannt und messbar.
  Neu als Entwurfsregel: Spiralschnitt auf der Eckskala (r₁ ~ C/n, kopfweise Konturen) und Linse auf der kritischen Höhe S₀.
- Muster (zum vierten Mal): Wo Beträge an der Grenze genau ausgeglichen sind, entscheidet Kohärenz. Das ist kein Zufall, sondern die Flachheit des Produkts.
