# Ideen für den Beweis (Phase 3)

Format je Idee: Nummer, Titel, Kern (2 Sätze), warum es gehen könnte, woran es scheitern könnte (Abgleich mit den bekannten
Sackgassen), Kosten, Bewertung K / V (Aussicht, Kosten, Risiko je 1–5).

## Ideen von Strang K

Stand 11.10. 20:30 (aus `date`). Bewertung K = (Aussicht, Kosten, Risiko), je 1–5, Kürzel der Sackgassen wie bei V.
Leitgedanke nach Landkarte A: offen ist eine Auslöschungsaussage in der Ecke; Ideen, die den Eingang nur über Normen sehen, sind tot (KA, GW, L2).

**K1. Kurvenkriterium: (H₀) ⟸ G_N ≤ CN^κ auf γ = Γout ∪ Linse(Γin) — Cauchy-Darstellung statt Eckformel.** [ausgeführt, k_kern/BERICHT §1]
Kern: Für jede Eingabe W, analytisch auf Int γ, ist J_N⁻¹T_NW = (1/2πi)∮_γ W(t)·J_N⁻¹T_Nk_t dt; mit W = J_∞k_η folgt |Π_Nk_η| ≤ C√N·max_γG_N/dist(η, K_γ),
K_γ = K' ∪ Linse, und Bernstein–Walsh gibt (H₀^{K_γ}). Es braucht kein D6, kein (E_sing), kein (a′)_r, keine Eckprojektoren — und es gilt für f_K.
Warum: G_N = sup_u|Y_u| ist Gs und Os Objekt; G10 liefert G_N aus |𝒟_n| ≤ C. Numerisch (k01): G_N auf Γout konstant, auf der Linse ≤ ~√N.
Scheitern: nur Umformulierung des Engpasses (beweist (H*_G) nicht); Abgleich: nicht LO (die Kurve ist global), nicht BM (G ist eine ausgewertete Grösse).
Kosten: erledigt bis auf Gegenlesen (V) und Kneser-Numerik. Bewertung K: (4 als Fundament, 1, 1).

**K2. Rechteck-Satz: vollständiger Konvergenzsatz für das Kleinste-Quadrate-Verfahren (N+1 Unbekannte, N+d Gleichungen) samt Zweig.**
Kern: Eckfreier Teil (X1 Schritt 4–5) + σ_min eines Rechteck-Abschnitts ist störungsstabil (X8 + X5 + (a″)_∞) ⇒ H1/H2 für das LS-Verfahren ohne
jede Eckaussage. Ergebnis: ein unbedingter Satz (bis (K2)/Gitterkonstanten) und die saubere Aussage „quadratisch − rechteckig = Rang-d-Eckteil“.
Warum: alle Bausteine liegen vor; die Familie ist im Rechteck-Fall eine Normstörung. Scheitern: der Zweig des LS-Verfahrens braucht einen eigenen
Startnachweis; LS-Normalgleichung mischt Tangente und Sekante (lösbar). Berührt den Engpass nicht. Kosten: 1 Tag (K schreibt, V prüft).
Bewertung K: (2 für die Vermutung / 4 als eigener Satz, 2, 2).

**K3. Kompositions-Vorkonditionierung auf der Gleichungsseite: T_N C_k mit k(0) = 0.** [Lemma in k_kern/BERICHT §0.3]
Kern: T_N C_k T_N = T_N C_k, also Π_N = (T_NC_kJ_∞|_{P_N})⁻¹T_NC_kJ_∞; für k = x·e^{εφ(x)} wirkt C_k auf Grad n wie e^{nε Re φ} — ein echtes
Phasenraumgewicht (Index × Winkel), das die verworfene beidseitige Agmon-Konjugation nicht liefern konnte. (H*_G) auf k(∂D) ⇔ ‖(T_NC_kJ_N)⁻¹‖ ≤ CN^κ.
Warum: verwandelt die Ecke in eine Stabilitätsaussage für quadratische Abschnitte eines neuen Operators δ∘(f∘k) + (m∘k)(δ∘k).
Scheitern: reine Akkretivität scheitert (die Eckmode transportiert Index 0.5N → 0.95N, ⟨C_kJu, u⟩ ≈ 0); nötig wäre ein Petrov–Galerkin-Testoperator,
den man nicht kennt. Die analoge MULTIPLIKATIVE Linksgewichtung ω = e^{Nφ} ist tot (verstärkt glatte Eingangsanteile ebenso — BRETT-Notiz 20:08 damit
überholt). Kosten: 1–2 Tage. Bewertung K: (2, 3, 4).

**K4. (H*_G) auf Γout und fern von L „weich“ beweisen, die Ecke isolieren.**
Kern: Auf Γout ist G_N numerisch exakt konstant (keine N-Abhängigkeit bis 1200): die Antwort auf k_t, t ∈ Γout, wird vom Kompositionskanal mit dem
nach innen kontrahierenden Zweig getragen. Ziel: Satz „G_N(t) ≤ C(t) für t ∈ Γout ∪ (Linse fern von L)“ mit einem Argument, das nur die Analytizität
der exakten Lösung J_∞⁻¹k_t ausserhalb einer Eckumgebung benutzt; übrig bliebe (H*_G) auf einem Kreisbogen |t − L| ≤ r₀ (= Os (KB) + Ws Eckzone).
Scheitern: J_∞⁻¹k_t hat bei L selbst eine schwache Singularität (Exponent s₀, Re 5.447), deren Abschneiden die Eckmoden anregt — der „weiche“ Teil
erbt die Ecke (Selbstähnlichkeit: (H*_G) bei t ist (H₀) für die Eingabe J_∞⁻¹k_t). Abgleich LO: genau das Risiko. Kosten: 1 Tag.
Bewertung K: (2, 2, 4).

**K5. Spielzeugmodell mit Randfixpunkt, Typ A: radiales Quadrat (M), Pseudodifferentialstruktur.**
Kern: J̃_N = P_N[Op(s(φ)^k) + M_m]P_N ohne Winkeltransport, Einträge zugeordnete Legendre-Funktionen mit gleichmässiger Asymptotik (Olver);
Ecke numerisch identisch mit Kneser (c = 0.0119, gleiches Profil). Ziel: (H*_G) im Modell über eine WKB-Konstruktion der adjungierten Gleichung.
Warum: bestes Modell der richtigen Ecke; reell. Scheitern: M §5 — nicht exakt lösbar; die Rate ist global (θ-abhängig) ⇒ man braucht trotzdem das
freie Randwertproblem. Kosten: Wochen. Bewertung K: (3, 5, 4).

**K6. Spielzeugmodell, Typ B: loxodromisches Möbius-Modell — Kreisbogengeometrie und Krawtchouk.** (= V3, Ergänzung)
Kern: Bei Möbius-f sind Γin/Γout-Grenzen, f(D̄), der A-Bogen und die Bucht Kreisbögen; das freie Randwertproblem (W §15/§20) und die Schweissung
(W §27.1) könnten explizit lösbar sein. J_N = R(S) + (dI + cS)^{−N}π_N(g) (alles Funktionen des einen Shifts S plus GL(2)-Darstellung).
Erst prüfen: universelle Ecke (κ₁, Ψ, κ_E) und (H*_G) auf einer Kurve wie in K1. Scheitern: nur ein abstossender Fixpunkt; nicht reell.
Kosten: Test 1/2 Tag (V), Analyse Tage. Bewertung K: (3, 3, 3).

**K7. c₂ als Fisher–Hartwig-Exponent der Ecke.** (= V2, Ergänzung)
Kern: log D_n = (n²/2)log λ + n log H − c₂ log n + …: der log n-Koeffizient −c₂ = 0.4273 (zwei Ecken, je 0.2136) ist der Kandidat für einen
lokalen FH-Exponenten. Gs Daten: dc₂/ds ≈ 0.4989 für μ → (1+s)μ, zweite Ableitung ≈ +2, Jet-Familie (gleiche L, μ) ändert c₂ nicht.
Kandidaten aus dem Malyuzhinets-Faktor Γ(½ ± iα(s−2)), α = log μ/2π, und β = log μ/2πi (κ₁ = 1/(1−β)) systematisch gegen Familien mit komplexen
μ-Richtungen und verschiedenen θ (radial) prüfen. Warum: eine exakte Formel bestätigt die lokale Natur der Ecke und gibt der Parametrix die Normierung.
Scheitern: hängt c₂ von globalen Daten ab, ist die Idee tot (Jet-Test spricht dagegen). Kosten: billig (V). Bewertung K: (2 für (H₀), 1, 2).

**K8. Abel-Koordinaten: die Ecke als flacher Zylinder ℂ/(γ_L − ½)ℤ, γ_L = 2πi/L.**
Kern: In Kneser-Koordinaten w ist J_∞ ≅ I + S_{1/2} (konstante Koeffizienten), und Ws Kegelvariable ist ζ^{κ₁} = exp(2πi(w − w₀)/(γ_L − ½)):
die Deckverschiebung 2πi − log μ ist die Gitterperiode γ_L − ½ („einmal um L, einmal f zurück“). Das Abschneiden ist eine Randschicht in Höhe
~log N. Lokale Parametrix als Wiener–Hopf-Problem für I + S_{1/2} auf einem Halbzylinder — mit eingebauter Kegelphase (die bei O fehlte).
Scheitern: WH (O §32.4) und vor allem W §27.5: die Amplitude ist global, ein lokales Problem liefert nur die Gestalt. Kosten: 2 Tage.
Bewertung K: (2, 3, 4).

**K9. Intuitives slog (Robbins) als lineares Zwillingsproblem.**
Kern: T_N[A∘exp − A − 1] = 0 ist derselbe Abschnittstyp (C_exp − I, Fixpunkt L auf |x| = ρ, Multiplikator L, κ₁ = 1.265 − 0.081i). Test der
Universalität (Rate, Profil, Satz K1 wörtlich übertragbar) und ein bedingter Satz „intuitives slog → Kneser-slog unter (H*_G)“ als Nebenprodukt.
Warum: lineares Problem ohne Familienfrage; offene Frage der Tetrationsliteratur. Scheitern: lenkt vom Engpass ab; gleiche Ecke, gleich schwer.
Kosten: Numerik 1/2 Tag. Bewertung K: (1 für die Vermutung, 2, 2).

**K10. Generator-/Levinson-Rekursion in n.** (= V1, Bewertung)
Kern und Prüfung wie V1. Meine Analyse: für die halbunendliche Matrix gilt exakt 𝐉S − L(F)𝐉 = L(w), w = M(x − F) (das ist G40 ohne Rand), also
𝒰S − (𝓛⁻¹L(F)𝓛)𝒰 = 𝓛⁻¹L(w) für 𝐉 = 𝓛𝒰; der Term 𝓛⁻¹L(F)𝓛 ist voll unterdreieckig ⇒ die Zeilen von 𝒰 erfüllen eine Rekursion mit VOLLER Historie.
Eine kurze Rekursion erwarte ich nicht (auch nicht für kubisches F mit bandbeschränktem L(F), L(w)). Billiger numerischer Test lohnt trotzdem.
Bewertung K: (1, 2, 4).

**K11. Transferoperator (Ruelle) auf der adjungierten Seite.**
Kern: J_∞ᵀ ≅ M_m + 𝓛 mit dem Perron–Frobenius-Operator 𝓛V(t) = Σ_{F(z)=t}V(z)/F'(z) (W §5); das Abschneiden ist dort „modulo Polynome ⊕ O(t^{−N−2})“.
Für das radiale Modell (konjugiert zu z² + c, c ausserhalb der Mandelbrotmenge) ist 𝓛 auf holomorphen Funktionen nahe der Julia-Menge nuklear mit
Spektrallücke; Lasota–Yorke-artige a-priori-Ungleichungen könnten |Y_u| auf der Linse kontrollieren. Scheitern: die Ecke L liegt nicht auf der
Julia-Menge-Struktur, die Ruelle braucht; Kneser ist nicht hyperbolisch-polynomial. Kosten: Tage. Bewertung K: (2, 4, 4).

**K12. Zwei-Variablen-Kern K_N(z, w) = e(z)ᵀJ_N⁻¹e(w) und Pluripotentialtheorie.**
Kern: G_N(t) = ‖K_N(·, ρ/t)‖, und K_N erfüllt T_N^{(z)}[K_N(F(z), w) + M(z)K_N(z, w)] = T_N[1/(1 − wz)] — ein Polynom vom Bigrad (N, N) mit
einer Funktionalgleichung in z. Siciak-/Bernstein–Walsh-Schranken in ℂ² auf einer Produktmenge, deren Gestalt die Funktionalgleichung erzwingt.
Scheitern: unklar, ob die Gleichung in ℂ² mehr Struktur hat als in ℂ (vermutlich NW-Niveau). Kosten: Tage. Bewertung K: (1, 4, 5).

**K13. Rechnergestützte starke Asymptotik (a-posteriori-Parametrix) ab N₁.** (≈ V6/V12)
Kern: Sobald eine äussere Parametrix (W-Modell) und eine innere (O) mit Fehlerschranken stehen, das Residuum in der Norm von G43 für alle N ≥ N₁
explizit abschätzen; endlicher Teil N < N₁ per Intervall-Bordering (V). Mit K1 genügt dafür die Zielgrösse G_N auf γ (keine Eckprojektoren).
Scheitern: die Parametrix existiert noch nicht. Kosten: Wochen. Bewertung K: (3 langfristig, 5, 4).

### Bewertung V zu Ks Ideen (11.10. 20:52; (Aussicht, Kosten, Risiko))
- **K1 (5 als Fundament, 1, 1)** — gegengelesen und **verifiziert** (v_pruefstand/BERICHT §1.3): Cauchy-Darstellung + Rationalität + Bernstein–Walsh,
  gilt für f_K. Ersetzt W1''/D6/(E_sing)/Q_N. Korrektur zur Stärke: mit Familie über X4a braucht es κ < 1.79 (nicht 2.79); gemessen 0.45.
- **K2 (3 als eigener Satz / 1 für die Vermutung, 2, 2)** — die Bausteine sind B (Ms (a″)_∞, X4a/X5/X8′) und hängen an (K2); genau die prüfe ich
  ohnehin. Eigener Startnachweis für den LS-Zweig nötig (das LS-System hat keinen σ-Homotopie-Zweig wie H2).
- **K3 (1, 3, 4)** — Π_N-Invarianz unter T_NC_k ist richtig (C_k unteres Dreieck), aber das ist eine Umformulierung ohne Hebel, solange kein
  Petrov–Galerkin-Test bekannt ist (Ks eigene Einschätzung).
- **K4 (2, 2, 4)** — die exakte N-Unabhängigkeit von G_N auf Γout (k01) ist der auffälligste Befund des Tages; ein Satz „r_t → s_t auf Γout“ wäre
  ein echter Teilbeweis, erbt aber die Ecke (Ks Selbstähnlichkeit). Wert: auch für die Familie (Überkonvergenz auf f(Γout)).
- **K5 (3, 5, 4)**, **K6 (3, 2, 3)** — K6 zuerst testen (V, billig); K5 nur, wenn K6 nicht die richtige Ecke hat.
- **K7 (2, 1, 2)** — Ks Zusatz dc₂/ds ≈ 0.4989 ist ein scharfer Prüfstein; ich rechne Kandidaten gegen Gs Familien (billig).
- **K8 (2, 3, 4)**, **K11 (2, 4, 4)**, **K12 (1, 4, 5)** — nach den Sackgassen (WH, lokale Amplitude, NW-Niveau) wenig Aussicht.
- **K9 (1 für das Ziel, 2, 2)** — als Nebenprodukt für den Tetrations-Graphen reizvoll (intuitives slog = Kneser ist eine offene Frage),
  K1 überträgt sich wörtlich. Nicht jetzt.
- **K10 (1, 2, 4)** — Ks Analyse (volle Historie) überzeugt mich; V1 damit zurückgestellt.
- **K13 (3 langfristig, 5, 4)** — = V6/V12; der endliche Teil (Intervall-Bordering von G_N auf γ bis N₁) ist sofort machbar (V).

## Ideen von Strang V

Stand 11.10. 20:08 (aus `date`). Bewertung V = (Aussicht, Kosten, Risiko), je 1–5 (5 = hohe Aussicht / hohe Kosten / hohes Risiko).
Sackgassen-Kürzel: KA = Kantorovich/Kugeln, LO = rein lokal (W §22–24), L2 = ℓ²-Extremalität, GR = Gonchar–Rakhmanov in der Ebene,
NW = n-te Wurzel allein, GW = andere Gewichtsradien, WH = lokales Wiener–Hopf, KI = Koeffizienten-Induktion, BM = Betragsmajoranten.

**V1. Verschiebungsrang 2 ⇒ Generator-Schur-Rekursion (Kailath) für 𝒟_n.**
Kern: Gs Lemma G40 sagt rank(J_nS² − L_nJ_nS + B_nJ_n) ≤ 2. Für Matrizen mit kleinem Verschiebungsrang erben Schur-Komplemente den Rang;
die Inverse wird dann durch wenige Generatorvektoren bestimmt (Gohberg–Semencul-Typ), und die Gradrekursion läuft über eine kleine
Transfermatrix statt über die ganze Vorgeschichte (G22b ist eine Volterra-Rekursion mit voller Historie).
Warum: würde aus der Volterra-Rekursion eine Rekursion fester Länge machen; All-n-Schranken wären dann Produkte kleiner Matrizen (wie bei Szegő).
Scheitern: L_n = T_F + S ist kein Shift, sondern operatorwertig; der Generator wächst mit n (G40 Grenze). Wenn die Generator-Rekursion nur
„O(n) Zahlen pro Schritt“ statt „O(1)“ wird, ist es KI in neuer Verpackung. Abgleich: nicht KI, solange nur ausgewertete Generatorfunktionen
majorisiert werden.
Kosten: 1–2 Tage Analyse; Test billig (V): Generatoren numerisch ausrechnen, prüfen, ob sie eine kurze Rekursion in n erfüllen.
Bewertung V: (2, 3, 4).

**V2. Universelle Pivot-Korrektur c₂ als Prüfstein und Normierung der lokalen Parametrix.**
Kern: G misst log(h_n/h_{n−1}) = log λ + c₂/n² mit c₂ = −0.427262 in kubisch UND Kneser (6 Stellen), während λ global ist (Grad-4-Störung
bei festem L, μ ändert λ). Also ist c₂ vermutlich eine Funktion nur von μ = f'(L) (und arg L) — eine rein lokale Grösse der Ecke.
Warum: Eine Formel c₂(μ) aus dem Mellin/Malyuzhinets-Eckmodell (O §33.4) wäre der erste exakt nachprüfbare Zusammenhang zwischen lokaler
Parametrix und globalen Daten; sie legt die Normierung fest, die Os Grenzproblem (P1)–(P3) heute fehlt.
Scheitern: c₂ könnte auch von globalen Daten abhängen (dann ist die Gleichheit Zufall der gleichen Fixpunktdaten); dann ist die Idee tot.
Kosten (V): billig — Gs neun Modellfamilien (family_*.json) enthalten Pivots bis n = 400 mit μ-, τ-, Jet-Variationen. Bewertung V: (3, 1, 2).

**V3. Loxodromisches Möbius-Modell als exakt strukturiertes Spielzeug mit Randfixpunkt.**
Kern: f = M⁻¹(μM), M = (x−L)/(x−L̄), hat eine Ecke bei L (Ds d24: Eckmode, Rate 0.077) und zugleich die exakte Faktorisierung
J_N = L_N((cx+d)^{−N})[π_N(g) + L_N(m/r)] (Xs Satz 1.5: Toeplitz × GL(2)-Darstellung, Krawtchouk-Matrixelemente).
Warum: Krawtchouk-Asymptotik (Baik–Kriecherbauer–McLaughlin–Miller) ist eine fertige RH-Theorie mit Band/Sättigungsbereichen; ein (H₀)-Beweis
im Modell wäre die Schablone für die Ecke. Erst prüfen (V, billig): hat das Modell dieselbe universelle Ecke wie kubisch/Kneser
(Os Profil Ψ, κ_E = ½, c₂)? Wenn nein, ist es kein Modell der richtigen Ecke.
Scheitern: Möbius hat nur einen abstossenden Fixpunkt (L̄ anziehend), nicht reell; die Ecke könnte anders aussehen. M §5: „exakte Lösbarkeit
und Eckmode schliessen sich aus“ — für dieses Modell nicht gezeigt. Kosten: Test 1 Tag (V); Beweis offen. Bewertung V: (3, 3, 3).

**V4. Kleiner Parameter ε = ln|μ| — Homotopie von der Möbius-Drehung (|μ| = 1, keine Ecke, (H*) mit C = 1.01) zu |μ| = 1.17.**
Kern: Die Eckrate c ist ungefähr proportional zu ln|μ| (D: ≈ 0.074 ln|μ|); bei ε = 0 gilt (H*) global (E R3.5). Man misst in einer
Modellfamilie mit festem arg μ, wie C_E(ε), κ_E(ε) von ε abhängen; ist κ_E stetig und klein, ist ein Störungssatz in ε mit expliziter Reichweite denkbar.
Warum: Die Ecke entsteht bei ε = 0 „aus dem Nichts“ — ein singulärer Störungsübergang, dessen Lage (ab welchem N die Ecke sichtbar wird, ~1/ε)
die Struktur verrät. Scheitern: genau dieser Übergang ist nicht gleichmässig in N (Grenzschicht der Breite 1/(εN)); dann liefert die Idee nur Einsicht.
Kosten: Messung 1 Tag (V), Beweis offen. Bewertung V: (2, 2, 4).

**V5. (H₀) gleichmässig über eine Modellumgebung ⇒ Familie (H) als Korollar.**
Kern: Statt (H₀) nur an f_K zu beweisen, formuliert man es für alle Keime f in einer C^{1,α}-Umgebung von f_K mit Fixpunkt L_f nahe dem
Kreis |x| = ρ (|L_f| − ρ = O(N^{−4})) und f'(L_f) nahe √L. Die Sekantenoperatoren 𝒥_P der Fast-Lösungen liegen in dieser Umgebung
(M_P statt f'∘f), also folgt (H-ii) für P aus derselben Aussage.
Warum: Xs Bogenmodell zeigt, dass die Ecke nur geglättete Fixpunktdaten sieht; ein Beweis von (H₀), der nur (L, μ, Univalenz, Randkurve)
benutzt, ist automatisch gleichmässig. Scheitern: der Fixpunkt neben dem Kreis ändert die Geometrie (Fixpunkt innen/aussen) — bei Abstand ≫ 1/N
könnte die Ecke qualitativ anders werden; dann braucht man die Abstandsskala ≪ 1/N explizit. Test (V): kubische Modelle mit |L| = ρ(1 ± δ),
δ = 1e-3…1e-8, Eckkonstanten messen. Kosten: Test 1 Tag. Bewertung V: (3, 2, 3).

**V6. Rechnergestützter All-N-Schluss in zwei Teilen: Intervall-Bordering bis N₁ + asymptotischer Satz ab N₁ mit expliziter Fehlerkonstante.**
Kern: V zertifiziert die Kreis-Majorante bzw. G_N auf γ* lückenlos für jedes N ≤ N₁ (Gs O(n²)-Bordering in Intervallarithmetik, Kneser bis 639 mit
propagierten Koeffizientenradien, kubisch bis 2000+ auf aicgolling); K liefert (G24): ‖χ_n^{−1}√n F_n − P_∞‖ ≤ A n^{−δ} ab n₁ mit explizitem A, δ.
Warum: trennt die Aufgabe sauber; der endliche Teil ist sicher machbar und senkt die Anforderung an den asymptotischen Teil (n₁ darf gross sein).
Scheitern: ohne asymptotischen Satz ist der endliche Teil wertlos (keine Aussage für alle N). Abgleich: kein BM, solange ausgewertete Funktionen
zertifiziert werden. Kosten: endlich ≈ Tage Rechenzeit Server; asymptotisch = der Engpass. Bewertung V: (3, 3, 4).

**V7. Schwächerer, jetzt erreichbarer Satz: „endlicher Satz bis N = 400 + bedingter Hauptsatz mit genau einer Hypothese“.**
Kern: Mit Zweig-Zertifikat bis 400, (K2) und Kopplung p31 gilt rigoros sup_{D̄}|p_N − f_K| ≤ ε_N für alle N ≤ 400 (explizite Tabelle),
und p_N → f_K für alle N, *falls* (H) mit Konstanten κ_E ≤ 2, C_E ≤ 1 gilt. Das ist ehrlich, prüfbar und als Zwischenergebnis für Janis vorzeigbar.
Warum: alle Teile sind endlich oder bereits verifiziert; der Rest ist sauber benannt. Scheitern: kein Beweis der Vermutung. Kosten: (K2) + Kette
+ Prüfungen (V, ~1–2 Tage). Bewertung V: (5 als Teilergebnis, 2, 1).

**V8. Grosses d statt Kante: den X3/X5-Weg (d = 4N bzw. 5N) als einzige Eckformulierung fixieren.**
Kern: Die offene „Kante (a″)_{αN} für kleines α“ wird nur gebraucht, wenn man den Eckteil klein halten will (d ~ r(N)); mit d = 4N folgt (a″)_{4N}
aus (a″)_∞ + X5, und die Eckbedingung wird die einteilige Randbedingung Θ^c ≤ 1 − ε für den Rang-4N-Eckteil (X3; Kneser Θ^c < 0 bei d = 2N, N ≤ 300).
Warum: streicht einen offenen Punkt ganz und reduziert die Eckaussagen auf EINE. Scheitern: **X hat das schon gemessen (X §3.5c):** Θ^{c,(2N)}
steigt linear (+0.028 je 100 Grade) und strebt nach Ws Profil gegen ≈ 0.82 (mit Aufschlag 0.99) — Reserve 0.18 bzw. 0.01, also grenzwertig, weil
der grosse Eckteil auf f_K(Γout) nicht abklingt. Rettung nur, wenn man für den grossen Eckteil einen besseren Rahmen als K' findet (Ecke wie W̄,
ohne (E_A)-Problem). Kosten: Buchhaltung + Rahmenfrage. Bewertung V (nach dem Befund): (2, 2, 3).

**V9. (Fam-Eck) scharf testen und als Ableitung der Eckrate nach den Fixpunktdaten beweisen.**
Kern: Xs offener Test: zwei Störungen mit gleichen geglätteten Fixpunktdaten (L_ε, μ_ε) und verschiedener Feinstruktur müssen dieselbe Eckantwort geben.
Gilt das, ist (Fam-Eck) = Glattheit der Ecken-Asymptotik in (L, μ), also ein Korollar eines (H₀)-Beweises, der glatt von den Daten abhängt.
Scheitern: Feinstruktur nahe L (auf Skala 1/N) zählt doch — dann ist die Familie ein eigenes Problem. Kosten: Test 1/2 Tag (V). Bewertung V: (3, 1, 2).

**V10. Falsifikationstest bei grossen N.**
Kern: Schur-Bordering kubisch bis n = 3000 auf aicgolling (Kreis-Majorante, γ*-Proben, Pivots), Kneser bis 639 mit zertifizierten Koeffizienten;
Ziel: Drift von κ_E, C_E und von √n|𝒟_n| sehen, bevor jemand einen Satz mit falscher Konstante beweist. E §5 warnt: Verfahren können über 100 Grade
konvergent aussehen und dann kippen. Scheitern: liefert nur Evidenz. Kosten: Server-Stunden, wenig Personenzeit. Bewertung V: (2, 2, 1).

**V11. Lyapunov–Schmidt auf eine explizite Eckbasis statt SVD.**
Kern: Die Vorwärts-Eckmoden sehen aus wie normierte abgeschnittene Singulärfunktionen bei L, Koeffizienten ~ k^s(L/ρ)^{−k} (Breite 1/N, Höhe √N;
Os u(L) → 0.42√N passt zu s ≈ −0.48, nahe Ws Kandidat Re α = −0.476). Man spaltet P_N = E_N ⊕ E_N^⊥ mit einer EXPLIZITEN Eckbasis E_N
(r(N)/2 komplexe Moden), löst auf E_N^⊥ mit (a″) (bewiesen bis (K2)) und reduziert auf eine kleine Gleichung für die Eckamplituden.
Warum: macht die Ecke zu einem endlichdimensionalen (r(N) ~ N/250) Problem mit expliziten Matrixelementen. Scheitern: r(N) wächst; die
Kopplung E_N ↔ E_N^⊥ ist genau die globale Amplitude (O: „Eckwert nicht lokal bestimmt“). Abgleich: LO, wenn die Eckbasis nur lokal gebaut wird.
Kosten: Test (Überlapp SVD-Moden mit Ansatz) 1/2 Tag (V). Bewertung V: (2, 2, 4).

**V12. Dual-Zertifikat im Small-Norm-Sinn (G43) zuerst numerisch prüfen.**
Kern: (H₀) ⇔ Schranke für die adjungierte Lösung y_N = J̃_N^{−T}u. Baut man ỹ_N explizit (Ws Kegel aussen + Os Parametrix innen), so muss das
Residuum J̃_N^T ỹ_N − u in einer Konturnorm klein sein, in der G43 greift. V misst für die beste verfügbare Näherung, ob das Residuum mit N fällt
(N^{−a}) — das entscheidet, ob sich die RH-Konstruktion lohnt, bevor K sie beweist.
Scheitern: ohne geschlossene Näherung (Os Grenzproblem nicht gelöst) gibt es nichts zu testen; dann zuerst V2/V3. Kosten: mittel. Bewertung V: (3, 3, 3).

**V13. N₀ senken: A_θ → A* ≈ 0.63 und zertifizierte Randkonstanten.**
Kern: Mit beweisbarem A_θ = 13.4 ist N₀ ≈ 318; mit dem gemessenen Optimum A* ≈ 0.63 wäre N₀ ≈ 93, und das Zertifikat bis 400 würde sogar
κ_E = 2.5, C_E = 1 (N₀ ≈ 222) abdecken. A* ist ein verallgemeinertes Eigenwertproblem (max ‖δ'∘f‖/(N‖J_∞δ‖)) — für festes N endlich zertifizierbar,
für alle N über eine L²-Bernstein-Ungleichung auf f(∂D) mit expliziter Konstante.
Warum: macht (A) robust gegen schlechtere Eck-Konstanten. Scheitern: die All-N-Bernstein-Konstante braucht konforme Abbildung mit Fehlerschranke.
Kosten: 1–3 Tage (V). Bewertung V: (3, 3, 2).

**V14. Kneser-Fassung der Eckreduktion nur mit Koeffizientenidentitäten.**
Kern: W1'' (Eckteil ⇐ (H*) auf γ*) ist nur für polynomiales f bewiesen (Eckformel D6). Statt einer analytischen D6 für f_K: den Eckanteil
T_i(ζ) = ⟨u_i, J_N^{−1}T_NJ_∞[k_ζ]_{>N}⟩ direkt über die Koeffizienten von f_K abschätzen (b_k bis 640 zertifiziert + (K2)), d. h. f_K =
Polynom vom Grad K + Schwanz mit |b_k| ≤ C₆k^{−6}; das Polynom erfüllt W1'' wörtlich, der Schwanz ist eine kleine Störung von J_∞ in einer Norm,
in der G10 nur endliche Identitäten braucht. Warum: schliesst die Lücke „Kneser ≠ kubisch“ (Landkarte B6.2) ohne (E_sing).
Scheitern: die Störung durch den Schwanz wirkt auf Werte ausserhalb D̄ (f(D̄) ⊄ D̄), dort ist k^{−6} nicht klein genug (Komposition verliert e^{Ng}).
Kosten: 1 Tag Analyse + Test. Bewertung V: (3, 2, 3).

**K14. (S_γ): gleichmässige Stabilität der konformen Abschnitte Ĵ_N = T_N^w C_k J_∞|_{P_N} als neues Hauptziel.** [11.10. 21:13, nach k04]
Kern: In der konformen Gleichungskoordinate w von D̄ ∪ Linse ist ‖Ĵ_N⁻¹‖ = 1.251 konstant (N = 100…800, kubisch), d. h. die Eckmoden verschwinden
ganz; mit Lemma K1″ folgt (H₀) mit Eck-Polynomfaktor 3/2. Ĵ_N ist der Block einer festen unendlichen Matrix
Ĵ_∞[i][j] = [w^i](Φ^j + M̂Ψ^j), Φ = (f∘k)/ρ, Ψ = k/ρ, M̂ = m∘k — dieselbe Zwei-Kanal-Bauart wie J̃ (dort Ψ = id).
Vermutetes Kriterium (zu prüfen, Kontrolle mit dünnerer Linse läuft): der jeweils dominante Kanal muss auf ∂𝔻 STRIKT expandieren
(|Φ| > 1 auf dem Γout-Bogen, |Ψ| > 1 auf dem Linsenbogen), Gleichheit nur in den Ecken; im Original ist |Ψ| = 1 auf ganz Γin (marginal) ⇒ Randzustände.
Beweiswege: (a) Abschnittsmethode für Ĵ_∞ (RSS: starker Limes 𝒥⁻¹ = J_∞⁻¹C_k⁻¹ beschränkt nach O7 + Kanten-Limes + Eckmodell);
(b) Parametrix: lokale Inverse durch Komposition mit den kontrahierenden Umkehrzweigen (Φ⁻¹ auf Γout, Ψ⁻¹ auf der Linse), deren
Abschneidefehler wegen strikter Expansion fern von L exponentiell klein sind; an L Eckmodell. Risiko: die Kanten-Limes-Operatoren haben
N-abhängige exponentielle Gewichte (keine feste RSS-Algebra); ob die Linse über dem Eckprofil Λ liegen muss (dann steckt die globale Eckrate
als Ungleichung c < 0.95·s(π) in (S_γ)), entscheidet der Kontrolllauf. Bewertung K: (4, 4, 3).

### Bewertung der V-Ideen durch K (11.10. 20:31)
| Idee | K: (Aussicht, Kosten, Risiko) | Begründung |
|---|---|---|
| V1 Generator-Schur | (1, 2, 4) | siehe K10: halbunendlich 𝐉S − L(F)𝐉 = L(w) exakt, LU-Konjugation macht die Historie voll; Test billig, Erwartung negativ |
| V2 c₂ universell | (2, 1, 2) | sehr gutes billiges Experiment (= K7); stützt Lokalität, beweist aber nichts an (H₀) |
| V3 Möbius lox. | (3, 3, 3) | einverstanden (= K6); Test zuerst: (H*_G) auf einer Kurve Γout ∪ Linse wie in K1, plus κ₁, Ψ |
| V4 kleines ln\|μ\| | (2, 2, 4) | O §11: nicht gleichmässig; als Messung nützlich. Zusatztest: reelles μ (arg μ = 0) hat Re κ₁ < 1 — W §19(b) sagt dann Scheitern voraus |
| V5 gleichmässig über Modellumgebung | (3, 2, 3) | wichtig für (H); mit K1 heisst die Familie: G^P_N ≤ CN^κ auf γ für alle P ∈ 𝒮_N ∩ B_N |
| V6 zweiteiliger Rechnerschluss | (3, 3, 4) | mit K1 genügt im endlichen Teil eine skalare Kurvenmajorante für G_N bzw. 𝒟_n auf γ (G11-Technik) |
| V7 endlich bis 400 + bedingt | (5 als Teilergebnis, 2, 1) | einverstanden; mit K1 ist die eine lineare Hypothese (H*_G) (+ Familie, (H-iii)) |
| V8 grosses d | (1, 2, 3) | durch K1 überholt (kein Rechteck-Abschnitt nötig) |
| V9 (Fam-Eck)-Test | (3, 1, 2) | billig, entscheidet, ob die Familie ein Korollar wird |
| V10 Falsifikation grosse N | (2, 2, 1) | sinnvoll; Messgrösse jetzt G_N auf γ (billiger als SVD, k01 bis 1200 in ~5 min) |
| V11 explizite Eckbasis | (2, 2, 4) | LO-Risiko; die Kopplung ist die globale Amplitude |
| V12 Dual-Zertifikat G43 | (3, 3, 3) | richtig gedacht; braucht erst eine Parametrix |
| V13 N₀ senken | (3, 3, 2) | Vs Bereich, wichtig für (A) |
| V14 Kneser-Eckreduktion | — | durch K1 erledigt (Cauchy-Darstellung braucht kein D6), falls Vs Prüfung von K1 hält |

## Auswahl

**Vorschlag K (11.10. 20:31), zur Abstimmung mit V:**
1. **A1 — Fundament K1** (K schreibt, V prüft den Beweis; Kneser-Numerik G_N auf γ mit `gemeinsam/kneser_halbexp_koeff_zert_640.json`).
   Danach ist die Zielaussage für den ganzen Eckteil (H*_G) (statt (Ẽ)/(H*)/W1'').
2. **A2 — (H*_G) im Modell mit Randfixpunkt** (K: Struktur, V: Messung). Erst loxodromisches Möbius-Modell (K6/V3), bei Misserfolg radiales Quadrat (K5).
3. **A3 — c₂ und Familie als Prüfsteine** (V: V2/K7, V9, V10 mit G_N auf γ).
4. **A4 — Rechteck-Satz** (K2) als gesicherter Zwischensatz, niedrigere Priorität.


---

## Runde 2 — nur Ecke (Auftrag Leitung 12.10.; K eröffnet 12.10. 00:34, V ergänzt und bewertet)

**Ziel der Ecke (Kopfform, BERICHT §11–12):** für die Köpfe 𝒟_n = ‖Q_n‖𝒜_n/h_n nahe L (Kegelskala d_n = n^{−a}, a = 1/κ₁)
(E₀) exponentiell: G_N = e^{o(N)} auf γ ∩ U(L) und Antwort auf glatte Daten mit Rate o(1); (E₁) polynomial: Amplitude B_n = √n|𝒟_n(L)| beschränkt,
Profilrest in W_{1.05} beschränkt (dann G37 ⇒ |𝒟_n| ≤ C/√n), für (S_γ) zusätzlich Fast-Orthogonalität der Köpfe.
**Ausgeschlossen (schon da gewesen):** Eckformel D6/(E_L)-Profil, Gonchar–Rakhmanov-Orthogonalität (W §25.2), tropische Bilanz/W12 (K3), lokale
Mellin-Parametrix allein (O, W §19), Fisher–Hartwig-Exponent c₂ (K7/V2), Abel-Zylinder (K8), Levinson-/Generator-Schur-Rekursion (K10/V1),
Möbius-/radiales Spielzeug (K5/K6/V3), Homotopie in ln|μ| (V4), a-posteriori-Parametrix/Dual-Zertifikat (K13/V6/V12, G43), explizite Eckbasis
(V11), Störung um ein gelöstes Modell (G-quadratisch), eigenständiges affines Eckmodell (§9.1(ii)), Spezialfunktions-/Gleichungssuche (§12.3).

Bewertung K: (Aussicht, Kosten, Risiko), je 1–5 (Aussicht 5 = gross; Kosten 5 = teuer; Risiko 5 = wahrscheinlich Sackgasse).

**R1. Aufblasen der Symbol-Singularität (algebraische Geometrie).** Gs gemischtes Symbol σ(z,s) = 1/(1−sF) + M/(1−sz) hat in (z_L, 1/z_L) einen
Punkt, an dem Zähler und beide Nenner verschwinden; lokal σ ≈ −[1/(σ' + με) + μ/(σ' + ε)] (homogen vom Grad −1, σ' = s/s_L − 1, ε = z/z_L − 1).
Aufblasen des Punktes ersetzt ihn durch eine projektive Gerade der Richtungen (σ':ε): Eckoperator = Mellin in der Radialvariablen ⊗ Operator auf ℙ¹.
Der Gradabschnitt (j, k ≤ n) wird auf der aufgeblasenen Fläche eine Kante in Richtung (1:1). *Neu:* Zwei-Variablen-Symbol statt Ein-Variablen-Mellin.
*Erster Test:* das ℙ¹-Problem explizit aufstellen und prüfen, ob seine Eigenrichtungen κ₁ und die Kegelskala a liefern. (3, 2, 3)

**R2. Renormierungsoperator auf den Kopfprofilen.** Der Schritt n → λn wirkt auf das normierte Profil S_n näherungsweise als Abbildung ℛ_λ
(Feigenbaum-Prinzip); S_∞ ist ihr Fixpunkt, die Amplitude B_n folgt aus dem Multiplikator von ℛ in Richtung der Konstanten. Beweisziel:
ℛ ist nahe S_∞ in W_{1.05} kontrahierend ⇒ Konvergenz von S_n MIT Rate, und B_n aus einem Eigenwert. *Neu:* Fixpunkt eines Skalierungsoperators
statt Asymptotik. *Erster Test:* aus Gs Profilen (n und 2n) die lineare Abbildung S_n ↦ S_{2n} (Taylor-Koeffizienten) schätzen, Spektrum ansehen. Hinweis: Gs Korrekturordnungen n^{−Re a}, n^{−1}, n^{−2Re a} (ETAPPE3 §4) sind genau die zu erwartenden
Eigenwerte λ^{−a}, λ^{−1}, λ^{−2a} von ℛ_λ — die Struktur ist also schon sichtbar. (3, 2, 3)

**R3. Mellin-Transformation im GRAD.** Dirichlet-Reihe H(s, ζ) := Σ_n n^{−s}√n𝒟_n(w_L + d_nζ); Selbstähnlichkeit in n wird zu Polen in s, die
Amplitude zum Residuum. Aus der exakten Volterra-Rekursion G22 bzw. der Ränderung eine Gleichung für H herleiten (Faltung in n ↦ Produkt in s).
*Neu:* Mellin in n, nicht im Ort. *Erster Test:* numerische Dirichlet-Reihe aus Gs Daten 50…1000, Polstruktur nahe s = 1 sichtbar? (2, 2, 4)

**R4. Kontinuierlicher Grad und Fluss (isomonodrome Struktur).** Statt hartem Abschnitt k ≤ n ein glattes Gewicht (1 − k/ν)_+^θ mit
kontinuierlichem ν; die Köpfe werden glatt in ν, ihre ν-Ableitung erfüllt eine Lax-/Toda-artige Gleichung (Verschiebungsrang 2, G40, spricht für
Integrabilität). Selbstähnliche Reduktion an der Ecke ⇒ eine gewöhnliche Differentialgleichung (Painlevé-Typ) für das Profil, deren Lösung über
Monodromiedaten beschrieben ist. *Neu:* Integrabilität/Isomonodromie. *Erster Test:* prüfen, ob die geglättete Matrix ebenfalls Verschiebungsrang ≤ 2
hat. (3, 4, 4)

**R5. Integrabler Operator nach Its–Izergin–Korepin–Slavnov.** Aus Verschiebungsrang 2 (G40) eine 2×2-Riemann–Hilbert-Aufgabe für die Resolvente
ablesen, deren Sprung direkt aus den Generatoren U_n, V_n gebaut ist; die lokale Parametrix an L ergibt sich aus dem Sprungverhalten, nicht aus einer
geratenen Mellin-Form. *Abgrenzung:* V1/K10 = schnelle Rekursion, G43 = globale Parametrix ohne diese Struktur. *Erster Test:* Generatoren U_n, V_n
numerisch für n = 200…800 auswerten — haben sie selbst ein skalenfreies Eckprofil? (3, 3, 3)

**R6. Γ-Konvergenz des Kopf-Variationsproblems.** Das Paar (Q_n, A_n) ist ein Sattelpunkt der Bilinearform ⟨⟨A, JQ⟩⟩ unter den Normierungen; auf der
Kegelskala zum Grenzfunktional übergehen (Γ-Konvergenz, angepasst an Sattelpunkte, z. B. über die Hermitisierung [[0, J],[J*, 0]]); Konvergenz der
Minimierer liefert S_n → S_∞ und über die Energie die Amplitude. *Neu:* Variationsmethode. *Erster Test:* Hermitisierte Form auf der Kegelskala
numerisch — ist der kleinste Singulärwert dort isoliert (Lücke)? (2, 3, 4)

**R7. Übergang Fixpunkt innen → Rand.** Eine Familie f_δ, deren Fixpunkt L_δ für δ > 0 im Inneren von D liegt (dort gelten Gs Beweistechniken mit
Kontraktionsreserve, vgl. quadratische Familie) und für δ → 0 auf den Rand wandert. Zweiparametrige Asymptotik in (n, δ): die Ecke ist die Grenzschicht
nδ^{…} ~ 1; gematchte Asymptotik verbindet den bewiesenen Innenbereich mit dem Randfall. *Neu:* Ecke als Grenzschicht eines Übergangs. *Erster Test:*
Köpfe für einige δ rechnen und sehen, ob das Profil von einer Variablen n·δ^{p} abhängt. (3, 3, 3)

**R8. Wendepunkt-Theorie im Grad.** Die Ränderung (n → n+1) als diskretes Übertragungssystem mit langsam veränderlichen Koeffizienten lesen;
die beiden WKB-Zweige (kompositionsdominiert, identitätsdominiert) verschmelzen an der Ecke — ein Wendepunkt in n mit komplexem Exponenten
a = 1/κ₁ statt Airy (2/3). Klassische Wendepunkt-Theorie (Olver) liefert lokale Lösungen und Anschlussformeln mit Fehlerschranken. *Neu:* WKB in n,
nicht im Phasenraum. *Erster Test:* lokale Zweige aus der Ränderung numerisch trennen (Wronski-artige Grösse in n). (2, 3, 4)

**R9. Eckangepasster Abschnitt auf der UNBEKANNTEN-Seite plus Vergleich.** Statt x^k-Abschnitt eine an die Kegelskala angepasste Basis (wie
gradierte Gitter bei Randintegralgleichungen mit Ecken) — dort ist Stabilität per Konstruktion beweisbar; danach p_N mit der angepassten Näherung
vergleichen: beide sind Polynome vom Grad ≤ N, die die Gleichung bis O(x^{N+1}) erfüllen, also reicht eine Abschätzung ihrer DIFFERENZ in einem
Raum, in dem die Ecke schon aufgelöst ist. *Neu:* Wechsel des Verfahrens + Vergleich statt direkter Stabilität. *Erster Test:* ob die Differenz
p_N − p_N^{angepasst} auf D̄ numerisch klein und eckfrei ist. (2, 3, 3)

**R10. Krein-Raum aus der reellen Struktur.** J̃ ist reell, L und L̄ sind vertauscht; die Bilinearform ⟨⟨A, JQ⟩⟩ zusammen mit der Konjugation
könnte J in einem indefiniten Skalarprodukt (Pontrjagin-Raum mit endlichem negativem Index) selbstadjungiert machen. Dann geben Pontrjagin-Sätze
Schranken für biorthogonale Köpfe über den negativen Index (endlich viele Eckmoden!). *Neu:* indefinite Metrik. *Erster Test:* numerisch nach einer
Gram-Matrix G mit J̃ᵀG = GJ̃ suchen und ihre Trägheit ansehen. (2, 2, 4)

**R11. Komplexe Interpolation in einem Modellparameter.** Holomorphe Familie f_τ (z. B. μ als Parameter bei festem L-Abstand), Köpfe 𝒟_n(L; τ)
holomorph in τ; mit Schranken auf zwei Randkurven des Parameterbereichs (|μ| = 1: Drehung, beschränkt; und ein zweiter lösbarer Rand) folgt eine
Schranke im Inneren aus Hadamard/Phragmén–Lindelöf. *Abgrenzung zu V4:* nicht störend, sondern komplex-analytisch in τ. *Erster Test:* gibt es einen
zweiten Rand mit bekannter Schranke (z. B. μ reell, wo das Symbol nicht spiralt)? (2, 2, 4)

**R12. Die künftigen Köpfe als Randwertaufgabe in (Index, Grad).** D_N = −Σ_{n>N} projizierte Köpfe (§11.3) ist universell; als Funktion zweier
skalierter Variablen (ξ = Index/N, λ = n/N) ist die Kopf-Familie eine Funktion auf einem Viertelraum, die eine lineare Gleichung in beiden
Richtungen erfüllt (Ränderung in λ, Matrixgleichung in ξ). Die Ecke wird dann ein ZWEIdimensionales Wiener–Hopf-Problem auf einem Quadranten
(bekannte Theorie: Malyuzhinets/Sommerfeld-Keil, Quadranten-Wiener–Hopf). *Neu:* Quadranten-Problem statt Halbachse. *Erster Test:* die
Kopfmatrix [Komponente ξ von c_n] über λ ∈ [1, 3] bei N = 300 aus Inversen ablesen und auf Skalenfreiheit in (ξ, λ) prüfen. (3, 3, 3)

**Vorläufige Wahl K (vor Vs Ideen):** R12 (macht aus meiner D_N-Formel ein lösbares Quadranten-Problem mit bekannter Theorie), R2 (billiger Test mit
Gs Daten, beweisnahe Struktur) und R5 (nutzt Gs Rang-2-Struktur direkt). Zurückgestellt: R4 (teuer), R3/R6/R8/R10/R11 (hohes Risiko).
→ V: bitte mindestens 5 eigene Ideen ergänzen, die sich von allen oben und von der Ausschlussliste unterscheiden, und R1–R12 bewerten.

### Erste Tests K zu R2, R5, R12 (12.10. 00:47; BERICHT §13)
- **R2 weiterführen:** linearisierter Skalierungsoperator aus Gs Profilen hat Eigenwerte 0.2890, 0.4955, 0.5451−0.0114i ≈ 2^{−2a}, 2^{−1}, 2^{−a}.
- **R5 verworfen (einfache Form):** J_n = Krylov(T_F) + T(M) hat kein Sylvester-Paar mit kleinem Rang; G40 ist zweiter Ordnung. Nebenbefund: φ = M(F − z)
  verschwindet in den Fixpunkten (für R1).
- **R12 zurückgestellt:** Köpfe skalenfrei in Norm und Lage, aber D_N entsteht durch mit N wachsende Auslöschung (oszillierender Kern) ⇒ erst Phasenreduktion.

### V ergänzt (12.10. 00:49) — sechs eigene Ideen zur Ecke (von R1–R12 und der Ausschlussliste verschieden)
Bewertung wie oben (Aussicht, Kosten, Risiko), dazu „Test“: wie schnell numerisch entscheidbar.

**VE1. Keil-Linse statt tangentialer Linse.** Bisher hat die Linse an L die Dicke O(τ^p), p > 1 (Dini-glatt, damit k' stetig). Die Geometrie erlaubt
aber eine Linse mit **Ecke** an L: f(ℬ) ⊂ D gilt, solange die Steigung unter U(τ)/τ ≈ 0.79 bleibt (V §3.2). Mit s(τ) ≈ c·τ ist die Dämpfung in der
Eckzone (τ ~ N^{−a}) e^{−cNτ} = e^{−c·N^{1−a}} → 0 — die Ecke würde mit N in den Steg wandern. Für die Kette von Satz V2 braucht man nur K1″(a)
(keine Glattheit von k'), die Eckwinkel-Regel für g_γ (Aussenwinkel von K_γ an L wird kleiner, Rate 5.447α bleibt > 3/2 mit grosser Reserve).
*Erster Test (gerechnet, 00:40):* kubisch, Steigung 0.5: ‖D_N‖ = 0.290 → 0.215 (N = 100 → 200, S = 0.0166) statt 0.477 → 0.512 (tangential);
bei S = 0.05: 0.184 → 0.195 statt 0.291 → 0.337. ‖Ĵ_N⁻¹‖ = 1.2505 bzw. 1.1089 (kleiner als tangential).
*Nachtrag 00:55 (N = 300/400):* Keil, S = 0.0166: ‖D_N‖ = 0.290, 0.215, 0.214, 0.216 (N = 100…400) — **sättigt**, fällt nicht weiter; S = 0.05: 0.184,
0.195, 0.200. Tangential dagegen langsam wachsend: 0.477, 0.492, 0.512, 0.542, 0.566, 0.584 (N = 100…500). ⇒ Starke Form (Ecke verschwindet)
**abgebrochen**; Nebenergebnis: mit Keil-Linse sättigt die Eck-Korrektur auf ein festes Paar (0.2164/0.2155 bei N = 400) — ein sauberes
Grenzobjekt für R12/R2, ohne das langsame Wachstum der tangentialen Linse. (4 → 2, 1, 3)

**VE2. Doppelkeil: Einbuchtung über Γout, so dass f(γ) ⊂ D̄.** Die Fehlskalierung (‖Ĵ_N‖ ~ e^{0.47N}, Ks §8.2) kommt allein daher, dass f(Γout)
ausserhalb D̄ liegt. K1″(a) verlangt aber nur W analytisch in Int γ — nicht D̄ ⊂ Int γ. Wählt man γ über Γout leicht INNERHALB von D, mit f(γ_out) ⊂ D̄
(möglich: {|f| ≤ ρ} enthält einen Keil an L, Steigung aus dem lokalen Bild L + μ(x − L)), so hat Ĵ beschränkte Einträge, und an L treffen sich zwei
Keile (Linse aussen, Einbuchtung innen). Dann greifen klassische Abschnittssätze vielleicht doch (beschränkte Folge), und die Eckgeometrie ist eine
andere. *Erster Test:* v12 mit negativem Profil über Γout; ‖Ĵ_N‖ und ‖Ĵ_N⁻¹‖ für N bis 800. (3, 2, 3), Test: Stunden.
*Ergebnis 00:58 (v36, kubisch, Einbuchtung bis Radius ≈ 0.55ρ, f(γ_out) ⊂ D̄ geprüft):* ‖Ĵ_N‖ = 2.5 / 2.8 / 12 (N = 50/100/200) — die Fehlskalierung
ist weg —, aber ‖Ĵ_N⁻¹‖ = 2.5e8 / 9e17 / 2e37 (Rate ≈ 0.42): die Gleichungen auf einer Kurve INNERHALB von D kontrollieren die hohen Koeffizienten
nicht mehr. Die Skalenfehlanpassung wandert vom Operator in die Inverse. **Abgebrochen** (spiegelbildlich zu Ks §8.2).

**VE3. Muckenhoupt-Gewicht an L auf der Gleichungsseite.** Statt Ĵ_N die Matrix T_N(ω)Ĵ_N mit ω(w) = ((w − w_L)(w − w_L̄))^β (analytisch in 𝔻, untere
Dreiecksmatrix, also mit den Abschnitten verträglich: T_N^w[ωg] = T_N(ω)T_N^w[g]). In K1″(a) wird dann sup_γ|ωW| statt sup_γ|W| gebraucht — das
Gewicht dämpft den 1/dist-Pol von W = J_∞k_η an L, während ‖(T_N(ω)Ĵ_N)⁻¹‖ ≤ ‖Ĵ_N⁻¹‖·‖T_N(ω)⁻¹‖ höchstens polynomial wächst. Ziel: β so, dass die
Eckmoden im gewichteten Raum genau auf der kritischen Mellin-Linie liegen (vgl. O7a, Re s = −½) — dann könnte die gewichtete Eck-Korrektur
gleichmässig klein sein und der Polynomfaktor κ_E sinken. *Erster Test:* ‖Ĵ_N⁻¹T_N(ω)⁻¹‖ und das gewichtete D_N für β ∈ {±¼, ±½}. (2, 1, 4), Test: Stunden.

**VE4. Faktorisierter Schur-Test mit optimiertem Diagonalgewicht (aus Ks Satz K5(5)).** (S_γ) folgt aus der Beschränktheit der beiden UNENDLICHEN
Dreiecksoperatoren U⁻¹W⁻¹ und WL̂⁻¹ für eine Diagonalmatrix W. Idee: W_n numerisch so wählen, dass beide Faktoren ausgeglichen sind (W_n ≈
(‖r̂_n‖/‖c_n‖)^{1/2}·Korrektur), dann für jeden Faktor einzeln einen Schur-Test (Zeilen-/Spaltensummen) beweisen — Dreiecksoperatoren mit
abklingendem Profil sind dafür viel zugänglicher als Ĵ_N⁻¹ selbst. *Abgrenzung:* nicht die Schur-REKURSION (K10/V1), sondern der Schur-TEST auf den
LU-Faktoren. *Erster Test:* LU von Ĵ_N (N = 200…600, arb), W optimieren, ‖U⁻¹W⁻¹‖, ‖WL̂⁻¹‖ und ihre Schur-Summen gegen N. (3, 2, 3), Test: 1 Tag.
*Erster Test gerechnet (v37, kubisch, Ks Linse, LU ohne Pivot in arb mit 5N + 800 Bit):* min_W ‖U⁻¹W⁻¹‖·‖WL̂⁻¹‖ ≈ 1.50 / 1.55 / 1.56 / 1.61 für
N = 60 / 100 / 150 / 200 (‖Ĵ⁻¹‖ = 1.2512), ohne Gewicht 50 / 790 / 2.5e4 / 7.9e5. Schur-Summen beim ausgeglichenen W (nicht auf Schur optimiert):
2.2 / 3.5 bei N = 200. ⇒ Mit einem Diagonalgewicht sind beide Dreiecksfaktoren O(1); (S_γ) zerfällt in zwei gewichtete Dreiecks-Schranken.
Server (Mittelpunkt-LU, 01:24): N = 300 / 400: 1.53 / 1.65; Keil-Linse N = 200 / 300: 1.56 / 1.61. Schur-Summen hängen vom (nicht Schur-optimierten)
Ausgleich ab (Produkt 7.5 / 10.9 / 9.8 bei N = 200 / 300 / 400). Nächster Schritt: Gestalt von w_n und Schur-optimiertes W bis N = 600. (3 → 4, 2, 3)

**VE5. Komplexe Interpolation in der Linsenhöhe.** F_S hängt analytisch von S ab; mit komplexem S wird Ĵ_N(S) eine holomorphe Matrixfamilie und
log‖Ĵ_N(S)⁻¹‖ subharmonisch. Für grosses Re S ist der Steg stark gedämpft, für die Ecke zählt nur die Geometrie an L. Hadamard-artige Schranken
in S (Dreigeraden-Satz) könnten eine Schranke bei der gewünschten Linse aus zwei leichteren Rändern (z. B. sehr dicke Linse, und Im S gross, wo die
Linse „rotiert“) gewinnen. *Abgrenzung zu R11:* Parameter ist die Koordinate (Linse), nicht das Modell (μ). (2, 3, 4), Test: 1–2 Tage.

**VE6. Diskrete Phasenraum-Lokalisierung in der w-Koordinate (IMS auf den Abschnitten).** Ms IMS-Zerlegung funktionierte auf Operatorebene; in w ist
die globale Schleife (Ks §8.2) durch die Linse gebrochen, also könnte jetzt ein lokales Prinzip auf den ABSCHNITTEN greifen: Ĵ_N = Σ χ_jĴ_Nχ_j +
Kommutatoren, χ_j = Toeplitz-Abschneidungen im Winkel (Ecke L, Ecke L̄, Steg) mal glatte Gewichte im Index; Eckstücke über die universelle
Eck-Korrektur (VE1/R12), Steg über die Linsendämpfung, Kommutatoren O(1/Breite). *Erster Test:* ‖[Ĵ_N⁻¹, χ]‖ und ‖Ĵ_N⁻¹ − Σχ_jĴ_N⁻¹χ_j‖ für
Abschneidungen der Breite 10°…40°, N = 100…400 (aus v26 schon nahe: fern von L fällt D_N wie e^{−SN}). (3, 2, 3), Test: Stunden.

### Bewertung V zu R1–R12 (Aussicht, Kosten, Risiko; Test = schnelle numerische Entscheidbarkeit)
| | V-Bewertung | Test | Kommentar |
|---|---|---|---|
| R1 Aufblasen | (3, 3, 3) | mittel | elegant; das ℙ¹-Problem ist aber im Gradabschnitt nicht kompakt (Kante (1:1)) — gleiche Schwierigkeit wie die Ecke selbst |
| R2 Renormierung | (3, 2, 3) | **schnell** | Gs Daten liegen vor; Spektrum von S_n ↦ S_{2n} ist in Stunden gemessen; Gs Korrekturordnungen passen |
| R3 Dirichlet in n | (2, 2, 4) | schnell, aber verrauscht | Polstruktur aus 50…1000 kaum auflösbar |
| R4 glatter Abschnitt | (2, 4, 4) | langsam | viel Theorie vor dem ersten Test |
| R5 IIKS | (3, 3, 3) | mittel | Rang 2 ist da; ob die Generatoren skalenfrei sind, ist schnell geprüft |
| R6 Γ-Konvergenz | (2, 3, 4) | mittel | Sattelpunkte schwer, Lücke unklar |
| R7 Fixpunkt innen → Rand | (3, 2, 3) | **schnell** | kubisches Modell mit L_δ = (1 − δ)L direkt in v12; Skalenvariable nδ^p in Stunden prüfbar |
| R8 WKB in n | (2, 3, 4) | mittel | komplexer Wendepunkt ohne Normalform |
| R9 angepasste Basis | (2, 3, 4) | mittel | Vergleich nur, wenn beide Verfahren dieselbe Gleichung (Taylor-Kollokation) lösen — sonst kein kleiner Unterschied |
| R10 Krein | (2, 2, 4) | schnell | J̃ hat keine erkennbare Symmetrie, die eine indefinite Metrik erzwingt |
| R11 Interpolation in μ | (2, 2, 4) | mittel | zweiter Rand mit Schranke fehlt |
| R12 Quadranten-WH | (3, 3, 3) | mittel | die beste Theorie-Anbindung; Test (Kopfmatrix über λ ∈ [1, 3]) aus vorhandenen Inversen |
**Vorschlag V für die Auswahl (mit K abzustimmen):** VE2 (Doppelkeil, Fehlskalierung weg) und VE4 (faktorisierter Schur-Test) von V; Ks R12 und R2,
beide mit der Keil-Linse (dort sättigt D_N, VE1-Nachtrag). R5 und R7 zurückstellen. (Ursprünglich VE1 vorne — starke Form 00:55 abgebrochen.)

### Runde 2, gezielt auf (H1) (K, 12.10. 01:22; nach KERNPROBLEM.md §3; Vs Vorschlag VE4 + R2 + Keil-Linse berücksichtigt)
**H-A. Ecke in Abel/Koenigs-Koordinaten: Kopf-Gleichgewicht als Erhaltungssatz.** In s = α(x) ist J_∞ ≅ Mult·(I + S_{½})·Mult (exakt, KERNPROBLEM §3.1);
die Ecke ist der Torus ℂ/(½ℤ + (2πi/L)ℤ), Gradskalierung ist Translation um −a·log λ/L. Für Lösungen von u(s+½) = −u(s), v(s+½) = −v(s) ist u·v eine
Funktion auf dem kompakten Torus — holomorph ⇒ KONSTANT (Liouville). Sind primaler und adjungierter Kopf nahe L genähert solche Eckmoden, so ist ihr
Produkt im Eckbereich (in der richtigen Paarung) nahezu konstant und τ_n ∝ 1/|Konstante| — (H1) ohne Amplitudenkonvergenz. *Abgrenzung:* K8 war nur die
Koordinate; neu sind die Gradskalierung als Translation und die Erhaltung des Produkts. *Erster Test:* Koenigs-Funktion des kubischen Modells an L
(über den inversen Zweig, kontrahierend mit 1/μ) ausrechnen, primalen Kopf Q_n und adjungierten Kopf 𝒜_n nahe L in s darstellen, Produkt (mit dem
Paarungsgewicht dx/x) auf einem Fundamentalbereich: ½-periodisch bis O(1/n)? nahezu konstant? Mittel ∝ 1/τ_n? (4, 2, 3)
**Auswahl (K, vorbehaltlich Vs Zustimmung; Vs Vorschlag aufgenommen):** H-A (K, Hauptstrang), VE4 (V, gewichtete Dreiecksfaktoren), R2 nur noch als
Teil von H-A (Renormierung = Translation) und mit Vs Keil-Linse als sauberem Grenzobjekt; R12 zurückgestellt (oszillierender Kern, BERICHT §13.3).


## Runde 3 — Kern mit den neuen Werkzeugen (Auftrag Leitung 12.10. ≈ 13:40; K eröffnet 12.10. 13:42; V ergänzt und bewertet)

**Neue Werkzeuge** (gab es in Runde 1/2 nicht):
- (W1) (CD1) (z − t)K_n = q_{n+1}(z)λ_n(t) − tν_n(z) + Ξ_n, Ξ_n = ∮K_n(z,s)(f(s) − s)K_n(f(s),t)ds/s. Bewiesen; V hat es bestätigt (kubisch und Kneser).
- (W2) explizite Eckprofile F_β (K10/K11, Gammaintegral).
- (W3) K19 mit (S_E): die Datenseite ist klassische Approximation auf E.
- (W4) Satz A (Budget γ = 0).
- (W5) Vs CF-Befunde: Repräsentant + Adjungierte ≈ −0.0024, Linsenfenster S ≈ 0.013–0.015.
- (W6) das feste Steg-Produkt ΔK_n(t,t)·m → κ(t) und B(t).
- (W7) Rang-1-Ränderung und die Rang-2-Verschiebung G40.

**Bekannte Fehlerart (KERNPROBLEM §2, §4):** Normverlust durch Kompositionstransport. Jede Abschätzung, die eine primale Grösse an f(Γout) (ausserhalb D) oder eine
adjungierte innerhalb D im Betrag braucht, verliert exponentiell (e^{0.47N}, ½·BW an f(Γout), §38).
**Pflicht je Idee:** welches Werkzeug, und warum diese Fehlerart nicht greift. Bewertung (Aussicht, Kosten, Risiko), je 1–5.

**K3-1. RH-Problem für (q_n, λ_n) aus (CD1) + (CD2) (Kandidat der Leitung). Bewertung K: global NEIN, lokal → K3-2.**
- (CD2) (Multiplikationskanal, (z − f(t))K = …) ist kein unabhängiger Satz. Die Differenz (CD1) − (CD2) ist (f(t) − t)K = −tν + Σ_{k<3}w_k(t)ν_k + J⁻¹T_n[(f − x)·J r_t].
  Wegen J r_t = T_n(Daten) + O(x^{n+1}) reduziert sich das auf eine Identität vom Rang 3 [Skizze, eine Rechnung].
  Zusammen liefern beide nur X = J⁻¹ (§41.4). Ein 3×3-RH-Problem aus (CD1) + (CD2) gibt es in dieser Form nicht.
- **R5 neu bewertet:** Global bestätigt sich der Einwand. Der Verschiebungsterm E = J S − S J hat vollen Rang (Kompositionskanal), und die Inverse ist
  „IIKS-Kern vom Rang 2 (Generatoren (q_{n+1}, ν_n) und (λ_n, −t)) plus Dyson-Term Ξ“. **Neu:** Der Dyson-Term trägt den Faktor f(s) − s, der in L, L̄ verschwindet.
  An der Ecke könnte der Kern also LOKAL integrabel sein (→ K3-2). Global bleibt ein RH-Problem mit Verschiebung (→ K3-4).

**K3-2. Integrabler Eckkern: (CD1) auf der Eckskala (lokales IIKS mit expliziten Generatoren).** (W1, W2, W7)
- Für z = L + d_nζ, t = L + d_nτ ist der Dyson-Integrand (f(s) − s)K(f(s),t) ≈ (μ − 1)(s − L)K(L + μ(s − L), t).
  Er ist auf der Eckskala um einen Faktor |s − L| ~ d_n kleiner als die Generatorterme, solange s im Eckbereich bleibt.
- Ist Ξ dort vernachlässigbar, dann ist der Eckkern ein integrabler Kern wie der Bessel-Kern an einer harten Kante:
  K_n ≈ [C·F_P(λ̃_Pζ)·F_A(τ) − t·ν(z)]/(z − t), mit den expliziten Profilen aus K10/K11.
- Das wäre die lokale Parametrix des RH-Ansatzes, und zugleich eine geschlossene Eckformel für K_n(z,t), nicht nur für P_n.
- *Warum die Fehlerart nicht greift:* Auf der Eckskala ist f eine Streckung um μ (Koenigs). f bildet den Eckbereich in einen Eckbereich ab und nicht nach f(Γout).
  Der Kompositionstransport bleibt lokal.
- *Risiko:* Das Integral ∮ ds läuft über eine ganze Kontur, nicht nur über den Eckbereich. Ob der Anteil fern von L auf der Eckskala klein ist, ist genau die Frage.
- *Erster Test (k60, ≈ 1–2 h):* Ξ_n(z,t) für z, t mit |ζ|, |τ| ∈ [0.5, 5] im Verhältnis zu (z − t)K_n und zu q_{n+1}(z)λ_n(t), bei n = 150/300/600.
  Fällt das Verhältnis wie d_n (bzw. n^{−1})? Dazu der Vergleich von q_{n+1}λ_n mit dem Produkt der Profile.
- Bewertung (4, 2, 3).

**K3-3. Bootstrap für (F3′) aus der Produkt-Reserve (mit V).** (W3, W4, W5, W1 nur über K13 kopfweise)
- Annahme entlang der σ-Homotopie (wie bei (EB_N)): P̃ ist auf der Linse δ₀-nah an f_K.
- Dann sind die Daten W = (f − q)∘P̃ auf der Linse klassisch (K19). Die Paarungen ⟨Λ_n, W⟩ werden auf der Linse abgeschätzt, wo die Adjungierte Reserve hat.
- Je Kopf n genügt der Keim mod x^{n+1}.
- Der neue CF-Repräsentant für N+1 ist dann wieder δ₀-nah, mit der Marge 0.0024 je Grad (Vs Produkt-Reserve).
- *Warum die Fehlerart nicht greift:* Alle Auswertungen liegen auf der Linse (Int γ). P̃ wird nur dort gebraucht, die Adjungierte nur dort (Reserve). Die primalen Köpfe gehen nur über
  ihre w-Koeffizienten im Stegband ein (CF-Hankel), nie über Werte an f(Γout).
- *Risiko:* Die Marge ist klein (0.0024), und die Beziehung „Bandkoeffizienten ↔ Paarungen“ muss ohne Betragsverlust gehen.
- *Erster Test (K + V, ≈ 1.5 h):* K schreibt das Gerüst und benennt die eine Ungleichung, die geschlossen werden muss. V misst sie mit dem gewichteten CF.
- Bewertung (3, 3, 3).

**K3-4. RH-Problem mit Verschiebung (Carleman/Haseman) und Verschweissung in Abel-Koordinaten.** (W1, K8)
- Die Kopfbedingungen [x^k](δ∘f + mδ) = 0 (k < n) sind Orthogonalität auf zwei Konturen: ∮_{f(γ₀)} δ(y)g(y)^{−k−1}g′(y)dy + ∮_{γ₀} δ(x)m(x)x^{−k−1}dx = 0.
  Dabei ist γ₀ ein kleiner Kreis um 0 und g = f⁻¹ nahe f(0). Das ist eine Zweipunkt-Hermite-Padé-Aufgabe (Punkte 0 und f(0)), gekoppelt durch f.
- Die Verschweissung über die Abel-Koordinate macht die Komposition zur Translation. Es entsteht ein Differenzen-RH-Problem (Borodin-Typ).
- *Warum die Fehlerart nicht greift:* Die Orthogonalität lebt auf kleinen Konturen um 0 und f(0); f wird nur lokal gebraucht.
- *Risiko:* Die Steilster-Abstieg-Verformung muss die Konturen bis |x| = ρ aufziehen. Dort trifft f(γ₀) auf f(Γout), also auf dieselbe Stelle wie vorher.
- *Langfristig.* Bewertung (3, 5, 4).

**K3-5. Szegő-artige Rekursion für das Kopfpaar aus Ränderung + (CD1) bei n und n+1.** (W1, W7)
- Die Differenz von (CD1) bei n und n+1 gibt (z − t)ΔK_{n+1} = q_{n+2}λ_{n+1} − q_{n+1}λ_n − tΔν + ΔΞ, mit ΔK_{n+1} = c·q_{n+1}λ_{n+1}.
- Wäre ΔΞ separiert, folgte eine Rekursion q_{n+2} = (cz − β_n)q_{n+1} + … (wie die Szegő-Rekursion bei OPUC).
- *Fehlerart:* ΔΞ ist der Kompositionskanal, also kommt die Fehlerart zurück. Bewertung (2, 3, 4).

**K3-6. B(t) exakt.** Grenzformel B = m(t)·(s₁(t) − Ξ_∞(t,t)/t) mit s₁ = J_∞⁻¹·1 [Skizze]. Geschlossene Form nicht in Sicht. Zurückgestellt (Leitung: höchstens 30 min). (1, 2, 4)

**Auswahl K (vorbehaltlich V):**
1. **K3-2 (Hauptstrang K), Zeitfenster 2 h, Abbruch:** Ist |Ξ_n|/|q_{n+1}λ_n| auf der Eckskala (|ζ|, |τ| ≤ 5) bei n = 600 nicht ≤ 0.1, und fällt es nicht mit n, ist der Eckkern nicht lokal integrabel → aufgeben.
   Fällt es wie n^{−c}: weiter zur expliziten Eckformel und ihrem Vergleich mit K10/K11.
2. **K3-3 (K Gerüst + V Messung), Zeitfenster 1.5 h K-Anteil, Abbruch:** Braucht die Schliessung der Ungleichung Werte der primalen Köpfe oder von P̃ ausserhalb der Linse
   im Betrag, ist es wieder die Fehlerart → aufgeben.

**Ergebnis K3-2 (K, 12.10. 14:01; BERICHT §42, k60):** gescheitert, in der strengen und in der schwachen Form. |Ξ|/|q_{n+1}λ_n| liegt auf der Eckskala bei 0.3 … 3.7 und wächst langsam
mit n. Auf der Diagonale gleicht Ξ das Eckprofil des Generators genau aus. Grund: f(s) − s ~ d_n und K ~ 1/d_n, das Produkt ist homogen vom Grad 0. Die Ecke ist echt zweikanalig.
Weiter mit K3-3.

**Ergebnis K3-3 (K, 12.10. 14:06; BERICHT §43):** Das Gerüst steht. (F3′) folgt aus EINER Ungleichung (BU) für die Bandkoeffizienten Σ_n [w^m](δ_n∘k)π_n/ℓ_nn, ohne Kompositionstransport.
Termweise schliesst sie nur auf der kritischen Linse S₀ ≈ 0.0125 mit Marge 0, und nur mit Spiralschnitt auf der Eckskala. Vs Marge 0.0024 ist Kohärenz der Paarungs-PHASEN;
Grössenschranken können sie nicht sehen. → V: R_P, R_A auf S₀ bis 3000; gewichtetes CF als direkter Test.
