# Prüfbericht Tetrations-Wissensgraph

425 Aussagen · 769 Kanten · 57 Experimente

**Effektiver Status:** bewiesen 155, numerisch 69, bewiesen_literatur 60, offen 34, bedingt 30, vermutung 29, definition 15, gescheitert 11, skizze 11, widerlegt 4, annahme 3, ueberholt 3, entlastet 1

**Zitate nicht am Original geprüft:** 104


## FEHLER (0)


## WARNUNG (0)


## HINWEIS (39)


**abschwaechbar** (22)

- `ar.kor.benachbart` — Voraussetzung unbenutzt: T = Kneser (A1, Korollar 1)
- `ar.lem.nullstellen` — Voraussetzung unbenutzt: T = Kneser
- `ar.satz.bpw` — Voraussetzung unbenutzt: T = Kneser
- `ar.satz.funktional_ks` — Voraussetzung unbenutzt: T = Kneser (A1, Satz 6)
- `ar.satz.orbitschranke` — Voraussetzung unbenutzt: T Kneser
- `br.satz.nicht_definierbar` — Voraussetzung unbenutzt: E monoton oder analytisch
- `br.satz.young_kreis` — Voraussetzung unbenutzt: H̃ differenzierbar
- `bw.prop.antisym_c1` — Voraussetzung unbenutzt: H = id + Φ orientierungserhaltender C¹-Kreisdiffeomorphismus
- `bw.prop.gruppoid` — Voraussetzung unbenutzt: reguläre Kneser-Zweige
- `bw.satz.imaginaere_richtung` — Voraussetzung unbenutzt: Kneser
- `fg.satz.allgemeine_loesung` — Voraussetzung unbenutzt: F(0)=G(0)=1
- `fp.satz.kreismittel` — Voraussetzung unbenutzt: Radius r = 1 (Shell–Thron-Rand)
- `fp.satz.rechteckwelle` — Voraussetzung unbenutzt: T analytisch oder regulär
- `fp.satz.spirale` — Voraussetzung unbenutzt: λ nicht reell
- `gr.satz.young` — Voraussetzung unbenutzt: T differenzierbar
- `it.satz.vertauschbar` — Voraussetzung unbenutzt: f^t reell auf der reellen Achse
- `kl.satz.basiswechsel_gleiches_L` — Voraussetzung unbenutzt: b, d reell
- `kl.satz.bruchsumme` — Voraussetzung unbenutzt: T regulär
- `kl.satz.grenzwertformel` — Voraussetzung unbenutzt: b reell
- `kn.satz.abel_bruchsumme` — Voraussetzung unbenutzt: A = Kneser
- `kn.satz.hexation_keine_halbebene` — Voraussetzung unbenutzt: H reell auf der Achse oder H(0) = 1
- `kn.satz.singulaer_minus2` — Voraussetzung unbenutzt: T differenzierbar in einer Umgebung von −1 (Paper I)

**fassung** (3)

- `kn.gescheitert.polynom_kantorovich` — wird von kn.beob.polynom_eckmode korrigiert; Status 'gescheitert' prüfen
- `kn.gescheitert.polynom_rekonstruktion` — wird von kn.satz.halbrang_quadratisch korrigiert; Status 'gescheitert' prüfen
- `kn.verm.tpid3` — wird von kn.verm.tpid3_bedingung korrigiert; Status 'vermutung' prüfen

**theorie-fortschritt** (5)

- `ar` — 29 von 31 Aussagen mit Theorie-Block
- `br` — 28 von 32 Aussagen mit Theorie-Block
- `gr` — 56 von 57 Aussagen mit Theorie-Block
- `it` — 19 von 20 Aussagen mit Theorie-Block
- `kn` — 116 von 126 Aussagen mit Theorie-Block

**uebersetzung-fortschritt** (6)

- `ar` — 31 von 31 Aussagen übersetzt (0 nicht aktuell, 0 gegengelesen), 29 von 31 Theorie-Blöcken
- `br` — 32 von 32 Aussagen übersetzt (0 nicht aktuell, 1 gegengelesen), 28 von 32 Theorie-Blöcken
- `gr` — 57 von 57 Aussagen übersetzt (0 nicht aktuell, 0 gegengelesen), 56 von 57 Theorie-Blöcken
- `it` — 20 von 20 Aussagen übersetzt (0 nicht aktuell, 0 gegengelesen), 19 von 20 Theorie-Blöcken
- `kn` — 126 von 126 Aussagen übersetzt (0 nicht aktuell, 1 gegengelesen), 116 von 126 Theorie-Blöcken
- `quelle` — 217 von 301 quelle-Gegenständen übersetzt

**voraussetzung-offen** (3)

- `it.satz.szekeres` — Benutzung unklar: f reell-analytisch (bzw. hinreichend glatt) nahe 0 mit 0 < f(x) < x
- `kn.satz.paulsen_komplex` — Benutzung unklar: Basis im Gültigkeitsbereich von Paulsen 2019 (nicht am Original geprüft)
- `kn.satz.schwellen_monoton` — Benutzung unklar: Diese Monotonie gilt für alle gebauten höheren Stufen

## Beweislücken

- `gr.verm.turmintegral_rate` — Skalierungsfunktion im Fenster (Normalform $w\mapsto-w+aw^3$ mit Parameter $x-e^{-e}$) nicht ausgerechnet; Konstante C offen.
- `kn.satz.hexation_logpunkt` — Dass θ bis zur Achse stetig ist (Carathéodory), wird vorausgesetzt und nicht aus der Streifen-Bedingung hergeleitet.
- `kn.satz.hexation_logpunkt` — Dass F_q nur an Log-Punkten der Form log₂ε + h(ε) oder verschachtelt gegen ∞ strebt, ist aus dem Bau von P geschlossen (P → −∞ nur bei w_b + jΠ_T), nicht vollständig bewiesen.
- `kn.satz.hexation_logpunkt` — Teil (iv): Ausser p sind weitere singuläre Werte von F_q ausserhalb der Achse nicht ausgeschlossen.
- `kn.konstr.zwei_sinh` — C^∞ nur sinngemäss übertragen, nicht ausgeschrieben.
- `kn.konstr.zwei_sinh` — Nirgends analytisch ist eine Vermutung (Forum).
- `kn.satz.grenzwert_trick` — Die genauen Bedingungen an δ (insbesondere abseits der reellen Achse) sind nur für den Basiswechsel ausgeschrieben (Paper I).
- `kn.satz.polynom_halbexp_hauptsatz` — Annahme (H) ist offen (kn.verm.polynom_halbexp_eckmoden).
- `kn.satz.polynom_halbexp_hauptsatz` — Phase 3: (H) in der Form bis an L ist für astronomisch grosse N vermutlich falsch (kn.beob.polynom_rand_plateau); tragfähige Fassung: kn.satz.polynom_hauptsatz_kr0.
- `kn.satz.polynom_halbexp_hauptsatz` — Die Einzelbeweise stehen in forschung/polynom_halbexp/ARCHITEKTUR.md §4 und den Strang-Berichten; die Leitung hat sie nicht Zeile für Zeile nachgeprüft.
- `kn.satz.polynom_halbexp_hauptsatz` — Eine Prüfung (f_K(e^z) ∉ (−∞,0] auf dem linken Bogen) ist nur numerisch.
- `kn.satz.kneser_e_zertifiziert` — Drei formale Lücken ohne Wirkung auf die Zahlen (Fehlerbudget in Prüfung B, eine Lipschitz-Schranke, Zweigtreue der DFT-Knoten), benannt in phase3/v_pruefstand/BERICHT §1.1.
- `kn.satz.polynom_eckprofil` — Lemma K12 (Übersetzung der Mellin-Randbedingung in Sektorbedingungen) nur skizziert.
- `kn.satz.polynom_eckprofil` — Dass die echten Köpfe gegen die Grenzlösung konvergieren, ist gemessen, nicht bewiesen; ebenso die Formel für C (Anpassung an die Aussenlösung).
- `kn.satz.polynom_approximation_ecke` — Dass E eine Quasischeibe ist (Jordan-Gebiet, Dreipunkt-Bedingung, Öffnung), ist nur numerisch auf einem Raster geprüft (kubisch und f_K).
- `kn.satz.polynom_hauptsatz_kr0` — Die drei Aussagen des Kerns (kn.verm.polynom_halbexp_kern).
- `kn.satz.polynom_hauptsatz_kr0` — Die Umstellung von H1/H2 auf K_{r₀} und K19 sind von Strang V gegengelesen, nicht von der Leitung Zeile für Zeile.
- `bw.satz.phasengesetz_gestoppt` — Die wörtliche Übertragung der Lemmas 5.2–5.4 (Paper I) auf $j\ge J_0$ ist behauptet, nicht ausgeschrieben.
- `bw.satz.phasengesetz_gestoppt` — (R3) — lokale Lipschitz-Stetigkeit von $A_b$ — folgt aus reeller Analytizität von $\mathrm{slog}_b^{\mathrm K}$; im Text nicht explizit gemacht.
- `bw.lem.randgesetz` — Kontinuumsnäherung statt rigoroser Durchgangszeit (Douady–Lavaurs-Theorie).
- `bw.lem.randgesetz` — Behauptete Korrektur $O(\ln(b-\eta))$ nicht hergeleitet.
- `bw.satz.ursache_kozykel` — Herleitung im Forschungslog, nicht als Beweis ausgeschrieben.
- `bw.satz.ursache_kozykel` — Zusammenhang mit der gemessenen Lage (Positionsgesetz) nur numerisch.
- `ar.lem.nullstellen` — Aus $\log((d+1)H)$ entsteht ein Term $d^2\log^2(d+1)$, der in der Schranke nicht explizit steht (für $H\ge d$ absorbierbar).
- `ar.lem.nullstellen` — Konstanten nicht ausgeführt; Fall $P(1,Y)\equiv0$ nur skizziert.
- `ar.satz.masser_typ` — Buchhaltung von κ nicht ausgeführt.
- `ar.satz.masser_typ` — Höhenkontrolle der Hilfspolynome ($\log H(P)\le c''d\log H$) übernommen, nicht nachgeprüft.
- `fp.verm.knospen_drehzahl` — Nahe der Wurzel der Knospe folgt die Drehzahl aus der Bifurkation ($f\approx L+\lambda(z-L)$ mit $\lambda\approx e^{2\pi ip/q}$); für die ganze Knospe fehlt ein Argument (Mandelbrot-Analogon: Douady–Hubbard über äussere Strahlen).
- `fp.satz.siegel_mittel` — Nach Herman (1985) ist $\Delta$ für diophantisches $\alpha$ unbeschränkt; Carathéodory greift dann nicht. Zu zeigen wäre, dass $\varphi$ Randwerte in $L^1(d\theta)$ hat (Hardy-Raum $H^1$), damit der Mittelwertsatz gilt.
- `fp.satz.siegel_mittel` — Ohne Jordan-Rand ist die Randdynamik nicht automatisch zu einer Drehung konjugiert; die Gleichverteilung des Turms bezüglich $\omega_L$ ist offen.
- `kl.satz.invariante_daempfung` — Dass die klassische periodische Phase (bw.satz.phasengesetz) als Grenzwert |λ| → 1 einer Familie gedämpfter Übergänge entsteht, ist eine Deutung, kein Grenzwertsatz.
- `kn.satz.theta_singularitaet` — Zweigwahl: bewiesen ist, dass die Zweige überall übereinstimmen, sobald sie an einem Punkt übereinstimmen; an z₀ = 0.05i ist der Abstand 4·10⁻⁴⁸ (nötig < 0.12), aber nur in Gleitkomma gerechnet, nicht mit Intervallen, und nur für Basis e.
- `br.satz.vergleichsfunktion` — Konvergenz für die Exponentialfakultät nur gemessen (doppelt exponentiell, ganze Werte auf 10⁻¹²), nicht bewiesen.
- `br.satz.vergleichsfunktion` — Für Tetration mit Produkt-Tetration als Vergleich (tid 1816) weder Konvergenz noch Glattheit geprüft.
- `gr.korr.grenzidentitaet` — Gleichmässige Abschätzung der unteren Faktoren, die den Grenzwert nicht ändern.

## Bedingte Aussagen (Wurzelbedingungen = unbewiesene Enden der Beweiskette)

- `kl.satz.folgen_vollmonoton` — kl.verm.vollmonoton
- `kl.satz.newton_reihe` — — · zusätzlich: Die Poincaré-Reihe Ψ am anziehenden Fixpunkt hat Konvergenzradius grösser als |s₀| (für b = √2 nachgerechnet: 1.316 > 0.632); allgemein aus der Holomorphie der regulären Tetration auf Re h > −2 zu zeigen, die bisher nur für √2 belegt ist (kl.satz.vier_sqrt2)
- `kn.konstr.pentation` — — · zusätzlich: λ_b = T_b′(L_b) > 1 und L_b einziger Fixpunkt von T_b in (−2, −1): gemessen (e 6.4607, 2 5.1718), bewiesen nur λ_b ≥ 1 am kleinsten Fixpunkt; gleichwertig zur Log-Konkavität von T_b auf (−1, 0), die für b > b₇ bricht (gr.gegen.tpid7_kneser)
- `kn.satz.hexation_struktur` — —
- `kn.konstr.hexation_regulaer` — —
- `kn.satz.hexation_keine_halbebene` — kn.beob.hexation_verzweigt
- `kn.satz.hexation_eindeutig_g0` — kn.beob.hexation_verzweigt
- `kn.satz.hexation_gebiete` — kn.beob.hexation_verzweigt
- `kn.satz.hexation_logpunkt` — kn.beob.hexation_verzweigt
- `kn.satz.halbexp_beschraenkt` — kn.beob.slog_scheibe
- `kn.satz.polynom_halbexp_reduktion` — kn.beob.slog_scheibe, kn.verm.polynom_halbexp_schranke
- `kn.satz.polynom_halbexp_hauptsatz` — kn.verm.polynom_halbexp_eckmoden
- `kn.satz.polynom_hauptsatz_kr0` — kn.satz.polynom_approximation_ecke, kn.verm.polynom_halbexp_kern
- `bw.satz.schwanzkonvergenz` — bw.ann.hyp22
- `bw.satz.phasengesetz` — bw.ann.hyp22
- `bw.lem.offset` — bw.satz.phasengesetz_gestoppt
- `bw.prop.kozykel` — bw.ann.hyp22
- `bw.satz.homoeo` — bw.ann.hyp22
- `bw.satz.antisymmetrie` — bw.ann.hyp22
- `bw.satz.c1_produkt` — bw.ann.hyp22 · zusätzlich: (S) differenzierbare strenge Monotonie und (R′) quantitative Kanalpositivität
- `bw.satz.antwortkern` — bw.ann.basis_diffbar, bw.ann.hyp22
- `bw.prop.hub` — bw.ann.hyp22
- `bw.satz.kern_eindeutig` — bw.ann.K1K2
- `bw.satz.kern_existenz` — bw.ann.K1K2
- `bw.satz.anker_rekonstruktion` — bw.ann.K1K2
- `bw.zert.8r1` — — · zusätzlich: korrekte Implementierung der Intervallarithmetik und die drei offenen Punkte des Zertifikats-README: Mikrosegment-Induktion, Hilbert-Röhren-Supremum, Γ-Steigungsbeweis
- `ar.satz.schanuel_orbit` — ar.verm.schanuel
- `gr.satz.ein_wendepunkt_e2` — — · zusätzlich: Annahme H: Motorwerte an den Kreisknoten auf 10⁻⁵ genau und |sexp| ≤ 1000 nahe der Achse; dass die Werte zu Kneser gehören, ist nur numerisch belegt
- `kn.satz.theta_singularitaet` — — · zusätzlich: Zweigwahl an einem Punkt (z₀ = 0.05i, Abstand 4·10⁻⁴⁸ bei nötigen < 0.12) nur in Gleitkomma geprüft, nicht mit Intervallen; nur Basis e
- `kn.satz.theta_koeffizienten` — — · zusätzlich: dieselbe Zweigprüfung wie bei kn.satz.theta_singularitaet

## Experimente

- `exp.fixpunkte` — aktuell, 42/42 Prüfungen. Fixpunkte über W_0/W_{−1} mit Residuum < 1e-50; für b > η alle geprüften |λ| > 1; Minimumsformel exakt bestätigt.
- `exp.turmgebiet` — aktuell, 2/2 Prüfungen. 109518 Basen mit |t|<0.97: alle konvergent (0 Gegenbeispiele). Knospen p/q, q ≤ 7: Periode q in 34/34 Fällen.
- `exp.sqrt2_abel` — aktuell, 14/14 Prüfungen. |p̂(k)| ≈ C k^γ exp(−2π²k/|ln λ_a|): freie Steigung trifft 2π²/|ln λ_a| auf -2.6e-04, -2.0e-04, -1.5e-04, -1.2e-04, 2.0e-05 (relativ); normierte Verhältnisse r_k über die Basen gleich bis 0.00043; r_k·k^(5/4) = 1 ± 0.0064 (k ≤ 10). Konkurrenz 2π²/ln λ_A um ≥ 20 % daneben. Basis √2: |p̂(1)| = 3.6260878e-25.
- `exp.carleman` — aktuell, 3/3 Prüfungen. Innere Punkte: Fehler 1.7e-5 (n=16) → 9.2e-9 (n=128), ≈ 8.04 Stellen; Abfall algebraisch ≈ n^-3.81, nicht geometrisch. Punkt 1.25 nahe |L_e|: 6.3e-9. Konvergenz gestützt, aber langsam; vorab erwartete 10 Stellen nicht erreicht.
- `exp.nichteindeutigkeit` — aktuell, 9/9 Prüfungen. Reguläre Tetration F(1/2) = 1.24362162767; zwei weitere reell-analytische, streng wachsende Lösungen mit T̃(1/2) = 1.15 bzw. 1.30 erfüllen die Gleichung auf 1e-40.
- `exp.singulaer` — aktuell, 4/4 Prüfungen. T'(−1) = 1.69973893; Restterm fällt wie ε (bei ε = 1e-12: -2.0e-12). Lineare Tetration: Konstante exakt 0.
- `exp.eta_baker` — aktuell, 4/4 Prüfungen. Konjugation exakt (Residuum 9.0e-60). Halbe Iterierte: c_2 = 1/4, c_3 = 1/48, …; |c_n|^(1/n) wächst von 0.196 (n=10) auf 1.630 (n=60) — divergent, wie Baker.
- `exp.kongruenz` — aktuell, 3/3 Prüfungen. V(a) für 359 Basen bestimmt (40 Stellen). Formel v₂(a²−1)−1 für a ≡ 5 (mod 10): 40/40 Treffer.
- `exp.randgesetz` — aktuell, 3/3 Prüfungen. Iteration: A = 3.23872 vs π√(2η/e) = 3.23893. Hub-Daten: A = 3.3822 mit ln-Term, 3.0155 ohne.
- `exp.young` — aktuell, 3/3 Prüfungen. Identität für 200 Zufallslifts (max. Abw. 6.0e-10), einen nicht-glatten Lift und 2 echte Phasen (max. 2.2e-12) bestätigt: Paper II Thm. B braucht nur Monotonie + Inverse.
- `exp.moden` — aktuell, 3/3 Prüfungen. Profilform über 2 Basen nahezu gleich (Spannweite ln ≤ 0.03 für k ≤ 6). m(2): phi_modes weicht 8.5e-37, mu_hub 1.8e-31 vom 180-Stellen-Wert ab.
- `exp.ableitung` — aktuell, 20/20 Prüfungen. Alle Ableitungsidentitäten an Kneser e auf ≤ 1e-20 bestätigt. T_e'(0) = 1.091767351258321, T_e''(0) = 0.54296642580339. Wendepunkte (Kneser): b=1.5: ['3.986806975'], b=2: ['-0.03024454275'], b=e: ['-0.5217647926'], b=3: ['-0.6038637144'], b=5: ['-0.8469203702'], b=10: ['-0.9967056974'], b=100: ['-1.181880344']. Regulär (1<b<η): b=1.2: keiner, b=1.3: keiner, b=1.41421: keiner, b=1.44: keiner.
- `exp.integrale` — aktuell, 11/11 Prüfungen. Kneser e: ∫₀¹T = 1.714521958890620488165975, ∫₁ᵉ slog = 1.0037598695684247472, ∫₋₁⁰T = 0.49921818282999133973; alle Identitäten auf ≤ 1e-28 (Produktform 1e-25). Regulär √2: ∫₀¹T = 1.2315849732804049508. PSLQ ohne Relation.
- `exp.basisableitung` — aktuell, 5/5 Prüfungen. Ableitungsrekursion und h'(a)-Formel auf ≤ 1e-35 bestätigt; (e − h)/√(2e³(η−a)/η) = 0.999998383423 bei η−a = 1e-12; nächste Ordnung ≈ -1.61658·√(η−a).
- `exp.slog` — aktuell, 5/5 Prüfungen. ⌈slog_e⌉ = log* an allen Punkten bis 10^1000; slog_e'(0) = 0.9159460564995334 = 1/T_e'(0); Stellenformel exakt; 2↑↑5 hat 19729, 3↑↑4 hat 3638334640025 Stellen.
- `exp.parameterraum` — aktuell, 8/8 Prüfungen. Periode 1 genau für |t| < 1 (0 Gegenbeispiele in beide Richtungen auf 520×520 Gitter); reelle Achse 2/1/Flucht; Basis i wie in der Notiz; Knospen-Drehzahl p/q in 17/17; parabolische Exponenten q=2: 0.501, q=3: 0.324.
- `exp.mittelwerte` — aktuell, 7/7 Prüfungen. Satz 1 (ganze Höhen) und Korollar 2 bestätigt. Für t = 1/2 weicht das Kreismittel von 1 ab und hängt vom Zweig von λ^t ab: r=0.5, Hauptzweig (−π, π]: (1.06186579 + 0.005130607658j); r=0.5, Zweig [0, 2π): (0.9919212687 - 0.08162539321j); r=0.8, Hauptzweig (−π, π]: (1.10575579 + 0.006400529086j); r=0.8, Zweig [0, 2π): (0.9861848451 - 0.1602791183j). Siegel-Fall: Cesàro-Fehler 2.8e-06 bei N = 2·10⁵.
- `exp.bruchsummen` — aktuell, 7/7 Prüfungen. Grenzwertformel: Fehler 5.2e-21 bei n = 120, Rate ≈ λ; danach wächst die Rundung (Faktor 1/λ je Rückschritt). Bruchsumme S(1/2) = 0.4277749202048676, S(x+1) − S(x) = T(x) exakt. Bruchprodukt-Formel regulär bis 3.8e-36, für T∘(id+θ) verletzt um 0.0185…0.2.
- `exp.koenigs_basiswechsel` — aktuell, 6/6 Prüfungen. √2 und d = 1.346032931 (λ_d = (ln 2)²): Ψ konjugiert f_d zu f_√2∘f_√2; Koenigs-Höhe exakt 2x; beliebiges Verhältnis 1.7 ebenfalls linear. Komplexe Spiralen: Ψ bildet den d-Orbit auf jeden zweiten b-Punkt ab.
- `exp.turmintegrale` — aktuell, 7/7 Prüfungen. Sophomore bestätigt; gerade 0.6767142737718128, ungerade 0.6313093708019882, Lücke 0.0454049029698246 (= Bild 2 auf 15 Stellen). ∫h = 0.62544581058879. I_n konvergiert wie n^(−1.502) (Vorhersage 3/2).
- `exp.gleiches_L` — aktuell, 21/21 Prüfungen. Basen mit gemeinsamem L bilden Familien u_k = (Log L + 2πik)/L. Reell: e^{π/2+2πk} (L = i), −exp((π/2+2πk)π/W(π)) (L = i·W(π)/π) und 8 gemischte Paare. Am gemeinsamen L ist α_b(T_d(x)) = C1 + C2 x + O(λ_d^x), C2 = log λ_d/log λ_b — keine bleibende Welle; mit verschiedenen L wird α_b(T_d) konstant.
- `exp.koenigs_kalkuel` — aktuell, 12/12 Prüfungen. Spektralreihe, Ableitungs-, Integral- und Shift-Kalkül und Lambert-Bruchsumme für √2 bestätigt. K = T_√2∘A_1.3 ist in Koenigs-Koordinaten die Potenz ρ = 0.38495759 (Verzweigungsexponent gemessen 0.38495759); Jacobi-Formel über Bruchsummen bestätigt. Quasikonformer Basiswechsel konjugiert, ist aber nur bei gleichem Multiplikator holomorph. Normalform e^w + κ bestätigt.
- `exp.eta_digamma` — aktuell, 12/12 Prüfungen. Basis η: C_T = 2.7982481542313876625 (aus ganzen Höhen), T_η(1/2) = 1.257153075054173, Σ_k[T(k) − e + 2e/k] = 10.79756227493931, Σ T′(k) = 1.4502290173934; Digamma-Zerlegung bestätigt. C_2 = 0.8125152587890625, C_e = 0.4338677391894938.
- `exp.lambert_zeta` — aktuell, 8/8 Prüfungen. Möbius-Modell: Lambert-Bruchsumme = q-Digamma (q = λ), Grenzwert ψ. Echte reguläre Tetration: T_b(1/2) → T_η(1/2) (Ordnung 2.0 in 1 − λ), S_b(1/2) → S_η(1/2). Kneser von oben: η+1e-1: 0.0421, η+1e-2: 0.00437, η+1e-3: 0.000438.
- `exp.abel_bruchsumme` — aktuell, 9/9 Prüfungen. Kompensierte Bruchsumme von slog_e konvergiert (A'' = O(1/y²)); S(1/2) = -0.115097609083, γ_A = 0.4951728249; Differenzengleichung und S'' = −Σ A'' bestätigt.
- `exp.funde` — aktuell, 14/14 Prüfungen. Original-Zahlen unabhängig bestätigt: −19/25 ↔ d = (1.425186736 - 0.04620336949j), σ = (-0.2004608936 + 0.003486040177j), Residuumsquotient → λ_d; A_e(L_2), K₁, A_2(L_e) (bis auf Zweig); b_links = (-0.08052815762548356 - 9.423207371442147e-52j); Spitzenkonstanten; Rechteckwelle |c_k|πk = |zA − zB|.
- `exp.superwurzel` — aktuell, 10/10 Prüfungen. ssrt(2) = 1.55961046946, ssrt(10) = 2.50618414559; Reihe um 1 bestätigt; Lösungszahl 0/2/1; dritte Superwurzel von 2: 1.47668433736; unendliche Superwurzel = y^(1/y) genau auf [1/e, e].
- `exp.tpid6` — aktuell, 4/4 Prüfungen. x_n − η ≈ 2π²η/(e n²): angepasst K = 10.49089 gegen 10.49067 (β = -3.77); bei n = 3·10⁵: n²(x_n − η) = 10.49057.
- `exp.qbinomial` — aktuell, 10/10 Prüfungen. q-Binomial-Reihe = reguläre Tetration (√2 und 1.3, reell und komplex), Fehler bei N = 80 unter 1e-30, geometrisch mit Faktor ≈ 0.395 je Glied (= q^(z+2), Singularität bei z = −2); nötig sind etwa 0.08 N² Stellen. Koenigs-Koeffizienten c_1…c_60 negativ für vier Basen. Möbius: reguläre Iterierte an beiden Fixpunkten gleich.
- `exp.cheta` — aktuell, 13/13 Prüfungen. usexp_√2(1/2) = 5.27364736917811 (ganz, Periode 19.2361i); cheta(1/2) = 9.68935067535323, C_ch = -1.25696115887251; cheta(x) − e ≈ 2e/|x| für x → −∞. Szekeres: Abel-Funktion von e^u − 1 total monoton bis Ordnung 8, Störung 0.001·sin 2πA verletzt das bei u = 0.0025.
- `exp.pentation` — aktuell, 9/9 Prüfungen. L = -1.8503545290271814185, λ = 6.46067129568184; e↑↑↑(1/2) = 1.63232474043606, e↑↑↑(−1/2) = 0.491054338635648, e↑↑↑1.5 = 6.48573606423933; Fox' Werte (2007) auf ≥ 24 Stellen bestätigt; Basis 2: 2↑↑↑1.5 = 2.702890571.
- `exp.tetra_euler` — aktuell, 10/10 Prüfungen. η₄ = 1.635324496715276399345345, e₄ = 3.088532271806717654482181 (Levenstein 2010 bestätigt); Fixpunktzahl 3 → 1 bei η₄; 1.6↑↑↑∞ = 2.26438440348, 1.7↑↑↑n divergiert; Selbst-Tetrawurzel maximal bei e₄.
- `exp.vermutungen` — aktuell, 9/9 Prüfungen. TPID 7 hält regulär und für Kneser bis b₇ = 10.2210634158 (T″(0) = T′(0)²), darüber nicht (Basis 100: 1.316); Re T ≥ 1 auf der rechten Halbebene (vier Basen); ∫_{−1}^0 T_b = 1/2 bei b = 2.69370284153056; Carleman-Eigenwerte → λ^k (1.3: Abstand 3.0e-16 bei N = 32).
- `exp.shell_thron` — aktuell, 6/6 Prüfungen. Fläche = 8.0509332809395096234; grösster Imaginärteil bei b = (0.815202742507 + 2.01664021991j); linkester Punkt (-0.380979728648 + 0.899700195546j), rechtester (2.04779052746 - 0.842045503531j) (beide abseits der reellen Achse).
- `exp.ableitungen_vorzeichen` — aktuell, 6/6 Prüfungen. Ungerade Ableitungen bis 25 überall positiv (20 Stellen); gerade wechseln genau einmal: x₂ = -0.52176479, x₁₀ = -0.185914, x₂₂ = 0.00779022, x₂₆ = 0.0465848; log(x+2)-Gesetz bestätigt.
- `exp.stufe5` — aktuell, 7/7 Prüfungen. η₅ = 1.7373559251516875731, e₅ = 3.3014072785233674934 (Pentation regulär am negativen Fixpunkt von T_b); Etas 1.44467, 1.63532, 1.73736; Eulers 2.71828, 3.08853, 3.30141.
- `exp.hexation` — aktuell, 10/10 Prüfungen. Reguläre Hexation unter η₅: 1.5↑↑↑↑½ = 1.32331315362, 1.5↑↑↑↑∞ = 1.63659373391, 1.6↑↑↑↑½ = 1.36258513424, 1.6↑↑↑↑∞ = 1.88635097523, 1.7↑↑↑↑½ = 1.39744726135, 1.7↑↑↑↑∞ = 2.3876372817; H(−3) = −2, Singularität bei s ∈ (−4, −3); über η₅ kein reeller Fixpunkt der Pentation; Basis 2: p anziehend, q abstossend.
- `exp.hexation_kneser` — aktuell, 11/11 Prüfungen. 2↑↑↑↑1.5 = 2.6728939593426194892 (Kneser an q, 31.9 Stellen); η₆ = 1.7991212361288126867, e₆ = 3.44223464771689; Etas 1.44467, 1.63532, 1.73736, 1.79912; Eulers 2.71828, 3.08853, 3.30141, 3.44223.
- `exp.ln_basis` — aktuell, 9/9 Prüfungen. ln b = T′(0)/T′(−1) = slog′(0)/slog′(1) und die Formen zweiter/dritter Ordnung an Kneser b = 2, e, 10 bestätigt, auch an einer gestörten Lösung; b₇ ⇔ Wendepunkt bei −1; Basisableitung: Summenformel ganzer Höhen und die lineare Differenzengleichung bei gebrochener Höhe bestätigt.
- `exp.turm_pde` — aktuell, 4/4 Prüfungen. Kobis Gleichung gilt für ganze Höhen exakt; für Kneser ist p = ∂_b slog − S·slog′ periodisch mit Mittel -3.15733299588 (b = 2) und Spanne 0.00057 — Kneser erfüllt sie bei gebrochener Höhe also nur fast.
- `exp.forum_formeln` — aktuell, 5/5 Prüfungen. Fox' rationale Fixpunktsummen folgen aus (e^z−1)/(e^z−z); die Newton-Reihe der ganzen Iterierten trifft die reguläre Tetration; x↑↑n stimmt bis (ln x)^n mit dem unendlichen Turm überein; i^i^i^… = 0.43828 + 0.36059i.
- `exp.halbe_exp` — aktuell, 5/5 Prüfungen. f = exp^[1/2] über Kneser: f(f(x)) = eˣ auf 1e-28, f(0) = 0.498563287941, f(1) = 1.64635423375, Nullstelle -0.696024740886.
- `exp.zwei_sinh` — aktuell, 4/4 Prüfungen. 2sinh-Tetration: T(½) = 1.6466507303037 gegen Kneser 1.64635423375119; die Höhenverschiebung θ wackelt mit Spanne 0.0002074.
- `exp.forum_rand` — aktuell, 7/7 Prüfungen. Ackermann A(4,n) = 2↑↑(n+3) − 3, f₃ ≥ 2↑↑n, Σ1/EF(n) = 1.6111149258…, Grenzfälle bei 0, quadratische Näherung C¹ nicht C², Nachbar-Turm = ln x/W(ln x).
- `exp.forum_iteration` — aktuell, 4/4 Prüfungen. Neun elementare Superfunktionen bestätigt; Möbius-Familie analytisch durch den parabolischen Fall; TPID 9 konvergiert; TPID 13 numerisch plausibel.
- `exp.turmintegral_periodisch` — aktuell, 7/7 Prüfungen. I(2m + r) → 2-periodische Grenzkurve I_∞(r) mit festen Ecken lim I₂ₙ, lim I₂ₙ₊₁, Mittel M = 0.654011822287, nur ungeraden Harmonischen und monotonen Flanken; die Form ist die Verteilung der Kreuzungshöhen und hängt von der Interpolation ab.
- `exp.grundlagen_identitaeten` — aktuell, 11/11 Prüfungen. Tauschregel, Taylor-Leiter (Bell-Polynome), Stufenleiter ln ln a = m-te Differenz von ln H′, Frei-und-fest unter T∘(id+p) und Exponentialität über slog an Kneser-Werten bestätigt.
- `exp.polynom_halbexp` — aktuell, 7/7 Prüfungen. Das Abschneiden p_N∘p_N = exp + O(x^(N+1)) trifft Kneser: |p_N(1) − T_e(1/2)| = 1.0e-16 bei N = 100, Rate ≈ 0.924 pro Grad; jeder Koeffizient und der Randwert p_N(L) → L konvergieren.
- `exp.halbexp_korridor` — aktuell, 7/7 Prüfungen. Ω = slog(D) ist waagrecht konvex, Re ≤ 0.31218, Breite ≥ 2.067 (oben → 2.3495), linker Rand mit Steigung −0.238; der Korridor Ω ∩ (Ω − ½) ist überall breiter als 1, und die Polynome p_N bilden ihn in D ab.
- `exp.polynom_schranke` — aktuell, 5/5 Prüfungen. Streng (Krawczyk) für N ≤ 120: eindeutige Lösung und |p_N| ≤ 2.2405 auf |x| ≤ |L|. Für alle N scheitert Kantorovich: |J⁻¹| wächst um ≈ 2.57 pro Grad, der Rest steht fast senkrecht auf den fast-singulären Richtungen bei L und L̄.
- `exp.polynom_team` — aktuell, 5/5 Prüfungen. Team-Daten bestätigt: max|Π_N f_K| = 2.232651 bis N = 1000, Randfehler ~N^-4.85, Innenfehler ~N^-7.91 (Soll −7.894), echte p_N quadratisch nah, Kneser-Koeffizienten = Darboux-Leitterm an L.
- `exp.polynom_phase3` — aktuell, 5/5 Prüfungen. ‖Ĵ_N⁻¹‖ = 1.251155 (kubisch, N = 100–200 konstant), Kneser 1.250337; ohne Linse wächst sie (1.366 → 1.852); (K2) mit den zertifizierten Koeffizienten verträglich; Eckprofil-Integral bestätigt.
- `exp.tpid10_abstaende` — aktuell, 8/8 Prüfungen. 2 − ηₙ fällt wie n^(−1.47); Vorhersage η₇ ≈ 1.8399 gegen Levenstein 2017 ≈ 1.839. TPID 11 (1) gilt für n = 4…6; die Quotienten 1.594, 0.9208 entscheiden (2) nicht.
- `exp.tpid3` — aktuell, 5/5 Prüfungen. TPID 3: s=3/4, p=2/3 (anziehend, x0 im Kreis): 6.1e-11 bei N = 32; s=3/4, p=3/2 (abstossend, x0 im Kreis): 0.00034 bei N = 32; s=3/10, p=2 (Polynom, x0 = 0.7): 1.0 bei N = 32; s=1/2, p=1/2 (x0 auf dem Rand): 0.23 bei N = 32; s=3/10, p=1/2 (x0 ausserhalb): 0.42 bei N = 32
- `exp.superwurzel_omega` — aktuell, 6/6 Prüfungen. Superwurzel-Formel exp(ω_{n−1}(ln Y)), h(x) = 1/ssrt(1/x), geschlossene Basisableitung und Ableitungen von ε_n, ω_n auf ≥ 20 Stellen bestätigt; die Grenzidentität gilt nur, wenn n → ∞ nach dem Einsetzen von ω_n geht.
- `exp.stoerungsreihe` — aktuell, 7/7 Prüfungen. Forumspolynom (tid 1827) als Störungsreihe von eˣ − 1 hergeleitet und exakt nachgebaut; ∫|Fehler| bei N = 8 = 2.41e-7, beste Ordnung N = 6, danach Divergenz. Näherungsleiter für Desmos gemessen.
- `exp.verwandte_funktionen` — aktuell, 11/11 Prüfungen. Produkt-Tetration: Basiswechsel als Konjugation, parabolischer Fixpunkt 0 (nicht formal konjugiert zu eˣ − 1), Fatou-Wert prodtet_e(½) = 1.5134281 (Forum 1.5134). Exponentielle Fakultät: D = -1.055214008713 (Forumswert −0.0613 nicht reproduziert), Fortsetzung konvergiert, 0.5!₂ = 0.7575097219, x!₂ → 0.573049 für x ↓ 0 (γ = 0.57722).
