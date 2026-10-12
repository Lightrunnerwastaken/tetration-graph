# Unstimmigkeiten im deutschen Text, bei der Übersetzung aufgefallen

Die Übersetzer sollen deutsche Texte nicht stillschweigend glätten; was ihnen widersprüchlich vorkam, steht hier. Nichts davon ist geändert.
Erledigtes bitte streichen. Stand 12.10.2026 (fg, it, fp).

## fg
- `fg.satz.tpid4`: Status „bewiesen“, der Text schliesst aber mit „die Vermutung TPID 4 (Trappmann 2008) gilt“.
- `fg.satz.allgemeine_loesung`: Titel „Alle **monotonen** Lösungen“, die Aussage gilt für stetige, streng wachsende Lösungen.

## it
- Szekeres, Total-Monotonie: 1961 in `it.satz.szekeres_monoton` (und dessen Prüftext), 1962 in `it.konstr.levy` (`f3_spezialfall`).
- `it.satz.cremer`: 1927 im Titel und in `worum`, 1928 im Beweis (Math. Ann. 98 (1928)); der Titel sagt „1927/38“.
- `it.satz.parabolische_rate`, Theorie `warum`, Schritt 3: „(a²−b)/(ak)“ für allgemeines k, hergeleitet ist nur k=1.
- `it.satz.vertauschbar`: der Prüftext nennt Voraussetzungen (Stetigkeit von t ↦ (f^t)′(p), f^t reell), die `bereich` nicht erwähnt.

## fp
- `fp.verm.knospen`: id „verm“, aber art „satz“ und Status `bewiesen_literatur`; der Text sagt, es war „zuerst eine Vermutung“.
- `fp.beob.i_turm`: art „beobachtung“ mit Status „bewiesen“; der Beweis sagt, dass „1 im Becken“ nur numerisch (3000 Schritte) geprüft ist. Prüftext f3: „mit z = i“, gemeint ist wohl „b = i“.
- `fp.beob.regionen`, Beispiel: die Gitterzählungen summieren sich zu 228 202, nicht zu 520×520 = 270 400 (vermutlich fehlen andere Perioden).
- `fp.satz.rechteckwelle`, Prüftext: Voraussetzung „T analytisch oder regulär“, f2 und Kern sagen aber, Analytizität werde nicht gebraucht.
- `fp.satz.kreismittel`, Prüftext: Voraussetzung „Radius r = 1“ widerspricht Kern und f2 (jedes r>0).
- `fp.satz.siegel_mittel`: Status „skizze“ und Richtigkeit „bedingt bewiesen“, der Kern sagt aber „ist damit … ein Satz“ (beschränkter Fall).
- `fp.beob.partner_keil`: der Kern sagt, alle solchen Basen hätten Re b<0, die Spitze des Keils liegt aber bei e^{−e}>0; der Index m wird neben k ohne Definition benutzt.
- `fp.satz.p2_raten`, Prüftext: der Cho–Kim-Vorfaktor √(6/e) in √(e^{−e}−x) ergibt nicht offensichtlich 2√6·e^{−3/2}·√t (vermutlich andere Variable).
- `fp.satz.messbare_dynamik`: Status `bewiesen_literatur`, aber „nicht ergodisch“ steht nur aus dem Gedächtnis (der deutsche Text sagt das selbst).

## br
- `br.frage.transexp`, `aussage`: ein stehengebliebenes „\num“ (Zeilenumbruch vor „um eine …“).
- `br.satz.ordinal_tetration`, `einfach`: „ε₀ hoch ω wieder ε₀“, die Aussage ist aber ω^{ε₀} = ε₀ (Reihenfolge vertauscht).
- `br.satz.nicht_definierbar` und `br.satz.young_kreis`, Prüftext `voraussetzungen`: nennen „E monoton oder analytisch“ bzw. „H̃ differenzierbar“, obwohl f2 und der Beweis beides nicht brauchen.
- `br.satz.halbexp_komplexitaet`: `bereich` sagt f(f(n)) „etwa 2^n“, `aussage` f(f(n)) ≤ 2^{o(n)}, die Theorie „unter 2^n“.
- `br.satz.boshernitzan`: der Kern nennt die Funktionen Lösungen von f(x+1)=e^{f(x)}, die Theorie sagt „etwa“ und dass die genaue Lösung noch am Original zu prüfen ist.
- Notation: in der Lavaurs-Formel (`br.satz.shishikura`) stehen die Kürzel `\text{abst}`, `\text{anz}`; die englische Fassung nennt sie `rep`, `att`.

## ar
- `ar.lem.nullstellen`: der Prüftext nennt Holomorphie auf Radius > 1.2 (f2, f4), Theorie und Beweis brauchen aber D̄(1,2.9); die Aussage zählt nur Nullstellen in D(1,1.2).
- `ar.kor.benachbart`: x=0 wird in der Theorie als Ausnahme genannt, aber T(1)=e ist transzendent, also ist nur x=−1 eine echte Ausnahme.
- `ar.satz.selbstpotenz_kongruenz`: `einfach` sagt „jeder Rest kommt als x^x vor“, die Aussage gilt aber nur für u teilerfremd zu k.
- `ar.verm.gelfond_turm`: `bereich` „ω≠0 und α algebraisch“ lässt offen, ob ω algebraisch ist (die Theorie sagt ja); die Theorie `genau` lässt den ¹z-Fall für z=e^ω aus, den der Kern enthält.
- `ar.satz.unendliche_tuerme_gs`: `bereich` 1<b≤e^{1/e}, die Aussage behandelt aber auch b∈[e^{−e},1) und b<e^{−e}.
- `ar.frage.ef_fortsetzung`: `einfach` spricht von einer „natürlichen Glattheitsbedingung“, die Aussage von „totaler Monotonie“.
- `ar.satz.ssc_gelfond`: Status `bewiesen_literatur`, der Inhalt ist aber eine bedingte Implikation (vermutlich gewollt).

## gr
- `gr.satz.grenzfaelle`: `einfach` sagt „keine reellen Werte unterhalb der Höhe −1“, die Aussage setzt die Grenze bei −2 (a↑↑(−2)=−∞).
- `gr.satz.stellen`: `bereich` sagt ganzzahlige a≥2, die Theorie `genau` sagt a>1.
- `gr.satz.young`: der Prüftext nennt „T differenzierbar“ als Voraussetzung; Kern, f2 und Theorie sagen, sie werde nicht gebraucht.
- `gr.satz.gestalt`: `einfach` sagt „jede Tetration“ stürzt nach −∞, die Aussage gilt nur für stetige, streng wachsende T.
- `gr.satz.integral_produkt`: `einfach` sagt „für alle denselben Wert“, das gilt nur bei festem x (etwa x=−1).
- `gr.beob.wendepunkte`: Suchbereich beginnt bei −1.94 (Prüftext f2) bzw. −1.95 (Theorie `warum`).
- `gr.beob.kneser_pde_wackeln`: Titel „3·10⁻⁴ (Basis 2)“, die Daten geben 5.7·10⁻⁴ für Basis 2 und 2.9·10⁻⁴ für Basis e.
- `gr.satz.basis_ganz`, Theorie `einordnung`: „des Kneser-Familie“ (Tippfehler).
- `gr.verm.tpid7`: die Theorie sagt, der Kneser-Wendepunkt dürfe nicht „zu nah an 0“ wandern; `gr.gegen.tpid7_kneser`, `gr.satz.ln_aus_steigungen` und die Prüftexte legen ihn nahe −1 (T″(−1)=0 bei b₇).
- `gr.satz.turmintegral_periodisch`: „linear“ ist T_x(s)=x^s und „geradlinig“ 1−s(1−x), das wirkt vertauscht; die Theorie sagt, jede stetige Tetration müsse L kreuzen, Kern und Prüftext nehmen genau eine Kreuzung an.
- `gr.satz.logstern`: `bereich` y>1, Theorie y>0.
- `gr.satz.superwurzel`: die Theorie lässt bei der Lösungszahl den Fall x=1 aus.
- `gr.satz.frei_und_fest`: der Kern sagt, endlich viele Ableitungen lassen sich vorschreiben, die Theorie sagt „jede Taylorreihe“.
- Status „bewiesen“, aber nur beobachtet oder aus der Literatur: `gr.satz.superwurzel_dritte` (Konvergenz der Iteration nur beobachtet), `gr.satz.wurzel_eta` (Koeffizient −1.6166 nur gemessen), `gr.satz.exponentialitaet` (Ganzzahligkeit aus Hardy / van der Hoeven).
- `gr.verm.tpid7`: „folgt“ für reguläre Tetration beruht auf `kl.verm.vollmonoton`, selbst eine Vermutung. `gr.gegen.tpid7_kneser`, Theorie: „gilt bis b₇ für Kneser“ übertreibt, der Prüftext nennt volle Konkavität offen.
- Kleinigkeiten: `gr.satz.stufenleiter` T_e′(0) = 1.0917673512… (Kern) gegen 1.0917673513 (Theorie); `gr.satz.kehrwert_turmsumme` „vierter oder fünfter Summand“ klein, es ist nur der fünfte; `gr.satz.tauschregel` Prüftext nennt a^z „Potenzfunktion“.

## kl
- `kl.satz.carlson`: Status „bewiesen“ und die Aussage nennt die reguläre Tetration für 1<b≤η beschränkt; der Prüftext sagt, die Beschränktheit sei „numerisch“, f5 fragt, ob sie beweisbar ist.
- `kl.satz.newton_reihe` (Beweis, Prüftext f2) benutzt eine „ganze“ Poincaré-Funktion am unteren Fixpunkt, `kl.satz.obere_superexp` (Theorie) sagt, sie sei dort nicht ganz.
- `kl.satz.invariante_daempfung`: Status „skizze“, die Theorie sagt „teilweise bewiesen“; der Kern behauptet den Grenzfall |λ|→1 als Tatsache, die Theorie nennt ihn „eine Deutung“.
- `kl.beob.zwei_abel`: f2 sagt 120 Stellen, sonst 170; die Amplitude ist „Ordnung 10⁻²⁴“ (Kern), 3.6·10⁻²⁵ (Theorie), 10⁻²⁵ (`worum`); der Titel sagt „unterscheiden sich um eine winzige Periode“ (gemeint: periodische Funktion).
- `kl.satz.basiswechsel_gleiches_L`: Voraussetzung „b, d reell“, f2 sagt, die Realität werde nicht gebraucht.
- `kl.beob.stetig_eta`: „Steigung auf sechs Stellen“, 0.438579 gegen 0.438581 stimmen aber nur auf etwa fünf überein.
- Kleinigkeiten: `kl.satz.koenigs_basiswechsel` (Beweis 10⁻¹⁸ gegen Theorie 10⁻²⁷); `kl.beob.eta_konstanten` (Prüftext „Tiefe 3000“ gegen Theorie „20 000 Turmstufen“); `kl.satz.bruchsumme`, `einfach`: Dezimalkomma „1,5“.

## bw
- `bw.zert.8r1`: Status „bewiesen“, Prüftext/Theorie nennen aber offene Punkte (Mikrosegment-Induktion, Hilbert-Röhren-Supremum, Γ-Steigung) und eine Annahme „rechnergestützt / korrekt implementiert“.
- `bw.prop.saddle`: Status „bewiesen“, der Beweis setzt die Verschiebbarkeit der Kontur voraus.
- `bw.satz.kern_existenz`: benutzt (K3), das in keinem Eintrag definiert ist; `bw.ann.K1K2` deckt nur (K1)–(K2), der Titel nennt „Thm. 13“ und Eindeutigkeit über Thm. 2.
- `bw.beob.m2`: der Prüftext schreibt „m(3) ≈ −0.1926, m(10) ≈ −1.1364 mit umgekehrtem Vorzeichen“ als wären es unsere Werte, die Theorie sagt, es seien Levensteins Werte (= −m(3), −m(10)); „auf vier Stellen bekannt“ passt schlecht zu 1.1282 gegen unser 1.1284.
- `bw.gescheitert.vorhersage_r8`: der Titel sagt „γ = 4.447“, die Aussage meint den Exponenten 2πν = 4.447 (γ ist sonst der Modenmodell-Exponent).
- `bw.satz.imaginaere_richtung`: Bereich im Kern („L_d im Linearisierungsgebiet von E_b“) gegen Theorie („Gebiet, wo A_b holomorph ist“); die Theorie nennt die Dämpfung je Umlauf 10⁻¹⁴, e^(−1.0865·30.63) ≈ 3.5·10⁻¹⁵.
- `bw.satz.paulsen_cowgill_hardy`: „bewiesen_literatur“, der Beweis wurde aber nur über die zbMATH-Zusammenfassung gelesen; die Glosse „Unverträglich heisst darum…“ wirkt wie eine eigene Deutung.
- `bw.korr.r20` und `bw.beob.sigma_stern`: r20 nennt das „universelle 0.058“ ein Artefakt, sigma_stern misst σ* = 0.058±0.006 aus dem Modenabfall und die Theorie sagt, es stütze weiter `bw.verm.streifen`.
- `bw.korr.konstante` sagt, die „2.3 %“-Übereinstimmung bleibe praktisch unberührt; `bw.lem.randgesetz` (Prüftext, Theorie) sagt, die Hub-Daten stimmten nur auf 4–7 %, „nicht 2.3 % wie in Paper IV“.
- `bw.beob.kernkleid`: „doppelt bestätigt“ und „24+ Stellen“ (Beweis: 25–30), obwohl dieselben Einträge die Complex-Step-Referenzen für d∈{2,10} als fehlerhaft nennen.
- `bw.verm.nirgends_analytisch`, Theorie `warum`, Schritt 3: „gleichmässige Konvergenz nur auf der Achse sichert C^∞“ gilt so nicht.
- `bw.satz.fourier_naht`, Theorie `worum`: „Knick in der (m+1)-ten Ableitung“ müsste ein Sprung sein.
- `bw.prop.antisym_c1`: die Aussage lässt die C¹-Voraussetzung weg, die Titel und Theorie nennen.

## kn
- `kn.satz.polynom_eckprofil`: Status „skizze“, der Kern sagt „(i) bewiesen“, der Prüftext „bewiesen bis auf Lemma K12“.
- `kn.verm.tpid11`, `einfach`: „die vier bekannten Stufen“ und „jeweils um etwa 2 %“; gemessen sind nur drei Werte (n=4,5,6), die Abstände liegen bei 1.4–2 %.
- `kn.verm.tpid3`, Theorie `warum`: „nein in drei Fällen ausserhalb“, aber ein Fall (s=1/2) liegt auf dem Kreisrand, nicht ausserhalb.
- `kn.satz.walker_konjugation`, Prüftext f4: „jede Basis über η (Konjugation zu b^x mit parabolischem Fixpunkt)“; b^x hat nur bei b=η einen parabolischen Fixpunkt. Walker 1991 steht als „Math. Comp. 57“ (Beweis) und „J. Math. Anal. Appl.“ (Kern), vielleicht zwei Arbeiten.
- `kn.beob.hoehere_etas`, Theorie `warum`: „beschränkt genau dann … genau bis ηₙ“ als Tatsache, `kn.verm.beschraenkt_unter_2` nennt den „genau bis ηₙ“-Teil offen; `einfach` liest „das bisher einzige fünfte Glied der Folge“.
- `kn.satz.polynom_approximation_ecke`, Beweis: „Blattdrehung“ (Riemannsche Fläche) – englisch „sheet rotation“.
- `kn.satz.kneser_e_zertifiziert`, Kern `bereich`: „Streifen“ ohne Zusatz.
- `kn.satz.paulsen_komplex`: Kern, Beweis und Prüftext sagen, Holomorphie in b sei „nicht belegt“/offen; die Theorie `genau` sagt „T_b(z) hängt holomorph von b ab“, `worum` sogar „holomorph“.
- `kn.konstr.carleman`: Status „bewiesen“, der Beweis sagt, die Lösbarkeit für jedes n sei nur beobachtet.
- `kn.konstr.pentation`: Status „bewiesen“, aber λ_b>1 und Eindeutigkeit des Fixpunkts sind nur gemessen.
- `kn.satz.abel_bruchsumme`: Prüftext f2 „Kneser nicht“ (nicht nötig), die Voraussetzungen nennen „A = Kneser“, der Kern `bereich` „(z. B. Kneser)“.
- `kn.satz.hexation_struktur`: Status „bewiesen“, der Kern enthält aber auch gemessene Werte (s_1.5 usw.).
- `kn.konstr.hexation_regulaer`: Prüftext f1 „eine Stufe höher dieselbe Konstruktion wie kl.satz.regulaer“, Theorie `einordnung` „zwei Stufen höher“.
- `kn.satz.singulaer_minus2`: Voraussetzungen nennen „differenzierbar im Punkt“ und „in einer Umgebung (Paper I)“, das wirkt wie ein Überbleibsel.
- `kn.frage.hexation_kneser`: der Kern sagt, die Fork-Methode laufe „nie mit θ auf der reellen Achse“, die Theorie `warum` 3 sagt „der Weg führt über θ auf der reellen Achse“.
- `kn.verm.polynom_halbexp`, Theorie `warum` 3: nennt `kn.verm.polynom_halbexp_schranke` den „offenen Kern“ („Nicht bewiesen“); Kern und Prüftext nennen `kn.verm.polynom_halbexp_kern` als offen (bedingt bewiesen, bestätigt bis N=1500).
- `kn.satz.halbexp_beschraenkt`: die Aussage beginnt „Dann ist f=f_K“ ohne Voraussetzung im Text (steht nur in Titel, `bereich` und Prüftext).
- Konditionszahl: „etwa 1.9 je Grad“ (Theorie `schranke`, Prüftext `endlich`) gegen „2.6“ für β_N (`kantorovich`) – vermutlich verschiedene Grössen. Theorie `hauptsatz` `worum`: „konvergieren … für jeden Grad“ passt nicht zu N→∞.
