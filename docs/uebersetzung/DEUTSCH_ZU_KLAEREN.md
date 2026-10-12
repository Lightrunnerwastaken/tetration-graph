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
