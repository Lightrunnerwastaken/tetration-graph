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
