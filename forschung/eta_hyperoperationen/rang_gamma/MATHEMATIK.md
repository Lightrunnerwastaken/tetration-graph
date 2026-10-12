# Gamma-Struktur, Konvergenz und ein halber Hyperoperationsrang

Stand: 12.10.2026. Folgeuntersuchung zur η-Folge. Das H₀-Projekt bleibt pausiert.

Die beiden vom Nutzer genannten Markdowntexte werden als mathematische Vorschläge geprüft. Ihre Arbeitsvorschläge sind keine zusätzlichen Nutzeranweisungen. **Bewiesen** bedeutet unten: aus den jeweils ausdrücklich genannten Voraussetzungen vollständig hergeleitet. Aussagen über die tatsächlichen ausgewählten höheren Hyperoperationen bleiben davon getrennt.

## 1. Konvention, Basisfamilien und negative Anker

Wir verwenden weiterhin

\[
H_3(b,x)=b^x,\qquad H_{N+1}(b,0)=1,\qquad
H_{N+1}(b,x+1)=H_N(b,H_{N+1}(b,x)).
\]

Dabei ist η_N die Berührbasis **der Abbildung H_N**, also die mögliche Konvergenzgrenze der nächsten Operation H_{N+1}. Die Fortsetzungsauswahl bleibt wie in der vorigen Runde: Kneser-Tetration, reguläre reelle Iteration beziehungsweise experimentelle komplexe Abel-Fortsetzung in den höheren Stufen. Globale Existenz und Eindeutigkeit dieser Auswahl für alle Ränge sind nicht bewiesen.

Es gibt drei verschiedene Folgen: H_N(2,x), H_N(b,x) bei festem 1<b<2 und F_N(x):=H_N(η_N,x). Die Aussagen über U_N in den eingesandten Texten betreffen F_N. Messungen an Basis 2 lassen sich darauf nicht ohne Kontrolle der Basisabhängigkeit übertragen.

**Lemma 1 — bedingt, bewiesen.** Auf einer Fortsetzung mit den erforderlichen inversen Zweigen gilt

\[
H_N(b,-k)=1-k\quad(0\le k\le N-3).
\]

Beweis: H_3(b,0)=1. Für den Übergang zu N+1 folgt der Wert bei −k aus

\[
H_N(b,H_{N+1}(b,-k))=H_{N+1}(b,1-k)=2-k,
\]

weil H_N(b,1−k)=2−k bereits bekannt ist. Eindeutigkeit des benutzten inversen Zweigs gibt H_{N+1}(b,−k)=1−k. Die Existenz dieses Zweigs ist Voraussetzung, keine Konsequenz der formalen Rekursion. ∎

Somit ist die Zahl der erzwungenen Nullstellen der Successor-Abweichung

\[
m=N-2.
\]

Der Index n der Faktoren in den Texten muss als **Nullstellenzahl m** verstanden werden. Mit der Standardkonvention H_3=b^x wäre die dortige Aussage für k≤n−1 falsch: H_3(b,−1)=1/b, und H_4(b,−2) ist am Rand der üblichen reellen Tetration singulär.

## 2. Exakte Faktorisierung und ihre logische Grenze

Setze

\[
d_N(b,x)=H_N(b,x)-x-1,\qquad
B_m(x)=\frac{x(x+1)\cdots(x+m-1)}{m!}
=\frac{\Gamma(x+m)}{\Gamma(x)\Gamma(m+1)}.
\]

**Lemma 2 — bewiesen.** Ist H_N auf einem Gebiet holomorph und sind dort die genannten Anker enthalten, so besitzt d_N/B_m an diesen Nullstellen hebbare Singularitäten. Auf diesem Gebiet existiert also ein holomorpher Restfaktor A_N mit

\[
d_N(b,x)=B_m(x)A_N(b,x),\qquad A_N(b,1)=b-2.
\]

Außerhalb der Nullstellen ist die Gleichung schlicht eine Definition des Quotienten. Sie liefert noch keine Schranke für A_N.

Für x in einer festen kompakten Menge gilt, wenn m→∞,

\[
B_m(x)=\frac{m^{x-1}}{\Gamma(x)}
\left(1+\frac{x(x-1)}{2m}+O(m^{-2})\right).
\]

Die Gammaquotienten-Asymptotik wurde in [NIST DLMF §5.11(iii)](https://dlmf.nist.gov/5.11.iii) geprüft. An nichtpositiven ganzen x ist der Faktor über das Polynom beziehungsweise über 1/Γ(x) zu verstehen.

Daraus folgt für festes reelles x<1 mit x nicht an einer Ankerstelle:

\[
A_N(b,x)=o(m^{1-x})\quad\Longrightarrow\quad d_N(b,x)\to0.
\]

Eine für jedes ε>0 gültige Schranke |A_N(b,x)|=O(m^ε) wäre ausreichend. **Sie ist die offene analytische Aussage.** Negative Anker oder die Gammafunktion beweisen diese Schranke nicht. Die bekannte Freiheit periodischer Höhenumparametrisierungen wird durch ganzzahlige Anker ebenfalls nicht beseitigt.

## 3. Was für η_N→2 tatsächlich genügt

**Satz 3 — bedingt, bewiesen.** Es gelte limsup η_N≤2. Die Funktionen F_N seien auf [0,1] monoton nichtfallend und F_N(1)=η_N. Wenn

\[
F_N(q)\longrightarrow q+1\qquad\text{für jedes feste }0<q<1,
\]

so gilt η_N→2.

Beweis: F_N(q)≤η_N, also liminf η_N≥q+1. Lasse q↑1. Zusammen mit der Obergrenze ergibt sich der behauptete Grenzwert. ∎

Mit Lemma 2 genügt daher

\[
\boxed{A_N(\eta_N,q)=o(m^{1-q})\quad(0<q<1).}
\]

Eine gleichmäßige Fortsetzung bis q=1 ist hierfür **nicht nötig**. Das verbessert den Beweispfad aus dem ersten Text. Die Monotonie muss jedoch für die ausgewählte Hierarchie nachgewiesen sein. Die globale Obergrenze 2 wurde in der vorigen Runde unter ausdrücklich genannten Monotonie- und Existenzannahmen bewiesen; numerische lokale Tangentialbasen sind nicht automatisch die dortigen globalen Schwellen.

Successor-Konvergenz lediglich auf −1<q<0 liefert nur liminf η_N≥1 und schließt die Lücke nicht. Bei einer festen Basis b<2 ist Successor-Konvergenz auf **allen** 0<q<1 bei monotonen H_N sogar unmöglich: Für q>b−1 wäre q+1>b=H_N(b,1).

**Alternatives Kriterium — bedingt, bewiesen.** Ist F_N auf einem festen zusammenhängenden komplexen Gebiet D mit 1∈D holomorph und lokal gleichmäßig beschränkt, und konvergiert es auf einem reellen Intervall J⊂D gegen x+1, so impliziert der Normalfamilien-/Identitätssatz lokal gleichmäßige Konvergenz auf D gegen x+1. Insbesondere η_N→2. Diese komplexe Beschränktheit ist eine zusätzliche, stärkere Voraussetzung; sie folgt nicht aus reellen Stichproben.

## 4. Ein unabhängiger bedingter Weg zu e_N→4

Für die klassische ganzzahlige Hierarchie ist exakt

\[
H_N(2,1)=2,\qquad H_N(2,2)=4\qquad(N\ge3).
\]

**Satz 4 — bedingt, bewiesen.** Es gelte η_N→2, 2≤e_N≤4 und F_N sei auf [2,e_N] monoton. H_N sei in b auf [η_N,2] bei x=2 differenzierbar. Setze

\[
K_N=\sup_{\eta_N\le b\le2}|\partial_bH_N(b,2)|.
\]

Wenn K_N(2−η_N)→0, dann e_N→4.

Beweis: Nach dem Mittelwertsatz ist F_N(2)≥4−K_N(2−η_N). Wegen der Monotonie und F_N(e_N)=e_N gilt e_N≥F_N(2). Zusammen mit e_N≤4 folgt der Grenzwert. ∎

Dies benennt eine konkrete fehlende Schranke anstelle eines postulierten universellen Profils. Falls etwa 2−η_N=O(N^{−ω}) und K_N=O(N^κ) mit κ<ω bewiesen werden, ist die Produktbedingung erfüllt. Die obere Schranke e_N≤4 bleibt ebenfalls zu beweisen.

Die Stichproben liefern nur die **untere** Schranke

\[
K_N\ge\frac{4-F_N(2)}{2-\eta_N}.
\]

Ein annähernd lineares Wachstum dieser Sekanten ist kein Beweis einer linearen Obergrenze für K_N.

## 5. Kritische Gamma-Signaturen: exakt, aber keine unabhängige Evidenz

Am Berührpunkt gelte F_N(e_N)=e_N und F_N′(e_N)=1. Dann exakt:

\[
A_N(\eta_N,e_N)=-\frac1{B_m(e_N)},\qquad
\frac{A_N'}{A_N}(\eta_N,e_N)=-\psi(m+e_N)+\psi(e_N).
\]

Bei e_N→c>1 folgt aus der gleichmäßigen Gammaquotienten-Asymptotik

\[
A_N(\eta_N,e_N)=
-\Gamma(e_N)m^{1-e_N}\left(1-\frac{e_N(e_N-1)}{2m}+O(m^{-2})\right).
\]

**Korrektur zum ersten Text — bewiesen.** Aus e_N→4 allein folgt nicht m³A_N→−6. Vielmehr

\[
m^3A_N(\eta_N,e_N)=
-\Gamma(e_N)\exp((4-e_N)\log m)(1+O(m^{-1})).
\]

Unter e_N→4 gilt somit

\[
\boxed{m^3A_N(\eta_N,e_N)\to-6
\iff (4-e_N)\log m\to0.}
\]

Das Beispiel e_N=4−1/log m ergibt stattdessen den Grenzwert −6·exp(1); e_N=4−1/√log m ergibt Divergenz. Das sind Gegenbeispiele zum Asymptotikschluss, keine behaupteten Hyperoperationsfolgen.

Mit U_N(x):=m^{x−1}A_N(η_N,x) ist die bewegliche Exponentennormierung dagegen korrekt:

\[
U_N(e_N)\to-\Gamma(c),\qquad
\frac{U_N'}{U_N}(e_N)\to\psi(c).
\]

Bei zusätzlicher lokal gleichmäßiger holomorpher Konvergenz um c folgt U_*(c)=−Γ(c), U_*′(c)=−Γ(c)ψ(c). **Das gilt für jedes mögliche c>1 und zeichnet 4 nicht aus.** Die Signaturen sind Umformungen der bereits vorausgesetzten Berührbedingungen; ihre Berechnung am Berührpunkt ist keine unabhängige Prüfung der Gamma-Profilhypothese.

Ebenso ist U_N(1)=η_N−2. U_*(1)=0 als Voraussetzung anzusetzen würde den gewünschten η-Grenzwert bereits voraussetzen.

## 6. Warum ein einziges holomorphes U_* nicht beide Bereiche tragen kann

**Satz 6 — bedingt, bewiesen.** Es gelte e_N→c>1. Auf einem festen zusammenhängenden Gebiet D, das c und ein reelles Intervall J⊂(−1,0) enthält, seien U_N holomorph. Auf J gelte F_N(x)→x+1. Dann können U_N nicht lokal gleichmäßig auf D gegen eine holomorphe Funktion U_* konvergieren.

Beweis: Für x∈J gilt

\[
U_N(x)=\frac{m^{x-1}}{B_m(x)}(F_N(x)-x-1)\to0,
\]

denn der erste Faktor konvergiert gegen Γ(x). Ein holomorpher Grenzwert müsste nach dem Identitätssatz auf ganz D null sein. Lokal gleichmäßige Konvergenz nahe c und e_N→c würden aber mit Abschnitt 5 U_*(c)=−Γ(c)≠0 liefern. Widerspruch. ∎

Dasselbe Argument gilt für J⊂(0,1). Der Satz widerlegt **unter der Successor-Annahme für die kritische Familie** ein ganzes oder gemeinsam holomorphes Universalprofil bis zum kritischen Punkt. Ein lokales Profil um c oder getrennte Profile mit einer singulären Übergangsschicht bleiben möglich. Falls die ursprüngliche Successor-Beobachtung nur Basis 2 betrifft, muss zunächst geklärt werden, ob sie für F_N überhaupt gilt.

## 7. Exakte Rangrekursion, vor jeder Linearisierung

Bei **fester Basis b** schreibe H_N(x)=x+1+d_N(x) und H_{N+1}(x)=x+1+d_{N+1}(x). Dann folgt exakt

\[
\boxed{d_{N+1}(x+1)-d_{N+1}(x)
=d_N(x+1+d_{N+1}(x)).}
\]

Dies ist die konkrete nichtlineare Gleichung hinter dem vorgeschlagenen Rang-Renormierungsoperator. Mit d_N=B_mA_N und d_{N+1}=B_{m+1}A_{N+1} wird sie, wo die Division zulässig ist,

\[
\frac{x+m+1}{m+1}A_{N+1}(x+1)
-\frac{x}{m+1}A_{N+1}(x)
=\frac{B_m(x+1+d_{N+1}(x))}{B_m(x+1)}
A_N(x+1+d_{N+1}(x)).
\]

Für die kritischen Profile ist rechts jedoch A_N(η_{N+1},·) erforderlich, nicht A_N(η_N,·). Deshalb ist eine Rekursion U_{N+1}=R_N(U_N) für diese Profile allein bisher **nicht geschlossen**. Man braucht die Basisfamilie oder kontrollierte Basissensitivitäten.

## 8. Lineares Modell: Gamma ist der richtige Faktor für Rangschritte

Setze formal f_ε(x)=x+1+εv(x)+O(ε²) und S_ε(x)=x+1+εw(x)+O(ε²), mit S_ε(x+1)=f_ε(S_ε(x)), S_ε(0)=1. Dann

\[
\boxed{w(x+1)-w(x)=v(x+1),\qquad w(0)=0.}
\]

Auf der Klasse der Polynome gibt es genau eine Lösung w: Zwei Lösungen unterscheiden sich durch ein periodisches Polynom, also eine Konstante, und die Normierung setzt diese auf null. Auf allgemeinen holomorphen Funktionen bleiben 1-periodische Summanden frei.

Für ν>0 definiere auf Re x>0

\[
B_\nu(x)=\frac{\Gamma(x+\nu)}{\Gamma(x)\Gamma(\nu+1)}.
\]

Die Gammafunktionalgleichung liefert exakt

\[
B_{\nu+1}(x+1)-B_{\nu+1}(x)=B_\nu(x+1).
\]

Somit sendet der normierte polynomielle lineare Rangschritt B_m nach B_{m+1}. Die bekannte Nullstellenstruktur wird hier aus der Rangrekursion erklärt; sie ist nicht bloß nachträglich herausdividiert.

**Satz 8 — im linearen Modell bewiesen.** Auf endlichen Linearkombinationen der B_ν, ν>0, definiert

\[
\mathcal L_sB_\nu=B_{\nu+s}\qquad(s\ge0)
\]

eine Halbgruppe mit L_s L_t=L_{s+t}. Insbesondere L_{1/2}²=L_1. Die Darstellung ist eindeutig: Für verschiedene ν_j sind die Funktionen B_{ν_j} linear unabhängig. Ihre Werte bei x=1,…,k sind Polynome in ν_j von aufsteigenden Graden 0,…,k−1; die zugehörige Vandermondematrix ist invertierbar.

Die Familie ist lokal analytisch in s und x, solange keine Gammapole betreten werden. Sie respektiert v(0)=0 und v(1)=0, weil B_ν(0)=0 und B_ν(1)=1. Bei nichtganzzahligem ν entstehen in der linken Halbebene neue Pole: Es handelt sich nicht um eine ganze Fortsetzung sämtlicher Höhen.

Beispiel:

\[
v(x)=x(x-1)=2B_2(x)-2B_1(x),
\]

\[
\mathcal L_{1/2}v(x)=2B_{5/2}(x)-2B_{3/2}(x),\qquad
\mathcal L_1v(x)=\frac{x(x-1)(x+1)}3.
\]

Hier sind zwei halbe Schritte mathematisch genau ein voller **linearer** Rangschritt. Es ist noch kein H_{N+1/2} der nichtlinearen, ausgewählten Hierarchie.

## 9. Derselbe halbe Schritt als fraktionelle Summation

Für Höhenfolgen v(0)=0 ist die Lösung bei ganzen k≥1

\[
(\mathcal L_1v)(k)=\sum_{j=1}^kv(j).
\]

Eine konkrete fraktionelle Fortsetzung ist

\[
(\mathcal L_sv)(k)=\sum_{j=1}^k
\frac{(s)_{k-j}}{(k-j)!}\,v(j).
\]

Die erzeugende Funktion des Kerns ist (1−z)^{−s}. Multiplikation der erzeugenden Funktionen beweist L_s L_t=L_{s+t}. Für v(k)=B_ν(k) lautet die erzeugende Funktion z(1−z)^{−ν−1}; die fraktionelle Summation ergibt z(1−z)^{−ν−s−1}, also genau B_{ν+s}(k).

Auf einer endlichen Höhengitterstrecke ist dies sogar eine exakte analytische Matrix-Halbgruppe. Mit dem rückwärts verschiebenden nilpotenten Matrixoperator E gilt

\[
\mathcal L_s=(I-E)^{-s}=\exp\left(s\sum_{j\ge1}\frac{E^j}{j}\right),
\]

wobei die Summe endlich ist. Diese Herleitung benötigt keine Literaturannahme über höherstufige Hyperoperationen.

## 10. Was ein echter halber Hyperoperationsrang zusätzlich braucht

Definiere nach Festlegung der analytischen Auswahl

\[
\mathcal T_b(f)=S,\qquad S(0)=1,\quad S(x+1)=f(S(x)).
\]

Die ganzzahlige Hierarchie erfüllt H_{N+1}(b,·)=T_b(H_N(b,·)). Ein konsistenter halber Rangschritt wäre eine Operator-Halbgruppe T_b^s mit T_b^1=T_b und

\[
\boxed{H_{N+s}(b,\cdot)=\mathcal T_b^s(H_N(b,\cdot)),\qquad
\mathcal T_b^{1/2}\circ\mathcal T_b^{1/2}=\mathcal T_b.}
\]

Dabei sind f und S **Funktionen als Zustände**, nicht einzelne Zahlenwerte. Die Formel σ^{-1}(λ^tσ(1)) für einen skalaren Fixpunkt konstruiert dagegen die fraktionelle **Höhe** H_{N+1}(b,t). Sie interpoliert den Rang nicht.

Wo S lokal invertierbar ist, kann die niedrigere Stufe explizit zurückgewonnen werden:

\[
\mathcal R(S)(y)=S(S^{-1}(y)+1)=f(y).
\]

Formal ist am Successor

\[
D\mathcal R[w](y)=w(y)-w(y-1)=\nabla w(y).
\]

Damit ist der lineare halbe **aufsteigende** Rangschritt genau eine ausgewählte fraktionelle inverse Rückwärtsdifferenz ∇^{−1/2}. Die periodische Auswahlfreiheit bleibt im nichtlinearen R erhalten: Für ψ(x+1)=ψ(x)+1 gilt R(S∘ψ)=R(S), soweit die inversen Zweige passen. Vor einer Einbettung muss deshalb eine geeignete Auswahlklasse beziehungsweise ein Quotient festgelegt werden.

**Hinreichendes Einbettungskriterium — bedingt, bewiesen.** Gibt es auf einer ausgewählten Zustandsklasse eine lokal biholomorphe Koordinate Φ, einen linearen Operator L mit gewähltem Logarithmus G und zulässigen Definitionsgebieten, so dass

\[
\Phi(\mathcal T_bf)=L\Phi(f),\qquad L=e^G,
\]

und e^{sG}Φ(f) im Koordinatengebiet bleibt, dann definiert

\[
\mathcal T_b^s(f)=\Phi^{-1}(e^{sG}\Phi(f))
\]

die gewünschte Halbgruppe. Beweis: Die Halbgruppeneigenschaft folgt unmittelbar aus e^{sG}e^{tG}=e^{(s+t)G}; bei s=1 gilt die gegebene Konjugation. ∎

**Offen sind hier:** die Zustandsräume und Domains der Komposition/Inversion; die periodische Auswahl; analytische Abhängigkeit von f und b; Kontrolle der nichtlinearen Terme; Existenz der Konjugation oder eines einbettenden Vektorfelds. Keiner dieser Punkte folgt aus einem Potenzfit von η_N.

Bereits die zweite Ordnung mischt die Newtonmoden. Für f_ε=x+1+εv ohne quadratischen Eingangsterm und S_ε=x+1+εw+ε²z+… gilt

\[
z(x+1)-z(x)=w(x)v'(x+1).
\]

Die Gamma-Verschiebung löst nur die erste Ordnung; dieser zusätzliche Term muss in einer nichtlinearen Konstruktion bearbeitet werden.

## 11. Warum das subkritische Tetrationsargument nicht direkt auf den Rang übergeht

Bei 1<b<exp(1/e) besitzt b^x einen reellen anziehenden Fixpunkt mit positivem Multiplikator kleiner 1. Die reguläre skalare Iteration gibt dadurch fraktionelle Höhen. Der Rangoperator T_b benötigt eine eigene Dynamik auf Funktionen.

**Lemma 11 — bewiesen.** Ein injektiver Fixpunkt f von T_b, dessen Definitionsbereich die benötigten Argumente enthält, muss f(x)=x+1 und b=2 erfüllen.

Beweis: Ein Fixpunkt erfüllt f(x+1)=f(f(x)). Injektivität gibt x+1=f(x). Dann f(1)=2=b. ∎

Für feste 1<b<2 gibt es deshalb keinen solchen regulären, injektiven Rangfixpunkt. Das schließt eine Ranginterpolation nicht aus, verhindert aber die unmittelbare Übernahme des einfachen anziehenden Fixpunktmechanismus der Tetration. Ein möglicher Ranggrenzwert kann auf größeren Höhenbereichen seine Analytizität oder Injektivität verlieren.

Die gelesenen Abschnitte in [Bracci–Contreras–Díaz-Madrigal, arXiv:1804.10465](https://arxiv.org/pdf/1804.10465), Einleitung und Satz 4.6 mit Beweis, verdeutlichen die skalare Analogie: Eine Koenigskoordinate erzeugt eine Halbgruppe erst zusammen mit passenden geometrischen Invarianzbedingungen ihres Bildgebietes. Das ist kein Existenzsatz für unseren Funktionaloperator.

## 12. Asymptotische Zeit und Übergangsschichten

Falls ein renormiertes Profil eine Entwicklung U_m=U_*+m^{−ω}V+… besitzt, wäre der formale Zwischenrang

\[
U_{m+s}^{\rm formal}=U_*+(m+s)^{-\omega}V+\cdots.
\]

Aus U_{m+1}−U_m≈m^{−1}G(U_m) entsteht die natürliche Zeit t=log m. Ein halber **Rang** entspricht daher

\[
\Delta t=\log\frac{m+1/2}{m}=\frac1{2m}+O(m^{-2}),
\]

und nicht einem festen halben Schritt in logarithmischer Zeit. Eine solche Näherung ist noch keine exakte Operator-Halbgruppe.

Der Gammafaktor selbst hat bei x=1 eine nachvollziehbare Übergangsskala:

\[
B_m\left(1+\frac{\xi}{\log m}\right)\longrightarrow e^\xi
\]

für beschränkte ξ. Im Modell eines langsam variierenden Restfaktors ist also 1/log m eine natürliche Höhenschichtbreite. **Das ist eine bewiesene Aussage über B_m, eine mögliche Messhypothese für H_N.** Ein tatsächlicher Restfaktor kann die Schicht verschieben oder ihre Skala verändern.

Schließlich ist eine gemeinsame Fehlerpotenz für η_N und e_N nicht automatisch. Hat lokal um c das Profil d_N die Form d_*+m^{−1}W_1+o(m^{−1}) in C² mit d_*(c)=−1, d_*′(c)=0 und d_*″(c)≠0, und gilt e_N→c, dann erzwingt die kritische Wertbedingung W_1(c)=0 und die Ableitungsgleichung

\[
e_N=c-\frac{W_1'(c)}{m\,d_*''(c)}+o(m^{-1}).
\]

Ein bei x=1 verschwindender 1/m-Modus kann somit η_N einen höheren Fehlerexponenten geben, während e_N einen 1/m-Term behält. Die entsprechenden Profilentwicklungen und Nichtentartung sind hier Voraussetzungen, nicht Ergebnisse.

## 13. Präzises Arbeitsergebnis

Bewiesen sind die korrigierte Ankerfaktorisierung unter Zweigvoraussetzungen, die exakte nichtlineare Ranggleichung, zwei hinreichende Konvergenzkriterien, die Einschränkung eines gemeinsamen holomorphen Profils und ein analytischer halber Rangschritt im ausgewählten linearen Modell.

**Nicht bewiesen:** η_N→2 und e_N→4 für die tatsächliche ausgewählte Hierarchie; die erforderlichen Restfaktor- und Basissensitivitätsschranken; ein nichtlinearer analytischer Hyperoperator H_{N+1/2}. Die Resultate reduzieren diese Aufgaben auf konkrete mathematische Aussagen und vermeiden, die Grenzwerte in ein universelles Profil hineinzudefinieren.
