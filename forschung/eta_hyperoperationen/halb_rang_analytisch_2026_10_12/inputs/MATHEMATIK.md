# Quadratischer halber Rang: bewiesene Teile und die Funktionsraumlücke

12.10.2026. **Kein Existenzbeweis für H₃.₅.** Gegenstand ist ein lokaler Rangschritt um den Successor und ein ausdrücklich deklariertes endliches Modell. Das Symbol K bezeichnet hier dessen Anzahl positiver Höhenknoten und hat nichts mit dem pausierten H₀-Projekt zu tun.

## 1. Ziel und Rangkonvention

Wir verwenden H₃(b,x)=bˣ, H₄=T_b, H₅=P_b und

\[
H_{r+1}(b,0)=1,\qquad H_{r+1}(b,x+1)=H_r(b,H_{r+1}(b,x)).
\]

Ein fraktionaler Rang soll über einen **ausgewählten Operator auf Funktionen** entstehen: 𝒯_s𝒯_t=𝒯_{s+t}, 𝒯₁(f) ist eine normierte Superfunktion von f. Erst dann wäre H₃.₅=𝒯₁/₂(H₃) begründet. Eine skalare halbe Iteration von bˣ interpoliert die Höhe und ergibt diese Rangoperation nicht.

Die Superfunktion ist ohne Auswahlregel nicht eindeutig. Deshalb gehört die Wahl der analytischen Klasse und der periodischen Freiheit zur Definition von 𝒯. Die folgenden endlichen Modelle treffen eine konkrete Wahl durch Polynominterpolation.

## 2. Exakte Ranggleichung und quadratischer Term — bewiesen

Setze f(x)=x+1+d(x), S(x)=x+1+e(x), S(0)=1. Dann

\[
e(x+1)-e(x)=d(x+1+e(x)).
\]

Für d=εv und e=εw+ε²z+O(ε³) folgt

\[
\Delta w(x)=v(x+1),\qquad
\Delta z(x)=w(x)v'(x+1),\qquad w(0)=z(0)=0.
\]

Bei polynomialem v und polynomialer Auswahl sind w,z eindeutig: Die Differenz zweier polynomialer Antidifferenzen ist konstant, die Normierung entfernt diese Konstante. Für allgemeine holomorphe Funktionen bleibt eine periodische Freiheit.

Für v(x)=x(x−1) erhält man genau

\[
w(x)=\frac{x(x-1)(x+1)}3,
\qquad
z(x)=\frac{x(x-1)(x-2)(x+1)(8x+1)}{60}.
\]

Einsetzen und Differenzieren verifiziert Δz=w(x)(2x+1); insbesondere z(3)=10.

## 3. Die endliche Gamma-Basis schliesst nicht — bewiesen

Sei B_ν(x)=Γ(x+ν)/(Γ(x)Γ(ν+1)), ν>0 nichtganzzahlig. Dann ΔB_{ν+1}(x)=B_ν(x+1). Für v=B_ν ist w=B_{ν+1} eine normierte lineare Lösung.

Bei x₀=−ν−1 besitzen w(x) und v(x+1) jeweils einen einfachen Pol mit nichtverschwindendem Residuum. Daher besitzt w(x)v′(x+1) bei x₀ einen **dreifachen Pol**. Jede endliche Linearkombination der B_μ und deren Differenzen besitzt höchstens einfache Pole. Ein z in einer endlichen solchen Gamma-Basis kann die quadratische Gleichung deshalb nicht erfüllen.

Der Vorschlag einer endlichen Gamma-Koeffizientenmatrix muss folglich eine Projektion oder eine erweiterte Basis enthalten. Parameterableitungen der Gamma-Moden und Produkte sind mögliche Erweiterungen; ihre genügende Abschluss- und Summierbarkeitseigenschaft ist hier nicht bewiesen.

## 4. Präzise Definition des endlichen Modells — bewiesen

Der Zustand ist c=(d(1),…,d(K)), mit d(0)=0. Er bestimmt das eindeutige Polynom p_c vom Grad höchstens K auf 0,…,K. Wir definieren

\[
f_c(x)=x+1+p_c(x),\quad S_c(0)=1,\quad S_c(k+1)=f_c(S_c(k)),
\]

und den vollen Modellschritt R_K(c)_k=S_c(k)−k−1, k=1,…,K. Dieser endliche Schritt ist exakt definiert und analytisch in c. Er bestimmt **keine** kanonische nichtganzzahlige Superfunktion; nur die integeren Höhenwerte werden hier erhoben.

Sei L die Matrix der Präfixsummen, W=L−I und D_K die Ableitungsmatrix der genannten Polynominterpolation. Dann

\[
R_K(c)=Lc+Q_1^K(c,c)+O(\|c\|^3),\qquad
Q_1^K(c,c)=L[(Wc)\odot(D_Kc)].
\]

Der Faktor W statt L ist wesentlich: In der Summe für z(k) steht w(j−1), nicht w(j). Polarisiert lautet die bilineare Form

\[
Q_1^K(u,v)=\tfrac12L[(Wu)\odot D_Kv+(Wv)\odot D_Ku].
\]

Der Rest ist für jedes feste K lokal kubisch. Eine K-unabhängige Restschranke wird nicht behauptet.

## 5. Expliziter halber Schritt bis zur zweiten Ordnung — bewiesen

Hebe c auf (c,q), q=c⊗c. Auf den linearen und quadratischen Koordinaten wirkt

\[
M_K=\begin{pmatrix}L&B_K\\0&L\otimes L\end{pmatrix},
\quad
B_K(q)=L\,\mathrm{diag}(WqD_K^T).
\]

Auf symmetrischen Tensoren entspricht dies genau dem quadratischen Jet des Modells. L−I ist nilpotent vom Index höchstens K; L⊗L−I hat Index höchstens 2K−1. Damit (M_K−I)^{3K−1}=0. Für beliebiges komplexes s ist daher

\[
M_K^s=\sum_{p=0}^{3K-2}{s\choose p}(M_K-I)^p
\]

eine endliche, wohldefinierte Summe. Sie erfüllt M_K^sM_K^t=M_K^{s+t}. Ihre obere Komponente definiert

\[
T_s^K(c)=L_sc+Q_s^K(c,c),
\quad
(L_sc)_k=\sum_{j=1}^k\frac{(s)_{k-j}}{(k-j)!}c_j.
\]

Somit, als Identität der Jets bis Grad zwei,

\[
T_s^K\circ T_t^K=T_{s+t}^K+O(\|c\|^3).
\]

Insbesondere zwei halbe Schritte ergeben den vollen Schritt modulo kubischer Terme. Für **festgelegtes** lineares H=L₁/₂ ist Q₁/₂ eindeutig: Auf bilinearen Formen hat die Abbildung Q↦HQ+Q(H·,H·) nur den Eigenwert 2 und ist invertierbar. Das ist eine endliche algebraische Aussage, keine uniforme Schranke für K→∞.

Die erste Zustandskoordinate bleibt erhalten: (T_s^Kc)₁=c₁. Da b=2+c₁ ist, bleibt die Basis in diesem Sinn erhalten. Auf c₁=0 bleibt auch c₂ erhalten. Diese Konsistenzbedingungen ersetzen keine analytische Auswahlregel.

## 6. Analytische lokale Modellhalbgruppe und formale volle Einbettung

**Bewiesen für jedes feste K:** G_K=log M_K ist eine endliche Summe. Auf Jets ist es eine Derivation: M_K^s respektiert die Multiplikation modulo Grad drei für ganzzahlige s; wegen der polynomialen Abhängigkeit gilt dies für alle s, und Differenzieren bei s=0 liefert die Leibnizregel.

Die obere Komponente von G_K entspricht deshalb einem polynomialen Vektorfeld

\[
V_K(c)=(\log L)c+V_{K,2}(c,c).
\]

Dessen lokale analytische Flüsse Φ_s^K erfüllen die Halbgruppengleichung exakt, und ihr quadratischer Jet ist T_s^K. Ihre Zeit-eins-Abbildung stimmt mit R_K **bis Grad zwei** überein. Der volle Modellschritt R_K wird damit noch nicht analytisch eingebettet.

**Literaturgestützte formale Aussage:** Weil DR_K(0)=L unipotent ist, liefert das formale Exponential-/Logarithmusprinzip einen eindeutigen formalen nilpotenten Generator des vollständigen R_K. Damit existieren konsistente formale Jets aller Ordnungen. Konvergenz dieser Reihe und ein Grenzübergang K→∞ folgen daraus nicht. Siehe Ribón, §2.1, Lemma 2.2; die Abgrenzung formal/analytisch steht bereits in seiner Einleitung.

## 7. Geschlossene quadratische Rangformel bei x=3 — bewiesen

Für v(x)=x(x−1) ist der lineare Schritt nach ganzzahligem n

\[
a_n(x)=L_nv(x)=2B_{n+2}(x)-2B_{n+1}(x)
=\frac{2(x-1)\Gamma(x+n+1)}{\Gamma(x)\Gamma(n+3)}.
\]

Es gelten a_n(1)=0, a_n(2)=2, a_n(3)=2(n+3) und

\[
a_n'(3)=2(n+3)(\mathsf H_{n+3}-1),
\]

wobei 𝖧_j die harmonische Zahl bezeichnet. Da die quadratischen Koordinaten bei Höhe 1 und 2 verschwinden, liefert die Jet-Komposition

\[
Q_n(v,v)(3)=4\sum_{j=0}^{n-1}(j+3)(\mathsf H_{j+3}-1).
\]

Die endliche Summe ist

\[
F(n)=2(n+2)(n+3)\mathsf H_{n+2}-(n+2)(3n+7)-4.
\]

Für komplexes s setzen wir die entsprechende explizite Funktion

\[
F(s)=2(s+2)(s+3)[\psi(s+3)+\gamma]-(s+2)(3s+7)-4.
\]

Diese Formel wird im nächsten Abschnitt als der Grenzwert **der hier festgelegten endlichen Modelle** hergeleitet. Sie soll keine Eindeutigkeit unter allen denkbaren analytischen Ranginterpolationen behaupten.

## 8. Exakte Newtonreihe und Modellgrenzen — bewiesen

Die Vorwärtsdifferenzen α_p=Δ^pF(0) lauten

\[
\alpha_0=0,\quad\alpha_1=10,\quad\alpha_2=22/3,\quad
\alpha_p=\frac{24(-1)^{p+1}}{(p-2)(p-1)p(p+1)}\quad(p\ge3).
\]

**Herleitung:** Schreibe 𝖧_{n+2}=∫₀¹(1−t^{n+2})/(1−t)dt. Für p≥3 verschwinden die polynomialen Beiträge der zweiten Differenzordnung. Wende ((t∂_t)²+5t∂_t+6) auf (t−1)^p an. Die drei verbleibenden Betaintegrale ergeben genau α_p.

Die Newtonreihe ist

\[
F(s)=\sum_{p=0}^{\infty}{s\choose p}\alpha_p.
\]

Zur Identifikation kann man den Rest ab p=3 auch als

\[
-4\int_0^1\frac{(1-t)^3}{t^3}
\left[(1-t)^s-1+st-\tfrac12s(s-1)t^2\right]dt
\]

schreiben. Explizit ist dieses Integral

\[
-4\lim_{\epsilon\to0}\left[
B(\epsilon-2,s+4)-B(\epsilon-2,4)
+sB(\epsilon-1,4)-\tfrac12s(s-1)B(\epsilon,4)\right].
\]

Die 1/ε-Terme heben sich auf. Nach Addition von 10s+(22/3)·binom(s,2) ergibt der konstante Term

\[
2(s+2)(s+3)[\psi(s+2)+\gamma]-3s^2-11s-12=F(s),
\]

wobei die letzte Gleichheit die Digamma-Verschiebung ψ(s+3)=ψ(s+2)+1/(s+2) benutzt. Am linken Integrationsrand entfernt die Klammer die ersten drei Taylorglieder. Für Re s>−4 ist das Integral regulär. Die Newtonreihe konvergiert dort lokal gleichmässig, da α_p=O(p⁻⁴) und binom(s,p)=O(p^{−Re s−1}) auf geeigneten kompakten Mengen; der scheinbare Digammapol bei s=−3 ist in F entfernbar.

**Warum das endliche Modell exakt die partielle Reihe liefert:** Für integeres k≥2 ist a_n(k)=2(k−1)(n+3)_{k−2}/Γ(k) polynomial in n vom Grad k−2. Daher ist D_Ka_n(3) polynomial in n vom Grad höchstens K−2. Der quadratische Koeffizient der n-ten Modelliteration bei Höhe 3 ist polynomial in n vom Grad höchstens K−1. Für n=0,…,K−1 stimmt er mit F(n) überein: In jedem benötigten vorherigen Schritt besitzt a_j als Höhenpolynom Grad j+2≤K und wird exakt differenziert. Folglich für alle s

\[
\boxed{Q_s^K(v,v)(3)=\sum_{p=0}^{K-1}{s\choose p}\alpha_p.}
\]

Die Carleman-Potenz ist ebenfalls polynomial in s und stimmt an allen integeren n mit diesem Modelljet überein; dies identifiziert beide Polynome. Der Grenzübergang ist hier somit analytisch begründet, nicht nur numerisch beobachtet.

## 9. Halbschritt, Generator und expliziter Fehler — bewiesen

Aus der geschlossenen Formel folgt

\[
\boxed{F(1/2)=341/12-35\log2
=4.15651534706858083706354241563048678\ldots}
\]

und der quadratische Generatorbeitrag bei dieser Eingaberichtung und Höhe

\[
\boxed{F'(0)=2\pi^2-13=6.739208802178717237\ldots .}
\]

Für s=1/2 sind alle Reihenterme ab p=3 positiv. Mit

\[
\left|{1/2\choose p}\right|
=\frac1{2p}\frac{{2p-2\choose p-1}}{4^{p-1}}
\le\frac1{2p\sqrt{\pi(p-1)}}
\]

und der Formel für α_p ergibt sich für K≥4

\[
\boxed{0<F(1/2)-Q_{1/2}^K(v,v)(3)
\le\frac8{3\sqrt\pi}(K-3)^{-9/2}.}
\]

Denn ein Reihenterm ist höchstens (12/√π)(p−2)^{−11/2}; die Summe ab p=K wird durch das Integral ab K−3 majorisiert. Auch die erste Asymptotik ist explizit:

\[
F(1/2)-Q_{1/2}^K(v,v)(3)
\sim\frac8{3\sqrt\pi}K^{-9/2}.
\]

Das ist ein kontrollierter **einzelner quadratischer Koeffizient**. Die Aussage darf nicht auf alle Höhen oder den vollen nichtlinearen Rangfluss ausgedehnt werden.

## 10. Drift bei nichtganzzahliger Höhe — numerisch

Nach dem Modellschritt wurden dessen K Werte wieder durch das Grad-K-Polynom auf 0,…,K interpoliert. Für dessen quadratischen Beitrag bei x=0.5:

| K | Q₁/₂, interpoliert bei x=0.5 |
|---:|---:|
| 16 | +0.007002656994 |
| 32 | −0.007670584698 |
| 48 | −0.018798624783 |
| 64 | −0.028087148559 |
| 80 | −0.036218850499 |

Auf K=16,…,80 beschreibt 0.04189316−0.00874521√K die Werte mit maximalem Rest 1.08·10⁻⁴. **Das ist kein bewiesenes Divergenzgesetz.** Bei x=1.5 und 2.5 ist die beobachtete Drift wesentlich schwächer. Der lineare Beitrag nähert sich dagegen an allen drei getesteten nichtintegeren Höhen der expliziten Gammaformel; bei x=0.5 ist sein Fehler für K=80 etwa −1.267·10⁻⁵.

Die Drift bei K=64 bleibt beim Wechsel 150→250 Dezimalstellen auf alle 65 gespeicherten signifikanten Stellen gleich. Es ist eine Modell-/Fortsetzungsfrage, keine erkennbare Rundungsfrage. Ein lokal gleichmässiger holomorpher Grenzwert dieses quadratischen Modells ist damit jedenfalls nicht durch die integeren Messungen belegt.

## 11. Ein rigoroses Problem der gewichteten Sequenzräume — bewiesen

Auf Zuständen sei ∥c∥_a=max_{1≤j≤K}|c_j|a^{−j}, a>1. Der lineare Halbschritt ist darin uniform beschränkt:

\[
\|L_{1/2}\|_a\le(1-a^{-1})^{-1/2},\qquad
\|L_{1/2}-I\|_a\le(1-a^{-1})^{-1/2}-1.
\]

Die Polynomableitung besitzt dagegen schon den Matrixeintrag |D_{1K}|=1/[K(K−1)]. Somit

\[
\|D_K\|_a\ge\frac{a^{K-1}}{K(K-1)}\longrightarrow\infty.
\]

Sogar auf dem Unterraum c₁=0, der die Basis 2 festhält, ist der quadratische volle Schritt unbeschränkt in K. Wähle u₂=a², v_K=a^K und alle anderen Komponenten null, K≥4; beide Eingaben haben Norm eins. Da |D_{3K}|=6/[K(K−1)(K−2)(K−3)], folgt

\[
\|Q_1^K\|_a\ge
\frac{3a^{K-1}}{K(K-1)(K-2)(K-3)}\longrightarrow\infty.
\]

Eine uniforme Halb-Rang-Konstruktion auf **allen** solchen gewichteten Sequenzen kann daher nicht durch diese Interpolationsmodelle begründet werden. Eine kleinere Klasse mit kontrollierter analytischer Fortsetzung könnte diese Gegenrichtungen ausschliessen. Ihre Konstruktion bleibt offen.

Für a=2 wird sogar der quadratische Halbschritt selbst unbeschränkt: Aus Q₁=HQ₁/₂+Q₁/₂(H·,H·) und ∥H∥₂≤√2 folgt ∥Q₁∥₂≤(√2+2)∥Q₁/₂∥₂. Somit

\[
\boxed{\|Q_{1/2}^K\|_2\ge
\frac{3\,2^{K-1}}{(2+\sqrt2)K(K-1)(K-2)(K-3)}\to\infty.}
\]

Das schliesst einen beschränkten Grenzoperator dieser Modelle in der genannten Operatornorm aus. Es ist kein Nichtexistenzbeweis auf einer engeren analytischen Funktionenklasse.

## 12. Bedingtes Funktionsraumlemma für die zweite Ordnung

**Bewiesen unter den genannten Annahmen:** Sei X ein Banachraum, H:X→X beschränkt und Q₁ eine beschränkte bilineare Abbildung X×X→X. Auf dem Banachraum der beschränkten bilinearen Abbildungen definiere

\[
\mathscr C Q=H\circ Q+Q\circ(H,H).
\]

Falls ∥𝒞−2I∥<2, existiert genau eine beschränkte Lösung Q₁/₂ von 𝒞Q₁/₂=Q₁, gegeben durch die Neumannreihe

\[
Q_{1/2}=\tfrac12\sum_{j\ge0}
[-(\mathscr C-2I)/2]^jQ_1.
\]

Insbesondere reicht 3α+α²<2 mit α=∥H−I∥. Dann ∥Q₁/₂∥≤∥Q₁∥/(2−3α−α²).

Für das gewichtete Sequenzmodell mit a=2 ist die lineare Bedingung erfüllt: α≤√2−1. Die quadratische Uniformität aus §11 scheitert trotzdem. Das zeigt genau, welche Voraussetzung eine gewählte analytische Klasse zusätzlich leisten muss. Dieses Lemma wurde nicht auf die tatsächlichen Hyperoperationen angewendet.

## 13. Warum hier keine kritische Basis des Halb-Rangs folgt — bewiesen

Eine kritische Tangenz des gesuchten F_b(x)=H₃.₅(b,x) würde

\[
F_\eta(e_*)=e_*,\qquad\partial_xF_\eta(e_*)=1
\]

erfüllen. Für d=F−x−1 gilt dort **d(e_*)=−1**. Eine uniformly kleine Successor-Störung mit |d|<1 auf dem betreffenden realen Gebiet hat überhaupt keinen Fixpunkt. Die lokale Entwicklung darf deshalb nicht ohne Restkontrolle mit ε=1 bis an die Tangenz fortgesetzt werden.

Insbesondere ist b=e keine kleine Exponential-Störung des Successors: d(1)=e−2≈0.7183 und d(2)=e²−3≈4.3891. Auch bei b=2 ist d(2)=1. Ein kleiner formaler Rangjet liefert H₃.₅,e daher nicht direkt.

Die Existenz eines Halb-Rangs ist ausserdem logisch von der Vermutung η_N→2 getrennt. In dieser Runde gibt es keinen neuen Beweis jener Grenzwertvermutung.

## 14. Der Fit und eine ausdrücklich ungeeignete Vergleichsinterpolation

Der früher ausschliesslich an den Rängen 3,…,6 angepasste Fit bleibt unverändert:

\[
\eta_{\rm fit}(r)=2-2.7849974381025855\,r^{-1.467173061332697}.
\]

Er gibt η_fit(3.5)=1.5568172767405143… . Das ist eine zu prüfende Vorhersage, keine Definition von H₃.₅.

Für F_θ=(1−θ)bˣ+θT_b(x) wurden lokale positive Tangenzpunkte numerisch gefunden:

| θ | Tangenzbasis | Tangenzhöhe |
|---:|---:|---:|
| 0.25 | 1.49124063769114 | 2.86463618944987 |
| 0.50 | 1.54254852987747 | 3.03219227221562 |
| 0.75 | 1.59296534842843 | 3.13914851739044 |

Für θ=0.5 unterscheiden sich 25/45-Digit-Läufe mit zugleich 48/80 Cauchypunkten und Radius 0.04/0.025 um 3.12·10⁻²⁹ in der Basis und 4.01·10⁻²⁸ in der Höhe. Das sind numerische Kontrollen, keine Intervallzertifikate. Die Differenz zur Fitbasis beträgt −0.01426874686, relativ −0.9165 %.

Diese Mischkurve ist **kein nachgewiesener fraktionaler Rang**. Die natürliche gesamte Mischfamilie verletzt sogar die erforderliche Rangrekursion: Für G_θ=(1−θ)T_b+θP_b ist bei b=2, θ=0.5, x=0.5

\[
G_\theta(x+1)-F_\theta(G_\theta(x))
\approx-0.0004591225955504457431.
\]

Bei b=e beträgt der entsprechende Rest −0.02346564845418… . Die beiden integeren Endpunkte θ=0,1 bestehen dieselbe Kontrolle mit Rest unter 2·10⁻³⁷. Damit wird konkret sichtbar, warum eine passende Zwischenkurve und eine plausible kritische Basis noch keine Rangdynamik konstruieren.

## 15. Präziser nächster Beweisbedarf

1. **Analytische Klasse und Auswahl:** Eine Klasse festlegen, auf der die normierte Superfunktion eindeutig ausgewählt wird. Die nötige Differentiation, Verschiebung und nichtlineare Komposition müssen kontrolliert werden. Eine endliche Gamma-Basis oder beliebige gewichtete Sequenzen reichen nach §3 und §11 nicht.
2. **Quadratische Fortsetzung als Funktion:** Q₁ und die kohomologische Gleichung aus §12 auf dieser Klasse lösen. Die bekannten integeren Koeffizienten sind Kontrollen. Zusätzlich auf einem festen komplexen Gebiet testen, das mindestens die relevanten nichtintegeren Höhen einschliesst; die Drift bei x=0.5 muss aufgelöst werden.
3. **Alle Ordnungen:** Einen konvergenten Generator mit quantitativen Koeffizienten-/Restschranken gewinnen. Formale Logarithmen allein reichen nicht. Eine Folge endlicher Jets muss sowohl die analytische Fortsetzung als auch den ausgewählten vollen Rangschritt treffen.
4. **Fortsetzung zum Exponentialprofil:** Einen kontrollierten Pfad vom lokalen Successor-Modell zum tatsächlichen f_b(x)=bˣ und speziell b=e etablieren. Auf diesem Pfad dürfen Existenz und Auswahl nicht stillschweigend wechseln.
5. **Erst danach kritische Basis:** Eine reelle Familie H₃.₅(b,x) definieren, Tangenz und Transversalität prüfen, dann η₃.₅ mit dem unveränderten Fit vergleichen. Für ein approximatives F braucht man C²-Kontrolle und eine positive untere Schranke für ∂_bF sowie ∂²_xF an der Tangenz.

**Gesamtstatus:** Quadratischer halber Schritt im endlichen Modell bewiesen; ein skalarer Grenzkoeffizient und dessen Fehlerrate bewiesen; Fortsetzung als Funktion und tatsächliches H₃.₅ offen. Die analytische Klasse ist jetzt ein konkreter Engpass, kein nur allgemein benannter Rest.
