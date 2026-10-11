# Runde 3: globaler Transport und ein Summenkriterium für (H₀)

**Strang G (GPT), 11.10.2026. Ausschließlich (H₀); K2 bleibt bei P.**

**Ergebnis:** Die kubische Vierzustandsrekursion ist auf geschnittenen festen Konturen exakt diagonalisiert; die notwendige zyklische Permutation der drei Umkehrzweige ist bestimmt. Für die inverse Momentaufgabe ist noch keine globale Parametrix mit polynomialem Rest konstruiert. Stattdessen ergibt die bestehende Biorthogonalität eine exakte, quellenhaltige Inversenentwicklung als Summe von Grad-Zuwächsen. Eine Schranke **C/√n für diese Zuwächse auf γ*** würde bereits **C′√N für die volle Majorante G_N** liefern. Die neuen Rechnungen stützen diese stärkere Aussage und zertifizieren die Zuwächse für sechs einzelne Grade auf vollständigen Eckkreisen. **Die All-n-Schranke, damit (H₀), und der Kneser-Übertragungssatz bleiben offen.**

Die Sätze G3–G7 und die ausführlichen Korrekturen I_f − B_y stehen im archivierten `RUNDE2_CAUCHY.md`. Runde 1 bleibt unverändert archiviert. Es werden keine zusätzlichen Freiheitsgrade, keine Grad-5N-Konstruktion und keine neuen K2-Experimente verwendet. Unabhängige P-Prüfung der neuen Herleitungen steht aus.

**Projektstatus:** „Vollständige eigene Herleitung“ bezeichnet die unten ausgeschriebenen Beweise, „Intervallprüfung“ die eigenen geprüften Rechnungen. Beides ist noch nicht von P unabhängig geprüft und wird nicht als P-verifiziertes Fundament des Gesamtbeweises ausgegeben.

## 1. Die genaue endliche Aufgabe

Wir untersuchen die quadratischen Jacobiabschnitte J_N des ursprünglichen Systems, mit N+1 Koeffizienten und N+1 Gleichungen. f ist dabei der feste Modellkeim; es wird keine neue nichtlineare Lösung p_N konstruiert. Wie zuvor ist

\[
\widetilde J_N=D_NJ_ND_N^{-1},\quad D_N=\operatorname{diag}(\rho^k),\quad
e_{N,x}=((\rho/x)^k)_{k=0}^N,\quad
v_N(x)=\widetilde J_N^{-1}e_{N,x}.
\]

Dabei ist **G_N(x)=‖v_N(x)‖₂** genau die SVD-freie Majorante von O §33.1. Für \(\widetilde J_N^Ty=u\) gilt \(Y_u(x)=u^Tv_N(x)\); folglich deckt G_N alle Einheits-u und alle Eckmoden zugleich ab. Die Transposition ist bilinear, die Norm hermitesch. Der Funktionenraum und die ρ-Norm werden nicht geändert.

Die Vorwärtsinterpretation bleibt

\[
T_N\{q_N\circ f+m q_N\}=T_N\frac{x}{x-z},\qquad
q_N(z)=\sum_{j=0}^N v_{N,j}(x)(z/\rho)^j.
\]

Insbesondere bleibt die Quelle vollständig im System. Die endliche Cauchy-Fassung aus Runde 2 ist \(\Pi_N^-\mathcal C_\Gamma Y=U\), mit \(\mathcal C_\Gamma Y=I_f-B_y\) auf dem inneren Gebiet. Die folgende Algebra ist eine Faktorisierung genau dieser Momentaufgabe; sie ersetzt weder den Projektor noch die Korrektur durch eine Näherung.

## 2. Was die vier Zustände global liefern

Für das kubische Modell \(f(z)=a+bz+cz^2+dz^3\) war in Runde 2 die exakte Transfermatrix

\[
\mathsf A(t)=
\begin{pmatrix}
b/(t-a)&c/(t-a)&d/(t-a)&0\\
1&0&0&0\\
0&1&0&0\\
0&0&0&1/t
\end{pmatrix}.
\]

Sie transportiert die **Testkerne**, nicht bereits die inverse Momentlösung. Ihre Ausgabereihe ist \((1,0,0,m(t))\). Das m(t) dieser Testkerne darf weiterhin nicht mit dem m(z) des gemeinsamen Cauchy-Integrals verwechselt werden.

### Satz G8 — exakte Diagonalisierung und Permutationssprung [vollständige eigene Herleitung; P-Prüfung offen]

Seien r₁(t), r₂(t), r₃(t) die drei einfachen Urbilder von t, \(\lambda_i=1/r_i\). Abseits von t=0,a und den kritischen Werten gilt auf jedem Gebiet mit einer analytischen Wahl der Zweige

\[
\mathsf A=V\operatorname{diag}(\lambda_1,\lambda_2,\lambda_3,1/t)V^{-1},
\qquad
V=\begin{pmatrix}
\lambda_1^2&\lambda_2^2&\lambda_3^2&0\\
\lambda_1&\lambda_2&\lambda_3&0\\
1&1&1&0\\
0&0&0&1
\end{pmatrix}.
\]

**Beweis:** \((t-a)\lambda_i^3-b\lambda_i^2-c\lambda_i-d=0\), weil f(1/λ_i)=t. Direktes Multiplizieren gibt die vier Eigenvektoren. Ferner
\(\det V=(\lambda_1-\lambda_2)(\lambda_1-\lambda_3)(\lambda_2-\lambda_3)\ne0\). ∎

Damit ist der **Transport** auf dem geschnittenen Gebiet exakt durch \(V\operatorname{diag}(\lambda_i^k,t^{-k})V^{-1}\) dargestellt, ohne asymptotischen Rest. Beim Weiterführen um eine Kontur, die beide kritischen Werte umschließt, permutieren die ersten drei Spalten zyklisch. Die vierte bleibt fest. In Matrixsprache: \(V_+=V_-P\), wobei P die entsprechende 3-Zykel-Permutation, erweitert um eine 1, ist. Ein globaler Ansatz muss diesen Sprung mitführen.

Für das **exakte** kubische Modell aus L=−W₋₁(−1) und μ=√L liegen die kritischen Werte rigoros in

\[
-0.399467528583751232394043810779\ \pm\ 3.5\,10^{-31},\qquad
-0.475120466990655645030804546668\ \pm\ 1.9\,10^{-31}.
\]

Beide haben Betrag <0,476<ρ. Auf dem festen Kreis |t|=ρ ergibt der Umlauf deshalb dieselbe Monodromie wie bei ∞: Für d z³≈t ist sie ein 3-Zykel. Die numerische Wurzelverfolgung mit 721 Punkten bestätigt (in ihrer Startlabelung) die Permutation (1,2,0). Die analytische Aussage beruht auf den zertifizierten kritischen Werten und der Deformation nach außen, nicht auf dem Raster.

Eine konkrete Schnittwahl ist der negative Halbstrahl bis zum größeren kritischen Wert, also (−∞,−0,39946752858…]. Auf seinem Komplement sind die drei Urbilder analytisch wählbar; ihre Kehrwerte sind in t=a meromorph. Der Kreis |t|=ρ kreuzt diesen Schnitt bei −ρ und trägt dort den beschriebenen 3-Zykel. Die Pole t=0,a bleiben in der Matrixdarstellung ausdrücklich erhalten.

**Die Wiederholung eines Eigenwerts bei L ist unschädlich für diese Basis.** Dort ist r₁(L)=L, also λ₁(L)=1/L, gleich dem vierten Eigenwert. Die anderen beiden λ bleiben verschieden; der kleinste Abstand der drei kubischen Eigenwerte ist rigoros >0,15516. V bleibt invertierbar. Die Matrix hat in L keinen dadurch erzeugten Jordanblock. Der beobachtete √N-Vorfaktor folgt somit nicht aus einer Singularität dieser diagonalisierenden Basis.

**Grenze:** Diese Darstellung führt die Zweige und ihre Sprünge mit. Sie bestimmt noch nicht die Amplituden, welche die N+1 Momentbedingungen erzwingen. Deshalb wäre es falsch, sie bereits als eine Parametrix für \(\widetilde J_N^{-1}\) oder als Beweis von (H₀) zu bezeichnen. Die freie Grenzgeometrie von W wird dadurch ebenfalls nicht zu einer exakten klassischen Orthogonalität in der Ebene.

Eine einfache quantitative Kontrolle zeigt dieselbe Grenze: Am inneren Kreispunkt \(t_n=L(1-n^{-\beta})\) besitzt \((\rho\mathsf A(t_n))^n\) den vierten Eigenwert \((\rho/t_n)^n\), dessen Betrag
\((1-n^{-\beta})^{-n}\ge\exp(n^{1-\beta})\) ist. Für β<1 übertrifft das jede feste Potenz von n. Eine Normabschätzung der reinen Transferpotenzen kann daher die beobachtete polynomiale Schranke nicht liefern. Es braucht die globalen Momentamplituden und ihre Auslöschung.

## 3. Exakte Inversenentwicklung über die Grad-Erweiterung

Hier wird W §25.1(c) ausdrücklich verwendet, nicht neu beansprucht. Neu untersucht wird die Norm des einzelnen Zuwachses auf γ* und auf den Eckkreisen.

### Satz G9 — Rang-eins-Zuwachs mit voller Quelle [vollständige eigene Herleitung; P-Prüfung offen]

Da f und m für alle Grade dieselben Keime sind, ist der linke obere Block von \(\widetilde J_n\) genau \(\widetilde J_{n-1}\). Schreibe

\[
\widetilde J_n=\begin{pmatrix}A&b\\r^T&d_n\end{pmatrix},\quad
A=\widetilde J_{n-1},\quad
s_n=d_n-r^TA^{-1}b,\quad
Q_n=\binom{-A^{-1}b}{1}.
\]

Falls A und \(\widetilde J_n\) invertierbar sind, gilt s_n≠0. Für jede verschachtelte rechte Seite \(g^{(n)}=(g_0,\dots,g_n)^T\) ist

\[
\widetilde J_n^{-1}g^{(n)}
=\iota_{n-1}\widetilde J_{n-1}^{-1}g^{(n-1)}
+Q_n\frac{g_n-r^TA^{-1}g^{(n-1)}}{s_n}.
\tag{G15}
\]

Dabei hängt der Zähler ausdrücklich von der gesamten alten Quelle ab. Er wird nicht zu g_n verkürzt. Die oberen n Gleichungen und die letzte Schur-Gleichung beweisen (G15) direkt. ∎

Die Verbindung zur Biorthogonalität lässt sich ohne neue Normierungskonvention ausdrücken:

\[
Q_n=\frac{\widetilde J_n^{-1}e_n}{p_n},\quad
A_n=\frac{\widetilde J_n^{-T}e_n}{p_n},\quad
p_n=(\widetilde J_n^{-1})_{nn}=s_n^{-1},\quad h_n=s_n.
\]

Beide Vektoren haben letzte Komponente 1. Sie erfüllen \((\widetilde J_nQ_n)_k=0\) für k<n, \((A_n^T\widetilde J_n)_j=0\) für j<n und \(A_n^T\widetilde J_nQ_n=h_n\). Das sind die monischen biorthogonalen Vektoren in der ρ-skalierten Basis. Für die Cauchy-Quelle ergibt sich exakt

\[
v_n(x)-\iota_{n-1}v_{n-1}(x)
=Q_n\frac{\mathcal A_n(x)}{h_n},\quad
\mathcal A_n(x)=\sum_{k=0}^n A_{n,k}(\rho/x)^k.
\tag{G16}
\]

Definiere die eine skalare Laurentfunktion

\[
\mathcal D_n(x):=\frac{\|Q_n\|_2}{h_n}\mathcal A_n(x),\qquad
\|v_n(x)-\iota_{n-1}v_{n-1}(x)\|_2=|\mathcal D_n(x)|.
\tag{G17}
\]

Sie ist keine zusätzliche Unbekannte des nichtlinearen Polynomproblems. Sie ist ein exaktes Objekt des ursprünglichen linearen Abschnitts. Ihre Messung braucht zwei Randobjekte und keine SVD: die letzte Spalte und die letzte Zeile der inversen Matrix.

Durch Iteration folgt die exakte finite Summe

\[
v_N(x)=\sum_{n=0}^N\iota_{n,N}Q_n\frac{\mathcal A_n(x)}{h_n}.
\tag{G18}
\]

Das ist W's biorthogonale Inversenentwicklung, hier mit der Cauchy-Quelle ausgewertet. Jede Stufe behält genau ihre Momentbedingungen. Die globale Amplitude steckt weiterhin in h_n und \(\mathcal A_n\); sie wird nicht als lokaler Eingang angenommen.

## 4. Eine hinreichende All-n-Abschätzung

### Satz G10 — Grad-Zuwächse umgehen eine vollständige RH-Parametrix [vollständige eigene Herleitung, bedingt; P-Prüfung offen]

Sei γ* eine feste kompakte Kurve ohne 0. Alle benötigten Hauptabschnitte seien invertierbar. Für ein n₀≥2 seien

\[
C_0=\sup_{x\in\gamma^*}G_{n_0-1}(x)<\infty,
\qquad
\sup_{x\in\gamma^*}|\mathcal D_n(x)|\le C_\Delta n^{-1/2}
\quad\text{für alle }n\ge n_0.
\tag{G19}
\]

Dann gilt für alle N≥n₀−1

\[
\sup_{\gamma^*}G_N
\le C_0+2C_\Delta(\sqrt N-\sqrt{n_0-1})
\le C_G\sqrt N,\quad
C_G=\max\{2C_\Delta,C_0/\sqrt{n_0-1}\}.
\tag{G20}
\]

**Beweis:** Dreiecksungleichung in (G18) und
\(\sum_{n=n_0}^Nn^{-1/2}\le\int_{n_0-1}^Nt^{-1/2}dt\). Die letzte Form folgt durch Division durch √N. ∎

Diese Hypothese ist **stärker** als das gewünschte G_N=O(√N): Aus einem G_N-Wert allein folgt kein n⁻¹ᐟ²-Zuwachs. Sie könnte auch durch einzelne Grad-Ausreißer scheitern. Die sechs unten geprüften Grade sind daher ein Test dieser Hypothese, kein Nachweis ihrer Gültigkeit für alle n.

Für die bloß polynomiale Version von (H₀) genügt schon die schwächere All-n-Hypothese \(\sup_{\gamma^*}|\mathcal D_n|\le C_\Delta\): Dann liefert dieselbe Summe \(\sup G_N\le C_0+C_\Delta(N-n_0+1)\), also κ=1. Der beobachtete n⁻¹ᐟ²-Abfall ist ein scharfes Prüfziel, keine notwendige Voraussetzung jedes Weges über die Grad-Erweiterung.

**Vollständige Eckkreise funktionieren im selben globalen Summenargument.** Sei r_n=ρ n⁻ᵝ mit β>0. Gilt \(|\mathcal D_n|\le C_\Delta/\sqrt n\) auf |x−L|=r_n, gilt es durch das Maximumprinzip auch in der zugehörigen Scheibe (sie enthält 0 nicht). Für x mit |x−L|≤r_N liegt x in allen Scheiben mit n≤N. Mit C₀ als Maximum von G_{n₀−1} auf der n₀-Scheibe beweist dieselbe Summe (G20) auf der ganzen N-Scheibe. Die globale Normierung der Interpolanten bleibt dabei Teil von \(\mathcal D_n\). Es ist kein quellenfreies lokales Argument.

**Explizite Konstanten:** Eine künftig bewiesene C_Δ=3,13 ab n₀=200 würde in diesem Kreisargument C_G=max{6,26,C₀/√199} ergeben. Das ist eine konditionale Konstante, nicht das heute erreichte 0,84 für G auf dem KB-Raster. Für die feste γ*-Kurve wäre entsprechend C_Δ≈1,96 ein numerisches Ziel, mit C_G=max{3,92,C₀/√199}. Die Kreisbedingung allein ersetzt nicht die Außenkontrolle auf der restlichen festen γ*.

Die Schnittstelle zu X bleibt W §23.6: Unter einer bewiesenen vollen γ*-Schranke und zertifizierter Geometrie würde
\(C_E=4,95 C_G+0,728\sqrt{1+1/N_0}\), κ_E=1/2 genügen. **4,95 und 0,728 sind dort numerische Geometriewerte**, keine neuen Zertifikate. Hieraus wird weder N₀ für den Gesamtbeweis noch Kneser-Auswahl behauptet.

## 5. Rigorose Kreis-Majorante für einen festen Grad

### Satz G11 — Kreisabbildung und endliche Taylor-Majorante [vollständige eigene Herleitung; P-Prüfung offen]

Setze w=ρ/x. Für r<|L|=ρ wird |x−L|=r auf den Kreis

\[
|w-c|=R,\qquad c=\frac{\rho\bar L}{\rho^2-r^2},\qquad
R=\frac{\rho r}{\rho^2-r^2}
\]

abgebildet. Sei \(\ell_{n,k}=(\widetilde J_n^{-1})_{nk}\) und
\(a_n(w)=\sum_{k=0}^n\ell_{n,k}w^k=\sum_{j=0}^n b_j(w-c)^j\). Dann ist

\[
\sup_{|x-L|=r}|\mathcal D_n(x)|
\le\|Q_n\|_2\sum_{j=0}^n|b_j|R^j.
\tag{G21}
\]

**Beweis:** Die Kreisformel folgt durch Quadrieren von |ρ−Lw|=r|w|. Ferner \(\mathcal D_n(x)=\|Q_n\|_2a_n(\rho/x)\). Es genügt die Dreiecksungleichung für das **vollständige endliche** Taylorpolynom. ∎

Es gibt hier keinen unkontrollierten analytischen Rest, keine Winkelquadratur und keine L²-Normänderung. Für das reelle Modell gilt derselbe Wert am konjugierten Eckkreis. Der Schritt nutzt die schon global bestimmten biorthogonalen Objekte; er widerspricht nicht dem Auslöschungsbefund des festen Cauchy-Integrals aus Runde 2.

**Vektorversion [vollständige eigene Herleitung; P-Prüfung offen]:** Schreibt man \(v_n(w)=\widetilde J_n^{-1}(1,w,\dots,w^n)^T=\sum_j b_j^{\rm vec}(w-c)^j\), liefert dieselbe Dreiecksungleichung \(\sup_{|w-c|=R}G_n\le\sum_j\|b_j^{\rm vec}\|_2R^j\). Diese zweite endliche Majorante wird für die Anfangsnorm und für eine zusätzliche direkte G-Prüfung verwendet.

## 6. Neue Messwerte und Zertifikate

| N | √N·abs(𝒟_N(L)) [num.] | √N·max γ*-Proben [num.] | √N·max 24 Kreiswinkel [num.] | √N·sup ganzer Kreis < [zert.] |
|---:|---:|---:|---:|---:|
| 200 | 0,863116 | 1,952311 | 3,122847 | 3,126349 |
| 300 | 0,747577 | 1,608543 | 2,701080 | 2,702708 |
| 400 | 0,696721 | 1,511118 | 2,515280 | 2,516129 |
| 600 | 0,648607 | 1,471419 | 2,339235 | 2,339471 |
| 800 | 0,623216 | 1,476676 | 2,246198 | 2,246259 |
| 1000 | 0,606643 | 1,519989 | 2,185399 | 2,185448 |


In der Tabelle bedeutet „γ*-Proben“ dieselben festen 20 Punkte wie in Runde 2, zusätzlich L; „Kreis“ ist **|x−L|=ρ N⁻ᵝ mit β=1/1,1182**. Die numerischen Kreiswerte verwenden 24 Winkel. Sie dienen nur dem Vergleich mit der vollständigen Intervallmajorante. Die Kurvenpunkte der Referenzlinse aus W sind weiterhin numerisch konstruiert.

**Zertifizierung:** Das kubische Modell wird für die Kreisprüfung direkt aus dem exakten Lambert-W-Fixpunkt mit acb berechnet. Die linearen Lösungen sind geprüfte arb-Enclosures, kein `algorithm="approx"`. Sämtliche n+1 Taylorterme in (G21) werden mit Intervallarithmetik summiert. Die skalare Schranke √n sup|𝒟_n|<3,13 ist nur für die aufgeführten sechs Grade zertifiziert. Das zertifiziert **nicht G_n**, nicht die ganze feste γ*-Kurve, nicht ausgelassene Grade und nicht eine asymptotische Aussage.

Die großen monischen Koeffizientennormen sind wesentlich: ‖Q_200‖₂≈3,9473·10⁷, ‖Q_1000‖₂≈2,4660·10⁴⁰. Trotzdem sind die Produkte (G17) auf den geprüften Kurven klein. Dies betrifft die **exakten linearen Interpolanten**, nicht die offen gelassene nichtlineare K2-Koeffizientenschranke. Separate grobe Schranken der beiden Faktoren würden diese Auslöschung verlieren.

Als Operator auf Koeffizientendaten hat der einzelne Rang-eins-Term Norm \(\|Q_n\|_2\|\ell_n\|_2\), also die ℓ²-Koeffizientennorm von 𝒟_n. Aus den **bereits gespeicherten** W-Inversen ergeben sich 1,7280 bei n=200 und 7795,68 bei n=1000. Der Quotient zur gespeicherten größten Inversensingulärzahl steigt von 0,933 auf 1,364. Die Daten tragen damit dieselbe Eckverstärkung. Fits \(\log\|\Delta_n\|=cn+b\log n+a\) liefern c=0,012136 auf 200…1000, c=0,011954 auf 400…1000 und c=0,011915 auf 600…1000. Das passt zur Referenzrate etwa 0,0118, bestimmt sie aber nicht rigoros; die logarithmischen Vorfaktoren driften weiterhin. Es wurde hierfür keine neue SVD gerechnet.

**Anfangsnorm [zert.]:** Für n₀=200 ergibt die Vektorversion von (G21) rigoros
C₀= sup_{|x−L|≤r₂₀₀}G₁₉₉(x)<12,5026,
also C₀/√199<0,886283<6,26.
Damit ist die Konstante **C_G=6,26** im Eckkreis-Summenkriterium vollständig explizit,
**falls** die noch offene Schranke |𝒟_n|≤3,13/√n für **alle** n≥200 bewiesen wird.

**Direkte G-Kreise [zert.]:** Mit denselben geprüften Inversen und der vollständigen Vektor-Taylormajorante gilt
G₂₀₀(x)<0,8893√200 auf |x−L|=r₂₀₀ und
G₁₀₀₀(x)<0,8041√1000 auf |x−L|=r₁₀₀₀.
Beide Aussagen gelten für alle Winkel und alle Einheits-u; durch reelle Koeffizienten ebenfalls am konjugierten Eckkreis.
Die vorherigen Rasterwerte waren 0,842994 bzw. 0,762158. Die finite Taylormajorante darf darüber liegen.
**Diese zwei G-Zertifikate sind Einzelgrad-Aussagen, keine All-N-(KB).**


### 6.1 Kneser-Gegenprobe [numerisch, keine Übertragung bewiesen]

| N | √N·max γ*-Proben kubisch | √N·max γ*-Proben Kneser | √N·max Kreiswinkel Kneser |
|---:|---:|---:|---:|
| 200 | 1,952311 | 1,947102 | 3,132060 |
| 300 | 1,608543 | 1,603481 | 2,706573 |
| 400 | 1,511118 | 1,506771 | 2,519085 |
| 600 | 1,471419 | 1,463445 | 2,341793 |

Eingabe sind die zertifizierten Kneser-Koeffizientenmittelpunkte bis 640; J_N benötigt a_{N+1}. Die Koeffizientenradien werden hier nicht propagiert. Es wurde daher kein Kneser-Intervallzertifikat erstellt.

Die Rang-eins-Identität und das Summenkriterium gelten auch für nicht-polynomiale Keime, weil nur die verschachtelten endlichen Matrizen verwendet werden. Die feste Cauchy-Korrektur bleibt dort I_F−B_y. Die kubische endliche Vierzustandsrekursion und die Rang-6-Vereinfachung werden ausdrücklich nicht übertragen. Auch eine Ähnlichkeit der neuen Zahlen beweist keine All-n-Stabilität des Produkts 𝒟_n unter der Kneser-Störung.

Die kubisch konstruierte Vergleichskurve ist außerdem kein Geometriezertifikat für Knesers γ*. Eine All-n-Übertragung muss die Wahl der festen Kontur und die globale Amplitude gemeinsam kontrollieren.

**Keine ungeprüfte Konstantenübernahme:** Bei Kneser liefert bereits das Kreisraster für N=200 den Wert √N·max|𝒟_N|≈3,132060, oberhalb der kubischen Schranke 3,13. 3,14 wäre hier nur ein numerischer Kneser-Kandidat, kein Zertifikat. Die kleine Modellabweichung ist bei dieser knappen Konstante relevant.

Eine eigenständige Kneser-Version von (G19) würde die zusätzliche Cauchy-Randspur-Lücke umgehen: (G15)–(G20) sind reine endliche Koeffizientenidentitäten und setzen keine holomorphe Fortsetzung von f_K über L voraus. Der Beweis einer solchen All-n-Version ist hier allerdings ebenso offen wie im kubischen Modell.

## 7. Prüfungen, Grenzen und konkrete offene Stelle

- Die Rang-eins-Identität mit voller Cauchy-Quelle, die monische biorthogonale Normierung und die Eigenbasis wurden unabhängig mit SymPy an exakten rationalen Parametern geprüft. Die allgemeinen Beweise stehen oben.
- Die kritischen Werte und die Separationswerte in L sind 512-Bit-Intervallrechnungen für den exakten kubischen Fixpunkt.
- Die Kreismajoranten sind echte finite Intervallprüfungen; die γ*-Tabellen und die Kneser-Gegenprobe sind hochpräzise numerische Rechnungen. Ihre Quell- und Modellkonventionen stehen in den JSON-Dateien.
- Präzisionsgegenprobe: Bei N=200 und 600 reproduziert ein Lauf mit **300 zusätzlichen Bit** alle veröffentlichten G- und Zuwachswerte; die maximale Abweichung in der Skala |neu−alt|/(1+|alt|) beträgt 0. Dies ist numerische Reproduzierbarkeit, keine Fehlerschranke der Approx-Solves. Zusätzlich wurden **alle Grade 200…210** gerechnet: Die Höchstwerte sind √n·|𝒟_n|=1.952311 an den festen γ*-Proben und 3.122847 auf dem Kreisraster. Es gab in diesem kurzen Fenster keinen größeren Ausreißer. Die übrigen Zwischen- und späteren Grade sind ungeprüft.
- Es wurden keine Daten fremder Stränge verändert, keine Nachrichten an andere Chats gesendet, keine K2-Läufe gestartet und keine SVD berechnet.

**Noch fehlend:** Invertibilität der nötigen Hauptabschnitte und (G19) **für alle n≥n₀**, mit einer expliziten C_Δ auf der festen γ* oder einer gewichteten Variante des Summenarguments. Die vierteilige Diagonalisierung bestimmt diese normierten Interpolantenamplituden noch nicht. Die neuen Zertifikate liefern überprüfbare Zielwerte für dieses verbleibende globale Problem; sie schließen es nicht.

Der gezielte nächste mathematische Schritt ist eine starke Abschätzung des **normierten adjungierten Interpolanten 𝒟_n**, mit seinem Schur-Pivot h_n, aus der gekoppelten Momentaufgabe. Die zyklischen Permutationssprünge aus G8 und die Quelle aus G15 müssen darin bleiben. Eine bloße n-te-Wurzel-Rate von Q_n oder A_n reicht weiterhin nicht. Eine vollständig quellenhaltige globale RH-Parametrix wäre ein möglicher Weg zu genau dieser Abschätzung; ihre Existenz wird hier nicht vorausgesetzt.

W §25.1(c) hat zudem unterschiedliche Nullstellenkurven von A_n und der Eckmode festgestellt. Deshalb darf Ws äußere Phase H der Eckmode nicht ungeprüft als äußere Phase von 𝒟_n eingesetzt werden. Das Summenkriterium verlangt eine neue starke Abschätzung der exakten Interpolanten, keine Wiederholung der bekannten n-te-Wurzel-Sackgasse.

## 8. Reproduzieren und Ablage

Eigener Forschungsordner: `forschung/polynom_halbexp/g_gpt_hilfssaetze/`.

Neue Dateien:

- `experimente/h0_border.py`: numerische Grad-Zuwächse auf den festen γ*-Proben und den N-Kreisen, kubisch und Kneser.
- `experimente/h0_border_circle_cert.py`: rigorose Majorante auf dem gesamten Kreis für einzelne kubische Grade.
- `experimente/h0_border_exact_checks.py`: exakte Algebra, kritische Werte, Separationswerte und als separate numerische Gegenprobe die Wurzelmonodromie.
- `experimente/h0_border_prefix_cert.py`: finite Anfangsnorm für das bedingte Summenargument.
- `experimente/h0_border_pair_norms.py`: Auswertung bereits vorhandener W-Inversen und Singularwerte, keine neue SVD.
- `experimente/h0_border_plot.py`: wissenschaftliche Abbildung aus den gespeicherten Daten.
- `ergebnisse/h0_border_*.json`: Rohdaten mit Statuskennzeichnung.
- `RUNDE2_CAUCHY.md`: unverändertes Archiv des vorigen Berichts.
- `PRUEFUNG_RUNDE3.md`, `MANIFEST.json`: Prüfliste und SHA256-Prüfsummen.

Interpreter: `C:/Users/janis/miniconda3/python.exe`. Beispiel im Forschungsordner:

```powershell
python -B experimente/h0_border.py --project-root .. --output-dir ergebnisse --degrees 200,300,400,600,800,1000
python -B experimente/h0_border.py --project-root .. --output-dir ergebnisse --model kneser_z --degrees 200,300,400,600
python -B experimente/h0_border_circle_cert.py --output-dir ergebnisse --degrees 200,300,400,600,800,1000
python -B experimente/h0_border_exact_checks.py --output-dir ergebnisse
python -B experimente/h0_border_prefix_cert.py --output-dir ergebnisse
python -B experimente/h0_border_circle_cert.py --output-dir ergebnisse --degrees 1000 --tag _G1000 --certify-G
```

Die Skripte erzeugen keine Bytecode-Dateien in den Ordnern anderer Stränge. Es gibt keine versteckten Kontur- oder Quadraturzertifikate: Der Status der vollen Kreismajorante folgt allein aus (G21) und den geprüften Matrixenclosures.

**Rechenzeit:** Der kombinierte N=1000-Lauf brauchte 23,1 Minuten auf einem Kern: 18,3 Minuten für die geprüfte Inverse, weitere 4,8 Minuten für die zusätzliche Vektormajorante. Damit lag dieser kombinierte Lauf etwas über dem Richtwert von ungefähr 20 Minuten. Die skalare Kreisprüfung wird im Skript jetzt vor der optionalen Vektorauswertung als Zwischenstand gesichert. Alle anderen neuen Läufe waren kürzer; es wird kein laufender Prozess zurückgelassen.
