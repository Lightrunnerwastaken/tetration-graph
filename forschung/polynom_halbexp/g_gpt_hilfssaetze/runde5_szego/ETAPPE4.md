# Runde 5 — Etappe 4: exakte Struktur und ein bewiesenes nichtlineares Modell

**Abgeschlossene Etappe 4, 11.10.2026. Nur (H₀); K2 bleibt bei P.**

## 1. Ergebnis und Grenze

Neu sind zwei exakte Strukturformeln für die endlichen Abschnitte: ein gemischtes Symbol in zwei Variablen und eine Verschiebungsidentität vom Rang höchstens zwei. Außerdem wird eine uniforme inverse Abschätzung für eine ganze quadratische Familie mit dem richtigen Gewicht **M=F′∘F** bewiesen. Sie liefert explizite Determinanten, Pivots, Q und exponentielle Schranken für die ausgewerteten A- und 𝒟-Funktionen. Das ist ein tatsächlich gelöstes nichtlineares Hilfsmodell.

Die strikte inverse Kontraktion dieses Hilfsmodells fehlt an den ursprünglichen Randfixpunkten. Der Beweis lässt sich deshalb nicht durch einen unkontrollierten Störungsschritt auf das kubische/Kneser-Modell übertragen. Dort bleiben die geschlossene Formel für λ und (H₀) offen. Der globale quellenhaltige Cauchy-Operator ist weiterhin die notwendige Aufgabe.

Alle Modelle und Identitäten stehen in der dimensionslosen Normierung der Etappe 1:

\[
p_j(z)=F(z)^j+M(z)z^j,\qquad
J_n=([z^k]p_j)_{0\le k,j\le n},\quad
Q_n=h_nJ_n^{-1}e_n,\quad A_n=h_nJ_n^{-T}e_n,
\]
\[
\mathcal D_n(w)=\|Q_n\|_2 e_n^TJ_n^{-1}e(w)
                =\|Q_n\|_2\,A_n(w)/h_n.
\]

Es wird kein f_K als Eingabe einer Konstruktion benutzt.

## 2. Das exakte Symbol besitzt zwei Kanäle

**Lemma G39 [bewiesen].** Als formale Reihe in s, beziehungsweise für genügend kleines |s|, gilt

\[
\sum_{j\ge0}p_j(z)s^j
=\frac1{1-sF(z)}+\frac{M(z)}{1-sz}
=\frac{1+M(z)-s(z+M(z)F(z))}
       {(1-sF(z))(1-sz)}.
\tag{G39}
\]

*Beweis:* Zwei geometrische Reihen und gemeinsamer Nenner. ∎

Das gemischte Symbol hat die beiden Faktoren 1−sF(z) und 1−sz. Für M≠0 und F≠z kann keiner der beiden Faktoren identisch mit dem Zähler gekürzt werden: Einsetzen von s=1/F ergibt 1−z/F, Einsetzen von s=1/z ergibt M(1−F/z). Bei M=0 wird der Identitätskanal hingegen exakt gekürzt.

Das erklärt algebraisch, warum der Fall τ=0 nicht als gleichmäßig kleine Änderung der positiven τ-Familie behandelt werden darf. Es **beweist nicht**, dass für jedes positive τ die Hauptordnung der Determinante gleich ist. Die entsprechende numerische Beobachtung bleibt eine Konjektur.

Ein festes eindimensionales Toeplitz-Symbol liegt hier nicht vor. Die in Etappe 1 gelesenen Szegő-Sätze können daher nicht direkt auf J_n angewendet werden. Ein Riemann–Hilbert-Ansatz müsste diese beiden Kanäle und die inverse Verzweigung von F abbilden.

## 3. Eine exakte Verschiebungsstruktur vom Rang zwei

Sei T_n die Jetprojektion bis Grad n, S_n die untere nilpotente Verschiebungsmatrix und T_{F,n} die untere Toeplitzmatrix der Multiplikation mit F, jeweils auf P_n. Setze

\[
L_n=T_{F,n}+S_n,\qquad B_n=S_nT_{F,n}.
\]

**Lemma G40 [bewiesen].** Für n≥1 gilt

\[
R_n:=J_nS_n^2-L_nJ_nS_n+B_nJ_n=U_nV_n^T,\qquad
V_n=(e_{n-1},e_n),
\]
\[
U_n=\bigl(-T_np_{n+1},\;T_n(zF p_n)\bigr).
\tag{G40}
\]

Insbesondere rank R_n≤2. Mit

\[
\mathsf C_n=\begin{pmatrix}L_n&-B_n\\I&0\end{pmatrix},\qquad
\mathsf H_n=\binom{J_nS_n}{J_n}
\]

folgt die Sylvester-Identität

\[
\mathsf C_n\mathsf H_n-\mathsf H_nS_n=\binom{-R_n}{0}.
\]

*Beweis:* Beide Summanden von p_j erfüllen punktweise
p_{j+2}=(F+z)p_{j+1}−zFp_j.
Multiplikation und Jetprojektion kommutieren bei den hier vorkommenden Multiplikationsoperatoren. Für die Spalten j≤n−2 ergibt das R_ne_j=0. Bei j=n−1 fehlt der abgeschnittene Term p_{n+1}, also R_ne_{n−1}=−T_np_{n+1}. Bei j=n verschwinden beide verschobenen Spalten, und R_ne_n=T_n(zFp_n). Blockmultiplikation ergibt die zweite Identität. ∎

Das ist eine exakte Struktur für **alle** holomorphen F,M, nicht eine numerische Annäherung. Für polynomiales F ist T_{F,n} zudem eine Bandmatrix. Für Kneser ist sie im Allgemeinen voll.

**Grenze:** Die Begleitmatrix hat Dimension 2(n+1), und die Generatoren U_n wachsen mit n. Daraus folgt noch keine geschlossene vierdimensionale Schur-Rekursion mit uniformen Konstanten. Vor allem dürfen die Koeffizientennormen der Generatoren nicht an die Stelle der ausgewerteten Buchtfunktion treten. Die Struktur ist ein Ansatz für einen gekoppelten Cauchy-/Resolventenbeweis, kein Beweis von (G19).

Exakte rationale Kontrollen bei n=2,…,8 bestätigen die beiden Matrixidentitäten. Die Kontrollen ergänzen den Beweis.

## 4. Ein vollständig explizites Modell mit nichtkonstantem M

**Satz G41 [bewiesen].** Sei F(z)=bz, M(z)=m+εz, b>1, m≥0 und a=|ε|<1+m. Setze d_j=b^j+m. Dann

\[
h_n=d_n,\quad Q_n=e_n,\quad
(J_n^{-1})_{nk}=\frac{(-\epsilon)^{n-k}}{\prod_{j=k}^n d_j},
\]
\[
A_{nk}=d_n\frac{(-\epsilon)^{n-k}}{\prod_{j=k}^n d_j},\qquad
\mathcal D_n(w)=\sum_{k=0}^n
 \frac{(-\epsilon)^{n-k}w^k}{\prod_{j=k}^n d_j}.
\]

Für n≥1 und |w|≤1 gilt die explizite Schranke

\[
|\mathcal D_n(w)|
\le b^{-n}\left(1+
 \frac{a}{1-a/(1+m)}\,b^{-(n-1)}\right).
\tag{G41a}
\]

Auf r≤|w|≤1 mit a<r(1+m) gilt ferner

\[
\mathcal D_n(w)=b^{-n}w^n(1+E_n(w)),\quad
|E_n(w)|\le mb^{-n}+
\frac{a}{r[1-a/(r(1+m))]}\,b^{-(n-1)}.
\tag{G41b}
\]

*Beweis:* J_n ist untere Bidiagonalmatrix mit Diagonale d_j und Unterdiagonale ε. Rücksubstitution ergibt die Inversenformel. Die letzte Spalte der Inversen ist e_n/d_n, also Q=e_n. Für n−k=ℓ≥1 behalte die Faktoren d_n≥b^n und d_{n−1}≥b^{n−1}; alle weiteren d_j sind ≥1+m. Summation der geometrischen Majorante ergibt G41a. Für G41b dividiere durch b^{−n}w^n und verwende |b^n/d_n−1|≤mb^{−n}; die übrigen Terme werden mit a/r summiert. ∎

Für b=2, m=2, ε=1/4, r=1/2 erhalten wir

\[
|\mathcal D_n(w)|\le \frac{14}{11}\,2^{-n}\quad(|w|\le1,n\ge1),
\qquad |E_n(w)|\le\frac{16}{5}\,2^{-n}.
\]

Hier sind die Q_n orthonormale Koordinatenvektoren. Mit der exakten Grad-Zuwachs-Identität ergibt sich auf |w|≤1 für alle N

\[
G_N(w)^2=\sum_{n=0}^N|\mathcal D_n(w)|^2
\le\frac19+\frac{196}{121}\sum_{n\ge1}4^{-n}
=\frac{709}{1089},
\qquad G_N(w)\le\frac{\sqrt{709}}{33}.
\]

Dies ist eine vollständige uniforme Schranke, ohne SVD. Das Gewicht in diesem Hilfsmodell ist **nicht** F′∘F, außer bei ε=0; der folgende Satz erfüllt diese zusätzliche Struktur.

## 5. Ein echtes nichtlineares Modell: uniforme quadratische Störung

**Satz G44 [bewiesen].** Für reelles |η|≤1/20 sei

\[
F_\eta(z)=2z+\eta z^2,\qquad
M_\eta(z)=F_\eta'(F_\eta(z))=2+2\eta F_\eta(z).
\]

Dann sind sämtliche J_n invertierbar. Für alle n≥0 gelten

\[
h_n=2^n+2,\qquad Q_n=e_n,\qquad
\det J_n=2^{n(n+1)/2} C(2,2)e^{-R_n},
\quad 0\le R_n\le 2^{1-n}.
\tag{G44a}
\]

Sei H_η der bei 0 fixierende inverse Zweig von F_η und
\(\mathsf T_\eta u(z)=F_\eta'(z)u(H_\eta(z))\).
Auf dem Hardyraum H²(D_{3/2}) ist I+𝕋_η invertierbar mit

\[
\|(I+\mathsf T_\eta)^{-1}\|<4.
\tag{G44b}
\]

Setze u_k=(I+𝕋_η)^{−1}H_η^k. Dann gelten **exakt**

\[
(J_n^{-1})_{nk}=[z^n]u_k,\quad
A_{nk}=(2^n+2)[z^n]u_k,\quad
\mathcal D_n(w)=\sum_{k=0}^n[z^n]u_k\,w^k,
\]

und die uniformen funktionalen Schranken

\[
|[z^n]u_k|\le4(2/3)^n(4/5)^k,\qquad
\sup_{|w|\le1}|\mathcal D_n(w)|\le20(2/3)^n,
\qquad \sup_{N,|w|\le1}G_N(w)\le60/\sqrt5.
\tag{G44c}
\]

### Beweis mit expliziten Konstanten

1. **Globaler inverser Zweig.** Setze R=3/2 und r=4/5. Auf |z|=r gilt für |ξ|≤R
   \[
   |2z|-|\eta z^2|-|\xi|
   \ge 2r-\tfrac1{20}r^2-R=\frac{17}{250}>0.
   \]
   Rouché liefert genau eine Lösung von F_η(z)=ξ in D_r. F′ hat dort keine Nullstelle. Dieser Zweig H_η ist holomorph auf einer Umgebung von \(\overline{D_R}\), H_η(0)=0 und |H_η|≤r.

2. **Subordination.** Für u∈H²(D_R) gilt \(\|u\circ H_\eta\|\le\|u\|\). Tatsächlich majorisiert die harmonische Poisson-Fortsetzung der Randwerte |u|² die Funktion |u|²; nach Komposition mit H_η bleibt die Majorante harmonisch. Ihr Kreismittel ist ihr Wert bei H_η(0)=0. Approximation durch Polynome liefert die Aussage für H². Ferner \(\|F_\eta'\|_{\infty,D_R}\le2+2|\eta|R=43/20\).

3. **Invarianter Schwanz.** V=z⁴H²(D_R) ist invariant unter 𝕋_η. Für u=z⁴v:
   \[
   \|\mathsf T_\eta u\|
   \le \frac{43}{20}r^4\|v\|
   =\frac{43}{20}(r/R)^4\|u\|
   =\frac{44032}{253125}\|u\|<\frac15\|u\|.
   \]
   Also besitzt I+𝕋_η auf V eine Inverse von Norm ≤5/4.

4. **Endlicher Kopf.** Auf den ersten vier Koeffizienten hat I+𝕋_η die Diagonale 3,2,3/2,5/4. Die einzigen Einträge unterhalb der Diagonale sind 2η,3η/4,−η²/8,η/4 an den Positionen (1,0),(2,1),(3,1),(3,2). In der **orthonormalen Hardybasis** z^j/R^j sind diese mit R^{k−j} zu multiplizieren. Die quadrierte Frobeniusnorm ist höchstens
   \[
   (2\eta R)^2+(3\eta R/4)^2+(\eta^2R^2/8)^2+(\eta R/4)^2
   \le\frac{4262481}{163840000}<\frac1{16}.
   \]
   Die kleinste Diagonale ist 5/4, also hat die inverse Kopfmatrix Norm ≤1.

5. **Blockinversion.** Bezüglich P₃⊕V ist I+𝕋_η untere Blockdreiecksmatrix. Die Kopplung Kopf→Schwanz hat Norm ≤43/20. Für f=f_P+f_V folgt
   \[
   \|u_P\|\le\|f_P\|,\qquad
   \|u_V\|\le\tfrac54(\|f_V\|+\tfrac{43}{20}\|f_P\|).
   \]
   Dreiecksungleichung und Cauchy–Schwarz liefern
   \[
   \|(I+\mathsf T_\eta)^{-1}\|
   \le\sqrt{(1+\tfrac54\tfrac{43}{20})^2+(\tfrac54)^2}
   =\frac{\sqrt{3881}}{16}<4.
   \]

6. **Exakte Quelle und endliche Abschnitte.** Die rückwärts geschriebene Operatorgleichung lautet
   \[
   u(z)+F_\eta'(z)u(H_\eta(z))=v(H_\eta(z)).
   \]
   Sie ist nahe 0 exakt äquivalent zu u∘F_η+M_ηu=v, weil M_η(H_η(z))=F_η′(z). Für v=z^k erhält man u_k und \(\|u_k\|_{H^2(D_R)}\le4r^k\). Der Hardy-Koeffizientenbound ergibt G44c. Da F_η(0)=0, hängen die ersten n+1 Koeffizienten der vorwärts geschriebenen Gleichung nur von den ersten n+1 Koeffizienten von u ab. Die endlichen Inversen sind deshalb **exakt** diese Jets, ohne Projektionsrest. Dreieckigkeit beweist G44a mittels G27; die Summation über k gibt 4/(1−r)=20. Schließlich ist
   \[
   \sum_{n\ge0}|\mathcal D_n(w)|^2
   \le400\sum_{n\ge0}(4/9)^n=720,
   \]
   also G_N≤√720=60/√5. ∎

### Warum dies die ursprüngliche Lücke nicht schließt

Der Satz benutzt einen **global** kontrollierten inversen Zweig auf einem festen größeren Kreis und eine strikt kontraktive Schwanzgleichung. Das unterscheidet ihn von einer isolierten lokalen Eckabschätzung. Die Quelle bleibt exakt erhalten.

Die strikte Inklusion H_η(D_R)⊂D_r mit r<1 schließt aber einen Fixpunkt von F_η auf |z|=1 aus. Bei einem solchen Randfixpunkt liegt auch H(L)=L auf |z|=1, und gerade die hier gewonnene Expansionsreserve ist nicht vorhanden. Hinzu kommt bei F(0)≠0, dass endliche Jets und inverse Komposition nicht so kommutieren. **Keine uniforme Störung bis zur ursprünglichen Geometrie wurde bewiesen.** Der dort beobachtete √n-Faktor wird durch G44 nicht erklärt.

## 6. Exakte Parameterableitung für λ-Tests

**Lemma G42 [bewiesen].** In einer holomorphen Parameterfamilie θ↦(F,M) mit invertierbaren Hauptabschnitten gilt lokal für n≥2

\[
\partial_\theta\log\frac{h_n}{h_{n-1}}
=\operatorname{tr}(J_n^{-1}\partial_\theta J_n)
-2\operatorname{tr}(J_{n-1}^{-1}\partial_\theta J_{n-1})
+\operatorname{tr}(J_{n-2}^{-1}\partial_\theta J_{n-2}),
\]
\[
(\partial_\theta J_n)_{kj}
=[z^k]\{jF^{j-1}\partial_\theta F+(\partial_\theta M)z^j\}.
\tag{G42}
\]

*Beweis:* Jacobi für jede nichtverschwindende Determinante; dann die zweite Differenz in n. ∎

Eine Vertauschung von ∂_θ mit n→∞ ist hierin **nicht** enthalten. Ein globaler λ-Kandidat müsste seine Variation mit dem Grenzwert dieser exakten zweiten Differenz in Einklang bringen. Das ist ein strengerer Test als die Übereinstimmung einer einzigen Kapazitätszahl.

## 7. Ein präzises globales Parametrixkriterium

**Lemma G43 [bewiesene bedingte Implikation; Hypothesen beim ursprünglichen Modell offen].** Für jedes n seien X_n,Y_n endlichdimensionale normierte Räume mit Isomorphismen zur Koeffizientenbeschreibung. Die X_n-Norm majorisiere die ausgewertete Funktion auf der verlangten Kontur sowie den gewichteten inneren Wienerbound. Sei ℒ_n:X_n→Y_n die exakte Darstellung von J_n^T; die Quelle e_n ist im entsprechenden Y_n dargestellt. Sei q_n^0>0 und \(\|Q_n\|_2\le C_Qq_n^0\).

Angenommen, eine rechte Parametrix \(\mathsf B_n:Y_n\to X_n\) erfüllt

\[
\mathcal L_n\mathsf B_n=I-\mathsf R_n,\qquad
\|\mathsf R_n\|_{Y_n\to Y_n}\le r_0<1
\]

und die **quellenbezogenen** Bounds für alle k≥0

\[
q_n^0\|\mathsf B_n\mathsf R_n^k e_n\|_{X_n}
\le C_B\,n^{-1/2}r_0^k .
\tag{G43}
\]

Dann ist ℒ_n invertierbar und die ausgewertete 𝒟-Funktion erfüllt

\[
\|\mathcal D_n\|_{X_n}\le
\frac{C_QC_B}{1-r_0}\,n^{-1/2}.
\]

*Beweis:* Neumann gibt (I−ℝ_n)^{-1}. Die rechte Inverse 𝔅_n(I−ℝ_n)^{-1} einer quadratischen endlichen Matrix ist ihre Inverse. Wende diese auf die exakte Quelle an, summiere die Bounds und multipliziere mit \(\|Q_n\|\). Die genaue Grad-Zuwachs-Identität wird benutzt, keine Koeffizientenmajorante vor der Auswertung. ∎

Ein stärkerer Operatorbound \(q_n^0\|\mathsf B_n\|_{Y_n\to X_n}\le C_B/\sqrt n\) bei \(\|e_n\|_{Y_n}\le1\) wäre hinreichend. Die Wahl von Y_n darf aber nicht ungeprüft auf die ungewichtete Koeffizienten-ℓ²-Norm beschränkt werden. Die schwierigste Quelle ist genau die Cauchy-Quelle der adjungierten Schurreihe.

G43 ist ein präzises Ziel für ein globales Small-Norm-Verfahren. Es behauptet **nicht**, eine solche Parametrix bereits für das kubische oder Kneser-Modell zu besitzen.

## 8. Beweisstrategie mit einzeln benannten Lücken

Die folgenden Schritte sind **Skizze**, nicht bewiesene Sätze für das ursprüngliche Modell.

1. **Exakte gekoppelte Form.** G39/G40 und die quellenhaltige Cauchy-Gleichung auf einer festen Kontur verwenden. Für kubisches F müssen dessen inverse Blätter und der Identitätskanal gemeinsam auftreten. Für Kneser muss der tatsächliche Korrekturterm **I_F−B_y** erhalten bleiben. B_y allein beschreibt die nichtpolynomiale Jetgleichung nicht.

2. **Globales Klebungsproblem.** Existenz, Eindeutigkeit und geeignete Uniformisierung der geklebten Fläche ℳ/Ω_y beweisen. Zwei intrinsische Normierungen an 0 und ∞ sowie ihre Beziehung zu den Schurquellen bestimmen. **Lücke L1:** Der numerische freie Rand ist kein bewiesenes Gleichgewichtsproblem mit passender Quelle.

3. **Äußere starke Profile.** Auf einer gemeinsamen festen Kontur quellenverträgliche Profile für JQ_n und A_n gewinnen. Nur dann liefert G31 den Kandidaten \(\lambda=\kappa_0/\kappa_\infty\). **Lücke L2:** Uniforme Profilreste und das Nichtverschwinden des Matchingintegrals H. Die n-te-Wurzel-Ebene liefert beides nicht.

4. **Innere Parametrix an beiden Ecken.** Die komplexe Skala \(n^{-a}\), a=1/κ, mit der äußeren Phase und der Quelle matchen. **Lücke L3:** Das ganze analytische Profil, seine Amplitude und sein Wienerrest auf Radius 1,05, einschließlich beider Seiten des sektoriellen Problems. Ein lokales Mellin-/Grenzproblem ohne globale Normierung reicht nicht.

5. **Uniforme Restinversion.** Die Parametrix im Sinne von G43 mit explizitem r₀<1, C_Q,C_B und n₁ abschließen. Alternativ zuerst nur die beiden G37-Bounds erreichen: \(B_n\le0,75\), \(\|S_n\|_{\mathcal W,1.05}\le4\) für n≥300 würden die **Eckkreise** mit \(C_\Delta=3\) kontrollieren. **Lücke L4:** Diese Bounds sind nur an gespeicherten Graden plausibel; auf der übrigen festen γ* fehlt eine uniforme funktionale Abschätzung ebenfalls.

6. **Endlicher Anfang und Störung.** Erst einen rigorosen endlichen Anfang zertifizieren, dann eine geometrisch kontrollierte Homotopie vom kubischen Modell zu Kneser mit dem vollständigen I_F−B_y. **Lücke L5:** Uniforme Konturen, Klebungstopologie und quellenbezogene Resolventenbounds entlang der Homotopie. Kleine Koeffizientenunterschiede beweisen keine uniforme Stabilität der hochgradigen Abschnitte.

G44 zeigt, dass ein solcher exakter inverser Cauchy-Weg im strikt expandierenden, bei 0 fixierenden Bereich funktioniert. Der Übergang zum Randfixpunkt ist die wesentliche neue Analyse, kein Routine-Störungslemma. Die reine Komposition G29 zeigt zusätzlich, dass der Vorfaktor sich zwischen Operatortypen ändern kann.

## 9. Status der Ergebnisse

| Aussage | Status |
|---|---|
| G25–G30: einfache Determinanten-/Schurmodelle und Zwei-Zweig-Integral | bewiesen |
| G31–G34: starke Matching- bzw. Summationsimplikationen | bewiesen unter präzisen, ursprünglicherweise offenen Hypothesen |
| G35, G37: komplexe Skala und uniforme Kreiseinbettung | bewiesen; skalare Zahlen intervallzertifiziert |
| G38: Konsequenz auf den alten Eckkreisen | bewiesen unter offenen Profil-/Amplitudenhypothesen |
| G39–G42, G44: Symbol, Rang-zwei-Struktur, explizite Modelle, Variation | bewiesen |
| G43: globale quellenbezogene Parametriximplikation | bewiesen unter offenen Hypothesen |
| Kubische/Kneser-Werte λ,c₂; natürliche normierte Profilform | hochpräzise endliche Numerik |
| λ als Quotient konformer Normierungen; ν=2 Im(1/κ) | Vermutungen |
| All-n-Amplitude, äußerer funktionaler Bound und Matching | offen |
| Geschlossener Szegő-Satz für das ursprüngliche Modell und (H₀) | nicht bewiesen |

**Literaturstatus unverändert zu Etappe 1:** Teile von Böttcher–Widom wurden gelesen; Cowen–MacCluer, Böttcher–Silbermann/Widom nur in der dort angegebenen Reichweite zitiert. Für G39–G44 wird kein ungelesenes Literaturtheorem als Beweis verwendet. Der Subordinationsschritt in G44 ist oben direkt bewiesen.

## 10. Nächster erfolgversprechender Schritt

Der engste sinnvolle Schritt ist die **globale inhomogene Parametrix für den exakten Zwei-Kanal-Cauchy-Operator**, zunächst kubisch. Ihre innere Norm soll das gesamte unnormierte Profil \(\sqrt n\,\mathcal D_n(w_L+n^{-a}\zeta)\) in \(\mathcal W_{1.05}\) erfassen. Das würde die Eckschranke liefern, ohne zuerst eine geschlossene Formel für S_∞ oder eine stabile Extrapolation von B_∞ zu verlangen. Die äußere Kontrolle auf der festen γ* muss im selben quellenhaltigen Verfahren mitgeführt werden.

Keine weitere Grad-5N-Konstruktion, keine Bearbeitung von K2, keine Wiederholung der lokalen oder ebenen Sackgassen.
