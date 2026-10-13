# Analytische Klassen und ein holomorpher quadratischer Halbrangjet

Stand: 2026-10-12. Konvention H₃(b,x)=bˣ, H₄=T_b. Rang und Höhe sind unterschiedliche Variablen. Der alte H₀-Beweis ist pausiert.

## 1. Was hier konstruiert wird

S ist ein normierter Rangschritt: S(f)(0)=1 und S(f)(x+1)=f(S(f)(x)). Am Fixpunkt f₀(x)=x+1 betrachten wir

\[
f_\varepsilon(x)=x+1+\varepsilon v(x),\qquad v(x)=x(x-1).
\]

Für einen ganzzahligen Rangzuwachs n betrachten wir zunächst die formalen Jets bis Ordnung ε², mit der Summationsauswahl aus §2:

\[
S^n(f_\varepsilon)(x)=x+1+\varepsilon a_n(x)+\varepsilon^2q_n(x)+O(\varepsilon^3).
\]

Die Koeffizienten werden hier analytisch in n fortgesetzt. **Bewiesen wird eine holomorphe Funktion q₁/₂(x) auf Re x>−1.** Sie setzt die intrinsischen quadratischen Halbrangwerte an allen ganzzahligen Höhen fort. Eine nichtlineare Ranggruppe auf einer offenen Funktionenmenge und ihre Anwendung auf bˣ werden damit noch nicht konstruiert. Insbesondere ist die Testfamilie fε keine Exponentialfamilie.

## 2. Eine Halbebenenskala für analytische Summation

**Bewiesen.** Für ρ>0, a≥0, p≥0 seien Ωρ={Re z>−ρ} und

\[
\|f\|_{\rho,a,p}=\sup_{z\in\Omega_\rho}
\frac{|f(z)|e^{-a(\Re z)_+}}{(1+|z|)^p}.
\]

Die Funktionen mit endlicher Norm bilden eine Banachklasse. Produkte addieren a und p. Ist 0<ρ′<ρ und 0<δ<ρ−ρ′, liefert Cauchy

\[
\|f'\|_{\rho',a,p}\le C_D\|f\|_{\rho,a,p},\quad
C_D=\delta^{-1}e^{a\delta}(1+\delta)^p.
\]

Polynome, Gammaquotienten Bν(z)=Γ(z+ν)/(Γ(z)Γ(ν+1)) und deren Ableitungen sind für geeignete ρ enthalten; außerdem bᶻ für a≥log b. Die Gammafunktion bei z=0 wird durch 1/Γ(z) fortgesetzt.

### Normierte Summation

Setze

\[
\mathfrak S f(z)=\int_0^z f(t)dt+\frac{f(z)-f(0)}2+
i\int_0^\infty\frac{f(it)-f(-it)-f(z+it)+f(z-it)}{e^{2\pi t}-1}\,dt.
\tag{A1}
\]

Der erste Weg ist das gerade Segment. Alle Punkte liegen in Ωρ. Die Differenzen im Zähler sind O(t) am Ursprung; am Unendlichen dominiert der Nenner jedes vertikale Polynom. Lokal gleichmäßige Konvergenz beweist Holomorphie. Es gilt

\[
\|\mathfrak S f\|_{\rho',a,p+1}\le C_S\|f\|_{\rho,a,p},\qquad
C_S=2+4C_D J_p,
\]

\[
J_p=\int_0^\infty\frac{t(1+t)^p}{e^{2\pi t}-1}dt.
\]

Für ganzzahliges p ist Jp die endliche Summe

\[
\sum_{j=0}^{p}\binom pj\frac{\Gamma(j+2)\zeta(j+2)}{(2\pi)^{j+2}}.
\]

Die endliche Abel–Plana-Formel gibt Sfrak f(k)=Σⱼ₌₁ᵏ f(j), insbesondere Sfrak f(0)=0. Quelle: Saharian, §2, (2.14), gelesen; siehe LITERATUR.md.

### Eindeutigkeit der ganzzahligen Abtastung

**Bewiesen.** Hat f in der Vereinigung dieser Klassen f(k)=0 für alle positiven ganzen k, dann f≡0. Auf Re z>0 ist e⁻ᵃᶻ f(z)/(1+z)ᵐ für m≥p beschränkt. Unter w=(z−1)/(z+1) bilden die Nullstellen k eine Folge mit divergierender Blaschke-Summe: Σ(1−|(k−1)/(k+1)|)=∞. Jensens Formel erzwingt die identisch verschwindende Funktion. Der Identitätssatz setzt dies auf Ωρ fort.

Folglich gilt analytisch

\[
\mathfrak S f(z+1)-\mathfrak S f(z)=f(z+1).
\tag{A2}
\]

Denn die Differenz beider Seiten gehört wieder zur Skala und verschwindet an allen nichtnegativen ganzen Zahlen. Auch die normierte Lösung dieser Differenzgleichung ist innerhalb der Skala eindeutig. Es wird keine unbeschränkte periodische Eichfunktion zugelassen.

### Der echte quadratische volle Rangschritt

Aus der Superfunktionsgleichung folgt

\[
a_1=\mathfrak S v,\qquad
q_1=\mathfrak S[(\mathfrak S v-v)v'].
\tag{A3}
\]

Die polarisierte Form ist ½ Sfrak[((Sfrak u−u)v′+(Sfrak v−v)u′)]. Produkte, Cauchy und (A1) geben explizite stetige Abbildungen zwischen Stufen der Skala. Im Allgemeinen erhöhen sich a zu 2a und p zu 2p+2. Dies ist **keine** selbstabbildende einzelne Banachklasse.

Für unser v:

\[
a_1(x)=x(x-1)(x+1)/3,\quad
q_1(x)=x(x-1)(x-2)(x+1)(8x+1)/60.
\]

Damit werden auch Produkte mit mehrfachen Polen außerhalb Ωρ behandelt; eine endliche Gamma-Modensumme ist nicht erforderlich.

## 3. Analytischer linearer Rangschritt

**Bewiesen.** Für s>0 definiere Kₛ(u)=Γ(u+s)/(Γ(s)Γ(u+1)). Ersetzt man in (A1) f(t) durch Kₛ(z−t)f(t), erhält man

\[
L_s f(z)=\int_0^z K_s(z-t)f(t)dt+\frac{f(z)-K_s(z)f(0)}2
\]
\[
\quad+i\int_0^\infty\frac{K_s(z-it)f(it)-K_s(z+it)f(-it)
-K_s(-it)f(z+it)+K_s(it)f(z-it)}{e^{2\pi t}-1}dt.
\tag{A4}
\]

Auf Ωρ′ mit ρ′<min(ρ,s) sind die Kernelpole ausgeschlossen. Der Gammaquotient wächst dort höchstens polynomial; dieselben Integralabschätzungen geben eine Stufe mit endlichem a und p. Am ganzzahligen k ist

\[
(L_sf)(k)=\sum_{j=1}^k\frac{(s)_{k-j}}{(k-j)!}f(j).
\]

Die Faltung dieser Gewichte, die Eindeutigkeit der Abtastung und f(0)=0 beweisen LₛLₜ=Lₛ₊ₜ und L₁=Sfrak. Vandermonde beweist

\[
L_s B_\nu=B_{\nu+s}.
\]

Für v gilt

\[
a_s(x)=2(B_{s+2}(x)-B_{s+1}(x)).
\tag{A5}
\]

Die Integrale wurden auf einem reellen und einem komplexen Höhenpunkt gegen die geschlossenen Gammaformeln geprüft, bei 50 und 75 Dezimalstellen. Das sind numerische Kontrollen, keine Intervallzertifikate.

## 4. Intrinsische quadratische Daten ohne endliche Polynomableitung

**Bewiesen.** Für k≥2 gilt

\[
a_n(k)=\frac{2(k-1)(n+3)_{k-2}}{\Gamma(k)},\quad
a_n'(k)=a_n(k)[H_{n+k}-H_{k-2}].
\]

An der Höhe m≥3 ist die Quelle (aₙ₊₁−aₙ)aₙ′ gleich

\[
P_m(n)[H_{n+m}-H_{m-2}],\quad
P_m(n)=\frac{4(m-1)(m-2)(n+4)_{m-3}(n+3)_{m-2}}
{\Gamma(m)\Gamma(m-1)}.
\]

Daraus folgt

\[
q_{n+1}(k)=\sum_{\ell=1}^k q_n(\ell)+
\sum_{m=3}^kP_m(n)[H_{n+m}-H_{m-2}].
\tag{A6}
\]

Es gibt rationale Polynome Uₖ,Vₖ mit

\[
q_n(k)=U_k(n)H_{n+2}+V_k(n).
\]

Man erhält sie durch rationale Antidifferenzen. Mit Rₖ=Σₘ₌₃ᵏPₘ und

\[
S_k=\sum_{m=3}^kP_m(n)\left[\sum_{h=3}^m\frac1{n+h}-H_{m-2}\right]
\]

lauten die Rekursionen

\[
\Delta U_k=\sum_{\ell<k}U_\ell+R_k,\quad U_k(-2)=0,
\]
\[
\Delta V_k=\sum_{\ell<k}V_\ell+S_k-\frac{U_k(n+1)}{n+3},\quad
V_k(0)=-U_k(0)H_2.
\tag{A7}
\]

Alle Brüche sind Polynome: Pₘ ist durch n+h für 3≤h≤m teilbar, und Uₖ(n+1) durch n+3. Die Rekursionen beweisen (A6) samt q₀=0 induktiv.

Die **gewählte kanonische Rangfortsetzung** ist

\[
q_s(k)=U_k(s)[\psi(s+3)+\gamma]+V_k(s).
\tag{A8}
\]

Für jedes feste k liegt sie in einer Klasse polynomialen Wachstums auf Re s>0; dort ist sie durch ihre Werte an ganzen Rängen eindeutig bestimmt. Diese Auswahl ist wesentlich. Es ist keine freie Kurvenanpassung in s.

Bei s=½ ist H₅/₂=46/15−2log 2; alle Werte haben daher die exakte Form rational+rational·log 2. Insbesondere q₁/₂(3)=341/12−35log 2.

## 5. Geschlossene erzeugende Struktur

**Bewiesen.** Uₛ(z)=ΣₖUₖ(s)zᵏ erfüllt

\[
U_s(z)=2z^2(1-z)^{-s-3}
\left[P_{s+2}\left(\frac{1+z}{1-z}\right)-1\right].
\tag{A9}
\]

Beweis: Die Quelle von (A7) hat die erzeugende Funktion
4(n+3)z³ ₂F₁(n+4,n+4;2;z)/(1−z). Pfaffs Transformation macht daraus Pₙ₊₂⁽¹,⁰⁾(t), t=(1+z)/(1−z). Nach Multiplikation mit (1−z)ⁿ teleskopiert die Rangsumme mittels

\[
(1-t)P_m^{(1,0)}(t)=P_m(t)-P_{m+1}(t).
\]

Die Bedingung U₋₂=0 bestimmt die Konstante. Für jedes feste Höhenkoeffizient ist die Fortsetzung in s durch die Abtastung eindeutig.

Sei αⱼ(s)=Δₖʲqₛ(0) der **Höhen-Newtonkoeffizient**. Die Binomialtransformation w↦w/(1+w) liefert

\[
\mathcal A_s(w):=\sum_{j\ge0}\alpha_j(s)w^j
=2w^2(1+w)^s\left[J_{s+2}(1+2w)-\tfrac32(1+2w)-W_s(w)\right].
\tag{A10}
\]

Hier wird J und W eindeutig durch die folgenden Potenzreihen definiert:

\[
J_\nu(1+2w)=\sum_{j\ge0}h_j(\nu)[H_\nu-H_j]w^j,
\quad h_j(\nu)=\frac{(-1)^j(-\nu)_j(\nu+1)_j}{(j!)^2},
\]

\[
W_s(w)=H_{s+2}-H_2+\sum_{j\ge1}w_j(s)w^j,
\]

\[
w_j(s)=(\mathfrak S_m p_j)(s+1)-p_j(1),\quad
p_j(m)=\frac{m(m-1)\cdots(m-j+1)(m+2)\cdots(m+j)}{(j!)^2}.
\tag{A11}
\]

Sfrakₘ auf einem Polynom ist die gewöhnliche normierte polynomial fortgesetzte Summe. Die leeren Produkte für j=1 geben p₁=m. Für s=½ ist jedes wⱼ rational.

Herleitung von (A10): Die Quelle in (A6), nach Pfaff, ist die Parameterableitung D=∂α−∂β von Pₙ₊₂⁽¹,⁰⁾. Die Jacobi-Verbindung bei α=ε, β=−ε gibt

\[
(1-t)D P_m^{(1,0)}=J_m-J_{m+1}+P_m/(m+1),
\]

wobei Jₘ=D Pₘ⁽⁰,⁰⁾. Die Summe m=2,…,n+1 teleskopiert. Durch Differentiation der hypergeometrischen Jacobi-Reihe entsteht die angegebene J-Reihe. Jeder feste w-Koeffizient ist ein Polynom in s mal Hₛ₊₂ plus ein Polynom; die ganzzahligen Rangwerte identifizieren ihn mit den aus (A8) berechneten Newtonkoeffizienten.

Eine zusätzliche Integral-Auswertung verwendet Jν(t)=Pν(t)·½log((t+1)/(t−1))−Qν(t), auf dem Zweig t>1, sowie die Abel–Plana-Rangsumme von Pₘ(t)/(m+1). Ihre numerische Kontrolle steht in generating.py. Die folgende Holomorphieaussage benötigt diese zusätzliche Integral-Auswertung nicht.

## 6. Holomorphie des quadratischen Halbrangjets

### Eine explizite schwache Schranke

**Bewiesen.** Für j≥3 und m in 0≤Re m≤3/2 gilt durch Spiegelung

\[
p_j(m)=\frac{(-1)^{j+1}\sin(\pi m)}{\pi(m+1)}
\frac{\Gamma(j-m)\Gamma(j+m+1)}{\Gamma(j+1)^2}.
\tag{A12}
\]

Auf der reellen Strecke [0,3/2] folgt |pⱼ(m)|≤315/(256j). Zum Beweis ist der Gammaquotient in Re m wachsend; für a=3/2 ist j·Γ(j−a)Γ(j+a+1)/Γ(j+1)² bei j=3 gleich 315π/256 und danach fallend.

Auf dem größeren Streifen −1/2≤Re m≤2 erhält man

\[
|p_j(m)|\le\frac{20}{\pi j}\cosh(\pi\Im m).
\]

Dabei wird |Γ(u+iy)|≤Γ(u), u>0, benutzt. Der maximale reelle Quotient ist bei a=2; j mal dieser Quotient ist höchstens 10. Cauchy mit Radius 1/2 auf den Vertikalen Re m=0,3/2 liefert

\[
|p_j'(m)|\le\frac{40e^{\pi/2}}{\pi j}\cosh(\pi\Im m).
\]

In (A1) begrenzt dies den Zähler durch 4t mal diese Ableitung. Da

\[
\int_0^\infty\frac{t\cosh(\pi t)}{e^{2\pi t}-1}dt
=\frac18-\frac1{2\pi^2},
\]

gilt für s=½

\[
|w_j|\le\frac{C_W}{j},\quad
C_W=\frac{315}{128}+\frac{160e^{\pi/2}}\pi
\left(\frac18-\frac1{2\pi^2}\right)<21.
\tag{A13}
\]

Für hⱼ(5/2) ergibt derselbe Quotientvergleich

\[
|h_j(5/2)|\le C_h/j,\qquad C_h=3465/256<14.
\]

Setze H=H₅/₂ und dⱼ=[wʲ](J₅/₂−3(1+2w)/2−W₁/₂). Dann d₀=0 und für j≥3

\[
|d_j|\le E(j)/j,\qquad E(j)=C_h(H+1+\log j)+C_W.
\]

Aus (A10), Σ|binom(½,k)|=2 und |binom(½,k)|≤½k⁻³/² folgt, mit m=j−2≥6,

\[
|\alpha_j|\le\frac{8E(m)}m+
\frac{2\sqrt2}{m^{3/2}}\left[|d_1|+|d_2|+E(m)(1+\log m)\right].
\tag{A14}
\]

Man trennt die Faltung bei k=m/2. Die wenigen kleinen Koeffizienten sind explizit. Insbesondere genügt die grobe globale Schranke

\[
|\alpha_j|\le1500(1+\log j)/j\quad(j\ge3).
\tag{A15}
\]

Zum Nachprüfen der Zahl genügen H<2, C_h<14, C_W<21, |d₁|+|d₂|<5 und (1+log m)/√m≤2/√e. Die kleinen j=3,…,7 sind rationale Zahlen plus rationale Vielfache von log 2; etwa 0.69<log 2<0.70 genügt zur Prüfung.

### Der analytische Kandidat

**Bewiesen.** Die Reihe

\[
\boxed{q_{1/2}(x)=\sum_{j=3}^\infty\alpha_j(1/2)\binom{x}{j}}
\tag{A16}
\]

konvergiert lokal gleichmäßig und absolut für **Re x>−1**. Denn auf jedem kompakten Teil dieser Halbebene gilt |binom(x,j)|≤C j⁻ᴿᵉˣ⁻¹, und (A15) ist summierbar. Sie definiert daher eine holomorphe Funktion und stimmt bei jedem nichtnegativen ganzen k mit (A8) überein. Insbesondere q₁/₂(0)=q₁/₂(1)=q₁/₂(2)=0.

Für x=½ gibt (A15) den **bewiesenen, sehr groben** Rest

\[
|R_J|\le750J^{-3/2}\left[\frac23(1+\log J)+\frac49\right].
\]

Diese Schranke beweist Konvergenz, zertifiziert aber die vielen numerisch stabilen Stellen nicht. Die im Bericht gezeigten Partialsummentabellen sind Messwerte, keine engen Fehlerintervalle für den Grenzwert.

## 7. Eine zweite analytische Klasse: Newtonkoeffizienten

**Bewiesen.** Koeffizienten |aⱼ|≤C(1+log(j+1))ᵖ/(j+1)ᵅ, α>0, definieren über Σaⱼ binom(x,j) eine analytische Klasse auf Re x>−α. Innerhalb dieser Klasse bestimmen die ganzzahligen Höhenwerte die Funktion, weil die Koeffizienten exakt deren Vorwärtsdifferenzen sind.

Differentiation besitzt die Formel

\[
a_k' =\sum_{m\ge1}\frac{(-1)^{m-1}}m a_{k+m}.
\tag{A17}
\]

Sie folgt aus ∂ₓbinom(x,j)=Σₘ₌₁ʲ(−1)ᵐ⁻¹binom(x,j−m)/m. Die Doppelsumme ist auf kompakten Teilmengen absolut konvergent. Die Ableitung behält α bei und verliert höchstens eine logarithmische Gewichtsstufe.

**Grenze:** Die ganze Klasse ist nicht unter Produkten abgeschlossen. Beispielsweise aⱼ=1/(j+1) erzeugt f(x)=(2ˣ⁺¹−1)/(x+1); sein Quadrat wächst an ganzen Höhen wie 4ᵏ/k² und besitzt keine nur polynomial beschränkten Newtonkoeffizienten. Die Halbebenenskala aus §2 behandelt Produkte; die Newtonklasse aus §7 behandelt die intrinsische Fortsetzung. Eine Zugehörigkeit der neuen q₁/₂-Funktion zur globalen Halbebenennormklasse aus §2 wird hier noch nicht bewiesen. Die volle Ranggruppe ist nicht durch bloßes Zusammensetzen dieser beiden Aussagen bewiesen.

## 8. Schärfere Struktur und offene Beweisfrage

**Skizze mit benannter Lücke.** Die Singularität von (A10) bei w=−1 sollte

\[
\mathcal A_{1/2}(w)=(1+w)^{1/2}
\left[\frac1\pi\log^2(1+w)+B\log(1+w)+C+\cdots\right]
\]

besitzen. Die Jacobi-Parameterableitung liefert den log²-Term. Für den log-Term ergibt die normierte Rangsumme formal

\[
B=\frac2\pi[H_{5/2}+\beta(7/2)],\qquad
\beta(z)=\tfrac12[\psi((z+1)/2)-\psi(z/2)].
\]

Damit wäre αⱼ=(−1)ʲj⁻³/² mal ein quadratisches Polynom in log j plus einen kleineren Rest. Der führende Koeffizient ist −1/(2π√π); der nächste ist 0.33751590910828857. **Noch nicht bewiesen ist hier eine gleichmäßige komplexe Randentwicklung der Rangsumme W und ein kontrollierter Transferrest.** Die bewiesene Holomorphie aus §6 benötigt diese Lücke nicht. Eine größere Halbebene Re x>−3/2 und schärfere Restintervalle dürfen erst danach als Satz bezeichnet werden.

## 9. Bezug zu H₃.₅ und die globale Wachstumsgrenze

Als analytischer Testjet steht jetzt auf Re x>−1

\[
x+1+\varepsilon\,2(B_{5/2}(x)-B_{3/2}(x))
+\varepsilon^2q_{1/2}(x)
\]

zur Verfügung. An ganzen Höhen stimmt er mit der kanonischen Rangfortsetzung des quadratischen Jets überein. Dies ist ein belastbarer Test für eine künftige Ranggruppe, keine vollständige Hyperoperation zur Basis b.

Auch die **lineare** Anwendung auf die echte Exponentialfamilie ist explizit: mit a=1−1/b gilt

\[
d_b(x)=b^x-x-1=\sum_{j\ge1}a^jB_j(x)-B_1(x),
\]

\[
L_s d_b(x)=B_{s+1}(x)\left[a\,{}_2F_1(1,x+s+1;s+2;a)-1\right].
\tag{A18}
\]

Bei s=½ ist das eine analytische erste Rangkorrektur; sie bewahrt die Werte 1 bei Höhe 0 und b bei Höhe 1. Die lokale gleichmäßige Konvergenz der Gamma-Reihe folgt aus dem Quotiententest. Die Euler-Integraldarstellung der hypergeometrischen Funktion mit erstem Parameter 1 gibt zudem endlichen horizontalen Exponentialtyp und polynomialen vertikalen Typ. Es fehlt die bilineare Fortsetzung Q₁/₂ für die Exponentialrichtung und die Konvergenz aller höheren Ordnungen. Eine kritische Basis dieses linearen Ausdrucks würde deshalb nicht η₃.₅ bestimmen.

**Bewiesen:** Eine globale Klasse endlichen horizontalen Exponentialtyps kann den ganzen Rangschritt bei b>exp(1/e) nicht aufnehmen. Die ganzzahlige Tetrationsfolge yₖ₊₁=bʸᵏ hat dann keinen positiven Fixpunkt, wächst gegen ∞ und erfüllt schließlich yₖ₊₁≥R yₖ für jedes feste R. Somit log yₖ/k→∞. Kein Gewicht eᵃᵏ(1+k)ᵖ reicht aus. Der Fitwert 1.5568172767405143 liegt bereits in diesem Bereich.

Der nächste vollständige Definitionsraum muss daher variable Wachstumsgewichte oder lokale holomorphe Bereiche entlang endlicher Höhen zulassen. Zusätzlich braucht er eine konkrete Auswahlregel für den Rangfluss. Die Eindeutigkeit nur anhand ganzzahliger Höhen folgt in einem beliebig großen Wachstumsraum nicht mehr. Die kritische Basis ist erst nach einer solchen Konstruktion ein sinnvoller Vergleich mit dem Fit.
