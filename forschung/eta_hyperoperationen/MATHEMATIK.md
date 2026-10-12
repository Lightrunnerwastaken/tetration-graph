# Mathematik zur η-Folge und zur gewählten Fortsetzung

Stand 12.10.2026. Keine Arbeit am pausierten H₀-Projekt. **Bewiesen** bedeutet hier eine vollständige Herleitung aus den ausdrücklich genannten Voraussetzungen. Die Erfüllung dieser Voraussetzungen durch sämtliche höheren analytischen Fortsetzungen ist damit nicht bewiesen.

## 1. Notation und Auswahl

Für b>1 sei H₃(b,x)=bˣ. Für r≥3 verlangen wir

\[
H_{r+1}(b,0)=1,\qquad
H_{r+1}(b,x+1)=H_r\bigl(b,H_{r+1}(b,x)\bigr).
\]

H₄ ist Tetration, H₅ Pentation, H₆ Hexation, H₇ Heptation, H₈ Octation, H₉ Nonation und H₁₀ Decation. Die **Berührbasis der Abbildung Hᵣ** ist ηᵣ; die zugehörige Höhe ist eᵣ. Die Grenze beschränkter Iteration betrifft die nächste Stufe Hᵣ₊₁. Wer nach der begrenzten Operation indizieren möchte, kann stattdessen Eₛ:=ηₛ₋₁ schreiben. Insbesondere E₄=η₃=e^(1/e).

Eine zusätzliche Auswahl \(\mathscr C\) der reell-analytischen Fortsetzungen wird vorausgesetzt. Ganz präzise wäre daher \(H_r^{\mathscr C}\), \(\eta_r^{\mathscr C}\). Im Bericht bleibt \(\mathscr C\) fest: Kneser-Tetration, reguläre Iteration für die ungeraden Stufen, anziehende reelle Iteration oder experimenteller Zwei-Bogen-Abel-Aufbau für die geraden Stufen.

## 2. Beschränktheit und Fixpunkte — bewiesen unter einfachen Voraussetzungen

**Satz M1.** Sei F:[1,∞)→[1,∞) stetig und monoton nichtfallend, F(1)=b>1. Die Folge y₀=1, yₖ₊₁=F(yₖ) ist genau dann beschränkt, wenn F einen Fixpunkt a≥1 besitzt. Ihr Grenzwert ist dann der kleinste solche Fixpunkt.

**Beweis.** y₁>y₀, und Monotonie liefert induktiv yₖ₊₁≥yₖ. Ein Fixpunkt a≥1 beschränkt die Folge von oben, denn yₖ≤a impliziert F(yₖ)≤F(a)=a. Eine beschränkte monotone Folge konvergiert gegen ℓ; Stetigkeit ergibt F(ℓ)=ℓ. Jeder Fixpunkt a≥1 liegt über allen yₖ, also über ℓ. Damit ist ℓ der kleinste. ∎

Für F=Hᵣ(b,·) ist yₖ=Hᵣ₊₁(b,k). Ein positiver lokaler Berührpunkt ist somit ein plausibler Schwellenkandidat. Dass er die **globale einzige** Grenze in b bildet, braucht weitere Monotonie- und Ausschlussargumente; eine lokale Newton-Rechnung allein liefert diese nicht.

## 3. Monotonie der Schwellen und die Obergrenze 2 — bedingter Satz, vollständig bewiesen

**Satz M2.** Angenommen, für alle betrachteten Ränge und alle b>1 existieren Hᵣ auf [1,∞), sind stetig und monoton in x sowie monoton nichtfallend in b, und erfüllen die obige Rekursion einschließlich Hᵣ(b,1)=b. Setze

\[
\mathcal E_r=\{b>1: H_r(b,\cdot)\text{ hat einen Fixpunkt }\ge1\}.
\]

Dann sind die nichtleeren Mengen \(\mathcal E_r\) Anfangsintervalle in (1,∞),

\[
\mathcal E_{r-1}\subseteq\mathcal E_r\subset(1,2)\qquad(r\ge4).
\]

Folglich besitzen ihre oberen Grenzen ηᵣ einen monotonen Grenzwert L mit η₃≤L≤2. **Die Gleichheit L=2 folgt daraus nicht.**

**Beweis des Einschlusses.** Hat Hᵣ₋₁(b,·) einen Fixpunkt a≥1, dann sind nach M1 alle ganzzahligen Werte Hᵣ(b,k)≤a. Wegen Monotonie in x gilt für jedes x≥1 auch Hᵣ(b,x)≤Hᵣ(b,⌈x⌉)≤a. Dagegen ist Hᵣ(b,1)−1=b−1>0. An einer Stelle x>a ist Hᵣ(b,x)−x<0. Der Zwischenwertsatz liefert einen Fixpunkt von Hᵣ.

**Beweis der Obergrenze.** Bei Basis 2 gilt H₃(2,k)=2ᵏ≥k+1 für alle ganzen k≥1. Induktiv in r bleiben alle ganzzahligen Iterationen ganzzahlig und erhöhen ihr Argument um mindestens 1. Damit gilt Hᵣ(2,k)≥k+1, und die Folge Hᵣ₊₁(2,k) ist unbeschränkt. Nach M1 hat Hᵣ(2,·) keinen Fixpunkt ≥1. Da Hᵣ(2,1)−1>0 und die Differenz stetig ist, bleibt Hᵣ(2,x)−x überall positiv. Monotonie in b schließt Fixpunkte auch für b>2 aus.

**Intervallform.** Für 1<b₁<b₂ und b₂∈Eᵣ liefert ein Fixpunkt a von Hᵣ(b₂,·) die Ungleichung Hᵣ(b₁,a)≤a. Zusammen mit Hᵣ(b₁,1)>1 folgt per Zwischenwertsatz ein Fixpunkt. ∎

Die Voraussetzungen sind für den Exponentialfall elementar. Für die gesamte hier experimentell gebaute Hierarchie sind sie **nicht bewiesen**. Insbesondere steht keine globale Monotonie in b aller geraden höheren Fortsetzungen zur Verfügung. Der Satz erklärt die erwartete Struktur, ersetzt aber keinen Existenz- und Auswahlsatz.

**Zusatz M2a, exakt für die normierte ganzzahlige Hierarchie bei Basis 2.** Für alle r≥3 und ganzen k≥1 gilt Hᵣ(2,k)≥2k. Beim Exponentialfall ist das 2ᵏ≥2k. Induktiv liefert die Höhenrekursion Hᵣ₊₁(2,k)≥2ᵏ≥2k. Außerdem

\[
H_{r+1}(2,3)=H_r(2,4)\ge2H_r(2,3),
\qquad H_r(2,3)\ge2^r.
\]

Die erste Gleichheit benutzt Hᵣ₊₁(2,2)=4. Für die Ungleichung bei r≥4 schreiben wir Hᵣ(2,4)=Hᵣ₋₁(2,Hᵣ(2,3)) und wenden die bewiesene ganzzahlige Schranke auf das ganzzahlige innere Argument an. Der Anfangsfall ist H₄(2,3)=16=2H₃(2,3). Somit divergieren bereits die Werte bei Höhe 3 mit dem Rang, unabhängig von der gewählten Fortsetzung zwischen ganzzahligen Höhen. Ein endlicher globaler Grenzgraph bei Basis 2 ist deshalb ausgeschlossen. Sinkende Werte bei Höhe 1.5 sagen nichts über die Werte ab Höhe 3.

## 4. Die erste kritische Basis — exakt bewiesen

Für H₃(b,x)=bˣ lautet die Fixpunktgleichung bˣ=x, also log b=(log x)/x. Die rechte Seite besitzt auf x≥1 ihr eindeutiges Maximum 1/e bei x=e. Daher η₃=e^(1/e), e₃=e. Der Multiplikator ist bei der Berührung x log b=1. Zusammen mit M1 ist dies genau die Schranke für beschränkte ganzzahlige Tetration.

## 5. Lokales Verhalten an einer höheren Berührung

Sei F(b,x)=Hᵣ(b,x), und gelte an (η,e)

\[
F(\eta,e)=e,\quad F_x(\eta,e)=1,\quad
\alpha:=F_b(\eta,e)>0,\quad \beta:=\tfrac12F_{xx}(\eta,e)>0.
\]

Unter C³-Regularität gilt lokal

\[
F(b,x)-x=\alpha(b-\eta)+\beta(x-e)^2
 +O\bigl((b-\eta)^2+|b-\eta||x-e|+|x-e|^3\bigr).
\]

**Bewiesen, lokaler Faltenmechanismus:** Zwei nahe Fixpunkte unterhalb der Schwelle haben führend Abstand \(\sqrt{\alpha(\eta-b)/\beta}\) zur Berührstelle. Oberhalb fehlen diese **lokalen** Fixpunkte. Andere Fixpunkte außerhalb der Umgebung werden hierdurch nicht ausgeschlossen.

**Satz M3, kritische Iterationsrate.** Konvergiert die echte Iteration bei b=η von unten gegen e, dann

\[
k\bigl(e-H_{r+1}(\eta,k)\bigr)\longrightarrow\beta^{-1}.
\]

**Beweis.** Mit tₖ=e−yₖ>0 gilt tₖ₊₁=tₖ−βtₖ²+O(tₖ³). Daher 1/tₖ₊₁−1/tₖ→β. Teleskopieren und Cesàro ergeben 1/(k tₖ)→β. ∎

Dies erklärt, weshalb eine endliche, langsam wachsende Iteration knapp oberhalb η keine verlässliche Bestimmung der Schwelle ist. Die häufig verwendete Engpasszeit π/√(αβ(b−η)) ist hier nur die lokale Differentialgleichungsheuristik; ein uniformer Fehlerbeweis dafür wurde nicht geführt.

## 6. Rekursion und ganzzahlige Höhen bestimmen die Fortsetzung nicht

**Satz M4, bewiesen.** Erfüllt S(x+1)=F(S(x)), S(0)=1, dann erfüllt auch

\[
\widetilde S_\varepsilon(x)=S\bigl(x+\varepsilon\sin(2\pi x)\bigr)
\]

diese Rekursion und dieselben Werte an sämtlichen ganzzahligen Höhen. Für |ε|<1/(2π) ist der neue Höhenparameter reell-analytisch und streng monoton, weil seine Ableitung 1+2πε cos(2πx)>0 ist. Das folgt unmittelbar aus ψ(x+1)=ψ(x)+1 und ψ(k)=k.

Dies gilt bei **festem unteren F**. Wird S danach selbst als untere Abbildung der nächsten Stufe verwendet, können deren Werte auch an ganzzahligen Höhen für **nichtganzzahlige Basen** wechseln: schon Hᵣ₊₁(b,2)=Hᵣ(b,b) benutzt die Höhe b. Bei ganzzahligen Basen bleiben die klassischen ganzzahligen Hyperoperationswerte erhalten. Es wird somit nicht behauptet, ein Wechsel des unteren S erhalte automatisch die gesamte alte Hierarchie bei reeller Basis.

Ist S(b,x)=Hᵣ(b,x) und besitzt seine Berührung Faltennichtentartung α≠0, Fₓₓ≠0, liefert der implizite Funktionensatz für die Berührbasis der veränderten Abbildung

\[
\left.\frac{d\eta_r(\varepsilon)}{d\varepsilon}\right|_{0}
=-\frac{\sin(2\pi e_r)}{\partial_bH_r(\eta_r,e_r)}.
\]

**Herleitung:** Die erste Berührgleichung nach ε ableiten. Der Koeffizient von e′ verschwindet, weil Fₓ−1=0. Übrig bleiben αη′+sin(2πe)=0. Die zweite Gleichung bestimmt e′. ∎

Für H₄ wurde α≈9.87526358547 numerisch gemessen, daraus η₄′(0)≈−0.0534687218276. Das ist eine konkrete Auswahlabhängigkeit. Eine solche Parameteränderung bewahrt im Allgemeinen zusätzliche globale Kneser-Auswahlbedingungen nicht. Deshalb sind die Schwellen für eine fixierte Fortsetzung durchaus sinnvolle Forschungsobjekte.

## 7. Konstruktion über lokale Konjugation — bedingtes analytisches Prinzip

Sei f=Hᵣ₋₁(b,·), p ein nichtneutraler Fixpunkt und λ=f′(p)>0, λ≠1. Ein lokaler Koenigs-Parameter σ mit σ(p)=0, σ′(p)=1 erfüllt σ(f(z))=λσ(z). Mit lokaler Umkehrung Ψ und geeigneter Normierung c entsteht

\[
H_r(b,x)=p+\Psi(c\lambda^x),\qquad H_r(b,0)=1.
\]

Die Formel muss außerhalb der lokalen Scheibe entlang einer festgelegten Vorwärts- oder Rückwärtskette analytisch fortgesetzt werden. An einem reellen abstoßenden p nimmt man lokal die anziehende Umkehrung von f; bei einem anziehenden positiven p nimmt man f selbst. Die Rekursion folgt dort exakt aus der Konjugationsgleichung. Die **global verwendbaren Zweige und Gebiete** sind eine zusätzliche Voraussetzung.

Das Programm bestimmt lokale Reihen durch Cauchy-Abtastung, löst die Schröder-Gleichung, normiert an 1 und vergleicht zwei zusätzliche Fortsetzungsschritte. Diese numerischen Kontrollen sind kein analytischer Konvergenzbeweis.

Für die unbeschränkte H₈ wird stattdessen um einen komplexen abstoßenden Fixpunkt q von H₇ eine lokale Abel-Funktion log σ/log λ gebildet und durch eine periodische Fourier-Korrektur mit einer reellen Abel-Reihe gekoppelt. Zwei Bogenlösungen werden verglichen. Die inverse normierte Abel-Funktion liefert H₈. Das ist eine **experimentelle Auswahl**, kein bewiesener allgemeiner Kneser-Satz für jede Stufe.

## 8. Ein einheitlicher Abel-Ausdruck für die Berührung — bewiesen lokal

Ist Bᵣ die lokale inverse Höhenfunktion von Hᵣ, so gilt Bᵣ(1)=0 und

\[
B_r\bigl(H_{r-1}(b,z)\bigr)=B_r(z)+1.
\]

Eine Berührung Hᵣ(b,e)=e, Hᵣ′(b,e)=1 ist genau gleichbedeutend mit Bᵣ(e)=e, Bᵣ′(e)=1. Das folgt aus dem lokalen Umkehrsatz und der Ableitungsregel für inverse Funktionen. Die kritischen Basen lassen sich deshalb einheitlich als Diagonalberührungen von **Hyperoperation oder ihrer Abel-Funktion** beschreiben.

Für die regulären Stufen mit Hᵣ=p+Ψ(cλˣ) und global entlang des gewählten Weges fortgesetztem Schröder-Parameter σ lauten die beiden Gleichungen äquivalent

\[
\log(\sigma(e)/c)=e\log\lambda,\qquad
\sigma'(e)=\sigma(e)\log\lambda.
\]

Die verwendeten Logarithmuszweige müssen zur reellen Normierung passen. Dieser Ausdruck enthält weiterhin globale Daten der vorherigen Stufe; er ist keine geschlossene Zahlenrekursion in ηᵣ. Er bietet aber einen möglichen numerischen Weg, extreme Vorwärtsiterationen zu vermeiden.

## 9. Reelle Koenigs-Approximation für eine teure nächste Stufe

Sei f=Hᵣ₋₁ auf dem betrachteten reellen Intervall analytisch, a ein anziehender Fixpunkt und 0<λ=f′(a)<1. Angenommen, die Iterationen von 1 und von allen verwendeten Anfangspunkten gelangen in dieselbe lokale Konjugationsumgebung. Dann gilt mit einem geeigneten Koenigs-Parameter κ

\[
f^k(z)=a+\lambda^k\kappa(z)+O(\lambda^{2k}).
\]

**Bewiesen unter diesen Voraussetzungen:** Die Lösung zₖ(t) für 0<t<1 der reellen Gleichung

\[
f^k(z_k(t))-a=\lambda^t\bigl(f^k(1)-a\bigr)
\]

konvergiert gegen die reguläre Fortsetzung Hᵣ(b,t), mit zₖ(t)=Hᵣ(b,t)+C(t)λᵏ+O(λ²ᵏ), sofern κ′ dort nicht verschwindet. Teilen durch λᵏ, Einsetzen der obigen Entwicklung und der lokale Umkehrsatz ergeben diese Aussage. Die Konstanten sind in dieser Untersuchung nicht uniform in r oder b kontrolliert.

Die extrapolierte Kombination (zₖ₊₂−λ²zₖ)/(1−λ²) hebt den führenden Fehler auf. In exakter Arithmetik ist ihr Fehler O(λ²ᵏ). Numerisch darf k nicht so groß werden, dass die Differenz fᵏ(1)−a im Rundungs- und Konstruktionsfehler verschwindet. Daher sind unterschiedliche Iterationslängen und ausreichend hohe Präzision erforderlich.

Dieser Ansatz wird für reelle Decation bei Basis 1.83 ausprobiert. Er vermeidet die komplexe Cauchy-Abtastung der verschachtelten Nonation. Er liefert **keinen Nachweis einer globalen komplexen Decation** und erschließt insbesondere nicht automatisch den negativen Fixpunkt, der für eine reguläre H₁₁ benötigt würde.

## 10. Was zur Vermutung noch fehlt

M2 liefert unter globalen Voraussetzungen nur die Existenz eines Grenzwerts L≤2. Für L=2 wäre etwa für **jede** feste Basis 1<b<2 zu beweisen, dass ab einem endlichen Rang Hᵣ(b,·) einen Fixpunkt ≥1 besitzt. Gleichwertig: solche Basen müssen schließlich in allen Eᵣ liegen. Weder endlich viele η-Werte noch ein Potenzfit leisten das.

Die Normierungen Hᵣ(2,1)=2 und Hᵣ(2,2)=4 für alle r sind exakt, sofern die jeweilige Hierarchie existiert. Sie machen 2 und 4 zu natürlichen Bezugspunkten, beweisen aber weder ηᵣ→2 noch eᵣ→4. Eine geschlossene Rekursion nur für (ηᵣ,eᵣ) wurde nicht gefunden; die Zwischenwerte der vorherigen Fortsetzung gehen wesentlich ein.
