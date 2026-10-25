# Bericht 29 — Transzendenz und Irrationalitätsmasse (Runde 25.10.2026)

## Kurz und einfach

Ob ein Kneser-Wert wie $T_e(\tfrac12)=1.6463\ldots$ ein Bruch ist, weiss niemand, und mit heutigen Mitteln kann es auch niemand zeigen. Diese Runde klärt, **was** man
trotzdem beweisen kann. (1) Die Summe $\tfrac12+\tfrac14+\tfrac1{16}+\tfrac1{65536}+\dots$ der Kehrwerte eines Turms ist **transzendent**, ihr Kettenbruch lässt sich ganz
hinschreiben; der Graph nannte das bisher „offen“. (2) Für die Ableitung $\mathrm{sexp}'(0)$ gibt es einfache, sichere Aussagen: $T'(1)=e\cdot T'(0)$, also sind
$T'(0)$ und $T'(1)$ nicht beide algebraisch. Bei ganzer Basis ist höchstens **ein** Wert der ganzen Folge $T'(-1),T'(0),T'(1),\dots$ algebraisch. (3) An
„fast allen“ Stellen sind die Werte transzendent und verhalten sich wie Zufallszahlen; das ist bewiesen, sagt aber nichts über eine bestimmte Stelle wie $\tfrac12$.
(4) Numerisch: 47 Kneser-Werte mit zusammen 21 485 gesicherten Kettenbruchgliedern verhalten sich in jeder Prüfung wie Zufallszahlen, ohne Spur einer Struktur.

Start 13:47:55, Ende etwa 14:15 (Uhr aus `date`).

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| B1 | Für ganzes $b\ge2$ ist $C_b=\sum_{k\ge1}1/(b\uparrow\uparrow k)$ eine **Liouville-Zahl** (Irrationalitätsexponent $\infty$, transzendent). Expliziter Kettenbruch durch Faltung: ist $S_N=[0;a_1..a_n]$ ($n$ gerade), $x_N=b^{e_{N+1}-2e_N}$, so gilt $S_{N+1}=[0;a_1..a_n,\,x_N-1,\,1,\,a_n-1,\,a_{n-1},..,a_1]$. Beispiel $C_2=[0;1,4,2,1,255,3,4,1,\,2^{65504}-1,\,5,3,255,1,2,4,1,\dots]$ | **bewiesen** | Faltungsregel exakt nachgerechnet: $b=2$ bis $N=5$, $b=3,4,5,7$ bis $N=3$, $b=10$ bis $N=2$ (`e2`) | Liouville elementar; Faltungslemma klassisch (van der Poorten–Shallit). **Korrigiert** `gr.satz.kehrwert_turmsumme` („ob $C_b$ transzendent ist, ist offen“, „Ist $C_2$ irrational?“) |
| B2 | Rang-Dichotomie der Kehrwertsummen, ganzes $b\ge2$: $\sum_k 1/H_r(b,k)$ divergiert für $r\le2$, ist **rational** $=1/(b-1)$ für $r=3$ und **Liouville** für jeden Rang $r\ge4$ | **bewiesen** | aus B1 (alle $H_r(b,k)$, $r\ge4$, sind Potenzen von $b$ mit Exponenten $E_{k+1}\ge 2E_k$, $E_{k+1}/E_k\to\infty$) | Formulierung vermutlich neu, Inhalt elementar |
| B3 | Jede differenzierbare reelle Lösung zur Basis $e$ mit $T(0)=1$: $T'(x+1)=T'(x)\,T(x+1)$, also $T'(-1)=T'(0)$, $T'(1)=e\,T'(0)$, $T'(n)=T'(0)\prod_{k=1}^nT(k)$. **$T'(0)$ und $T'(1)$ sind nicht beide algebraisch** | **bewiesen** | Hermite–Lindemann ($e$ transzendent) | Folklore, nicht im Graph |
| B4 | Ganze Basis $b\ge2$, jede differenzierbare Lösung mit $T(n)=b\uparrow\uparrow n$: $T'(n)=T'(-1)(\ln b)^{n+1}\prod_{k=0}^{n}b\uparrow\uparrow k$. Die ganze Folge $(T'(n))_{n\ge-1}$ enthält **höchstens eine algebraische Zahl** und liegt in $\mathbb Q(T'(-1),\ln b)$ (Transzendenzgrad $\le2$) | **bewiesen** | Quotienten sind $(\ln b)^m\cdot$(rational $\ne0$), $\ln b$ transzendent | vermutlich neu (elementar) |
| B5 | Algebraische Basis $b\notin\{0,1\}$: unter $T'(-1),T'(0),T'(1)$ ist höchstens eine algebraisch (Quotienten $\ln b$, $b\ln b$, $b(\ln b)^2$) | **bewiesen** | wie B4 | vermutlich neu (elementar) |
| B6 | Unter Schanuel (Basis $e$): $\{T'(n):n\ge0\}$ enthält höchstens eine algebraische Zahl; ist $T'(0)$ algebraisch, so sind $T'(1),\dots,T'(n)$ algebraisch unabhängig; allgemein $\operatorname{trdeg}\mathbb Q(T'(0),..,T'(n))\ge n$ | **bedingt** (vollständig unter SC) | Quotienten $T'(n)/T'(m)=\prod_{k=m+1}^nT(k)$, Schanuel-Orbit mit $\alpha=1$ | Folgerung aus `ar.satz.schanuel_orbit`, für Ableitungen neu im Graph |
| B7 | Jet-Identitäten (Basis $e$, jede analytische Lösung, $T(h)=\sum a_kh^k$): $\mathrm{slog}'(0)=1/a_1$, $\mathrm{slog}''(0)=(a_1^2-2a_2)/a_1^3$; allgemein liegen alle Taylor-Koeffizienten von $T$ bei $-1$ und von $\mathrm{slog}$ bei $0$ in $\mathbb Q(a_1,a_2,\dots)$, die bei $n\ge1$ in $\mathbb Q(e,e^e,..,T(n))[a_1,a_2,\dots]$ | **bewiesen**, numerisch bestätigt | Kneser-Werte auf 199 Stellen (`e4`) | Identitäten vermutlich bekannt; Folge: in der 50-Konstanten-Suche (`ar.beob.keine_relationen`) waren $\mathrm{slog}'(0)$, $\mathrm{slog}''(0)$ nicht unabhängig von $a_1,a_2$ |
| B8 | Keine arithmetische Information über $\mathrm{sexp}'(0)$ in der reellen Gleichung: für **jedes** $c>0$ gibt es eine reell-analytische streng wachsende Tetration zur Basis $e$ mit $\tilde T(0)=1$, $\tilde T'(0)=c$ | **bewiesen** | Poisson-Kern-Konstruktion (unten) | Gegenstück zu `ar.satz.keine_arith_info`, vermutlich neu formuliert |
| B8' | Gleichzeitig $\tilde T(\tfrac12)=r\in(1,e)$ und $\tilde T'(0)=c>0$ frei wählbar | hergeleitet mit Lücke | Konvexitätsargument skizziert, Positivität der gemischten Kerne nicht ausgeschrieben | — |
| B9 | Logische Form: „Schanuel ⇒ $T(\tfrac12)$ transzendent für **jede** analytische wachsende Lösung“ ist **gleichwertig zu „Schanuel ist falsch“** | **bewiesen** | exotische Lösung mit $\tilde T(\tfrac12)=\tfrac32$ existiert unbedingt | schärft `ar.satz.keine_arith_info` f5 zu einer Äquivalenz |
| B10 | Metrisch: für jede reell-analytische Lösung (jede Basis $b>1$) ist $\{x: x,T(x)\text{ algebraisch abhängig}\}$ abzählbar; für Lebesgue-fast alle $(x_1..x_n)$ sind die $2n$ Zahlen $x_i,T(x_i)$ algebraisch unabhängig (also (KS) mit $2n$ statt $n$); für fast alle $x$ hat $T(x)$ Irrationalitätsexponent 2 und die Chintschin-, Lévy- und Gauss–Kuzmin-Statistik | **bewiesen** | Nullstellenmengen analytischer Funktionen + $C^1$-Diffeomorphismus | Folklore; zeigt, dass Statistik nichts über einzelne Stellen beweist |
| B10' | Dasselbe für fast alle **Basen** $b>\eta$ bei festem $x$ (z. B. $T_b^{\mathrm K}(\tfrac12)$) | hergeleitet mit Lücke | Lücke: reelle Analytizität und Nichtkonstanz von $b\mapsto T^{\mathrm K}_b(x)$ (Paulsen–Cowgill), hier nicht nachgeprüft | — |
| B11 | Fixpunkte: $L=-W(-1)$ ist transzendent (Hermite–Lindemann), ebenso der Multiplikator $\lambda=L$; für algebraisches $b$ ist der nichtreelle Fixpunkt $L_b$ transzendent (Gelfond–Schneider). Der Multiplikator $\lambda_b=L_b\ln b$ ($b\ne e$) ist unter Schanuel transzendent, unbedingt offen | bewiesen / bedingt | $\lambda e^{-\lambda}=\ln b$; Schanuel mit $(\lambda,\ln b)$ | Folklore, nicht im Graph |
| B12 | **Population**: 47 Kneser-Werte ($T_e$ an 26 rationalen Stellen, $\mathrm{slog}_e$ an 12, $T_2$ an 5, $\mathrm{slog}_2$ an 4), 21 485 gesicherte Teilnenner: Gauss–Kuzmin $\chi^2=12.0$ (FG 8, $p=0.15$; Zufalls-Pools 2.1–23.2), Chintschin 2.664 (Pools 2.645–2.718, Soll 2.685), Lévy $1.180\pm0.006$ (Soll 1.1866, $z=-1.06$), grösster Teilnenner 154 488 ($p=0.18$), $\mu_{\max}\le2.231$ (Median 2.120 gegen 2.131 bei Zufallszahlen) | **numerisch** | Intervall-Kettenbruch mit Fehler $10^{-470}$ bei $\ge494$ wahren Stellen; Gegenprobe Funktionalgleichung auf $10^{-519}$; Kontrollen $e$ ($\chi^2=780$) und $\sqrt2$ (2995) werden erkannt | erweitert A1 (ein Wert, 1936 Glieder) auf eine Population |
| B13 | Ableitungs-Konstanten ($\mathrm{sexp}_e'(0)$, $a_2$, $T''(0)$, $\int T$, $0!_2$, $\mathrm{sexp}_2'(0)$) auf 150 Stellen: Kettenbrüche generisch, $\mu_{\max}\le2.15$; $\mathrm{sexp}'(0)=[1;10,1,8,1,2,1,1,2,1,6,3,\dots]$ | numerisch | Gegenprobe gegen 160er-Stufe: gleiche Glieder | — |
| B14 | Für **keine** nichtganze Stelle ist ein Kneser-Wert als irrational bekannt; an ganzen Stellen nur $T(-1)=0$, $T(0)=1$, $T(1)=e$ (mit $\mu(e)=2$); $T(2)=e^e$ offen. Ein unbedingter Irrationalitätsbeweis für $T_e(\tfrac12)$ bräuchte eine Kneser-spezifische Methode (B8, B9) | **negativ** | Literaturstand; B8/B9 zeigen, dass die reelle Gleichung allein nie reicht | bekannt |
| B15 | $C_e=e^{-1}+e^{-e}+e^{-e^e}+\dots$: Transzendenz offen. Der Liouville-Weg aus B1 bräuchte ein Mass der algebraischen Unabhängigkeit von $e,e^e,e^{e^e}$, das niemand hat | negativ / offen | Rest nach drei Gliedern $\approx e^{-3.8\cdot10^6}$ | — |

## Beweise und Details

### B1 — Kehrwertsumme eines Turms
Sei $b\ge2$ ganz, $e_1=1$, $e_{k+1}=b^{e_k}$; dann $b\uparrow\uparrow k=b^{e_k}$ und $C_b=\sum_k b^{-e_k}$. Für $S_N=\sum_{k\le N}b^{-e_k}$ ist
$b^{e_N}S_N=\sum_{k\le N}b^{e_N-e_k}\equiv1\pmod b$, also $S_N=p_N/q_N$ gekürzt mit $q_N=b^{e_N}$. Es gilt
$0<C_b-S_N<2b^{-e_{N+1}}=2q_N^{-e_{N+1}/e_N}$ und $e_{N+1}/e_N=b^{e_N}/e_N\to\infty$. Also ist $C_b$ eine Liouville-Zahl und transzendent. ∎

*Kettenbruch.* Faltungslemma: ist $p/q=[a_0;a_1..a_n]$ mit Nenner $q=q_n$, so ist $p/q+(-1)^n/(xq^2)=[a_0;a_1..a_n,x,-a_n,..,-a_1]$. Normalisieren mit
$[\dots,x,-y]=[\dots,x-1,1,y-1]$ ergibt für gerades $n$ die Regel in B1 (Nullen werden mit den Nachbarn zusammengefasst; bei ungeradem $n$ zuerst
$a_n\to a_n-1,1$). Hier ist $b^{-e_{N+1}}=1/(q_N^2x_N)$ mit $x_N=b^{e_{N+1}-2e_N}$, und $x_N$ ist ganz, weil $b^{e_N}\ge2e_N$. Folge: $C_b$ hat einen vollständig
expliziten Kettenbruch, und die riesigen Glieder $x_N-1$ zeigen den Liouville-Charakter direkt ($C_2$: Glied $2^{65504}-1$ nach $q_8=2^{16}$, Exponent $\approx4096$).
Für $b\ge3$ beginnt er $[0;\,b-1,\,1,\,b^{b-2}-1,\,b,\,b^{b^b-2b}-1,\,1,\,b-1,\,b^{b-2}-1,\,1,\,b-1,\dots]$ (nachgerechnet für $b=3,4,5,7$).

### B2 — Rang-Dichotomie
$r=1,2$: $\sum1/(b+k)$, $\sum1/(bk)$ divergieren. $r=3$: $\sum_{k\ge1}b^{-k}=1/(b-1)$. $r\ge4$: induktiv ist $H_r(b,k)=b\uparrow\uparrow M_k$ mit streng wachsendem
$M_k$, also $1/H_r(b,k)=b^{-E_k}$ mit $E_k=e_{M_k}$, $E_{k+1}\ge b^{E_k}$; derselbe Beweis wie B1 (auch die Faltung, da $E_{k+1}\ge2E_k$). ∎
Bemerkenswert: der Sprung von „rational“ zu „Liouville“ liegt genau zwischen Potenz und Tetration.

### B3–B6 — Ableitungen
Ableiten von $T(x+1)=b^{T(x)}$ gibt $T'(x+1)=\ln b\,T(x+1)\,T'(x)$. Aus $T(0)=1$ folgt (reell) $T(-1)=0$.
*Basis $e$:* $T'(0)=T'(-1)\cdot T(0)=T'(-1)$, $T'(1)=T(1)T'(0)=eT'(0)$. Wären $T'(0),T'(1)$ algebraisch (und $\ne0$), so wäre $e$ algebraisch. ∎ (B3)
*Ganze Basis:* $T'(n)/T'(m)=(\ln b)^{n-m}\prod_{k=m+1}^nb\uparrow\uparrow k$ mit rationalem, positivem Produkt; algebraisch wäre das nur, wenn $\ln b$ algebraisch
wäre — nach Hermite–Lindemann nicht. Also höchstens ein algebraisches Glied. ∎ (B4) *Algebraische Basis:* die drei Quotienten $\ln b$, $b\ln b$, $b(\ln b)^2$ sind
transzendent. ∎ (B5)
*Unter Schanuel (B6):* $T(k)=\exp^{\circ k}(1)$; nach `ar.satz.schanuel_orbit` (mit $\alpha=1$) sind $T(1),\dots,T(n)$ algebraisch unabhängig. Für $0\le m<n$ ist
$T'(n)/T'(m)=\prod_{k=m+1}^nT(k)$ ein Produkt algebraisch unabhängiger Zahlen, also transzendent. Ist $T'(0)$ algebraisch, so erzeugen $T'(1..n)$ über $\bar{\mathbb Q}$
denselben Körper wie $T(1..n)$ (Transzendenzgrad $n$), also sind sie algebraisch unabhängig. ∎
Gegensatz: bei ganzer Basis hat die ganze Ableitungsfolge Transzendenzgrad $\le2$, bei Basis $e$ wächst er (unter SC) unbeschränkt.

### B7 — Jet-Identitäten
Basis $e$: $T(-1+h)=\log T(h)=\log(1+a_1h+a_2h^2+\dots)$; die Koeffizienten sind Polynome in $a_k$ mit **rationalen** Koeffizienten
($T'(-1)=a_1$, $T''(-1)=2a_2-a_1^2$). $\mathrm{slog}$ ist die Umkehrung bei $-1$: $\mathrm{slog}'(0)=1/T'(-1)=1/a_1$,
$\mathrm{slog}''(0)=-T''(-1)/T'(-1)^3=(a_1^2-2a_2)/a_1^3$. Bei $n\ge1$: $T(n+h)=\exp(T(n-1+h))$ gibt $T(n)\cdot$Polynom. ∎
Numerisch (Kneser, 200 Stellen): alle drei Differenzen $<10^{-199}$ (`e4_jet_identitaeten.py`). Der Kettenbruch von $\mathrm{slog}'(0)$ ist der um eine Stelle
verschobene von $\mathrm{sexp}'(0)$ (`e3`).

### B8 — $\mathrm{sexp}'(0)$ ist durch die reelle Gleichung nicht festgelegt
$T$ Kneser, $\gamma=c/T'(0)>0$. Poisson-Kern $P_r(x)=(1-r^2)/(1-2r\cos2\pi x+r^2)$: analytisch, 1-periodisch, Mittel 1, Werte in
$[\tfrac{1-r}{1+r},\tfrac{1+r}{1-r}]$. Für $\gamma\ge1$ wähle $r$ mit $\tfrac{1+r}{1-r}\ge\gamma$ und $g=1+t(P_r-1)$, $t=(\gamma-1)/(P_r(0)-1)\in[0,1]$; für $\gamma<1$
$g=1+t(P_r(\cdot-\tfrac12)-1)$, $t=(1-\gamma)/(1-\tfrac{1-r}{1+r})\in(0,1]$. In beiden Fällen ist $g>0$ (Mischung von $1$ und einem positiven Kern), Mittel 1, $g(0)=\gamma$.
$\varphi(x)=\int_0^xg$ ist analytisch, streng wachsend, $\varphi(x+1)=\varphi(x)+1$, $\varphi(0)=0$. $\tilde T=T\circ\varphi$ erfüllt $\tilde T(x+1)=e^{\tilde T(x)}$, $\tilde T(0)=1$,
$\tilde T'(0)=T'(0)\gamma=c$. ∎

### B9 — Logische Form
Sei $P$: „Für jede reell-analytische streng wachsende Lösung $T$ mit $T(0)=1$ gilt: SC ⇒ $T(\tfrac12)$ transzendent.“ Nach `ar.satz.keine_arith_info` existiert
unbedingt eine solche Lösung mit $T(\tfrac12)=\tfrac32$. Also ist $P$ gleichwertig zu „SC ⇒ falsch“, d. h. zu ¬SC. Jede bedingte Transzendenzaussage über
$T^{\mathrm K}_e(\tfrac12)$ muss eine Kneser-Eigenschaft benutzen, die die exotischen Lösungen nicht haben. ∎

### B10 — metrische Aussagen
$T$ ist als Funktion nicht algebraisch (wächst schneller als jede algebraische Funktion). Für $P\ne0$ in $\mathbb Q[X_1..X_n,Y_1..Y_n]$ ist
$F(x)=P(x,T(x_1),..,T(x_n))$ analytisch und nicht identisch null (Induktion über $n$: Koeffizienten bei festen $x_1..x_{n-1}$ generisch $\ne0$, dann ist
$x_n\mapsto\sum c_{\alpha\beta}x_n^\alpha T(x_n)^\beta\not\equiv0$). Die Nullstellenmenge hat Mass 0 (für $n=1$: diskret); abzählbar viele $P$. Für $n=1$ ist die
Ausnahmemenge abzählbar. Da $T$ lokal ein $C^1$-Diffeomorphismus ist, hat das Urbild jeder Nullmenge (z. B. $\mu>2$, nicht Chintschin-typisch) Mass 0. ∎

### B12 — Population (Verfahren)
`e1_kf_statistik.py`: Kettenbruch des Intervalls $[x-10^{-470},x+10^{-470}]$ (nur gemeinsame Glieder beider Endpunkte zählen). Effektiver Exponent
$\mu_i=1+\ln q_{i+1}/\ln q_i$, Maximum ab $i=20$. Kontrollen: 200 Zufallszahlen (gleiche Länge, gleiche Pipeline) für $\mu_{\max}$ und Lévy, 40 Zufalls-Pools aus je
47 Zahlen für $\chi^2$ und Chintschin. Doppelte Stellen ($0.1=\tfrac1{10}$ usw.) entfernt. Gauss–Kuzmin-Klassen $1,2,3,4,5\text{–}6,7\text{–}10,11\text{–}20,21\text{–}50,>50$:
gezählt 8882, 3672, 2129, 1223, 1522, 1431, 1229, 814, 583; erwartet 8917, 3651, 2001, 1265, 1512, 1442, 1255, 840, 602. $T_e(\tfrac12)$ allein: 454 Glieder,
$\mu_{\max}=2.140$ (A1 fand mit 1936 Gliedern ebenfalls 2.14). Die Statistik **beweist nichts** (B10): sie zeigt nur, dass keine verborgene Struktur
(wie bei $e$ oder quadratischen Irrationalzahlen) in den ersten ~460 Gliedern steckt.

## Vorschläge für Grapheinträge

- `ar.satz.kehrwertsumme_liouville` — „$\sum1/(b\uparrow\uparrow k)$ ist für ganze $b\ge2$ Liouville, mit explizitem Kettenbruch“ — satz — bewiesen (dazu: f4/f5 von `gr.satz.kehrwert_turmsumme` korrigieren)
- `ar.satz.rang_dichotomie_kehrwerte` — „Kehrwertsummen: divergent (Rang ≤2), rational (3), Liouville (≥4)“ — satz — bewiesen
- `ar.satz.ableitungsorbit` — „$T'(1)=eT'(0)$, $T'(-1)=T'(0)$: $T'(0),T'(1)$ nicht beide algebraisch“ — satz — bewiesen
- `ar.satz.ableitungen_ganze_basis` — „Ganze Basis: höchstens ein algebraisches $T'(n)$; alles in $\mathbb Q(T'(-1),\ln b)$“ — satz — bewiesen
- `ar.satz.schanuel_ableitungen` — „Unter Schanuel höchstens ein algebraisches $T'(n)$, $n\ge0$“ — satz — bedingt
- `ar.satz.jet_identitaeten` — „$\mathrm{slog}'(0)=1/\mathrm{sexp}'(0)$, $\mathrm{slog}''(0)=(a_1^2-2a_2)/a_1^3$“ — satz — bewiesen
- `ar.satz.keine_arith_info_ableitung` — „$\mathrm{sexp}'(0)$ ist für exotische Lösungen beliebig“ — satz — bewiesen
- `ar.satz.schanuel_logik` — „Schanuel allein beweist $T(\tfrac12)$ transzendent nur, wenn Schanuel falsch ist“ — bemerkung — bewiesen
- `ar.satz.metrisch` — „Fast alle Stellen: $2n$ algebraisch unabhängige Zahlen, $\mu=2$“ — satz — bewiesen
- `ar.satz.fixpunkt_transzendent` — „$L$ und $L_b$ (algebraisches $b$) transzendent; $\lambda_b$ unter Schanuel“ — satz — bewiesen/bedingt
- `ar.beob.kf_population` — „47 Kneser-Werte, 21 485 Kettenbruchglieder: Gauss–Kuzmin, Chintschin, Lévy, $\mu\approx2$“ — beobachtung — numerisch

## Offene Fragen

1. Ist $C_e=\sum e^{-e\uparrow\uparrow(k-1)}$ transzendent? Ist $C_b$ für rationales, nicht ganzes $b$ irrational?
2. Ist der Multiplikator $\lambda_2=L_2\ln2$ unbedingt transzendent? (Gleichwertig: kein algebraisches $\lambda$ mit $\lambda e^{-\lambda}=\ln2$.)
3. Gibt es eine Kneser-spezifische Eigenschaft (z. B. das Verhalten $T(z)\to L$ für $\operatorname{Im}z\to\infty$), die in eine Transzendenzmethode passt? B9 zeigt, dass ohne sie nichts geht.
4. Ist $\mathrm{sexp}'(0)$ algebraisch unabhängig von $e,e^e,\dots$? (Unter SC folgt nur $\operatorname{trdeg}\ge n$ für $n+1$ Zahlen.)

## Dateien

- `e1_kf_statistik.py` — Population, Kontrollen; Ergebnis `cache/e1_ergebnis.json`
- `e2_kehrwertsumme_liouville.py` — exakte Kettenbrüche von $S_N$, Prüfung der Faltungsregel
- `e3_ableitungen_kf.py` — Kettenbrüche der Ableitungs-Konstanten (150 Stellen, Gegenprobe 160er-Stufe)
- `e4_jet_identitaeten.py` — Jet-Identitäten auf 199 Stellen
- `cache/` — Kopien der Eingabedaten (`te_values_520.json`, `te2_values_520.json` aus `Tetration/Research/arith-half/results`, `konstanten_200.json` aus `forschung/identitaeten_2026_10_14/daten`)
