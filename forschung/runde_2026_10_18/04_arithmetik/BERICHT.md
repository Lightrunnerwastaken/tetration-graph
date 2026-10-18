# Bericht 04 — Arithmetik der Türme und Hyperoperationen (Runde 18.10.2026)

## Kurz und einfach

Rechnet man Türme wie 7^7^7^… nicht mit gewöhnlichen Zahlen, sondern „nach Teilbarkeit durch eine Primzahl p“ (p-adisch), dann kommt der Turm zur Ruhe.
Diese Runde zeigt, **wohin** er kommt und **wie genau**:

1. **Formel für den Grenzwert.** Der p-adische Grenzwert ist dieselbe Reihe wie beim gewöhnlichen unendlichen Turm,
   $\sum (n+1)^{n-1}(\log a)^n/n!$ (Lambert-Reihe), nur mit dem p-adischen Logarithmus und einer Einheitswurzel als Drehung. Geprüft auf 60 Stellen für 668 Fälle.
2. **Genaue Zahl der „fertigen“ Ziffern.** Wie viele p-adische Endziffern von $a\uparrow\uparrow H$ schon endgültig sind, ist **exakt** $V\cdot(H-t)$ — zwei kleine Zahlen
   $V$ (Geschwindigkeit) und $t$ (Anlaufhöhe), beide elementar bestimmbar. Bewiesen; 49 652 Fälle ohne Abweichung.
3. **Pentation und höher.** Jede Pentation, Hexation … ist ein Turm mit riesiger Höhe; daher haben alle Ränge ab 4 denselben Grenzwert, und die Zahl fertiger Ziffern
   folgt sofort. Beispiel: $3\uparrow\uparrow\uparrow3$ hat genau **7 625 597 484 986** endgültige Dezimal-Endziffern.
4. **Die nächsten Ziffern folgen einem Gesetz.** Die erste noch nicht fertige Ziffernreihe von $a\uparrow\uparrow H$ ist $\kappa\cdot u^H$ — eine geometrische Folge mit einer
   neuen Konstante $\kappa_p(a)$ (p-adische Koenigs-Konstante). Bei der Pentation wird dieser Block sogar **konstant**: dort steht im Exponenten der Turmgrenzwert selbst.
5. **Keine p-adische Höhen-Interpolation ab Rang 4.** $a+n$, $a\cdot n$, $a^n$ lassen sich p-adisch stetig in $n$ fortsetzen, Tetration und alle höheren Ränge **nie** — weil sie konvergieren.
6. **Exponentielle Fakultät** $n^{(n-1)^{\cdots^1}}$: modulo $m$ schliesslich periodisch, und die kleinste Periode ist **genau** $\tilde\Lambda(m)$, das kgV der Carmichael-Kette
   von $m$ mit auf 2 gekapptem Zweieranteil. Bewiesen, an 850 Moduln bestätigt. Dazu eine profinite Fortsetzung $\hat E$ mit $\hat E(0)=0$, $\hat E(-1)=1$.

## Befunde

| # | Aussage | Status | Beleg / Genauigkeit | Neuheit-Vermutung |
|---|---|---|---|---|
| B1 | Exakte Bewertung: $p$ ungerade, $p\nmid a$: $v_p(a\uparrow\uparrow H-\alpha_p(a))=V_p(a)(H-t_p(a))$ für $H\ge t_p(a)$, mit $V_p(a)=v_p(a^{p-1}-1)$, $t_p(a)$ = Stabilisierungshöhe von $a\uparrow\uparrow n \bmod \operatorname{ord}_p(a)$ (Start $a\uparrow\uparrow(-1)=0$). $p=2$, $a$ ungerade: $v_2(\cdot)=v_2(a-1)+H(v_2(a^2-1)-1)$ für alle $H\ge0$. $p\mid a$: $v_p(a\uparrow\uparrow H)=v_p(a)\,a\uparrow\uparrow(H-1)$ | **bewiesen** | E1: 49 652 Fälle ($p\le13$, $a<400$, $H<30$, 40 Stellen), 0 Abw.; E3a: $H=100,300$ | Kern (konstante Geschwindigkeit $V$) bekannt (Ripà–Onnis); exakte Form mit Anlaufhöhe $t_p(a)$ als Ordnungs-Stabilisierung vermutlich neu |
| B2 | Für $H<t_p(a)$ ist **keine** Formel möglich: Summen von Einheiten können zufällig durch $p$ teilbar sein (506 Fälle in E1, z. B. $a=2,p=3,H=0$) | bewiesen (Gegenbeispiele) | E1 erste Fassung | — |
| B3 | Dezimal: Zahl endgültiger Endziffern $D(a\uparrow\uparrow M)=\min(v_2\text{-Formel},v_5\text{-Formel})$ | bewiesen (aus B1) | E3c: $a<60$, $M<40$, 0 Abw. | bekannt im Kern |
| B4 | Alle Ränge $r\ge4$ haben denselben p-adischen Grenzwert $\alpha_p(a)$; $H_r(a,n)=a\uparrow\uparrow M_r(a,n)$, also $v_p(H_r(a,n)-\alpha)=V(M-t)$. Pentation: $M=a\uparrow\uparrow\uparrow(n-1)$; z. B. $D(3\uparrow\uparrow\uparrow3)=7625597484986$, $D(2\uparrow\uparrow\uparrow4)=65534$ | bewiesen | Folgerung aus B1; Tabelle E3c | Formulierung für höhere Ränge vermutlich neu (elementar) |
| B5 | Lambert-Formel: $\alpha_p(a)=u\sum_{n\ge0}(n+1)^{n-1}(u\log_pa)^n/n!$, $u=\omega(a)^c$, $c$ = stabiler Rest des Turms mod $p-1$ ($p=2$: $u=\pm1\equiv a \bmod 4$) | **bewiesen** | E2a: 668 Fälle, 60 Stellen, $p\in\{2,3,5,7,11,13,101\}$; Gegenprobe 25 Stellen | Verfeinert Shapiro–Shapiro Prop. 5.4; explizite Reihe nicht gefunden → vermutlich neu |
| B6 | $a\mapsto\alpha_p(a)$ ist auf jeder Klasse (gleiches $a \bmod p$, gleiches $c$) lokal analytisch und eine **Isometrie**: $v_p(\alpha(a')-\alpha(a))=v_p(a'-a)$; insbesondere nirgends lokal konstant, injektiv auf Klassen | bewiesen | E2b: 2 992 Fälle, 0 Abw. | vermutlich neu |
| B7 | Multiplikator am Fixpunkt: $\lambda=\alpha\log_pa=\log_p\alpha$, $v_p(\lambda)=V_p(a)$ — die Kongruenzgeschwindigkeit **ist** die Bewertung des Multiplikators | bewiesen | E2c (Kontrolle beider Ausdrücke, 120 Stellen) | Deutung vermutlich neu |
| B8 | p-adische Koenigs-Konstante $\kappa_p(a)=\lim(a\uparrow\uparrow n-\alpha)/\lambda^n\neq0$ existiert; $v_p(\kappa)=-V t_p(a)$ ($p$ ungerade), $v_2(\kappa)=v_2(a-1)$; Konvergenz $v(\kappa_{n+1}-\kappa_n)=Vn+c$ | bewiesen (Koenigs-Satz im Nichtarchimedischen) + numerisch | E2c: 11 Paare, 120 Stellen | neue Invariante (vermutlich neu) |
| B9 | Ziffern-Gesetz: frischer Block $F_j(H)=((a\uparrow\uparrow H-\alpha)/p^{v})\bmod p^j\equiv\mathrm{unit}(\kappa)\,u^H$, $u=\mathrm{unit}(\lambda)$; also periodisch in $H$ mit Periode $\operatorname{ord}(u \bmod p^j)$ | bewiesen (aus B8) + numerisch | E3b, E6a: 7 Paare, je ≥67 Höhen, 0 Abw. | vermutlich neu |
| B10 | Pentation: der frische Block von $a\uparrow\uparrow\uparrow n$ ist für grosse $n$ **konstant** $=\mathrm{unit}(\kappa)\,u^{\hat\alpha(a)}$, $\hat\alpha(a)\in\hat{\mathbb Z}$ profiniter Turmgrenzwert (Tetration: geometrisch; Pentation: eingefroren) | bewiesen (B9 + Stabilisierung von $M_n$ mod jedem Modul) | E6b nennt $\beta_j$ (z. B. $p=5,a=2$: $\beta_4=357$); direkte Prüfung unmöglich (Höhen zu gross) | vermutlich neu |
| B11 | Keine p-adisch stetige Höhen-Interpolation für Rang $\ge4$, auf **keiner** arithmetischen Progression; Ränge 1–3 interpolieren ($a^n$ für $a\equiv1 \bmod p$, sonst klassenweise mod $p-1$) | bewiesen | elementar (2 Zeilen) | Dichotomie-Formulierung vermutlich neu, Inhalt elementar |
| B12 | Profiniter Turmgrenzwert $\hat\alpha:\hat{\mathbb Z}\to\hat{\mathbb Z}$ stetig; $\hat\alpha(0)=0$, $\hat\alpha(1)=1$, $\hat\alpha(-1)=-1$ | bewiesen | E6c: drei Moduln bis $2.5\cdot10^8$ | Periode in $a$ ist OEIS A127699 (bekannt); Werte trivial |
| B13 | **Exponentielle Fakultät:** $EF(n) \bmod m$ schliesslich periodisch mit kleinster Periode $P(m)=\tilde\Lambda(m)=\operatorname{kgV}_{j\ge0}\widetilde{\lambda^j(m)}$, $\tilde x$ = $x$ mit Zweieranteil auf höchstens $2^1$ gekappt; äquivalent $P(m)=\operatorname{kgV}(\tilde m,P(\lambda(m)))$ | **bewiesen** | E4b: alle 218 $m<260$ und 632 $m\in[260,1500)$ ($\tilde\Lambda\le3\cdot10^4$), 0 Abw. | in OEIS nicht gefunden; Gegenstück A127699 hat **keine** Kappung → vermutlich neu |
| B14 | $EF(n)\bmod 2^k=[n \text{ ungerade}]$ für grosse $n$ (Periode 2) | bewiesen | E4, E4c | elementar |
| B15 | Vorperiode: $EF(n)\bmod m$ periodisch ab $N(m)\le6$ für alle $m<400$; $N(2^k)=6$ für $21\le k$ (bis $\approx EF(5)$), also $N$ unbeschränkt, aber wie eine Umkehrung von $EF$ wachsend | bewiesen für $2^k$, numerisch sonst | E4b, E4c | — |
| B16 | Profinite exponentielle Fakultät $\hat E:\hat{\mathbb Z}\to\hat{\mathbb Z}$, stetig, $\hat E(x)=x\star\hat E(x-1)$; $\hat E(0)=0$, $\hat E(1)=1$, $\hat E(-1)=1$, $\hat E(2)=(0;2,2,2,\dots)$ (2-adisch 0, sonst 2), $\hat E_2(x)=[x\text{ ungerade}]$ | bewiesen | E4c-Tabelle für $2^{20},3^{12},5^8,7^6,10^6$ | vermutlich neu (profinite Antwort auf TPID 2 „EF(0)“: 0) |
| B17 | Stabilisierungshöhe $T(m)=\max_a t(a,m)$ ist meist die Carmichael-Kettenlänge $\ell(m)$: für $m\le2000$: $T=\ell$ in 1965, $T=\ell-1$ in 21, $T=\ell+1$ in 13 Fällen (alle $8\mid m$, kleine Basen) | numerisch | E5 ($a<2m+60$) | wohl bekannt; nicht weiter verfolgt |

## Beweise und Details

### B1 — exakte Bewertung
*Werkzeug (LTE).* $p$ ungerade, $p\nmid a$, $o=\operatorname{ord}_p(a)$, $x\ge0$ ganz: $v_p(a^x-1)=0$ für $o\nmid x$ und $=v_p(a^o-1)+v_p(x)$ für $o\mid x$;
$v_p(a^o-1)=v_p(a^{p-1}-1)=:V$, da $(p-1)/o$ zu $p$ teilerfremd ist.

Setze $D_H=a\uparrow\uparrow(H+1)-a\uparrow\uparrow H$ ($D_{-1}=1$) und $d_H=v_p(D_H)$. Wegen $D_H=a^{a\uparrow\uparrow(H-1)}(a^{D_{H-1}}-1)$ ist
$d_H=V+d_{H-1}$, falls $o\mid D_{H-1}$, sonst $d_H=0$. Sei $t$ die kleinste Höhe $\ge-1$, ab der $a\uparrow\uparrow n \bmod o$ konstant ist.
Für $H\ge t+1$ gilt $o\mid D_{H-1}$. Für $t\ge0$ ist $o\nmid D_{t-1}$ (Minimalität), also $d_t=0$; für $t=-1$ ($o=1$) ist $d_0=v_p(a-1)=V$.
In beiden Fällen $d_H=V(H-t)$ für $H\ge t$. Da $\alpha-a\uparrow\uparrow H=\sum_{j\ge H}D_j$ mit streng wachsenden Bewertungen ist, folgt $v_p(a\uparrow\uparrow H-\alpha)=d_H$. ∎

$p=2$, $a$ ungerade: $v_2(a^x-1)=v_2(a-1)$ für ungerades $x$, $=v_2(a^2-1)+v_2(x)-1$ für gerades $x$. $D_{-1}=1$ ist ungerade, alle $D_H$ ($H\ge0$) gerade; daraus $d_0=v_2(a-1)$, $d_{H}=d_{H-1}+v_2(a^2-1)-1$. ∎
Für $H<t$ gibt es keine Formel (B2): dort sind mehrere $D_j$ Einheiten, ihre Summe kann durch $p$ teilbar sein.

### B4 — höhere Ränge
$H_r(a,n)=H_{r-1}(a,H_r(a,n-1))$; induktiv ist jedes $H_r(a,n)$ ($r\ge4$, $n\ge1$) ein Turm $a\uparrow\uparrow M$ mit $M=M_r(a,n)\to\infty$. Damit gilt B1 mit $H=M$.
Beispiel $3\uparrow\uparrow\uparrow3=3\uparrow\uparrow7625597484987$: $v_2=1+2M$, $v_5=M-t_5(3)=M-1$ ($\operatorname{ord}_5 3=4$, $V_5(3)=v_5(80)=1$, Turm mod 4 stabil ab Höhe 1), also $D=M-1$.
Neue fertige Ziffern pro Pentationsstufe: $V_{10}(a)\,(a\uparrow\uparrow\uparrow n-a\uparrow\uparrow\uparrow(n-1))$ (für $10\nmid a$, sobald derselbe Primfaktor das Minimum stellt) — die „Kongruenzgeschwindigkeit“ der Pentation wächst tetrational.

### B5 — Lambert-Formel
$Y(z)=\sum_{n\ge0}(n+1)^{n-1}z^n/n!$ erfüllt formal $Y=e^{zY}$ und konvergiert für $v_p(z)>1/(p-1)$ (da $v_p(n!)\le(n-1)/(p-1)$); dort ist $Y\equiv1$, also $v_p(zY)=v_p(z)$, und die Komposition mit $\exp$ ist zulässig.
Ab der Stabilisierungshöhe mod $p-1$ liegen alle Exponenten in der Klasse $c$, also $a^x=\omega(a)^c\langle a\rangle^x=u\,e^{Lx}=:g(x)$ mit $L=\log_pa$, $v_p(L)\ge1$.
$g$ ist auf $\mathbb Z_p$ eine Kontraktion ($|g(x)-g(y)|=|L||x-y|$) mit einzigem Fixpunkt; der Turm ist ein $g$-Orbit, sein Grenzwert also dieser Fixpunkt. $x=uY(uL)$ erfüllt $x=u\,e^{Lx}$. ∎
($p=2$: $u=\pm1\equiv a\bmod4$, $v_2(L)\ge2$, $v_2(n!)\le n-1$; gleiches Argument.) Reelles Gegenstück: $h(a)=-W(-\ln a)/\ln a$ — dieselbe Reihe.

### B6 — Isometrie
Für $n\ge2$ hat der Koeffizient von $(z^n-z'^n)/(z-z')$ in $Y(z)-Y(z')$ Bewertung $\ge(n-1)-v_p(n!)>0$; also $v_p(Y(z)-Y(z'))=v_p(z-z')$. Gleiche Klasse ⇒ gleiches $u$ und
$L(a)-L(a')=\log_p\langle a/a'\rangle$ mit $v_p=v_p(a-a')$. ∎ Ableitung: $\alpha'(a)=\alpha^2/(a(1-\lambda))$, eine Einheit.

### B7–B10 — Multiplikator, Koenigs-Konstante, Ziffern-Gesetz
$g'(\alpha)=L\,u\,e^{L\alpha}=L\alpha=\log_p\alpha=:\lambda$ (aus $\alpha=a^\alpha$), $v_p(\lambda)=v_p(L)=V$. Da $0<|\lambda|<1$, gibt es (Koenigs, gilt wörtlich in nichtarchimedischen Körpern)
eine analytische Linearisierung $\sigma(g(x))=\lambda\sigma(x)$ nahe $\alpha$ mit $\sigma(x)=x-\alpha+O((x-\alpha)^2)$; daraus $(a\uparrow\uparrow n-\alpha)/\lambda^n\to\kappa\ne0$ mit Fehler $O(|\lambda|^n)$.
$v_p(\kappa)$ folgt aus B1. Der frische Block ist $\mathrm{unit}(\kappa)\,\mathrm{unit}(\lambda)^H \bmod p^j$, sobald der Fehler unter $p^{-j}$ liegt.
Pentation: $M_n=a\uparrow\uparrow\uparrow(n-1)$ ist ein Turm, also mod $\operatorname{ord}(u\bmod p^j)$ schliesslich konstant $=\hat\alpha(a)$ mod dieser Ordnung ⇒ Block konstant. ∎

Zahlenbeispiele (E2c, die letzten Ziffern des Einheitsteils von $\kappa$ mod $p^{12}$): $\kappa_5(2)$: 207392627, $v=-2$; $\kappa_3(2)$: 485993, $v=-1$; $\kappa_2(3)$: 583, $v=1$; $\kappa_7(3)$: 9302931928, $v=-1$.

### B11 — Interpolations-Dichotomie
Sei $r\ge4$, $S=n_0+d\mathbb N$, $F$ stetig auf dem p-adischen Abschluss von $S$ mit $F=H_r(a,\cdot)$ auf $S$. Dann $n_0+dp^k\to n_0$ p-adisch und $\to\infty$ reell, also
$F(n_0)=\lim H_r(a,n_0+dp^k)=\alpha$; ebenso $F(n_0+d)=\alpha$ — Widerspruch zu $H_r(a,n_0)<H_r(a,n_0+d)$. ∎ Für Rang 3 ist $n\mapsto a^n$ für $a\equiv1 \bmod p$ stetig (sogar analytisch), sonst auf Klassen mod $p-1$.
Folgerung: Eine p-adische Tetration mit Höhen in $\mathbb Z_p$ gibt es nicht; auch der Schröder-Weg scheitert, weil $\lambda^s=p^{Vs}u^s$ für $s\notin\mathbb Z$ nicht erklärt ist.

### B13 — Periode der exponentiellen Fakultät
*Obere Schranke.* CRT: $P(m)=\operatorname{kgV}_{p^k\| m}P(p^k)$. $p=2$: B14 ($n$ gerade ⇒ $n^{\text{riesig}}\equiv0$; $n$ ungerade ⇒ $\lambda(2^k)\mid 2^k\mid EF(n-1)$ ⇒ $\equiv1$), also $P(2^k)=2$.
$p$ ungerade: $EF(n)\bmod p^k$ hängt von $n\bmod p^k$ und $EF(n-1)\bmod\lambda(p^k)$ ab, also $P(p^k)\mid\operatorname{kgV}(p^k,P(\lambda(p^k)))$.

*Untere Schranke* ($p$ ungerade). Zuerst (ii), dann (i). (i) $P(\lambda(p^k))\mid P(p^k)$: Sei $d$ eine Periode; nach (ii) ist $p^k\mid d$, also $n+d\equiv n \pmod{p^k}$. Für $q^j\|p-1$ wähle $n$ in Klassen mod $p$ mit $q^j\mid\operatorname{ord}_p n$; aus
$n^{E(n)}\equiv n^{E(n+d)}$ folgt $EF(n-1)\equiv EF(n+d-1)\pmod{q^j}$ auf einer Vereinigung von Klassen mod $p$. Da $P(q^j)$ nur Primteiler $<p$ hat, trifft diese Menge jede Klasse mod $P(q^j)$ (CRT),
also ist $d$ Periode von $EF\bmod q^j$. Die Komponente $p^{k-1}$ ist trivial ($EF\bmod p^{k-1}$ ist Reduktion). (ii) $p^k\mid P(p^k)$ (benutzt nur die obere Schranke): Wäre der $p$-Anteil von $P$ kleiner, so teilte $P$ das unten definierte $D$. Induktiv hat $P(\lambda(p^k))$ genau den $p$-Anteil $p^{k-1}$.
$D:=\operatorname{kgV}(p^{k-1},P(\lambda(p^k)))$ ist keine Periode: für $k=1$ zerstört $p\nmid D$ das Nullstellenmuster ($EF(n)\equiv0\iff p\mid n$); für $k\ge2$ nimm $n\equiv2 \pmod p$,
dann $e=EF(n-1)\equiv1\pmod p$ und $E(n+D)=E(n)$, also $(n+D)^e/n^e=(1+D/n)^e$ mit $v_p((1+D/n)^e-1)=k-1<k$. ∎
Zusammen $P(p^k)=\operatorname{kgV}(p^k,P(\lambda(p^k)))$, und durch Abwickeln $P(m)=\tilde\Lambda(m)$.
Beispiele: $P(5)=10$, $P(23)=2530$, $P(16)=2$, $P(100)=50$. Unterschied zu A127699 (Turm $n\uparrow\uparrow\infty$ mod $m$, Periode $\operatorname{kgV}(m,a(\lambda m))$ ohne Kappung): bei $EF$ ist der Exponent $EF(n-1)$ für ungerade $n$ durch jede Zweierpotenz teilbar.

### B16 — profinite exponentielle Fakultät
$\hat E(x)\bmod m:=EF(n)\bmod m$ für $n\equiv x \pmod{\tilde\Lambda(m)}$, $n\ge N(m)$; verträglich, da $\tilde\Lambda(m)\mid\tilde\Lambda(m')$ für $m\mid m'$.
$\hat E(0)=0$ ($m\mid n$). $\hat E(-1)=1$: $n$ ungerade, $EF(n-1)$ gerade, $n\equiv-1$ mod $p^k$ ⇒ $(-1)^{\text{gerade}}=1$. $\hat E_p(x)=1$ für $x\equiv1 \pmod p$. $\hat E(2)$: Exponent $\equiv\hat E(1)=1 \pmod{\lambda(p^k)}$ ⇒ $2$ in $\mathbb Z_p$ ($p$ ungerade), 0 in $\mathbb Z_2$.
Im Sinn von TPID 2 (ar.frage.ef_fortsetzung): profinit ist „EF(0)“ $=0$, nicht 1 — die profinite Fortsetzung beantwortet die reelle Frage nicht, zeigt aber, dass $\hat E(-1)=\hat E(1)=1$.

## Reproduktion
Python: `C:/Users/janis/miniconda3/python.exe`, nur Ganzzahlen (kein gmpy2 vorhanden). Im Ordner ausführen.
- `ar_lib.py` — Türme mod $m$ (iterativ über die Carmichael-Kette), Stabilisierungshöhe, p-adischer Log, Lambert-Reihe, schnelle Turmfolge `tower_run`.
- `e1_bewertung.py` (≈1 min): Gegenproben `tower_mod` (7 871 + 2 999 Fälle, darunter $2\uparrow\uparrow5=2^{65536}$ exakt und eine unabhängige φ-Implementierung), B1.
- `e2_lambert_koenigs.py` (≈1 min): B5, B6, B7, B8.
- `e3_hoehere_raenge.py` (≈40 s): B1 bei grossen Höhen, B9, B3/B4-Tabelle.
- `e4_expfakultaet.py`, `e4b_perioden.py` (≈2 min), `e4c_vorperiode2.py`: B13–B16. `e5_stabhoehe.py` (≈2.5 min): B17. `e6_block_profinit.py` (3 s): B9, B10, B12.

## Vorschläge für Grapheinträge
- Exakte p-adische Bewertung von $a\uparrow\uparrow H-\alpha_p(a)$: $V_p(a)(H-t_p(a))$ — satz — bewiesen (Kanten: br.satz.ripa, ar.satz.padische_tuerme)
- Alle Hyperoperationen ab Rang 4 haben denselben p-adischen Grenzwert; Ziffernzahl der Pentation — korollar — bewiesen
- p-adischer Turmgrenzwert als Lambert-Reihe in $\log_p a$ — satz — bewiesen (verfeinert ar.satz.padische_tuerme)
- $a\mapsto\alpha_p(a)$ ist klassenweise eine analytische Isometrie — satz — bewiesen
- Kongruenzgeschwindigkeit = Bewertung des Multiplikators $\lambda=\log_p\alpha$ — satz — bewiesen
- p-adische Koenigs-Konstante $\kappa_p(a)$ und Ziffern-Gesetz der frischen Ziffern — satz — bewiesen (mit Koenigs-Satz), numerisch bestätigt
- Frischer Ziffernblock der Pentation ist schliesslich konstant ($\kappa u^{\hat\alpha(a)}$) — satz — bewiesen
- Keine p-adische Höhen-Interpolation ab Rang 4 (Dichotomie zu Rang ≤ 3) — satz — bewiesen
- Periode der exponentiellen Fakultät mod $m$ ist $\tilde\Lambda(m)$ — satz — bewiesen
- Profinite exponentielle Fakultät $\hat E$ mit $\hat E(0)=0$, $\hat E(\pm1)=1$ — konstruktion — bewiesen
- Stabilisierungshöhe meist = Carmichael-Kettenlänge — beobachtung — numerisch

## Offene Fragen
1. Ist $\kappa_p(a)$ irgendwie geschlossen ausdrückbar (etwa durch eine Reihe wie B5)? Ist $\kappa_p(a)$ als Funktion von $a$ klassenweise analytisch (zu erwarten) und gibt es Relationen zwischen $\kappa_p(a)$ für verschiedene $p$?
2. Genaue Formel für die Vorperiode $N(m)$ der exponentiellen Fakultät (Daten: $N\le6$ für $m<1500$; $N(2^k)$ springt bei $k=3,19,21$).
3. Lässt sich die Lambert-Formel auf gemischte Türme $b_1^{b_2^{\cdots}}$ (periodische Basen) übertragen — Fixpunkt eines Produkts von Kontraktionen?
4. Gibt es ein p-adisches Gegenstück zur Kneser-Frage mit Höhen in einem anderen Ring (z. B. $s\in\mathbb Z$ mit $p$-adischer Abhängigkeit von der Basis $a$)? B11 schliesst nur Höhen in $\mathbb Z_p$ aus.
5. Für $H<t_p(a)$: Verteilung der Bewertungen (B2) — zufällig wie $1/p$?
