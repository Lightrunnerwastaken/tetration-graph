# Runde 5 — Etappe 2: nächste Pivotordnung und Schurnormen

**Zwischenstand nach Etappe 2. Nur (H₀), K2 bleibt bei P.**

## 1. Numerischer Befund: 1/n², nicht 1/n

Es werden die lückenlosen hochpräzisen h_n aus Runde 4 benutzt; die Analyse bildet Quotienten neu aus den 60-stelligen h_n, nicht aus den veröffentlichten Gleitkommaquotienten. Fits bei 65–75 Dezimalstellen ändern die Rechengenauigkeit, machen aus einem endlichen Fit aber keinen asymptotischen Satz.

Numerischer Ansatz:
\[
\log(h_n/h_{n-1})=\log\lambda+c_2n^{-2}+c_3n^{-3}+c_4n^{-4}+\cdots.
\]

| Modell | spätes Fitfenster | λ, numerisch | c₂, numerisch | freier 1/n-Koeffizient | Fehler im getrennten Kontrollfenster |
|---|---|---:|---:|---:|---:|
| kubisch | 800…1000 | 1.17379722692 | -0.427262048 | 5.40e-10 | 6.544e-14 |
| Kneser | 450…600 | 1.17373123578 | -0.427263047 | 4.85e-09 | 5.083e-13 |

Der Kontrollfit trainiert kubisch auf 400…700 und wird auf 701…1000 geprüft; Kneser trainiert auf 250…450, Kontrolle 451…600. Der letzte Fehler ist der maximale Fehler von **log(h_n/h_{n−1})**, nicht ein Intervallradius. Alle Fits und die Fensterabhängigkeit stehen in `old_asymptotics.json` und `second_order.json`.

Ein fitfreier Test ist
p_n=log₂[(log r_{2n}−log r_n)/(log r_{4n}−log r_{2n})], r_n=h_n/h_{n−1}.
Kubisch erhält man p_100=2,03077, p_150=2,01990, p_200=2,01470 und p_250=2,01165. Kneser: p_100=2,03104, p_150=2,02009. Das passt zu einer 1/n²-Hauptkorrektur mit weiteren inversen Potenzen; eine exponentielle Hauptkorrektur passt dazu nicht. Ein winziger echter 1/n-Term außerhalb der geprüften Range ist trotzdem nicht ausgeschlossen.

Die späte kubische Schätzung lautet c₂≈−0,42726205; Kneser c₂≈−0,42726305. Die Differenz von ungefähr 10^−6 ist noch mit höheren Korrekturen und den kürzeren Kneser-Fenstern vereinbar. Bei den festen Eckdaten der jet_ε-Familie bleibt c₂ nahe diesem Wert, während λ deutlich wechselt. Variation von μ oder des M-Multiplikators verändert c₂. **Konjektur, nicht bewiesen:** λ ist global; der logarithmische Determinantenterm wird von lokaler Ecken- und Multiplikatorinformation bestimmt. Keine geschlossene lokale Formel für c₂ wurde gefunden.

## 2. Von einer bewiesenen Quotientenabschätzung zu starker Determinantenasymptotik

**Lemma G32 [bewiesene bedingte Implikation].** Sei h_n>0, λ>0, a∈ℝ, C≥0. Für alle n≥n₀ gelte
\[
\left|\log\frac{h_n}{h_{n-1}}-\log\lambda+\frac a{n^2}\right|\le\frac C{n^3}.
\]

Dann existiert H>0 und für n≥max(n₀,1)
\[
\left|\log\frac{h_n}{H\lambda^n}-\frac an\right|\le\frac{|a|+C/2}{n^2}.
\tag{G32}
\]

*Beweis:* Die Folge log h_n−n log λ konvergiert wegen der summierbaren Differenzen. Ihr Grenzwert ist log H. Subtraktion vom Grenzwert gibt −∑_{j>n}(log r_j−log λ). Es gilt |∑_{j>n}j^{−2}−1/n|≤1/n² und ∑_{j>n}j^{−3}≤1/(2n²). Einsetzen liefert G32. ∎

**Korollar G33 [bewiesene bedingte Implikation].** Unter denselben Hypothesen gilt für ein C_D∈ℝ
\[
\log D_n=\frac{n(n+1)}2\log\lambda+(n+1)\log H+a\log n+C_D+E_n,\qquad
|E_n|\le\frac{3|a|/2+C/2}{n}.
\]

*Beweis:* Summiere G32 mit dem endlichen Anfangsteil. Der Schwanz der O(1/j²)-Reste ist ≤(|a|+C/2)/n. Für H_n=∑_{j≤n}1/j gilt |H_n−log n−γ_E|≤1/(2n); kombiniere. ∎

Damit würde der gemessene Wert a≈0,427262 zu einem Faktor **n^{0,427262} in D_n**, aber zu **keiner Potenz n^a in h_n** führen: h_n≈Hλ^n exp(a/n). Dieser Unterschied ist wesentlich. Die aus den Fits gewonnenen Zahlen erfüllen die Hypothese nur auf einem endlichen Bereich; G32–G33 sind für das ursprüngliche Modell noch nicht anwendbar.

## 3. ‖Q_n‖/h_n: Exponentialteil gut sichtbar, Potenz noch nicht stabil

Die gespeicherten Normen legen kubisch eine Exponentialrate γ_Q≈0,0954 nahe. Da log λ≈0,160244 ist, fällt ‖Q_n‖/h_n ungefähr mit exp(−0,0648n), bis auf Potenzen und weitere Faktoren. Ein Fit von log(‖Q_n‖/h_n) mit n, log n, 1, n^{−1}, n^{−2}, n^{−3} auf 700…1000 ergibt etwa
\[
\log(\|Q_n\|/h_n)\approx-0,0648298\,n-0,34330\log n+\text{weitere Terme}.
\]

**Dieser log-n-Exponent ist nicht belastbar als Grenzexponent.** Auf 300…1000 ergibt derselbe Ansatz −0,4704, auf 500…1000 −0,3334; bei Kneser 400…600 etwa −0,5752. Mehr Fitglieder reduzieren den endlichen Rest, ohne die Potenz eindeutig zu bestimmen. Eine Potenz wurde deshalb nicht als Ergebnis festgelegt. Q wird nicht mit einem neuen SVD-Verfahren analysiert.

**Lemma G34 [bewiesene Implikation].** Wenn ‖Q_n‖=C_Qq^nn^β(1+o(1)) und h_n=Hλ^n(1+o(1)), dann
‖Q_n‖/h_n=(C_Q/H)(q/λ)^nn^β(1+o(1)).
Der √n-Vorfaktor von 𝒟 folgt daraus **nur zusätzlich**, wenn die ausgewertete A-Seite die komplementäre Rate (λ/q)^n n^{−β−1/2} trägt.

*Beweis:* Division und die exakte Identität 𝒟_n(w)=‖Q_n‖A_n(w)/h_n. ∎

Dies ist kein Versuch, die Bucht durch eine Koeffizientennorm zu kontrollieren. Die Aussagen über Normen sind diagnostisch; die nötige A-Abschätzung muss auf der ausgewerteten Kurve gelten.

## 4. Was im einfachen Modell bereits vollständig bewiesen ist

G27 aus Etappe 1 hat explizite exponentielle Quotientenfehler und einen expliziten Determinantenschwanz. G29 liefert für das affine reine Randmodell alle Schurgrößen sowie
log|𝒟_n(1)|=−¼log(πn)−1/(16n)+ε_n, |ε_n|≤1/(360n³).
Das ist eine starke Asymptotik mit Konstanten im einfachsten nichttrivialen Modell, aber mit einem anderen Vorfaktor. Eine Störung dieses Satzes zum kubischen Modell kann den √n-Faktor nicht bloß übernehmen.

## 5. Noch fehlender Schritt

Für das ursprüngliche Modell fehlen bewiesene λ, H, a, C und n₀ in G32, außerdem eine starke Q-Normasymptotik und eine starke **funktionale** A-Asymptotik. Numerik mit +300-Bit-Kontrollen beseitigt hier Rundungszweifel, nicht die All-n-Lücke. Etappe 3 untersucht deshalb die vollständigen komplexen Profile und die notwendige Phase.
