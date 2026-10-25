# Eintrag Bericht 15 „Monodromie des regulären Erzeugers J_L“ (25.10.2026)

Dateien: `daten/aussagen/r25_15_it.yaml`, `daten/aussagen/r25_15_kn.yaml`, `daten/theorie/r25_15_it.yaml`, `daten/theorie/r25_15_kn.yaml`.
`python baue.py --nur-pruefen`: 0 Fehler, keine Warnung aus diesen Dateien.

## Angelegte ids

| id | Titel kurz |
|---|---|
| it.satz.erzeuger_riemannflaeche | Riemannfläche von J_L = u-Ebene, Ψ Überlagerung über ℂ∖P (bewiesen, M1) |
| it.satz.erzeuger_kranzrekursion | Blätter ↔ ℤ^(∞), a verschiebt m₁, b_k verschiebt m_{k+1} nur bei m₁…m_k = 0 (bewiesen, M2, M9) |
| it.kor.erzeuger_verzweigung_logarithmisch | jedes exp^k(0) auf dem Hauptblatt logarithmisch, keine algebraischen Verzweigungen (bewiesen, M3) |
| it.satz.erzeuger_monodromiegruppe | G = IMG(exp) = ℤ≀ℤ≀…, nicht frei, nicht endlich erzeugt (bewiesen, M4, M4a, M4b) |
| it.kor.erzeuger_umlaufwerte | nach b_k^j: J_L = L^{k+1} J^{(j)} (bewiesen, Werte numerisch, M5) |
| it.kor.kneser_erzeuger_sigma | J_K = J_L · Q(σ^{2πi/L}) (bewiesen, M6) |
| it.satz.regulaere_slog_riemannflaeche | Riemannfläche der regulären slog, Monodromie ⊃ ⊕ℤ (Skizze, M8) |
| it.satz.erzeuger_monodromie_basen | dieselbe Gruppe für jede Basis b > η (Skizze, M10) |
| kn.kor.theta_verzweigung_unendlich | Knesers θ: Verzweigungspunkt unendlicher Ordnung an jeder ganzen Zahl (bewiesen, M7, M7a) |

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| M1 | Beweis gelesen: keine kritischen Punkte (kn.satz.theta_singularitaet (i)); Wegehebung für beliebige Keime über $\ell_n\circ\gamma$ (geprüft, dass dieser Weg $P$ meidet, weil $\exp^n(p_j)=p_{j+n}$); Überlagerungskriterium (Forster 4.19); Keim-Argument $\Psi(cu)=\Psi(u)\Rightarrow c=1$ | richtig, bewiesen. |
| M2 | Rechteck- und Halbebenenargument nachvollzogen (Abschnitte von $a$ umschliessen Rechtecke in Re ≤ −1; $b_k$ bei $m=0$ geradlinig homotop in der oberen Halbebene); Injektivität der Kodierung ergänzt ($u\mapsto(m_1,u')$ ist bijektiv, $0^\infty$ nur für $u=0$); `m1_monodromie.py A` neu gerechnet | Ausgabe identisch mit `m1_ausgabe_A.txt` (48/48, $s_1$, $J^{(1)}$). Bestanden. |
| M3 | Klassifikation über Bahnen von $b_k$, Liouville-Schritt geprüft | richtig. |
| M4, M4a, M4b | treue Kranzprodukt-Darstellung, $b_{k+1}\notin G_k$, Relationen über getrennte Träger, IMG-Argument geprüft; `m2_relationen.py` neu gerechnet | Ausgabe identisch mit `m2_ausgabe.txt`. Bestanden. |
| M5 | $\Psi'(L^is_j)=\Psi'(s_j)$ per Induktion nachgerechnet; $|L||s_1|=3.348$ nachgerechnet; $J^{(1)}$ gegen it.satz.regulaerer_erzeuger_fortsetzung | gleich auf alle 20 Stellen des Bestands. Bestanden. |
| M6 | Umformung geprüft | richtig; benutzt nur die unbedingt bewiesene Produktformel von it.satz.kneser_erzeuger (dessen Bedingung betrifft nur den Glattheitsgrad) — Kante deshalb `stuetzt_sich_auf` mit Notiz statt `folgt_unmittelbar_aus`. |
| M7 | Beweis geprüft (F_L lässt 0 aus ⇒ unbeschränkt; $\mathrm{sexp}_K$ biholomorph um −1; Konjugate von $a$ haben nur unendliche Bahnen) | richtig. Der Teil „unbeschränkt“ steht schon in kn.satz.theta_singularitaet (ii); neu ist die Verzweigungsordnung. „Logarithmisch“ ist topologisch gemeint (Blätterzahl) — im Eintrag klargestellt, damit kein Widerspruch zu (iii) „weder Potenz noch Logarithmus“ (Wachstum) entsteht. M7a: Re(Lθ) aus `m3_ausgabe.txt` nachgerechnet (1.0979 bei y = 0.3, 1.4069 bei y = 0.004); θ(0.01i) = 26.2051+5.2268i stimmt mit kn.satz.theta_singularitaet (26.205+5.227i). Skript nicht neu gerechnet (liest Daten aus anderen Ordnern). |
| M8 | Beweis gelesen | **herabgestuft** von „bewiesen (Überlagerung, Untergruppe, Quotient)“ auf **Skizze**: $H/N\cong G$ braucht, dass der Kern von $H\to G$ genau $N=\bigoplus\mathbb Z$ ist; für Elemente aus dem Kern von $\pi_1(\mathbb C\setminus P)\to G$ (z. B. $[b_1,ab_1a^{-1}]$) ist nicht gezeigt, dass ihre Verschiebungen von $k$ endlich getragen sind. Überlagerung und $N\cong\bigoplus\mathbb Z$ sind in Ordnung. |
| M9 | — | in it.satz.erzeuger_kranzrekursion eingearbeitet. |
| M10 | — | Skizze wie im Bericht (zwei Lücken benannt). |

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **it.satz.regulaerer_erzeuger_fortsetzung**, Feld `pruefung.f5_neue_frage`: „Welche Darstellung von π₁(ℂ∖P) ist die Monodromiegruppe? Ist jeder Punkt exp^k(0) ein echter Verzweigungspunkt?“ → „Beantwortet: Kranzrekursion, Monodromiegruppe IMG(exp) = ℤ≀ℤ≀… (it.satz.erzeuger_kranzrekursion, it.satz.erzeuger_monodromiegruppe); jedes exp^k(0) ist auf dem Hauptblatt logarithmisch (it.kor.erzeuger_verzweigung_logarithmisch). Offen: Monodromie der regulären slog (it.satz.regulaere_slog_riemannflaeche).“ Kante `motiviert` → it.satz.erzeuger_kranzrekursion.
2. **kn.satz.theta_singularitaet**, Feld `aussage` (iii): „Die Singularität ist weder eine Potenz noch ein Logarithmus.“ → „Die Singularität wächst weder wie eine Potenz noch wie ein Logarithmus; topologisch ist sie ein Verzweigungspunkt unendlicher Ordnung (kn.kor.theta_verzweigung_unendlich).“ Kante `verwandt` → kn.kor.theta_verzweigung_unendlich.
3. **it.satz.kneser_erzeuger**, Feld `pruefung.f5_neue_frage` ergänzen: „Knesers Erzeuger als Funktion von σ: it.kor.kneser_erzeuger_sigma.“
4. Nach dem Eintragen von Bericht 07 (anderer Bearbeiter, `r25_07_kn.yaml`): Kante zwischen it.satz.regulaere_slog_riemannflaeche und den dortigen slog-Monodromie-Aussagen (z. B. kn.kor.slog_monodromie_L) prüfen; deren Fragen nach der Monodromiegruppe berühren M4/M8.

## Neue Quellen (Vorschlag, nicht angelegt)

- O. Forster, *Riemannsche Flächen*, Springer 1977, Satz 4.19 (lokaler Homöomorphismus mit Wegehebung ist Überlagerung) — benutzt in it.satz.erzeuger_riemannflaeche.
- V. Nekrashevych, *Self-similar groups*, AMS 2005 (iterierte Monodromiegruppen, Kranzrekursion) — Sprache von it.satz.erzeuger_kranzrekursion / it.satz.erzeuger_monodromiegruppe; Literaturabgleich für IMG transzendenter Abbildungen steht aus.

## Stränge und Ideenkarte

- Bestehender Strang **it.erzeuger**: die acht it-ids aufnehmen. `verlauf` ergänzen: „Der reguläre Erzeuger lebt auf der u-Ebene; seine Monodromie ist ein Zählwerk (IMG(exp) = ℤ≀ℤ≀…), jede Verzweigung ist logarithmisch.“
- kn.kor.theta_verzweigung_unendlich in den bestehenden Strang **kn.theta**; `verlauf` ergänzen: „… an den ganzen Zahlen singulär, mit unendlich vielen Blättern …“.
- Ideenkarte: nichts Neues; über die Stränge automatisch in **dynamik.erzeuger** bzw. der Idee von kn.theta.
