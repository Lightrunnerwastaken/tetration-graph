# Eintrag Bericht 13 „Rangkoordinate Φ“ (25.10.2026)

Dateien: `daten/aussagen/r25_13_kn.yaml`, `daten/theorie/r25_13_kn.yaml`. `python baue.py --nur-pruefen`: 0 Fehler, keine Warnung aus diesen Dateien.

## Angelegte ids

| id | Titel kurz |
|---|---|
| kn.konstr.rangkoordinate | Rangkoordinate Φ aus ln(p_r − b); Φ(Tg) = Φ(g)+1, Φ(H_r) = r, falls der Grenzwert existiert (bedingt) |
| kn.beob.rangkoordinate_gradient | Gradient nach θ-Moden: −6.125, 8.060, −2.84 / −8.042, −11.258, 9.91; kein Frequenzabfall (numerisch) |
| kn.beob.rangkoordinate_nichtlinear | Φ differenzierbar, nicht linear; Q_sc ≈ 17, Q_cc ≈ 29 (numerisch) |
| kn.verm.rangkoordinate_rangwachstum | dieselbe Welle wirkt höher stärker, ≈ 1.6 ln r (Vermutung) |
| kn.beob.halbrang_niveauflaeche | gleiche Φ ⇒ gleiche Leiter darüber, aber verschiedenes H_3.5 (numerisch) |
| kn.verm.halbrangleiter_oben_kanonisch | Halbrang-Leiter oben kanonisch, unten nur bis auf das Blatt Φ = const (Vermutung) |
| kn.gescheitert.rangkoordinate_basis2 | bei b = 2 keine Fixpunkt-Messgrösse (gescheitert) |

## Gegenprobe

| Befund | geprüft wie | Ergebnis |
|---|---|---|
| B1 Abel-Gleichung | Beweis gelesen (Indexverschiebung $s_n(Tg)=s_{n+1}(g)$, Stützstellen geben $\Phi(H_r)=r$) | richtig; Existenz nur numerisch → Status **bedingt** (Bericht: „bewiesen“ für die Abel-Gleichung). Beispielwerte 3.969339 … aus `cache/r4_s1_0.005…` nachgelesen ($n=5\dots8$). |
| B2 Gradient | **frisch gerechnet** mit `gegenprobe.py` (benutzt `phi.py` unverändert, eigener Speicher `gegenprobe_cache/`, 28 s): Basisleiter, dann $H_4\circ(x\pm0.005\,\eta)$ für sin1 und (1−cos)1 | ∂Φ/∂sin1 = −6.125149, ∂Φ/∂(1−cos)1 = −8.0422898 ($n=8$) — identisch mit `daten/exp_…_grad.json`. Bestanden. |
| B3 Punktmodell | `punktmodell.py` neu ausgeführt | Reste bis 7.75 (Modell mit η, η′, η″), bis 11.4 mit weniger Ableitungen. Bestanden. |
| B4 Nichtlinearität | Zahlen in `daten/exp_…_add.json` nachgelesen, Q-Werte nachgerechnet | Q_sc = 17.25, Q_s1s2 = −1.57, Q_ss = −1.02 (symmetrische Konvention); zusätzlich Q_cc ≈ 29 aus c(1−cos, 0.02) abgeleitet, passt zu Runde 18. Bestanden. |
| B5 Rangwachstum | Folgen aus `cache/r4…r8_s1_±0.005` neu ausgewertet (eigenes Kurzskript, nur lesend) | −6.125149, −6.525833, H6 −6.821844, H7 −7.05593, H8 −7.236063, alle Folgen monoton fallend: die Schranken stimmen. Bestanden. |
| B6 Basis 1.5 | nicht neu gerechnet | aus dem Bericht übernommen (2–3 Stellen). |
| B7/B8 halber Rang | `half2.py` neu ausgeführt (`gegenprobe_half2_ausgabe.txt`) und beide Leitern $T^ng^*$ mit $t^*=0.62762$/$0.62300$ frisch gerechnet | Leiterwerte und Unterschiede reproduziert (2.82e−3, 9.26e−5, 9.26e−6, 1.68e−6, 5.45e−7, …, 5.97e−8 bei n = 10), $g^*(0.5)$, $g^*(1.5)$ gleich. **Neu bemerkt:** beide $g^*$ haben $\Phi\approx3.4991$ (6 Punkte; 4 Punkte 3.4989), nicht 3.5 — die quadratische Interpolation in $t$ durch drei Stützstellen trifft die Nullstelle nur auf ~$10^{-3}$. Der „Rest-Unterschied 5·10⁻⁵“ im Bericht ist $\Phi_{\rm lin}-\Phi_{\rm geo}$. Die Aussage (gleiche Φ ⇒ gleiche Zukunft, verschiedene Gegenwart) hält; im Eintrag präzisiert. |
| B9 Struktur | Herleitung gelesen | als **Vermutung** eingetragen (Existenz der langsamen Kurve/Blätterung offen). |
| B10 | Werte 3.49888/3.49913 in der half2-Ausgabe bestätigt (3.4988799, 3.499131) | in kn.konstr.rangkoordinate eingearbeitet. |
| B11 | — | als gescheitert eingetragen. |

## Änderungen am Bestand (Vorschläge, nicht ausgeführt)

1. **kn.beob.rangverschiebung**, Feld `aussage`: „Gradient $(dc/d\varepsilon_{\sin},dc/d\varepsilon_{\cos})\approx(-6.13,-7.74)$“ → „Gradient $\approx(-6.125,-8.04)$ (symmetrische Differenzen, kn.beob.rangkoordinate_gradient); $-7.74$ ist die Sekante bei $\varepsilon=0.01$, der quadratische Term ist in der Kosinusrichtung gross ($Q_{cc}\approx29$)“. Begründung: B2, beim Eintragen frisch nachgerechnet.
2. **kn.verm.rangkoordinate**, Feld `aussage`: „(bei $\sqrt2$ Gewichte $\approx(-6.13,-7.74)$ für $(\sin,1-\cos)$)“ → „(bei $\sqrt2$ Gewichte $(-6.125,-8.042)$, kn.beob.rangkoordinate_gradient)“; ergänzen: „Rechenbare Fassung: kn.konstr.rangkoordinate (Existenz numerisch auf 5–6 Stellen); die Zählung ist bei √2 bestätigt (kn.beob.halbrang_niveauflaeche).“ Status bleibt `vermutung`.
3. **kn.frage.halber_rang**, Feld `aussage`, ergänzen: „Bei $b=\sqrt2$ legt die Rangkoordinate nur ein Blatt fest: Kandidaten mit gleichem Φ haben $H_{3.5}(\sqrt2,0.5)\approx1.2228$–$1.2234$ je nach Familie (kn.beob.halbrang_niveauflaeche).“

## Stränge und Ideenkarte

- Bestehender Strang **kn.halbrang**: alle sieben ids als Mitglieder aufnehmen. `verlauf` ergänzen: „Die Rangkoordinate Φ ist rechenbar (Abel-Gleichung exakt, Existenz numerisch); Φ = 3.5 legt die Leiter darüber fest, nicht den halben Rang selbst.“
- Ideenkarte: nichts Neues; über den Strang automatisch in **stufen.halbrang**.
