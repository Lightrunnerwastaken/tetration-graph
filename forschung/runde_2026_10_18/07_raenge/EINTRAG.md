# Eintrag in den Graphen — Bericht 07 „Ränge“

Dateien: `daten/aussagen/r18_07_kn.yaml`, `daten/theorie/r18_07_kn.yaml` (8 Aussagen, je ein Theorie-Block).

## Angelegte ids
- kn.satz.koenigs_leiter_monoton | Koenigs-Stufen steigen gegen p_r, b < p_{r+1} < p_r (B1, bewiesen)
- kn.verm.grenzprofil_rampe | Grenzprofil H_r(b,x) → min(x+1, b) (B2, Vermutung)
- kn.satz.plateau_modell | (δ, λ)-Modellabbildung der Leiter (B3, Skizze mit Lücke O(Lδ))
- kn.verm.leiter_rate | r·μ_r → k = (2−b)/(b−1), δ_r ~ (b−1)k^k/(r^k ln r) (B4)
- kn.beob.rangverschiebung | θ-Störung = feste Rangverschiebung c(ε), Querabweichungen überexponentiell klein (B5+B6, numerisch)
- kn.verm.rangkoordinate | Abel-Funktional Φ(Tf) = Φ(f)+1; halber Rang = eine Gleichung (B7)
- kn.satz.konvex_eta_grenze | konvexe Ketten ⇒ η_r → 2, Schranke b + b(b−1)^{r−2} = 2 (B8, bedingt)
- kn.gescheitert.konvexitaetsweg | kritische F_N ab N = 5 nicht konvex (S-Form); dazu (4−e_r)/(2−η_r) ohne Grenzkonstante 3 (B9+B10)

Nicht eingetragen: B11 (Leiter auf negativen Höhen zu langsam, rein technisch).

## Vorschlag Stränge
- Bestehender Strang `kn.halbrang`: kn.beob.rangverschiebung, kn.verm.rangkoordinate.
  verlauf: „18.10.: Im beschränkten Bereich wird die θ-Freiheit hoch in der Leiter zu einer einzigen Rangverschiebung; ein halber Rang ist dort eine Gleichung für ein Abel-Funktional Φ.“
- Bestehender Strang `kn.hoehere_stufen`: kn.satz.koenigs_leiter_monoton, kn.verm.grenzprofil_rampe, kn.satz.plateau_modell, kn.verm.leiter_rate.
  verlauf: „18.10.: Beschränkte Leitern bis Rang 16 — Fixpunkte fallen beweisbar, das Profil geht vermutlich gegen die Rampe min(x+1, b), mit Rate k/r.“
- Bestehender Strang `kn.etas`: kn.satz.konvex_eta_grenze, kn.gescheitert.konvexitaetsweg.
  verlauf: „18.10.: Konvexität würde η_r → 2 beweisen, gilt aber an den kritischen Stufen ab N = 5 nicht.“

## Vorschlag Ideenkarte
- `stufen.halbrang` (Gebrochene Ränge): kn.verm.rangkoordinate als Kernbefund.
- `stufen.etas`: kn.satz.konvex_eta_grenze, kn.gescheitert.konvexitaetsweg.
- Optional neue Unteridee unter `stufen`: `stufen.grenzprofil`, Name „Die Leiter im Grenzfall“, kurz: „Für hohe Ränge wird jede beschränkte Rechenart zur Rampe min(x+1, b).“
