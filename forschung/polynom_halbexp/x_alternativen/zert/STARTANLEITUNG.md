# Zweig-Zertifikat 193 → 400 — Startanleitung und Wiederaufnahme (Strang X, 11.10.2026 04:50)

**Stand:** vorbereitet und getestet, der volle Lauf ist NICHT gestartet (Team angehalten). Auf aicgolling liegen `~/tetration-zert`
(Skripte + daten/) und das eigene venv `~/tetration-venv` (python-flint 0.9.0, mpmath 1.3.0); es laufen keine Prozesse von mir.

**Getestet (Server, Ergebnisse in `test_server/`):** Übergänge 193 → 194 und 194 → 195 ok (NC = 7, Reserve 0.4999998, t3 = 2e-7, je ≈ 27 s),
Kette 193 → 195 ok (Ps Box von p_193 ⊂ B_193, K_194 ⊂ B_194). Ausserdem 300 → 301 ok (NC = 10, t1 = 2.8e-72, t3 = 9.6e-7, 309 s).

**Was die Skripte tun:** `zert_kern.py` = Kopie von Ps p24 (p_pruefung/, nur gelesen) mit zwei Änderungen: NC als Parameter und volle
Newton-Schritte an den Tschebyschow-Knoten (Ps Sehnen-Newton bleibt bei ~1e-51 stehen; für M > 220 muss der Röhrenradius r wie 1.88^{−M}
fallen). `zert_job.py N`: Newton-Mitte von p_N (Start: Ps Box bei 193, Zweigdaten bis 300, sonst T_N f_K aus den zertifizierten Kneser-
Koeffizienten), Startbox B_N = Mitte ± 2^{−(N+150)}ρ^{−k} (1000-mal grösser als die Ausgabebox des Vorgängers), Krawczyk-Übergang
mit NC = 7 + ⌈(M−205)/30⌉, bei Misserfolg NC + 3, + 6. `zert_kette.py N1 N2`: prüft P-Box ⊂ B_193 und K_{N+1} ⊂ B_{N+1} für alle N.
**Lehren aus zwei abgebrochenen Läufen:** (1) ohne volle Newton-Schritte scheitern alle M ≳ 230 (t3 bis 1e43); (2) eine feste Startbox 1e-150
ist für M < 338 KLEINER als die Ausgabebox des Vorgängers (2^{−(M+160)}) — die Kette würde reissen; daher jetzt 2^{−(N+150)}.

**Start (zwei getrennte Aufrufe, Git Bash, `export PATH="/c/Users/janis/bin:$PATH"`):**
1. vorher `aic 'uptime; nproc'` — höchstens so viele Arbeiter, dass ≥ 16 Kerne frei bleiben (lauf.sh: 12 Arbeiter bei 28 Kernen);
   alte Ergebnisse wegschieben: `aic 'cd ~/tetration-zert && mkdir -p alt3 && mv log boxen erg lauf.log alt3/ 2>/dev/null'`
2. `aic 'cd ~/tetration-zert && setsid nohup timeout 43200 bash ~/tetration-zert/lauf.sh > ~/tetration-zert/lauf.log 2>&1 < /dev/null &'`
   (der aic-Aufruf kehrt nicht zurück → mit `timeout 30` aufrufen; der Lauf läuft weiter), danach
   `aic 'ps -u janis -o comm=,args= | awk "\$1==\"python\" && /zert_job/" | wc -l'` → 12.
3. Fortschritt: `aic 'cd ~/tetration-zert; grep -h ENDE log/job_*.log | wc -l; grep -l "ok = False" log/job_*.log'`.
   Erwartet: ≈ 25–30 CPU-Stunden, mit 12 Arbeitern ≈ 2.5 h; die längsten Übergänge (M ≈ 400, NC = 14) ≈ 15–20 min.
4. Ende: lauf.sh ruft zert_kette.py 193 400 auf → `log/kette.log`, `kette_193_400.json`. Ergebnisse zurückholen:
   `aic 'cd ~/tetration-zert && tar -cf - log boxen erg kette_193_400.json' | tar -C <x_alternativen/zert/lauf400> -xf -`.
5. Abbrechen nur über die Sitzung: SID des `timeout … lauf.sh`-Prozesses suchen und alle Prozesse dieser SID beenden
   (`kill -TERM $(ps -u janis -o pid=,sid= | awk -v s=SID "\$2==s {print \$1}")`), nie `pkill -f` mit einem Muster aus dem eigenen Befehl.
6. Danach: BRETT-Eintrag, P prüft die Kette (zert_kette.py ist kurz; die Strenge liegt in Ps uebergang()).
