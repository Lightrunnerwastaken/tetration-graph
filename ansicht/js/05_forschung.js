// ===================================================================== Forschung
// Brett statt Graph: Lücken, Wurzeln, Experimente, Korrekturen, nächste Fragen.
function zeichneForschung() {
  const f = $("#forschung");
  const eintrag = (id, text, meta) => `<button class="eintrag" data-id="${id}"><div class="eintrag-titel">${punkt(K.get(id)._status_eff)} ${mathe(K.get(id).titel)}</div>${text ? `<div class="eintrag-text">${mathe(text)}</div>` : ""}${meta ? `<div class="eintrag-meta">${meta}</div>` : ""}</button>`;
  const spalte = (titel, erkl, inhalt) => `<section class="spalte"><h4>${titel}</h4><p class="klein">${erkl}</p>${inhalt || '<p class="klein">Leer.</p>'}</section>`;
  const geb = (k) => `<span class="punkt" style="background:${GEB.get(k.gebiet).farbe};width:7px;height:7px;border-radius:50%;display:inline-block"></span>${esc(GEB.get(k.gebiet).name)}`;

  const neu = G.knoten.filter((k) => k._neu_pruefen?.length);
  const korr = G.knoten.filter((k) => ["korrektur", "gegenbeispiel", "gescheitert"].includes(k.art) || ["widerlegt", "ueberholt"].includes(k.status));
  const wurzeln = wurzelZaehlung();
  const maxW = Math.max(1, ...wurzeln.map((w) => w[1]));
  const luecken = G.knoten.filter((k) => k.beweis?.luecken?.length);
  const verm = G.knoten.filter((k) => k.art === "vermutung");
  const fragen = G.knoten.filter((k) => k.naechste_frage);
  const versionen = G.knoten.filter((k) => k.versionen?.length);
  const ungeprueft = {};
  G.knoten.forEach((k) => (k.quellen || []).forEach((q) => {
    const ref = typeof q === "string" ? q : q.ref;
    const gp = (typeof q === "object" && q.geprueft) || Q[ref]?.geprueft;
    if (gp === "nein") (ungeprueft[ref] = ungeprueft[ref] || []).push(k.id);
  }));

  f.innerHTML = [
    spalte("Neu prüfen", "Aussagen, deren Beweiskette eine widerlegte oder überholte Aussage enthält.",
      neu.map((k) => eintrag(k.id, "", k._neu_pruefen.map((p) => p.map((x) => esc(K.get(x).titel)).join(" → ")).join("<br>"))).join("")),
    spalte("Woran die bedingten Sätze hängen", "Unbewiesene Enden der Beweisketten, nach Zahl der abhängigen Sätze.",
      wurzeln.map(([w, c]) => `<button class="eintrag" data-id="${w}"><div class="eintrag-titel">${punkt(K.get(w)._status_eff)} ${mathe(K.get(w).titel)}</div><div class="eintrag-meta">${c} bedingte Sätze</div><div class="wurzel-balken" style="width:${(100 * c) / maxW}%"></div></button>`).join("")),
    spalte("Beweislücken", "Benannte offene Schritte in Skizzen und bedingten Beweisen.",
      luecken.map((k) => eintrag(k.id, k.beweis.luecken.map((l) => "• " + l).join("\n"), geb(k))).join("")),
    spalte("Experimente", "Jedes Experiment prüft sich selbst; „veraltet“ heisst: Skript geändert, neu rechnen.",
      G.experimente.map((e) => {
        const pr = e.ergebnis?.pruefungen || [];
        return `<button class="eintrag" data-exp="${e.id}"><div class="exp-mini"><div class="bild">${e.svg || ""}</div><div><div class="eintrag-titel">${esc(e.titel)}</div><div class="eintrag-meta">${pr.filter((p) => p.ok).length}/${pr.length} Prüfungen · ${esc(e.zustand)}</div></div></div></button>`;
      }).join("")),
    spalte("Vermutungen und Gegenbeispielsuche", "Was gesucht wurde und mit welchem Ergebnis.",
      verm.map((k) => eintrag(k.id, kurz(reinText(k.gegenbeispielsuche || "keine Suche dokumentiert"), 220), `${esc(S.status[k._status_eff].label)} · ${geb(k)}`)).join("")),
    spalte("Korrekturen, Gegenbeispiele, gescheiterte Ansätze", "Bleiben erhalten und sind mit den betroffenen Aussagen verknüpft.",
      korr.map((k) => eintrag(k.id, "", `${esc(S.arten[k.art])} · ${esc(S.status[k._status_eff].label)}`)).join("")),
    spalte("Fassungen und Versionskonflikte", "Aussagen mit dokumentierten Änderungen.",
      versionen.map((k) => { const v = k.versionen[k.versionen.length - 1]; return eintrag(k.id, v.aenderung, `${esc(v.datum)}${v.grund ? " · " + esc(reinText(v.grund)) : ""}`); }).join("")),
    spalte("Nächste prüfbare Fragen", "Konkrete nächste Schritte, je eine Zeile.",
      fragen.map((k) => eintrag(k.id, k.naechste_frage, geb(k))).join("")),
    spalte("Ungeprüfte Quellen", "Aus dem Gedächtnis zitiert; vor Weiterverwendung nachschlagen.",
      Object.entries(ungeprueft).sort((a, b) => b[1].length - a[1].length).map(([ref, ids]) => `<div class="eintrag"><div class="eintrag-titel">${esc(Q[ref].autor)}: ${esc(Q[ref].titel)}</div><div class="eintrag-meta">${ids.length}× zitiert · ${ids.map((i) => `<button class="verweis" data-id="${i}">${esc(i)}</button>`).join(", ")}</div></div>`).join("")),
  ].join("");
}
