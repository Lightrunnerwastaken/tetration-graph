// ===================================================================== Tafel
// Rechts einklappbar: erscheint erst, wenn etwas gewählt ist (oder auf Wunsch).
const tafel = $("#tafel");
const punkt = (st) => `<span class="punkt" style="background:${farbe(st)}"></span>`;
const knotenLink = (id) => {
  const k = K.get(id);
  if (!k) return esc(id);
  return `${punkt(k._status_eff)}<button class="verweis" data-id="${id}">${mathe(k.titel)}</button>`;
};
function tafelZeigen(html) {
  tafel.innerHTML = `<div class="tafel-leiste"><button class="tafel-zu" data-tafel="zu" aria-label="Tafel schliessen">×</button></div><div class="tafel-innen">${html}</div>`;
  tafel.scrollTop = 0;
  document.body.classList.add("tafel-offen");
  $("#tafel-knopf").setAttribute("aria-expanded", "true");
}
function tafelSchliessen() {
  document.body.classList.remove("tafel-offen");
  $("#tafel-knopf").setAttribute("aria-expanded", "false");
}
function statusband(liste) {
  const z = {};
  liste.forEach((k) => (z[k._status_eff] = (z[k._status_eff] || 0) + 1));
  return `<div class="statusband">${RINGFOLGE.filter((s) => z[s]).map((s) => `<span title="${esc(S.status[s].label)}: ${z[s]}" style="flex:${z[s]};background:${farbe(s)}"></span>`).join("")}</div>
    <div class="legende-reihe">${RINGFOLGE.filter((s) => z[s]).map((s) => `<span>${punkt(s)}${esc(S.status[s].label)} ${z[s]}</span>`).join("")}</div>`;
}
function wurzelZaehlung() {
  const z = new Map();
  for (const k of G.knoten) for (const w of k._wurzelbedingungen || []) if (k._status_eff === "bedingt") z.set(w, (z.get(w) || 0) + 1);
  return [...z.entries()].sort((a, b) => b[1] - a[1]);
}

function zeigeUebersicht() {
  const st = G.statistik;
  const wurzeln = wurzelZaehlung();
  const neu = G.knoten.filter((k) => k._neu_pruefen?.length);
  tafelZeigen(`
    <h2>Tetration und ihre Nachbarn</h2>
    <p class="klein">${G.knoten.length} Aussagen in ${G.gebiete.length} Gebieten und ${(G.straenge || []).length} Forschungssträngen, ${G.kanten.length} Beziehungen, ${G.experimente.length} Experimente. Gebaut ${esc(G.gebaut)}.</p>
    ${statusband(G.knoten)}
    <p class="klein" style="margin-top:10px">Mausrad oder die Stufenleiste links: Gebiete → Verläufe → Aussagen → Inhalt. Ein Klick auf ein Gebiet oder einen Strang zoomt hinein, ein Klick auf eine Aussage stellt sie mit ihren Nachbarn in den Mittelpunkt.</p>
    <h3>Hier anfangen: Grundlagen</h3>
    <p class="klein">Elementar bewiesen und an der Kneser-Tetration nachgerechnet; jede Aussage hat eine Erklärung in einfachen Worten. Mit den Knöpfen „Einstieg / Vertiefung / Forschung“ unten links lässt sich die Karte auf ein Niveau beschränken.</p>
    ${(G.straenge || []).filter((s) => s.gebiet === "gr").map((s) => `<button class="eintrag" data-strang="${s.id}"><div class="eintrag-titel">${esc(s.name)}</div><div class="eintrag-text">${mathe(s.verlauf)}</div></button>`).join("")}
    <h3>Woran die bedingten Sätze hängen</h3>
    <ul class="liste">${wurzeln.map(([w, c]) => `<li>${knotenLink(w)}<span class="klein" style="margin-left:auto">${c} Sätze</span></li>`).join("")}</ul>
    <h3>Neu prüfen</h3>
    ${neu.length ? `<ul class="liste">${neu.map((k) => `<li>${knotenLink(k.id)}</li>`).join("")}</ul>` : `<p class="klein">Nichts markiert.</p>`}
    <h3>Prüfstand</h3>
    <p class="klein">${st["befund:FEHLER"] || 0} Fehler, ${st["befund:WARNUNG"] || 0} Warnungen, ${st["befund:HINWEIS"] || 0} Hinweise · ${st.zitate_ungeprueft} Zitate noch nicht am Original geprüft.</p>
    <h3>Gebiete</h3>
    <ul class="liste">${G.gebiete.map((g) => `<li><span class="punkt" style="background:${g.farbe}"></span><button class="verweis" data-gebiet="${g.id}">${esc(g.name)}</button><span class="klein" style="margin-left:auto">${G.knoten.filter((k) => k.gebiet === g.id).length}</span></li>`).join("")}</ul>`);
}

function zeigeGebiet(gid) {
  const g = GEB.get(gid);
  const liste = G.knoten.filter((k) => k.gebiet === gid);
  const st = (G.straenge || []).filter((s) => s.gebiet === gid);
  tafelZeigen(`
    <div class="kopfzeile"><span class="marker"><span class="punkt" style="background:${g.farbe}"></span>Gebiet</span></div>
    <h2>${esc(g.name)}</h2>
    <p class="klein">${mathe(g.kurz)}</p>
    ${statusband(liste)}
    <h3>Forschungsstränge</h3>
    ${st.map((s) => `<button class="eintrag" data-strang="${s.id}"><div class="eintrag-titel">${esc(s.name)}</div><div class="eintrag-text">${mathe(s.verlauf)}</div><div class="eintrag-meta">${s.mitglieder.length} Aussagen</div></button>`).join("")}`);
}

function zeigeStrang(sid) {
  const s = STR.get(sid);
  if (!s) return;
  const g = GEB.get(s.gebiet);
  const glieder = s.mitglieder.map((i) => K.get(i)).filter(Boolean);
  // Reihenfolge entlang der Beweiskette: Grundlagen zuerst
  const inS = new Set(s.mitglieder);
  const tiefe = new Map();
  const t = (i, weg = new Set()) => {
    if (tiefe.has(i)) return tiefe.get(i);
    if (weg.has(i)) return 0;
    weg.add(i);
    const d = (K.get(i)._abhaengig || []).filter((x) => inS.has(x));
    const v = d.length ? 1 + Math.max(...d.map((x) => t(x, weg))) : 0;
    tiefe.set(i, v);
    return v;
  };
  glieder.forEach((k) => t(k.id));
  glieder.sort((a, b) => tiefe.get(a.id) - tiefe.get(b.id) || (b._tragweite || 0) - (a._tragweite || 0));
  tafelZeigen(`
    <div class="kopfzeile"><span class="marker"><span class="punkt" style="background:${g.farbe}"></span><button class="verweis" data-gebiet="${g.id}">${esc(g.name)}</button></span><span class="marker">Forschungsstrang</span></div>
    <h2>${esc(s.name)}</h2>
    <div class="aussage mathe" style="--status:${g.farbe}">${mathe(s.verlauf)}</div>
    ${statusband(glieder)}
    <h3>Verlauf, von den Grundlagen zu den Folgerungen</h3>
    <ol class="verlauf">${glieder.map((k) => `<li style="--tiefe:${tiefe.get(k.id)}">${knotenLink(k.id)}<span class="klein">${esc(S.arten[k.art])}${k._neu_pruefen?.length ? " · neu prüfen" : ""}</span></li>`).join("")}</ol>`);
}

function experimentKarte(eid, mitBild = true) {
  const e = EXP.get(eid);
  if (!e) return "";
  const erg = e.ergebnis || {};
  const pr = erg.pruefungen || [];
  const ok = pr.filter((p) => p.ok).length;
  return `<div class="exp"><div class="exp-kopf"><b>${esc(e.titel)}</b><span class="klein" style="margin-left:auto">${ok}/${pr.length} Prüfungen · ${esc(e.zustand || "fehlt")}</span></div>
    ${mitBild && e.svg ? `<div class="exp-bild">${e.svg}</div>` : ""}
    <div class="exp-text"><div class="klein">Frage: ${mathe(e.frage)}</div><p style="margin-top:6px">${mathe(erg.zusammenfassung || "noch nicht gelaufen")}</p>
    <details><summary class="klein">Einzelprüfungen, Skript</summary><ul class="exp-pruef">${pr.map((p) => `<li class="${p.ok ? "" : "fehl"}">${esc(p.name)}${p.wert != null ? ` <span class="klein">[${esc(typeof p.wert === "object" ? JSON.stringify(p.wert) : p.wert)}]</span>` : ""}</li>`).join("")}</ul>
    <p class="klein">${esc(e.skript)} · gelaufen ${esc(erg.gelaufen || "—")} (${esc(erg.laufzeit_s ?? "?")} s)</p></details></div></div>`;
}
function zeigeExperiment(eid) {
  const e = EXP.get(eid);
  tafelZeigen(`<div class="kopfzeile"><span class="marker">Experiment</span><span class="knoten-id">${esc(eid)}</span></div>
    <h2>${esc(e.titel)}</h2>${experimentKarte(eid)}
    <h3>Stützt</h3><ul class="liste">${(e.stuetzt || []).map((s) => `<li>${knotenLink(s)}</li>`).join("")}</ul>`);
}

function zeigeKnoten(id) {
  const k = K.get(id);
  const geb = GEB.get(k.gebiet);
  const st = k._status_eff;
  const strang = STR.get(k._strang);
  const T = [];
  T.push(`<div class="kopfzeile">
    <span class="marker"><span class="punkt" style="background:${geb.farbe}"></span><button class="verweis" data-gebiet="${geb.id}">${esc(geb.name)}</button></span>
    ${strang ? `<span class="marker"><button class="verweis" data-strang="${strang.id}">${esc(strang.name)}</button></span>` : ""}
    <span class="marker">${esc(S.arten[k.art])}</span>
    <span class="marker ebene" title="${esc(S.ebenen[k.ebene])}">Ebene ${k.ebene}</span>
    <span class="marker niveau niveau-${k._niveau}" title="${esc(S.niveau[k._niveau] || "")}">${NIVEAU_NAME[k._niveau] || ""}</span>
    <span class="knoten-id">${esc(k.id)}</span></div>`);
  T.push(`<h2>${titelHtml(k.titel)}</h2>`);
  T.push(`<div class="kopfzeile"><span class="marker">${punkt(st)}${esc(S.status[st]?.label)}</span>
    ${st !== k.status ? `<span class="klein">eingetragen: ${esc(S.status[k.status]?.label)} → abgeleitet: ${esc(S.status[st]?.label)}</span>` : ""}
    <span class="klein">${esc(S.ebenen[k.ebene])}</span></div>`);
  if (k._neu_pruefen?.length) T.push(`<div class="warnbox"><b>Neu prüfen.</b> Hängt an widerlegten oder überholten Aussagen:<br>${k._neu_pruefen.map((p) => p.map(knotenLink).join(" → ")).join("<br>")}</div>`);
  if (k._teilweise_entlastet_durch) T.push(`<div class="hinweisbox">Teilweise entlastet durch ${k._teilweise_entlastet_durch.map(knotenLink).join(", ")}.</div>`);
  if (k.einfach) T.push(`<div class="einfach"><div class="einfach-kopf">In einfachen Worten</div>${mathe(k.einfach)}</div>`);
  T.push(`<div class="aussage mathe" style="--status:${farbe(st)}">${mathe(k.aussage)}</div>`);
  if (k.theorie) T.push(`<div class="tafel-lehrbuch"><div class="tafel-lehrbuch-kopf"><span>Ausführlich erklärt</span>
    <button class="verweis" data-lehrbuch="${id}">im Lehrbuch lesen</button></div>${theorieHtml(k, true)}</div>`);
  else if (buch.kapitelVon(id)) T.push(`<p class="klein"><button class="verweis" data-lehrbuch="${id}">Im Lehrbuch an dieser Stelle lesen</button></p>`);
  if (k.formen?.length) {
    T.push(`<h3>Gleichwertige Formen</h3><ul class="formen">${k.formen.map((f) => {
      const tx = typeof f === "string" ? f : f.tex, no = typeof f === "string" ? "" : f.notiz;
      return `<li><div class="form-mathe">${tex(tx, true)}</div>${no ? `<div class="klein">${mathe(no)}</div>` : ""}</li>`;
    }).join("")}</ul>`);
  }
  if (k.bereich) T.push(`<p class="klein" style="margin-top:8px"><b>Bereich:</b> ${mathe(k.bereich)}</p>`);
  if (k.objekte?.length) T.push(`<p class="klein"><b>Objekte:</b> ${k.objekte.map((o) => G.objekte[o] ? tex(G.objekte[o].tex, false) + " " + esc(G.objekte[o].name) : esc(o)).join(" · ")}</p>`);
  if (k.voraussetzungen?.length) {
    T.push(`<h3>Voraussetzungen — im Beweis benutzt?</h3><ul class="voraussetzungen">${k.voraussetzungen.map((v) => `
      <li><span class="benutzt ${v.benutzt}" title="benutzt: ${v.benutzt}">${{ ja: "✓", nein: "✗", unklar: "?" }[v.benutzt] || "?"}</span>
      <span>${mathe(v.text)}${v.ref ? ` · ${knotenLink(v.ref)}` : ""}${v.notiz ? `<div class="klein">${mathe(v.notiz)}</div>` : ""}${v.benutzt === "nein" ? `<div class="klein" style="color:var(--s-annahme)">unbenutzt — kann gestrichen werden</div>` : ""}</span></li>`).join("")}</ul>`);
  }
  const bed = [...new Set([...(k.bedingungen || []), ...(k._bedingungen_auto || [])])];
  if (bed.length || k.bedingungen_text || (st === "bedingt" && k._wurzelbedingungen?.length)) {
    T.push(`<h3>Bedingungen</h3>`);
    if (bed.length) T.push(`<ul class="liste">${bed.map((b) => `<li>${knotenLink(b)}</li>`).join("")}</ul>`);
    if (k.bedingungen_text) T.push(`<p class="klein">${mathe(k.bedingungen_text)}</p>`);
    if (k._wurzelbedingungen?.length) T.push(`<p class="klein">Unbewiesene Enden der Beweiskette: ${k._wurzelbedingungen.map(knotenLink).join(", ")}</p>`);
  }
  if (k.beweis) {
    T.push(`<h3>Beweis</h3><p class="klein">${esc({ vollstaendig: "vollständig", literatur: "in der Literatur", skizze: "Skizze", bedingt: "bedingt", keiner: "keiner" }[k.beweis.art] || k.beweis.art)}</p>`);
    if (k.beweis.text) T.push(`<div class="mathe" style="font-size:15px">${mathe(k.beweis.text)}</div>`);
    if (k.beweis.luecken?.length) T.push(`<p class="klein" style="margin-top:8px"><b>Offene Beweisschritte</b></p><ul class="liste luecken">${k.beweis.luecken.map((l) => `<li>• ${mathe(l)}</li>`).join("")}</ul>`);
  }
  if (k._abhaengig?.length || k._traegt?.length) {
    T.push(`<h3>Beweisabhängigkeiten</h3>`);
    if (k._abhaengig?.length) T.push(`<p class="klein">stützt sich auf</p><ul class="liste">${k._abhaengig.map((d) => `<li>${knotenLink(d)}</li>`).join("")}</ul>`);
    if (k._traegt?.length) T.push(`<p class="klein" style="margin-top:6px">darauf bauen auf (${k._tragweite} insgesamt, transitiv)</p><ul class="liste">${k._traegt.map((d) => `<li>${knotenLink(d)}</li>`).join("")}</ul>`);
  }
  const weitere = G.kanten.filter((e) => (e.von === id || e.nach === id) && KANTENKLASSE[e.typ] !== "beweis");
  if (weitere.length) {
    T.push(`<h3>Weitere Beziehungen</h3><ul class="liste">${weitere.map((e) => e.von === id
      ? `<li><span class="klein">${esc(S.kanten[e.typ]?.label)}</span> ${knotenLink(e.nach)}</li>`
      : `<li>${knotenLink(e.von)} <span class="klein">${esc(S.kanten[e.typ]?.label)} diese Aussage</span></li>`).join("")}</ul>`);
  }
  if (k.gegenbeispielsuche) T.push(`<h3>Gegenbeispielsuche</h3><div class="mathe" style="font-size:15px">${mathe(k.gegenbeispielsuche)}</div>`);
  if (k.experimente?.length) { T.push(`<h3>Experimente</h3>`); for (const eid of k.experimente) T.push(experimentKarte(eid)); }
  if (k.quellen?.length) {
    T.push(`<h3>Quellen</h3>${k.quellen.map((q) => {
      const ref = typeof q === "string" ? q : q.ref;
      const d = Q[ref] || {};
      const gp = (typeof q === "object" && q.geprueft) || d.geprueft || "nein";
      return `<div class="quelle">${esc(d.autor || "")}: <i>${esc(d.titel || ref)}</i>. ${esc(d.ort || "")}${q.stelle ? ` — <b>${esc(q.stelle)}</b>` : ""}
        <span class="geprueft ${gp}">${esc(S.quellen_pruefung[gp])}</span>${d.notiz ? `<div class="klein">${esc(d.notiz)}</div>` : ""}</div>`;
    }).join("")}`);
  }
  if (k.bewertung) {
    const b = k.bewertung;
    T.push(`<h3>Bewertung</h3><dl class="bewertung">
      <dt>Richtigkeit</dt><dd>${mathe(b.richtigkeit)}</dd>
      <dt>Neuheit</dt><dd>${esc(S.neuheit[b.neuheit] || b.neuheit)}${b.recherche ? `<div class="klein">Recherche: ${mathe(b.recherche)}</div>` : ""}</dd>
      <dt>Nutzen</dt><dd>${b.nutzen_stufe ? `<b>${esc(b.nutzen_stufe)}</b> — ` : ""}${mathe(b.nutzen)}</dd></dl>`);
  }
  if (k.pruefung) T.push(`<h3>Die fünf Prüffragen</h3><dl class="fragen">${Object.entries(S.prueffragen).map(([f, t]) => k.pruefung[f] ? `<dt>${esc(t)}</dt><dd>${mathe(k.pruefung[f])}</dd>` : "").join("")}</dl>`);
  if (k.naechste_frage) T.push(`<h3>Nächste prüfbare Frage</h3><p>${mathe(k.naechste_frage)}</p>`);
  if (k.versionen?.length || k.nachgeprueft?.length) {
    T.push(`<h3>Fassungen und Nachprüfungen</h3><ul class="zeitleiste">${(k.versionen || []).map((v) => `<li><span class="datum">${esc(v.datum)}${v.fassung ? " · " + esc(v.fassung) : ""}</span><div>${mathe(v.aenderung)}${v.grund ? `<div class="klein">Grund: ${mathe(v.grund)}</div>` : ""}</div></li>`).join("")}
      ${(k.nachgeprueft || []).map((n) => `<li><span class="datum">${esc(n.datum)} · nachgeprüft gegen</span> ${knotenLink(n.gegen)}<div class="klein">${mathe(n.ergebnis)}</div></li>`).join("")}</ul>`);
  }
  if (k._verwandt?.length) {
    const max = Math.max(...k._verwandt.map((v) => v.score));
    T.push(`<h3>Stark verwandt</h3><ul class="liste verwandt">${k._verwandt.map((v) => `<li><span>${knotenLink(v.id)}</span><span class="balken" style="width:${(100 * v.score) / max}%"></span><span class="warum">${esc(v.warum.join(" · "))}</span></li>`).join("")}</ul>`);
  }
  if (k.befunde?.length) T.push(`<h3>Befunde des Prüfprogramms</h3>${k.befunde.map((b) => `<div class="befund ${b.stufe}">${esc(b.stufe)} · ${esc(b.regel)}: ${mathe(b.text)}</div>`).join("")}`);
  T.push(`<p class="klein" style="margin-top:22px">Datei: daten/aussagen/${esc(k._datei)}</p>`);
  tafelZeigen(T.join(""));
}
