// ===================================================================== Lehrbuch (Theorie-Reiter, Lesemodus)
// Kapitel aus daten/lehrbuch.yaml; Stränge werden zu ihren Aussagen aufgeklappt; jede Aussage zeigt ihren Theorie-Block.
const TEIL_NAMEN = (S.theorie_teile) || { worum: "Worum es geht", genau: "Genau", warum: "Warum es gilt", beispiel: "Beispiel", bilder: "Bilder", einordnung: "Einordnung" };

function theorieHtml(k, kompakt) {
  const t = k.theorie;
  if (!t) return "";
  const teile = [];
  for (const teil of ["worum", "genau", "warum", "beispiel", "bilder", "einordnung"]) {
    if (!t[teil]) continue;
    if (teil === "bilder") { teile.push((t.bilder || []).map((b) => bildHtml(b)).join("")); continue; }
    teile.push(`<div class="lb-teil lb-teil-${teil}"><h4>${esc(TEIL_NAMEN[teil])}</h4>${md(t[teil])}</div>`);
  }
  return `<div class="lb-theorie${kompakt ? " kompakt" : ""}">${teile.join("")}</div>`;
}

const buch = (() => {
  const wurzel = $("#buch");
  const kapitel = (G.lehrbuch && G.lehrbuch.kapitel) || [];
  const KAP_VON = new Map();          // Aussage -> Kapitel
  const IDS_VON = new Map();          // Kapitel -> Aussagen in Lesereihenfolge
  for (const kap of kapitel) {
    const ids = [];
    for (const t of kap.teile || []) {
      if (t.strang) (STR.get(t.strang)?.mitglieder || []).forEach((i) => ids.push(i));
      if (t.knoten) ids.push(t.knoten);
    }
    ids.forEach((i) => KAP_VON.set(i, kap.id));
    IDS_VON.set(kap.id, ids);
  }
  let aktuell = null;

  function fortschritt(kid) {
    const ids = IDS_VON.get(kid) || [];
    return { mit: ids.filter((i) => K.get(i)?.theorie).length, alle: ids.length };
  }

  function navHtml() {
    return `<div class="lb-nav-kopf">Inhalt</div><ol class="lb-nav-liste">${kapitel.map((kap, n) => {
      const f = fortschritt(kap.id);
      const offen = kap.id === aktuell;
      const unter = offen ? `<ol class="lb-nav-unter">${(kap.teile || []).filter((t) => t.strang || t.knoten || t.abschnitt).map((t) => {
        if (t.abschnitt) return `<li><button data-sprung="lb-abschnitt-${esc(t.anker || "")}">${esc(t.abschnitt)}</button></li>`;
        if (t.strang) return `<li><button data-sprung="lb-strang-${t.strang}">${esc(STR.get(t.strang)?.name || t.strang)}</button></li>`;
        return `<li><button data-sprung="lb-${t.knoten}">${esc(reinText(K.get(t.knoten)?.titel || t.knoten))}</button></li>`;
      }).join("")}</ol>` : "";
      return `<li class="${offen ? "offen" : ""}"><button data-kapitel="${kap.id}" aria-current="${offen ? "true" : "false"}">
        <span class="lb-nav-nr">${n}</span><span class="lb-nav-titel">${esc(kap.titel)}</span>
        ${f.alle ? `<span class="lb-nav-fort" title="${f.mit} von ${f.alle} Aussagen mit ausführlichem Text">${f.mit}/${f.alle}</span>` : ""}</button>${unter}</li>`;
    }).join("")}</ol>`;
  }

  function knotenBlock(id) {
    const k = K.get(id);
    if (!k) return "";
    const st = k._status_eff, geb = GEB.get(k.gebiet);
    const kopf = `<div class="lb-kopf">
        <span class="lb-art">${esc(S.arten[k.art])}</span>
        <span class="lb-status" style="--s:${farbe(st)}">${esc(S.status[st]?.label || st)}</span>
        <span class="lb-niveau niveau-${k._niveau}">${NIVEAU_NAME[k._niveau] || ""}</span>
        <button class="lb-details" data-details="${id}">Details und Quellen</button></div>`;
    let rumpf = "";
    if (k.theorie) {
      rumpf = `<div class="lb-aussage">${mathe(k.aussage)}</div>${theorieHtml(k)}`;
    } else {
      rumpf = `<div class="lb-aussage">${mathe(k.aussage)}</div>`
        + (k.einfach ? `<div class="lb-teil"><h4>In einfachen Worten</h4><p>${mathe(k.einfach)}</p></div>` : "")
        + (k.beweis?.text ? `<div class="lb-teil"><h4>Kurz begründet</h4><p>${mathe(k.beweis.text)}</p></div>` : "")
        + `<p class="lb-folgt">Der ausführliche Text zu dieser Aussage folgt in einer späteren Etappe.</p>`;
    }
    return `<section class="lb-knoten${k.theorie ? "" : " knapp"}" id="lb-${id}" style="--geb:${geb?.farbe || "var(--akzent)"}">
      ${kopf}<h3 class="lb-titel">${titelHtml(k.titel)}</h3>${rumpf}</section>`;
  }

  function kapitelHtml(kap, n) {
    const teile = [`<header class="lb-kapitel-kopf"><div class="lb-kapnr">Kapitel ${n}</div><h1>${esc(kap.titel)}</h1>
      ${kap.untertitel ? `<p class="lb-unter">${mathe(kap.untertitel)}</p>` : ""}</header>`];
    for (const t of kap.teile || []) {
      if (t.abschnitt) teile.push(`<h2 class="lb-abschnitt" id="lb-abschnitt-${esc(t.anker || "")}">${esc(t.abschnitt)}</h2>`);
      if (t.text) teile.push(`<div class="lb-fliess">${md(t.text)}</div>`);
      if (t.bild) teile.push(bildHtml(t.bild, t.unterschrift));
      if (t.strang) {
        const s = STR.get(t.strang);
        teile.push(`<h2 class="lb-strang" id="lb-strang-${t.strang}">${titelHtml(s?.name || t.strang)}</h2>${s?.verlauf ? `<p class="lb-verlauf">${mathe(s.verlauf)}</p>` : ""}`);
        if (t.text_vor) teile.push(`<div class="lb-fliess">${md(t.text_vor)}</div>`);
        (s?.mitglieder || []).forEach((i) => teile.push(knotenBlock(i)));
      }
      if (t.knoten) teile.push(knotenBlock(t.knoten));
    }
    const vor = kapitel[n - 1], nach = kapitel[n + 1];
    teile.push(`<nav class="lb-blaettern">${vor ? `<button data-kapitel="${vor.id}">← ${esc(vor.titel)}</button>` : "<span></span>"}
      ${nach ? `<button data-kapitel="${nach.id}">${esc(nach.titel)} →</button>` : ""}</nav>`);
    return teile.join("");
  }

  function zeigeKapitel(kid, ziel, sanft) {
    const n = kapitel.findIndex((k) => k.id === kid);
    if (n < 0) return;
    if (aktuell !== kid) {
      aktuell = kid;
      wurzel.innerHTML = `<nav class="lb-nav" aria-label="Inhalt des Lehrbuchs">${navHtml()}</nav>
        <article class="lb-seite" tabindex="-1">${kapitelHtml(kapitel[n], n)}</article>`;
    }
    const seite = $(".lb-seite", wurzel);
    if (ziel) {
      // erst nach dem Layout (das Detailfeld kann gerade aufgehen) rollen, sonst bricht der Browser das Rollen ab
      requestAnimationFrame(() => requestAnimationFrame(() => {
        const elz = document.getElementById(ziel);
        if (!elz) return;
        elz.scrollIntoView({ block: "start", behavior: sanft && BEWEGUNG ? "smooth" : "auto" });
        elz.classList.remove("lb-blitz"); void elz.offsetWidth; elz.classList.add("lb-blitz");
      }));
    } else seite.scrollTop = 0;
    zustand.kapitel = kid;
    hashSetzen();
  }

  wurzel.addEventListener("click", (e) => {
    const t = e.target.closest("[data-kapitel],[data-sprung],[data-details],[data-id]");
    if (!t) return;
    e.stopPropagation();
    const d = t.dataset;
    if (d.kapitel) zeigeKapitel(d.kapitel);
    else if (d.sprung) zeigeKapitel(aktuell, d.sprung, true);
    else if (d.details) { zustand.gewaehlt = d.details; zeigeKnoten(d.details); hervorheben(d.details); hashSetzen(); }
    else if (d.id) zuKnoten(d.id);
  });

  function hervorheben(id) {
    wurzel.querySelectorAll(".lb-knoten.gewaehlt").forEach((x) => x.classList.remove("gewaehlt"));
    if (id) document.getElementById("lb-" + id)?.classList.add("gewaehlt");
  }
  function zuKnoten(id) {
    const kid = KAP_VON.get(id);
    if (!kid) { zustand.gewaehlt = id; zeigeKnoten(id); return; }
    zeigeKapitel(kid, "lb-" + id);
    hervorheben(id);
  }
  return {
    zeige(kid) { wurzel.hidden = false; zeigeKapitel(kid || aktuell || zustand.kapitel || kapitel[0]?.id); },
    verberge() { wurzel.hidden = true; },
    zuKnoten, hervorheben,
    kapitelVon: (id) => KAP_VON.get(id),
    aktuell: () => aktuell,
  };
})();
