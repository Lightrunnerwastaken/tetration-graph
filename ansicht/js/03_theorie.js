// ===================================================================== Theorie (SVG)
// Definitionen, Sätze, Vermutungen in Beweisstufen: links Grundlagen, rechts Folgerungen.
const theorie = (() => {
  const svg = $("#leinwand");
  const defs = el("defs", {}, svg);
  for (const [id, v] of [["pfeil-beweis", "var(--kante)"], ["pfeil-rot", "var(--kante-rot)"], ["pfeil-gruen", "var(--kante-gruen)"],
    ["pfeil-evidenz", "var(--s-numerisch)"], ["pfeil-struktur", "var(--kante-leicht)"]]) {
    const m = el("marker", { id, viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" }, defs);
    el("path", { d: "M0,1 L10,5 L0,9 z", style: `fill:${v}` }, m);
  }
  const welt = el("g", {}, svg);
  let tf = { k: 1, x: 0, y: 0 };
  const setzeTf = () => welt.setAttribute("transform", `translate(${tf.x},${tf.y}) scale(${tf.k})`);
  function zoomUm(f, cx, cy) {
    const k = Math.max(0.15, Math.min(4, tf.k * f));
    tf.x = cx - ((cx - tf.x) * k) / tf.k; tf.y = cy - ((cy - tf.y) * k) / tf.k; tf.k = k;
    setzeTf();
  }
  svg.addEventListener("wheel", (e) => {
    e.preventDefault();
    const r = svg.getBoundingClientRect();
    zoomUm(Math.exp(-e.deltaY * 0.0015), e.clientX - r.left, e.clientY - r.top);
  }, { passive: false });
  let zug = null;
  svg.addEventListener("pointerdown", (e) => {
    if (e.target.closest(".knoten")) return;
    zug = { x: e.clientX, y: e.clientY, tx: tf.x, ty: tf.y, bewegt: false };
    svg.setPointerCapture(e.pointerId);
    svg.classList.add("zieht");
  });
  svg.addEventListener("pointermove", (e) => {
    if (!zug) return;
    const dx = e.clientX - zug.x, dy = e.clientY - zug.y;
    if (Math.abs(dx) + Math.abs(dy) > 3) zug.bewegt = true;
    tf.x = zug.tx + dx; tf.y = zug.ty + dy;
    setzeTf();
  });
  svg.addEventListener("pointerup", () => {
    if (zug && !zug.bewegt && zustand.gewaehlt) waehle(null);
    zug = null;
    svg.classList.remove("zieht");
  });

  let gez = { knoten: new Map(), kanten: [] };
  let lay = null;
  function layout() {
    if (lay) return lay;
    const alle = G.knoten.filter((k) => THEORIE_ARTEN.has(k.art)).map((k) => k.id);
    const alleSet = new Set(alle);
    const hatNachfolger = new Set();
    alle.forEach((i) => (K.get(i)._abhaengig || []).forEach((d) => alleSet.has(d) && hatNachfolger.add(d)));
    const isoliert = alle.filter((i) => !(K.get(i)._abhaengig || []).some((d) => alleSet.has(d)) && !hatNachfolger.has(i));
    const ids = alle.filter((i) => !isoliert.includes(i));
    const set = new Set(ids);
    const deps = new Map(ids.map((i) => [i, (K.get(i)._abhaengig || []).filter((d) => set.has(d))]));
    const stufe = new Map();
    const tiefe = (i, weg = new Set()) => {
      if (stufe.has(i)) return stufe.get(i);
      if (weg.has(i)) return 0;
      weg.add(i);
      const d = deps.get(i);
      const s = d.length ? 1 + Math.max(...d.map((x) => tiefe(x, weg))) : 0;
      stufe.set(i, s);
      return s;
    };
    ids.forEach((i) => tiefe(i));
    const maxS = Math.max(...stufe.values());
    const schichten = Array.from({ length: maxS + 1 }, () => []);
    const rang = new Map(G.gebiete.map((g, i) => [g.id, i]));
    ids.forEach((i) => schichten[stufe.get(i)].push(i));
    schichten.forEach((s) => s.sort((a, b) => rang.get(K.get(a).gebiet) - rang.get(K.get(b).gebiet) || a.localeCompare(b)));
    const ZEILEN = 13, SPALTE = 262, ZEILE = 38;
    const spaltenStart = [];
    let xLauf = 0;
    schichten.forEach((s, i) => { spaltenStart[i] = xLauf; xLauf += Math.max(1, Math.ceil(s.length / ZEILEN)) * SPALTE; });
    const pos = new Map();
    const setze = () => schichten.forEach((s, i) => {
      const sp = Math.max(1, Math.ceil(s.length / ZEILEN)), pro = Math.ceil(s.length / sp);
      s.forEach((id, j) => {
        const c = Math.floor(j / pro), rr = j % pro, n = Math.min(pro, s.length - c * pro);
        pos.set(id, { x: spaltenStart[i] + c * SPALTE, y: (rr - (n - 1) / 2) * ZEILE });
      });
    });
    setze();
    for (let runde = 0; runde < 6; runde++) {
      for (const s of schichten.slice(1)) {
        const m = new Map(s.map((i) => { const d = deps.get(i); return [i, d.length ? d.reduce((u, x) => u + pos.get(x).y, 0) / d.length : pos.get(i).y]; }));
        s.sort((a, b) => m.get(a) - m.get(b));
      }
      setze();
    }
    const maxY = Math.max(...[...pos.values()].map((p) => p.y));
    const rasterY = maxY + 150;
    isoliert.sort((a, b) => rang.get(K.get(a).gebiet) - rang.get(K.get(b).gebiet) || a.localeCompare(b));
    const proZeile = Math.max(3, Math.floor(xLauf / SPALTE));
    isoliert.forEach((id, j) => pos.set(id, { x: (j % proZeile) * SPALTE, y: rasterY + Math.floor(j / proZeile) * ZEILE }));
    lay = { ids: [...ids, ...isoliert], pos, schichten, spaltenStart, rasterY: isoliert.length ? rasterY : null };
    return lay;
  }

  function zeichneKnoten(id, p, parent, r) {
    const k = K.get(id);
    const g = el("g", { class: "knoten", transform: `translate(${p.x},${p.y})`, tabindex: 0, role: "button", "aria-label": k.titel }, parent);
    g.dataset.id = id;
    if (k._neu_pruefen?.length) el("circle", { r: r + 6, style: "fill:none;stroke:var(--s-widerlegt);stroke-width:1.6;stroke-dasharray:3 2.5" }, g);
    const f = form(k.art, r, g);
    const st = stilFuer(k);
    f.setAttribute("style", `fill:${st.fill};stroke:${GEB.get(k.gebiet)?.farbe};stroke-width:2.4;${st.strich ? "stroke-dasharray:" + st.strich : ""}`);
    const t = el("text", { x: r + 6, dy: "0.35em", class: "etikett" }, g);
    t.textContent = kurz(reinText(k.titel), 38);
    g.addEventListener("click", (e) => { e.stopPropagation(); waehle(id); });
    g.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); waehle(id); } });
    gez.knoten.set(id, g);
  }

  function zeichne() {
    welt.innerHTML = "";
    gez = { knoten: new Map(), kanten: [] };
    const { ids, pos, schichten, spaltenStart, rasterY } = layout();
    const rT = (i) => Math.min(14, 5 + 1.8 * Math.sqrt(K.get(i)._tragweite || 0));
    const oben = Math.min(...[...pos.values()].filter((p) => rasterY == null || p.y < rasterY - 50).map((p) => p.y));
    schichten.forEach((s, i) => {
      if (!s.length) return;
      el("text", { x: spaltenStart[i] - 6, y: oben - 34, class: "schicht-text" }, welt).textContent = i === 0 ? "Stufe 0: Grundlagen" : `Stufe ${i}`;
    });
    if (rasterY != null) el("text", { x: -6, y: rasterY - 34, class: "schicht-text" }, welt).textContent = "Ohne Beweisverbindung zu anderen Theorie-Aussagen";
    const gk = el("g", {}, welt);
    const set = new Set(ids);
    for (const e of G.kanten) {
      if (!set.has(e.von) || !set.has(e.nach)) continue;
      const kl = KANTENKLASSE[e.typ] || "leicht";
      if (!["beweis", "gruen", "rot"].includes(kl)) continue;
      const [a, b] = kl === "beweis" ? [e.nach, e.von] : [e.von, e.nach];
      const pa = pos.get(a), pb = pos.get(b), x1 = pa.x + rT(a) + 3, x2 = pb.x - rT(b) - 5, mx = (x1 + x2) / 2;
      const p = el("path", { d: `M${x1},${pa.y} C${mx},${pa.y} ${mx},${pb.y} ${x2},${pb.y}`, class: `kante ${kl}`, "marker-end": `url(#pfeil-${kl})` }, gk);
      p.dataset.von = e.von; p.dataset.nach = e.nach;
      gez.kanten.push(p);
    }
    const gn = el("g", {}, welt);
    ids.forEach((i) => zeichneKnoten(i, pos.get(i), gn, rT(i)));
    filter();
    const bb = welt.getBBox(), W = svg.clientWidth, H = svg.clientHeight;
    const k = Math.max(0.62, Math.min((W - 80) / bb.width, (H - 130) / bb.height, 1.4));
    tf = k <= 0.62 ? { k, x: 24 - k * bb.x, y: 56 - k * bb.y } : { k, x: W / 2 - k * (bb.x + bb.width / 2), y: (H + 30) / 2 - k * (bb.y + bb.height / 2) };
    setzeTf();
  }
  function filter() {
    for (const [id, g] of gez.knoten) g.style.display = sichtbar(K.get(id)) ? "" : "none";
    for (const p of gez.kanten) {
      const a = gez.knoten.get(p.dataset.von), b = gez.knoten.get(p.dataset.nach);
      p.style.display = a?.style.display === "none" || b?.style.display === "none" ? "none" : "";
    }
  }
  function hervorheben(id) {
    for (const g of gez.knoten.values()) g.classList.remove("gedimmt", "gewaehlt");
    for (const p of gez.kanten) p.classList.remove("gedimmt", "hervor");
    if (!id || !gez.knoten.has(id)) return;
    const nah = new Set([id]);
    for (const p of gez.kanten) {
      if (p.dataset.von === id || p.dataset.nach === id) { nah.add(p.dataset.von); nah.add(p.dataset.nach); p.classList.add("hervor"); }
      else p.classList.add("gedimmt");
    }
    for (const [i, g] of gez.knoten) if (!nah.has(i)) g.classList.add("gedimmt");
    const g = gez.knoten.get(id);
    g.classList.add("gewaehlt");
    const m = /translate\(([-\d.]+),([-\d.]+)\)/.exec(g.getAttribute("transform"));
    const x = +m[1] * tf.k + tf.x, y = +m[2] * tf.k + tf.y, W = svg.clientWidth, H = svg.clientHeight;
    if (x < 60 || x > W - 60 || y < 60 || y > H - 60) { tf.x += W / 2 - x; tf.y += H / 2 - y; setzeTf(); }
  }
  return {
    zeige() { svg.hidden = false; zeichne(); hervorheben(zustand.gewaehlt); },
    verberge() { svg.hidden = true; },
    hervorheben, filter,
    zoom(f) { zoomUm(f, svg.clientWidth / 2, svg.clientHeight / 2); },
    passt: (id) => THEORIE_ARTEN.has(K.get(id)?.art),
  };
})();
