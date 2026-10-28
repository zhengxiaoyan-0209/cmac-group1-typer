/* =======================================================================
   app.js — UI wiring
   Markers arrive through <script src="assets/markers.js"> (window.CMAC_MARKERS)
   rather than a fetch call, so opening index.html straight off disk (file://) works.

   One file = one sample. A whole folder of assemblies can be dropped in and each
   one is typed separately, then collected into a summary table.
   ======================================================================= */
(function () {
  "use strict";
  const db = window.CMAC_MARKERS;
  const core = window.CMACCore;
  const $ = function (id) { return document.getElementById(id); };

  let samples = [];      /* {name, nContigs, len, contigs} */
  let results = [];      /* {name, res, bp, n} */
  let sel = 0;           /* index of the sample currently shown */

  /* ------------------------------------------------------ header meta */
  function fmt(n) { return n.toLocaleString("en-US"); }

  $("refs").innerHTML =
    "Reference <code>" + db.meta.reference.split(" ")[0] + "</code> &middot; " +
    fmt(db.meta.ref_length) + " bp<br>" +
    "Panel: " + db.meta.n_markers + " group-specific SNPs &middot; " +
    "cohort n=" + (db.meta.n_group_members["1"] + db.meta.n_group_members["2"] + db.meta.n_group_members["3"]);

  $("methodMeta").innerHTML =
    "Reference genome: " + db.meta.reference + " (" + fmt(db.meta.ref_length) + " bp).<br>" +
    "The 85 Group 1&ndash;specific SNPs were derived from an alignment of " +
    (db.meta.n_group_members["1"] + db.meta.n_group_members["2"] + db.meta.n_group_members["3"]) +
    " genomes (Group 1 " + db.meta.n_group_members["1"] +
    " / Group 2 " + db.meta.n_group_members["2"] +
    " / Group 3 " + db.meta.n_group_members["3"] + "), " +
    "selected for presence in 100% of Group 1 genomes and 0% of all others. " +
    "Group 2 and Group 3 yield no usable locus under that criterion, so this panel only calls Group 1.";

  $("footMeta").innerHTML =
    db.meta.note + " &nbsp;|&nbsp; Sequences never leave your browser.";

  /* ------------------------------------------------------ file input */
  const drop = $("drop"), fileInput = $("file");

  drop.addEventListener("click", function () { fileInput.click(); });
  ["dragenter", "dragover"].forEach(function (e) {
    drop.addEventListener(e, function (ev) { ev.preventDefault(); drop.classList.add("over"); });
  });
  ["dragleave", "drop"].forEach(function (e) {
    drop.addEventListener(e, function (ev) { ev.preventDefault(); drop.classList.remove("over"); });
  });
  drop.addEventListener("drop", function (ev) {
    const items = ev.dataTransfer.items || [];
    for (let i = 0; i < items.length; i++) {
      const it = items[i].webkitGetAsEntry && items[i].webkitGetAsEntry();
      if (it && it.isDirectory) walkDir(it);
    }
    handleFiles(ev.dataTransfer.files);
  });
  fileInput.addEventListener("change", function () {
    handleFiles(fileInput.files);
    fileInput.value = "";
  });

  function walkDir(entry) {
    entry.createReader().readEntries(function (es) {
      es.forEach(function (en) {
        if (en.isDirectory) walkDir(en);
        else en.file(function (f) { handleFiles([f]); });
      });
    });
  }

  async function readText(file) {
    if (/\.gz$/i.test(file.name) && typeof DecompressionStream === "function") {
      const ds = new DecompressionStream("gzip");
      const stream = file.stream().pipeThrough(ds);
      return new Response(stream).text();
    }
    return file.text();
  }

  async function handleFiles(list) {
    for (let i = 0; i < list.length; i++) {
      const f = list[i];
      try {
        const text = await readText(f);
        addSample(f.name, text);
      } catch (e) {
        samples.push({ name: f.name, error: true });
        renderSamples();
      }
    }
    renderSamples();
  }

  function addSample(name, text) {
    const parsed = core.parseFasta(text);
    if (!parsed.contigs.length) return;
    samples.push({
      name: name,
      nContigs: parsed.contigs.length,
      len: parsed.contigs.reduce(function (a, c) { return a + c.seq.length; }, 0),
      contigs: parsed.contigs
    });
  }

  $("addPaste").addEventListener("click", function () {
    const t = $("paste").value.trim();
    if (!t) return;
    addSample("Pasted sequence " + (samples.length + 1), t);
    $("paste").value = "";
    renderSamples();
  });

  function totalBp() {
    return samples.reduce(function (a, s) { return a + (s.len || 0); }, 0);
  }

  function renderSamples() {
    const el = $("filelist");
    if (!samples.length) { el.hidden = true; $("run").disabled = true; return; }
    el.hidden = false;
    el.innerHTML = samples.map(function (s, i) {
      return '<div class="row"><div>' +
        '<div class="nm">' + s.name + "</div>" +
        '<div class="mt">' + (s.error ? "Could not be read (unsupported gzip?)"
          : fmt(s.nContigs) + " sequences &middot; " + fmt(s.len) + " bp") +
        "</div></div>" +
        '<button class="x" data-i="' + i + '" title="Remove">&times;</button></div>';
    }).join("") +
      '<div class="row" style="background:var(--panel2)"><div class="nm">' +
      samples.length + " samples in total</div>" +
      '<div class="mt">' + fmt(totalBp()) + " bp</div></div>";
    el.querySelectorAll(".x").forEach(function (b) {
      b.addEventListener("click", function () {
        samples.splice(+b.dataset.i, 1); renderSamples();
      });
    });
    $("run").disabled = totalBp() === 0;
  }

  $("clear").addEventListener("click", function () {
    samples = []; results = []; sel = 0;
    renderSamples();
    $("resultCard").hidden = true;
    $("detailCard").hidden = true;
  });

  /* ----------------------------------------------------- main routine */
  const yieldUI = function () {
    return new Promise(function (r) { setTimeout(r, 0); });
  };

  function opts() {
    return {
      minCov: Math.max(0, Math.min(100, +$("minCov").value || 0)) / 100,
      hi: Math.max(0, Math.min(100, +$("hi").value || 0)) / 100,
      lo: Math.max(0, Math.min(100, +$("lo").value || 0)) / 100
    };
  }

  async function runAnalysis() {
    const o = opts();
    if (!samples.length) return;

    $("run").disabled = true;
    $("progress").hidden = false;
    results = []; sel = 0;

    for (let si = 0; si < samples.length; si++) {
      const s = samples[si];
      $("progtxt").textContent = "Sample " + (si + 1) + "/" + samples.length + " · " + s.name;
      $("bar").style.width = (si / samples.length * 100).toFixed(1) + "%";
      await yieldUI();

      const ctx = core.prepareQuery(s.contigs, db, o);
      const rows = [];
      const CH = 8;
      for (let i = 0; i < db.snps.length; i += CH) {
        const end = Math.min(db.snps.length, i + CH);
        for (let j = i; j < end; j++) {
          rows.push(core.classifyRow(db.snps[j], core.callLocus(db.snps[j], ctx)));
        }
        await yieldUI();
      }
      results.push({ name: s.name, res: core.summarize(rows, db, o), bp: s.len, n: s.nContigs });
    }

    $("bar").style.width = "100%";
    $("progress").hidden = true;
    $("run").disabled = false;
    renderAll();
  }

  $("run").addEventListener("click", runAnalysis);

  /* ---------------------------------------------------------- render */
  const ICON = { pos: "✓", neg: "✕", warn: "!" };
  const LABEL = { GROUP1: "Group 1 allele", NON: "Reference allele", OTHER: "Third allele", NA: "Undetermined" };
  const LEVEL = { pos: "Group 1", neg: "Not Group 1", warn: "Inconclusive" };

  function renderAll() {
    $("resultCard").hidden = false;
    $("detailCard").hidden = false;
    if (sel >= results.length) sel = 0;

    const many = results.length > 1;
    $("summaryWrap").hidden = !many;
    if (many) renderSummaryTable();

    renderVerdict();
    renderTabs();
  }

  function renderSummaryTable() {
    const tb = $("tblSummary").querySelector("tbody");
    tb.innerHTML = results.map(function (r, i) {
      const s = r.res.stats;
      return '<tr class="' + (i === sel ? "sel" : "") + '" data-i="' + i + '" style="cursor:pointer">' +
        "<td><b>" + r.name + "</b></td>" +
        '<td><span class="tag ' + (r.res.level === "pos" ? "GROUP1" : r.res.level === "neg" ? "NON" : "NA") + '">' +
        LEVEL[r.res.level] + "</span></td>" +
        '<td class="num">' + (s.score * 100).toFixed(1) + "%</td>" +
        '<td class="num">' + s.called + " / " + s.n_total + "</td>" +
        '<td class="num">' + (r.bp / 1e6).toFixed(2) + " Mb</td>" +
        "</tr>";
    }).join("");
    tb.querySelectorAll("tr").forEach(function (tr) {
      tr.addEventListener("click", function () {
        sel = +tr.dataset.i;
        renderAll();
      });
    });
  }

  function renderVerdict() {
    const r = results[sel], s = r.res.stats;
    const v = $("verdict");
    v.className = "verdict " + r.res.level;
    $("vicon").textContent = ICON[r.res.level];
    $("vmain").textContent = r.res.verdict;
    $("vsub").textContent =
      r.name + " · match score " + (s.score * 100).toFixed(1) + "% · " +
      s.called + "/" + s.n_total + " loci readable · " +
      "thresholds ≥" + (r.res.thresholds.hi * 100).toFixed(0) + "% / ≤" +
      (r.res.thresholds.lo * 100).toFixed(0) + "%";

    $("mScore").textContent = (s.score * 100).toFixed(1) + "%";
    $("mScoreBar").style.width = (s.score * 100).toFixed(1) + "%";
    $("mCov").textContent = s.called + " / " + s.n_total;
    $("mCovNote").textContent = "coverage " + (s.coverage * 100).toFixed(1) +
      "% (minimum " + (r.res.thresholds.minCov * 100).toFixed(0) + "%)";
    $("mSize").textContent = (r.bp / 1e6).toFixed(2) + " Mb";
    $("mSizeNote").textContent = r.n + " sequences";

    $("lg1").textContent = s.n_g1;
    $("lnon").textContent = s.n_non;
    $("lother").textContent = s.n_other;
    $("lna").textContent = s.n_na;
  }

  function renderTabs() {
    renderTable();
    renderGeneSummary();
  }

  function renderTable() {
    const r = results[sel].res;
    const q = ($("q").value || "").toLowerCase();
    const f = $("fv").value;
    const tb = $("tbl").querySelector("tbody");
    let n = 0;
    tb.innerHTML = r.rows.map(function (x) {
      if (f && x.verdict !== f) return "";
      if (q && (String(x.pos) + x.gene + x.locus_tag + x.product).toLowerCase().indexOf(q) < 0) return "";
      n++;
      return "<tr>" +
        '<td class="num">' + fmt(x.pos) + "</td>" +
        '<td><span class="tag ' + x.verdict + '">' + LABEL[x.verdict] + "</span></td>" +
        '<td class="mono">' + x.ref + "</td>" +
        '<td class="mono"><b>' + x.alt + "</b></td>" +
        '<td class="mono">' + x.obs + "</td>" +
        "<td>" + x.gene + "</td>" +
        '<td class="mono">' + x.locus_tag + "</td>" +
        '<td class="prod">' + x.product + "</td>" +
        '<td class="mono">' + (x.verdict === "NA" ? (x.note || "-") : (x.side + x.anchor_len)) + "</td>" +
        "</tr>";
    }).join("");
    $("tblcount").textContent = "Showing " + n + " of " + r.rows.length + " loci";
  }

  function renderGeneSummary() {
    const tb = $("tblGeneSum").querySelector("tbody");
    tb.innerHTML = db.gene_summary.map(function (g) {
      return "<tr><td><b>" + g.locus + "</b></td><td class='prod'>" + g.product + "</td>" +
        '<td class="num">' + g.n_snp + "</td>" +
        '<td class="mono">' + g.span + "</td><td>" + g.category + "</td></tr>";
    }).join("");
  }

  /* ---------------------------------------------------- interaction */
  $("q").addEventListener("input", function () { if (results.length) renderTable(); });
  $("fv").addEventListener("change", function () { if (results.length) renderTable(); });
  ["minCov", "hi", "lo"].forEach(function (id) {
    $(id).addEventListener("change", function () { if (samples.length) runAnalysis(); });
  });

  $("tabs").addEventListener("click", function (e) {
    const b = e.target.closest(".tab");
    if (!b) return;
    $("tabs").querySelectorAll(".tab").forEach(function (t) { t.classList.remove("active"); });
    b.classList.add("active");
    ["snp", "geneSum", "method"].forEach(function (k) {
      $("pane-" + k).hidden = k !== b.dataset.tab;
    });
  });

  /* --------------------------------------------------------- export */
  function download(text, name) {
    const blob = new Blob(["\uFEFF" + text], { type: "text/csv;charset=utf-8" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  $("csvSummary").addEventListener("click", function () {
    if (!results.length) return;
    const lines = [["sample", "result", "match_score_pct", "loci_typed", "group1_allele",
      "reference_allele", "third_allele", "undetermined", "sequences", "total_bp"].join(",")];
    results.forEach(function (r) {
      const s = r.res.stats;
      lines.push([r.name, LEVEL[r.res.level], (s.score * 100).toFixed(2),
        s.called + "/" + s.n_total, s.n_g1, s.n_non, s.n_other, s.n_na,
        r.n, r.bp].join(","));
    });
    lines.push("", "# reference," + db.meta.reference,
      "# thresholds,Group1>=" + (results[0].res.thresholds.hi * 100).toFixed(0) +
      "% / non-Group1<=" + (results[0].res.thresholds.lo * 100).toFixed(0) +
      "% / min coverage " + (results[0].res.thresholds.minCov * 100).toFixed(0) + "%");
    download(lines.join("\n"), "group1_typing_summary.csv");
  });

  $("csv").addEventListener("click", function () {
    if (!results.length) return;
    const r = results[sel];
    const head = ["sample", "position", "call", "ref_allele", "group1_allele", "observed",
      "gene", "locus_tag", "product", "feature", "region_kb", "note"];
    const lines = [head.join(",")];
    r.res.rows.forEach(function (x) {
      lines.push([r.name, x.pos, LABEL[x.verdict], x.ref, x.alt, x.obs,
        x.gene, x.locus_tag, '"' + x.product + '"', x.feature, x.region_kb,
        '"' + (x.note || "") + '"'].join(","));
    });
    lines.push("", "# result," + LEVEL[r.res.level] + " - " + r.res.verdict,
      "# match_score," + (r.res.stats.score * 100).toFixed(2) + "%",
      "# loci_typed," + r.res.stats.called + "/" + r.res.stats.n_total,
      "# reference," + db.meta.reference);
    download(lines.join("\n"), "group1_typing_" + r.name.replace(/[^\w.-]+/g, "_") + ".csv");
  });
})();
