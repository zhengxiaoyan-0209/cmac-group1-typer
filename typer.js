#!/usr/bin/env node
/* =======================================================================
   typer.js — command-line batch version (shares assets/core.js with the web page)

      node scripts/typer.js genome1.fna genome2.fa.gz -o out.csv

   Options
     --min-cov 70   minimum locus coverage (%)
     --hi 90        Group 1 call threshold (%)
     --lo 10        non-Group 1 call threshold (%)
     --detail       also write a per-locus CSV (out.detail.csv)
     -o FILE        write summary to FILE (default: stdout)
   ======================================================================= */
"use strict";
const fs = require("fs");
const path = require("path");
const zlib = require("zlib");

const BASE = path.join(__dirname, "..");
const core = require(path.join(BASE, "assets", "core.js"));
const db = JSON.parse(fs.readFileSync(path.join(BASE, "data", "markers.json"), "utf8"));

/* -------------------------------- arguments -------------------------- */
const argv = process.argv.slice(2);
const files = [];
let out = null, detail = false;
const opts = { minCov: 0.70, hi: 0.90, lo: 0.10 };

for (let i = 0; i < argv.length; i++) {
  const a = argv[i];
  if (a === "-o") { out = argv[++i]; }
  else if (a === "--detail") { detail = true; }
  else if (a === "--min-cov") { opts.minCov = +argv[++i] / 100; }
  else if (a === "--hi") { opts.hi = +argv[++i] / 100; }
  else if (a === "--lo") { opts.lo = +argv[++i] / 100; }
  else if (a.charAt(0) === "-") { console.error("Unknown option: " + a); process.exit(2); }
  else files.push(a);
}

if (!files.length) {
  console.error("Usage: node scripts/typer.js <genome.fa> [genome2.fna.gz ...] [-o out.csv]");
  process.exit(2);
}

/* -------------------------------- read file -------------------------- */
function readFasta(file) {
  let raw;
  if (/\.gz$/i.test(file)) {
    raw = zlib.gunzipSync(fs.readFileSync(file));
  } else {
    raw = fs.readFileSync(file);
  }
  return raw.toString("utf8");
}

/* -------------------------------- main loop -------------------------- */
const summary = [];
const detailRows = [];
const LABEL = { GROUP1: "Group 1 allele", NON: "Reference allele", OTHER: "Third allele", NA: "Undetermined" };

files.forEach(function (f) {
  let parsed;
  try {
    parsed = core.parseFasta(readFasta(f));
  } catch (e) {
    console.error("Could not read " + f + ": " + e.message);
    return;
  }
  if (!parsed.contigs.length) {
    console.error("No sequence found in " + f);
    return;
  }
  const bp = parsed.contigs.reduce(function (a, c) { return a + c.seq.length; }, 0);
  const t0 = Date.now();
  const r = core.typeGenome(parsed.contigs, db, opts);
  const ms = Date.now() - t0;

  summary.push({
    file: path.basename(f),
    verdict: r.verdict,
    level: r.level,
    score: (r.stats.score * 100).toFixed(2),
    called: r.stats.called + "/" + r.stats.n_total,
    g1: r.stats.n_g1, non: r.stats.n_non,
    other: r.stats.n_other, na: r.stats.n_na,
    contigs: parsed.contigs.length,
    bp: bp,
    ms: ms
  });

  if (detail) {
    r.rows.forEach(function (x) {
      detailRows.push([
        path.basename(f), x.pos, LABEL[x.verdict], x.ref, x.alt, x.obs,
        x.gene, x.locus_tag, '"' + x.product + '"'
      ].join(","));
    });
  }

  process.stderr.write(
    "  " + path.basename(f).padEnd(28) +
    String(r.verdict).padEnd(30) +
    "score " + (r.stats.score * 100).toFixed(1) + "%  " +
    "typed " + r.stats.called + "/" + r.stats.n_total +
    "  (" + ms + " ms)\n"
  );
});

/* -------------------------------- output ----------------------------- */
const head = ["sample", "result", "match_score_pct", "loci_typed", "group1_allele",
  "reference_allele", "third_allele", "undetermined", "sequences", "total_bp", "ms"];
const lines = ["\uFEFF" + head.join(",")];
summary.forEach(function (s) {
  lines.push([s.file, s.verdict, s.score, s.called, s.g1, s.non, s.other, s.na,
    s.contigs, s.bp, s.ms].join(","));
});
const text = lines.join("\n");
if (out) {
  fs.writeFileSync(out, text, "utf8");
  console.error("\nSummary -> " + out);
} else {
  console.log(text);
}

if (detail && detailRows.length) {
  const dpath = out ? out.replace(/\.csv$/, ".detail.csv") : "typer_detail.csv";
  const dh = ["sample", "position", "call", "ref_allele", "group1_allele", "observed",
    "gene", "locus_tag", "product"].join(",");
  fs.writeFileSync(dpath, "\uFEFF" + dh + "\n" + detailRows.join("\n"), "utf8");
  console.error("Locus detail -> " + dpath);
}
