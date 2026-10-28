/* =======================================================================
   validate.js - offline self-consistency checks for the typing engine

     node scripts/validate.js

   Checks
   ------
   T1  the reference itself (FDAARGOS_1195 belongs to Group 2) -> must be "not Group 1"
   T2  the 85 loci replaced by their Group 1 allele in silico -> "Group 1", 85/85
   T3  T2 plus random background differences -> recall must stay high
   T4  the same genome broken into contigs, orientations shuffled -> unaffected
   T5  only half of the loci converted -> must land in the inconclusive band
   T6  the whole reference fed in as its reverse complement -> unaffected
   ======================================================================= */
"use strict";
const fs = require("fs");
const path = require("path");
const crypto = require("crypto");

const BASE = path.join(__dirname, "..");
const core = require(path.join(BASE, "assets", "core.js"));
const db = JSON.parse(fs.readFileSync(path.join(BASE, "data", "markers.json"), "utf8"));

const refPath = path.join(BASE, "data", "reference_CP069516.1.fna");
const refSeq = fs
  .readFileSync(refPath, "utf8")
  .split(/\r?\n/)
  .filter((l) => !l.startsWith(">"))
  .join("")
  .toUpperCase();

const OPTS = { minCov: 0.70, hi: 0.90, lo: 0.10 };

/* --------------------------- reproducible pseudo-random numbers -------- */
function makeRng(seed) {
  let h = crypto.createHash("sha256").update(String(seed)).digest();
  let i = 0;
  return function () {
    if (i >= h.length - 4) {
      h = crypto.createHash("sha256").update(h).digest();
      i = 0;
    }
    const v = h.readUInt32BE(i);
    i += 4;
    return v / 4294967296;
  };
}

/* ------------------------- building modified genomes ------------------- */
function applyAlt(seq, positions) {
  const arr = seq.split("");
  positions.forEach(function (p) { arr[p - 1] = p.__a; });
  return arr.join("");
}

function buildAltGenome(fraction) {
  const arr = refSeq.split("");
  db.snps.forEach(function (m, i) {
    if (i < Math.round(db.snps.length * fraction)) arr[m.pos - 1] = m.alt;
  });
  return arr.join("");
}

function addNoise(seq, rate, seed) {
  const rng = makeRng(seed);
  const B = ["A", "C", "G", "T"];
  const arr = seq.split("");
  let n = 0;
  for (let i = 0; i < arr.length; i++) {
    if (rng() < rate) {
      const b = B[Math.floor(rng() * 4)];
      if (b !== arr[i]) { arr[i] = b; n++; }
    }
  }
  return { seq: arr.join(""), n: n };
}

function shard(seq, meanLen, seed) {
  const rng = makeRng(seed);
  const contigs = [];
  let i = 0, k = 0;
  while (i < seq.length) {
    const L = Math.max(2000, Math.round(meanLen * (0.4 + 1.2 * rng())));
    const piece = seq.slice(i, i + L);
    const flip = rng() < 0.5;
    contigs.push({
      name: "contig_" + (++k) + (flip ? "_rc" : ""),
      seq: flip ? core.revcomp(piece) : piece
    });
    i += L;
  }
  return contigs;
}

/* -------------------------------- timing and assertions ---------------- */
let pass = 0, fail = 0;
function check(name, cond, detail) {
  if (cond) { pass++; console.log("  [PASS] " + name + (detail ? "  — " + detail : "")); }
  else { fail++; console.log("  [FAIL] " + name + (detail ? "  — " + detail : "")); }
}

function run(label, seq, extra) {
  const t0 = Date.now();
  const r = core.typeGenome([{ name: label, seq: seq }], db, OPTS);
  const ms = Date.now() - t0;
  console.log(
    "\n" + label +
    "\n   call=" + r.verdict +
    "  | typed " + r.stats.called + "/" + r.stats.n_total +
    "  | score " + (r.stats.score * 100).toFixed(1) + "%" +
    "  | G1=" + r.stats.n_g1 + " NON=" + r.stats.n_non +
    " OTHER=" + r.stats.n_other + " NA=" + r.stats.n_na +
    "  | " + ms + " ms" + (extra ? "\n   " + extra : "")
  );
  return r;
}

console.log("=".repeat(74));
console.log("C. macginleyi Group 1 typing engine - offline validation");
console.log("Reference: " + db.meta.reference + "  " + db.meta.ref_length + " bp");
console.log("Loci in panel: " + db.meta.n_markers);
console.log("=".repeat(74));

/* ---- T1: the reference itself belongs to Group 2 ---- */
const t1 = run("T1 Reference itself (FDAARGOS_1195, Group 2)", refSeq);
check("called not Group 1", t1.level === "neg", t1.verdict);
check("match score = 0", t1.stats.score === 0, (t1.stats.score * 100).toFixed(1) + "%");
check("all 85 loci readable", t1.stats.called === 85, t1.stats.called + "/85");
check("no unexpected third alleles", t1.stats.n_other === 0, "third=" + t1.stats.n_other);

/* ---- T2: Group 1 assembled in silico ---- */
const alt = buildAltGenome(1.0);
const t2 = run("T2 Group 1 in silico (all 85 loci set to the Group 1 allele)", alt);
check("called Group 1", t2.level === "pos", t2.verdict);
check("match score = 100%", t2.stats.score === 1, (t2.stats.score * 100).toFixed(1) + "%");
check("85/85 carry the Group 1 allele", t2.stats.n_g1 === 85, t2.stats.n_g1 + "/85");

/* ---- T3: random background differences ---- */
[0.001, 0.005, 0.01].forEach(function (rate) {
  const nz = addNoise(alt, rate, "noise" + rate);
  const t = run(
    "T3 Group 1 + " + (rate * 100).toFixed(1) + "% random noise (" + nz.n + " extra mutations)",
    nz.seq
  );
  check("called Group 1", t.level === "pos", t.verdict);
  check("at least 90% of loci readable", t.stats.called >= 77, t.stats.called + "/85");
  /* 1% divergence (one change every ~132 bp) is far beyond what is seen within
     C. macginleyi; here we only require the call itself to stay correct. */
  check("match score >= 95%", t.stats.score >= 0.95, (t.stats.score * 100).toFixed(1) + "%");
});

/* ---- T4: broken into contigs, orientations shuffled ---- */
[30000, 8000].forEach(function (ml) {
  const contigs = shard(alt, ml, "shard" + ml);
  const t0 = Date.now();
  const r = core.typeGenome(contigs, db, OPTS);
  const ms = Date.now() - t0;
  console.log(
    "\nT4 broken into contigs (mean " + ml + " bp, " + contigs.length + " pieces, orientations shuffled)" +
    "\n   call=" + r.verdict +
    "  | typed " + r.stats.called + "/" + r.stats.n_total +
    "  | score " + (r.stats.score * 100).toFixed(1) + "%" +
    "  | G1=" + r.stats.n_g1 + " NA=" + r.stats.n_na + "  | " + ms + " ms"
  );
  check("called Group 1", r.level === "pos", r.verdict);
  check("match score = 100%", r.stats.score === 1, (r.stats.score * 100).toFixed(1) + "%");
  check("at least 95% of loci readable", r.stats.called >= 81, r.stats.called + "/85");
});

/* ---- T5: only half of the loci converted ---- */
const half = buildAltGenome(0.5);
const t5 = run("T5 half of the loci set to the Group 1 allele (recombination / chimera)", half);
check("must not be called Group 1", t5.level !== "pos", t5.verdict);
check("must not be called not-Group 1", t5.level !== "neg", t5.verdict);
check("lands in the inconclusive band", t5.level === "warn", t5.verdict);

/* ---- T6: whole genome fed in as its reverse complement ---- */
console.log("\nT6 Reverse-complemented input");
const flipped = core.parseFasta(
  ">flip\n" + core.revcomp(refSeq)
).contigs;
const t6 = core.typeGenome(flipped, db, OPTS);
console.log("  whole reference reversed -> " + t6.verdict +
  "  | typed " + t6.stats.called + "/" + t6.stats.n_total +
  "  | score " + (t6.stats.score * 100).toFixed(1) + "%");
check("still called not Group 1", t6.level === "neg", t6.verdict);
check("all loci still readable", t6.stats.called === 85, t6.stats.called + "/85");

console.log("\n" + "=".repeat(74));
console.log(pass + " passed, " + fail + " failed");
console.log("=".repeat(74));
process.exit(fail ? 1 : 0);
