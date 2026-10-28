/* =======================================================================
   core.js — C. macginleyi Group 1 typing engine
   Loadable both in the browser (window.CMACCore) and in Node (module.exports),
   so the web page and the offline checks run exactly the same code.

   How it works
   ------------
   The published workflow aligns each genome against CP069516.1, reads the alleles
   at 85 Group 1 specific positions and judges the "match ratio". This engine
   reaches the same answer without an alignment, by anchoring:

     1. take the 25 bp immediately beside the SNP as an anchor;
     2. locate that anchor in the query and in its reverse complement;
     3. once anchored, read the single adjacent base - the observed allele;
     4. if it is absent, fall back to a 20 bp and then a 16 bp anchor;
     5. if the anchor matches in several places and different bases come back,
        the locus is scored NA (a missed call beats a wrong call).

   Requiring one whole 51-mer probe to match exactly would lose every locus that
   carries any additional variant nearby. Anchoring tolerates other variation in
   the flank, so recall is much higher, while a 25 bp string is expected to occur
   far less than once in 2.5 Mb, so specificity still holds.

   Haplotype-aware flanks
   ----------------------
   The 85 loci cluster: the densest window holds six positions within 25 bp
   (around 2,390 kb). A flank taken from the reference alone would cover a
   neighbouring mutated base, so a genuine Group 1 strain fails to anchor there.
   Every flank is therefore also stored as a "Group 1 haplotype" version, and both
   are searched, which keeps recall without favouring either allele.
   ======================================================================= */
(function (factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else window.CMACCore = factory();
})(function () {
  "use strict";

  /* ------------------------------------------------------------ helpers */
  function revcomp(s) {
    const K = 65536, n = s.length, parts = [];
    for (let st = 0; st < n; st += K) {
      const e = Math.min(n, st + K), m = e - st, a = new Array(m);
      for (let i = st; i < e; i++) {
        const c = s.charCodeAt(i);
        a[e - 1 - i] = c === 65 ? 84 : c === 84 ? 65 : c === 67 ? 71 : c === 71 ? 67 : 78;
      }
      parts.push(String.fromCharCode.apply(null, a));
    }
    return parts.join("");
  }

  function sanitize(seq) {
    return seq.toUpperCase().replace(/[^ACGTN]/g, "N");
  }

  function findAll(hay, needle, limit) {
    const out = [];
    let i = hay.indexOf(needle);
    while (i !== -1) {
      out.push(i);
      if (out.length >= limit) break;
      i = hay.indexOf(needle, i + 1);
    }
    return out;
  }

  const COMP = { A: "T", T: "A", C: "G", G: "C", N: "N" };

  function parseFasta(text, maxContigs) {
    const contigs = [];
    let name = null, buf = [];
    const lines = text.split(/\r?\n/);
    for (let i = 0; i < lines.length; i++) {
      const ln = lines[i];
      if (ln.charCodeAt(0) === 62) {            // '>'
        if (name !== null) contigs.push({ name: name, seq: buf.join("") });
        if (maxContigs && contigs.length >= maxContigs) {
          return { contigs: contigs, truncated: true };
        }
        name = ln.slice(1).trim();
        buf = [];
      } else if (name !== null) {
        buf.push(ln);
      }
    }
    if (name !== null) contigs.push({ name: name, seq: buf.join("") });
    return { contigs: contigs, truncated: false };
  }

  /* ------------------------------------- assembly as one pseudo-genome */
  const SPACER = "N".repeat(200);   /* keeps an anchor from matching across contigs */

  function prepareQuery(contigs, db, opts) {
    opts = opts || {};
    const fwd = contigs.map(function (c) { return sanitize(c.seq); }).join(SPACER);
    return {
      fwd: fwd,
      rev: revcomp(fwd),
      fallback: opts.fallback || db.meta.anchor_fallback || [25, 20, 16],
      limit: opts.limit || 12
    };
  }

  /* -------------------------------------------------- typing one locus */
  /* Search order: anchor length (25 -> 20 -> 16) first, then left/right flank,
     then the two haplotype versions. Length has to be the outer loop - otherwise
     the reference version would be shortened all the way to 16 bp before the
     other haplotype is tried, and a 16 bp anchor easily hits elsewhere in the
     genome, returning an unrelated base. */
  function callLocus(marker, ctx) {
    const sides = ["L", "R"];
    for (let li = 0; li < ctx.fallback.length; li++) {
      const len = ctx.fallback[li];

      for (let oi = 0; oi < sides.length; oi++) {
        const side = sides[oi];
        const refver = marker[side];
        const altver = marker[side + "a"];
        if (!refver || len > refver.length) continue;

        const versions = [refver];
        if (altver && altver !== refver) versions.push(altver);

        const needles = [];
        for (let v = 0; v < versions.length; v++) {
          const src = versions[v];
          const n = side === "L" ? src.slice(src.length - len) : src.slice(0, len);
          if (n.length === len && n.indexOf("N") < 0 && needles.indexOf(n) < 0) {
            needles.push(n);
          }
        }

        const seen = Object.create(null);
        let hits = 0;
        for (let s = 0; s < 2; s++) {
          const X = s === 0 ? ctx.fwd : ctx.rev;
          for (let ni = 0; ni < needles.length; ni++) {
            const occ = findAll(X, needles[ni], ctx.limit);
            for (let k = 0; k < occ.length; k++) {
              const bIdx = side === "L" ? occ[k] + len : occ[k] - 1;
              if (bIdx < 0 || bIdx >= X.length) continue;
              seen[X[bIdx]] = true;
              hits++;
            }
          }
        }
        if (hits === 0) continue;              /* try the other side, or a shorter anchor */
        const alleles = Object.keys(seen);
        if (alleles.length > 1) {
          return { status: "ambiguous", side: side, len: len, hits: hits };
        }
        return { status: "ok", base: alleles[0], side: side, len: len, hits: hits };
      }
    }
    return { status: "missing" };
  }

  /* ---------------------------------------------------- row formatting */
  function classifyRow(m, call) {
    let verdict = "NA", note = "";
    if (call.status === "ok") {
      if (call.base === m.alt) verdict = "GROUP1";
      else if (call.base === m.ref) verdict = "NON";
      else { verdict = "OTHER"; note = "third allele " + call.base; }
    } else if (call.status === "ambiguous") {
      note = "anchor matched repeatedly (repeat region)";
    } else {
      note = "flank not anchored (variant / deletion / across contig boundary)";
    }
    return {
      pos: m.pos, ref: m.ref, alt: m.alt,
      gene: m.gene, locus_tag: m.locus_tag, product: m.product,
      feature: m.feature, region_kb: m.region_kb,
      obs: call.base || "-", verdict: verdict, note: note,
      side: call.side || "", anchor_len: call.len || 0, hits: call.hits || 0
    };
  }

  /* ---------------------------------------------------------- summary */
  function summarize(rows, db, opts) {
    opts = opts || {};
    const minCov = opts.minCov != null ? opts.minCov : 0.70;
    const hi = opts.hi != null ? opts.hi : 0.90;
    const lo = opts.lo != null ? opts.lo : 0.10;

    let nG1 = 0, nNon = 0, nOther = 0, nNA = 0;
    for (let i = 0; i < rows.length; i++) {
      const v = rows[i].verdict;
      if (v === "GROUP1") nG1++;
      else if (v === "NON") nNon++;
      else if (v === "OTHER") nOther++;
      else nNA++;
    }
    const called = nG1 + nNon + nOther;
    const total = rows.length || db.snps.length;
    const coverage = called / total;
    const score = called ? nG1 / called : 0;

    let verdict, level;
    if (called === 0) {
      verdict = "No locus could be read — check that this really is a C. macginleyi genome";
      level = "warn";
    } else if (coverage < minCov) {
      verdict = "Too few loci readable (" + called + "/" + total + ") for a call";
      level = "warn";
    } else if (score >= hi) {
      verdict = "Group 1";
      level = "pos";
    } else if (score <= lo) {
      verdict = "Not Group 1";
      level = "neg";
    } else {
      verdict = "Inconclusive (match score " + (score * 100).toFixed(1) + "%, between the two thresholds)";
      level = "warn";
    }

    return {
      rows: rows, db: db,
      stats: {
        n_total: total, n_g1: nG1, n_non: nNon,
        n_other: nOther, n_na: nNA, called: called,
        coverage: coverage, score: score
      },
      thresholds: { minCov: minCov, hi: hi, lo: lo },
      verdict: verdict, level: level
    };
  }

  /* ------------------------------------- one-shot entry point for Node */
  function typeGenome(contigs, db, opts, onProgress) {
    opts = opts || {};
    const ctx = prepareQuery(contigs, db, opts);
    const rows = [];
    for (let i = 0; i < db.snps.length; i++) {
      rows.push(classifyRow(db.snps[i], callLocus(db.snps[i], ctx)));
      if (onProgress && i % 10 === 0) onProgress(i + 1, db.snps.length);
    }
    return summarize(rows, db, opts);
  }

  return {
    revcomp: revcomp, sanitize: sanitize, parseFasta: parseFasta,
    prepareQuery: prepareQuery, callLocus: callLocus, classifyRow: classifyRow,
    summarize: summarize, typeGenome: typeGenome,
    COMP: COMP
  };
});
