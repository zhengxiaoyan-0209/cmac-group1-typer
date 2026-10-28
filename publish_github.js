#!/usr/bin/env node
/* =======================================================================
   publish_github.js — push this repository to GitHub and switch on Pages.

   Why not plain "git push": on some networks github.com:443 is unreachable
   while api.github.com works. This script therefore talks to the REST API
   over HTTPS and never touches the git protocol.

       GITHUB_TOKEN=ghp_xxx node scripts/publish_github.js
       GITHUB_TOKEN=ghp_xxx REPO_NAME=cmac-group1-typer node scripts/publish_github.js

   The token needs the classic "repo" scope (repo creation + contents +
   Pages). It is read from the environment only, never written to disk.

   Safe to re-run: if the repository already exists it is reused, and the
   working tree is committed as a single new commit on top.
   ======================================================================= */

const fs = require("fs");
const path = require("path");
const { execFileSync } = require("child_process");

const ROOT = path.resolve(__dirname, "..");
const REPO = process.env.REPO_NAME || "cmac-group1-typer";
const TOKEN = process.env.GITHUB_TOKEN;
const DESCRIPTION =
  "Typing tool for Corynebacterium macginleyi Group 1 based on 85 " +
  "group-specific SNP loci. Runs entirely in the browser.";

if (!TOKEN) {
  console.error("GITHUB_TOKEN is not set.");
  console.error("  export GITHUB_TOKEN=ghp_xxx   (Linux / macOS / git bash)");
  console.error("  set GITHUB_TOKEN=ghp_xxx     (Windows cmd)");
  process.exit(1);
}

const API = "https://api.github.com";

async function api(method, endpoint, body, tries) {
  tries = tries || 3;
  const headers = {
    Authorization: "Bearer " + TOKEN,
    Accept: "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
    "User-Agent": "cmac-group1-typer-publisher",
  };
  if (body) headers["Content-Type"] = "application/json";

  for (let i = 0; i < tries; i++) {
    const res = await fetch(API + endpoint, {
      method: method,
      headers: headers,
      body: body ? JSON.stringify(body) : undefined,
    });
    const text = await res.text();
    let json = null;
    try { json = JSON.parse(text); } catch (e) { /* non-JSON body */ }
    if (res.ok) return json;
    if (res.status >= 500 && i < tries - 1) {
      await new Promise((r) => setTimeout(r, 1500));
      continue;
    }
    const msg = (json && json.message) || text.slice(0, 300);
    const err = new Error(method + " " + endpoint + " -> " + res.status + " " + msg);
    err.status = res.status;
    throw err;
  }
}

/* --------------------------------------------------- files to publish */
function listFiles() {
  const out = execFileSync("git", ["ls-files"], { cwd: ROOT, encoding: "utf8" });
  return out.split(/\r?\n/).filter(Boolean);
}

async function main() {
  const me = await api("GET", "/user");
  const owner = me.login;
  console.log("GitHub account : " + owner + (me.name ? " (" + me.name + ")" : ""));
  console.log("Repository     : " + owner + "/" + REPO);

  /* 1. create the repository if it does not exist yet */
  let repo;
  try {
    repo = await api("GET", "/repos/" + owner + "/" + REPO);
    console.log("Repository     : already there, reusing it");
  } catch (e) {
    if (e.status !== 404) throw e;
    repo = await api("POST", "/user/repos", {
      name: REPO,
      description: DESCRIPTION,
      homepage: "",
      private: false,
      auto_init: false,          /* keep it empty, we supply our own history */
      has_issues: false,
      has_wiki: false,
      has_projects: false,
    });
    console.log("Repository     : created");
  }

  /* 2. upload every tracked file as a blob */
  const files = listFiles();
  console.log("\nUploading " + files.length + " files:");
  const tree = [];
  for (const rel of files) {
    const abs = path.join(ROOT, rel);
    const buf = fs.readFileSync(abs);
    const blob = await api("POST", "/repos/" + owner + "/" + REPO + "/git/blobs", {
      content: buf.toString("base64"),
      encoding: "base64",
    });
    tree.push({
      path: rel.split(path.sep).join("/"),
      mode: "100644",
      type: "blob",
      sha: blob.sha,
    });
    console.log("   " + rel.padEnd(34) + (buf.length / 1024).toFixed(1) + " KB");
  }

  /* 3. one commit on top of whatever is already there */
  const newTree = await api("POST", "/repos/" + owner + "/" + REPO + "/git/trees", {
    tree: tree,
  });

  let parents = [];
  try {
    const ref = await api(
      "GET",
      "/repos/" + owner + "/" + REPO + "/git/ref/heads/main"
    );
    parents = [ref.object.sha];
  } catch (e) {
    if (e.status !== 404 && e.status !== 409) throw e;
    /* empty repository: no parent */
  }

  const commit = await api("POST", "/repos/" + owner + "/" + REPO + "/git/commits", {
    message:
      (parents.length ? "Update" : "Add") +
      " C. macginleyi Group 1 typing tool\n\n" +
      "85 Group 1 specific SNP loci typed by flank anchoring, in the browser,\n" +
      "no upload and no server.",
    tree: newTree.sha,
    parents: parents,
  });
  console.log("\nCommit         : " + commit.sha.slice(0, 7));

  if (parents.length) {
    await api("PATCH", "/repos/" + owner + "/" + REPO + "/git/refs/heads/main", {
      sha: commit.sha,
      force: false,
    });
  } else {
    await api("POST", "/repos/" + owner + "/" + REPO + "/git/refs", {
      ref: "refs/heads/main",
      sha: commit.sha,
    });
  }

  /* 4. make sure main is the default branch */
  await api("PATCH", "/repos/" + owner + "/" + REPO, {
    default_branch: "main",
    description: DESCRIPTION,
  });

  /* 5. switch on GitHub Pages */
  const url = "https://" + owner.toLowerCase() + ".github.io/" + REPO + "/";
  let pages = null;
  try {
    pages = await api("POST", "/repos/" + owner + "/" + REPO + "/pages", {
      source: { branch: "main", path: "/" },
    });
    console.log("Pages          : enabled");
  } catch (e) {
    if (e.status === 409) {
      pages = await api("PUT", "/repos/" + owner + "/" + REPO + "/pages", {
        source: { branch: "main", path: "/" },
      });
      console.log("Pages          : already enabled, source reset to main /");
    } else {
      throw e;
    }
  }

  await api("PATCH", "/repos/" + owner + "/" + REPO, { homepage: url });

  console.log("\n" + "-".repeat(66));
  console.log("URL            : " + url);
  console.log("-".repeat(66));
  console.log(
    "The first build takes roughly 1-3 minutes. Check progress with:\n" +
      "  https://github.com/" + owner + "/" + REPO + "/deployments\n" +
    "\nIf it still 404s after five minutes, open Settings -> Pages in the\n" +
    "repository and confirm Source is \"Deploy from a branch\" / main / (root)."
  );
}

main().catch(function (e) {
  console.error("\nFAILED: " + e.message);
  process.exit(1);
});
