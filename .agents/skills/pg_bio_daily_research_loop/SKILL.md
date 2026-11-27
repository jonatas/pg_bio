---
name: pg_bio_daily_research_loop
description: >-
  Autonomous daily research, worktree prototyping, benchmarking, PR creation, blog publishing,
  and video backlog pipeline for pg_bio. Scans daily literature (arXiv, bioRxiv, OpenAlex, PubMed)
  in computational protein design, structural bioinformatics, protein language models, and molecular
  dynamics alongside PostgreSQL engine/extension advancements, tracks attempts in local storage
  (SQL ledger in bio_demo + data/research_ledger.json), proves theories in isolated git worktrees,
  opens benchmarked PRs, publishes companion posts on ../ideia.me, and queues video episodes.
---

# `pg_bio` Autonomous Daily Research & Worktree Pipeline (`pg_bio_daily_research_loop`)

Use this skill for the recurring daily research loop (and whenever the user asks to scan new papers/PostgreSQL optimizations, track feature attempts, spin up worktrees, benchmark against the local `bio_demo` research database, open PRs, publish companion `ideia.me` blog posts, and enqueue video episodes).

---

## 1. Context Efficiency & Sub-Agent Delegation

To keep the main coordinator conversation lean across recurring runs:
1. **Always check Local Storage FIRST** (`pg_bio_research_attempts` & `pg_bio_video_backlog` in `postgresql://jonatas@localhost:28818/bio_demo` and [`data/research_ledger.json`](../../data/research_ledger.json)) so you never re-research or re-implement an already-attempted paper or PostgreSQL technique.
2. **Delegate Literature & PostgreSQL Discovery to Sub-Agents (`invoke_subagent`)**:
   * **Sub-Agent A (`Bioinformatics & PLM Literature Researcher`)**: Searches the past 24–72 hours of papers, preprints, and datasets across computational protein design, structural bioinformatics, protein language models (PLMs), enzymology, and molecular dynamics (using `literature-search-arxiv`, `literature-search-biorxiv`, `literature-search-openalex`, `pubmed-database`, or web search). Returns a concise structured digest of 2–3 high-leverage candidates with exact math/algorithms.
   * **Sub-Agent B (`PostgreSQL & Vector/Spatial Systems Researcher`)**: Researches recent PostgreSQL 18+ improvements, index access methods (GiST, SP-GiST, BRIN, HNSW, PRISM/RaBitQ, SIMD fastscan), varlena/TOAST optimizations, and extension techniques, benchmarking how they can boost `pg_bio` against our `2.34M`-protein `bio_demo` core research database (`localhost:28818`).

---

## 2. Local Storage & Progress Ledger (SQL + JSON Dual Ledger)

Every research scan, hypothesis test, worktree prototype, PR, blog post, and video backlog item MUST be persisted in both:
1. **Core Research DB (`postgresql://jonatas@localhost:28818/bio_demo`)**:
   * Table `public.pg_bio_research_attempts`:
     * `attempt_id SERIAL PRIMARY KEY`
     * `attempt_code VARCHAR(64) UNIQUE NOT NULL` (e.g., `'ATT-2026-11-26-ENZREF-PSLIM'`)
     * `run_date DATE NOT NULL DEFAULT CURRENT_DATE`
     * `category VARCHAR(40) NOT NULL` (`'PLM_Kinetics'`, `'Structural_Bio'`, `'Molecular_Dynamics'`, `'Postgres_Systems'`, `'Cheminformatics'`)
     * `paper_title TEXT NOT NULL`
     * `paper_ref TEXT NOT NULL` (arXiv ID, DOI, or repo URL)
     * `pg_strategy TEXT NOT NULL` (how it translates to `pg_bio` Rust/SIMD/index primitives)
     * `worktree_path TEXT` (e.g., `'.worktrees/enzref-pslim'`)
     * `branch_name VARCHAR(120)` (e.g., `'feat/enzref-role-contrast-pslim-tethers'`)
     * `status VARCHAR(32) NOT NULL` (`'Proposed'`, `'Theory_Testing'`, `'Worktree_Active'`, `'Benchmarked'`, `'PR_Opened'`, `'Merged'`, `'Published'`, `'Rejected'`)
     * `theory_test_name TEXT`
     * `benchmark_summary JSONB` (speedup vs Python/baseline, throughput, recall/accuracy, memory/TOAST footprint)
     * `pr_number INTEGER`
     * `pr_url TEXT`
     * `blog_post_slug TEXT`
     * `video_backlog_id INTEGER`
     * `notes TEXT`
     * `created_at TIMESTAMPTZ DEFAULT now()`, `updated_at TIMESTAMPTZ DEFAULT now()`
   * Table `public.pg_bio_video_backlog`:
     * `video_id SERIAL PRIMARY KEY`
     * `episode_slug VARCHAR(100) UNIQUE NOT NULL`
     * `title TEXT NOT NULL`
     * `layer_name VARCHAR(80) NOT NULL`
     * `attempt_code VARCHAR(64) REFERENCES pg_bio_research_attempts(attempt_code)`
     * `pr_number INTEGER`
     * `blog_post_slug TEXT`
     * `bio_problem TEXT NOT NULL`
     * `cs_mechanism TEXT NOT NULL`
     * `speedup_metric TEXT NOT NULL`
     * `status VARCHAR(32) NOT NULL DEFAULT 'Queued'` (`'Queued'`, `'Scaffolded'`, `'Rendered'`, `'Published_IdeiaMe'`)
     * `created_at TIMESTAMPTZ DEFAULT now()`, `updated_at TIMESTAMPTZ DEFAULT now()`
2. **Git-Tracked JSON Snapshot ([`data/research_ledger.json`](../../data/research_ledger.json))**:
   * Synced via `uv run python .agents/skills/pg_bio_daily_research_loop/scripts/manage_research_ledger.py sync` so worktrees and offline sessions always have the latest attempts, benchmarks, and video queue.

---

## 3. End-to-End 7-Stage Execution Lifecycle

### Stage 1: Check Ledger & Inspect Active Worktrees
Run:
```bash
uv run python /Users/jonatas/code/pg_bio/.agents/skills/pg_bio_daily_research_loop/scripts/manage_research_ledger.py status
```
Check `git -C /Users/jonatas/code/pg_bio worktree list` to see if an active worktree already has a feature in flight or ready for benchmarking/PR.

### Stage 2: Parallel Literature & PostgreSQL Discovery (Sub-Agents)
Launch two sub-agents in parallel to scan:
1. **Bioinformatics / Protein Design / PLM / MD papers & datasets from the past day** (or recent high-impact preprints not yet in `pg_bio_research_attempts`).
2. **PostgreSQL & Index/SIMD/Storage engine optimizations** tested against the local `bio_demo` database (`2,340,766` proteins, `scope40_esm2_650m_ref`, `bio_substrates`, `swarm_pocket_ranges`, `orphan_discoveries`).

### Stage 3: Score Candidates & Prove Theory with Small Tests
Rank candidate ideas by:
* **Scientific Impact**: Does it eliminate wet-lab or discovery false positives (e.g., dead binders, split hydrophobic cores, steric clashes, anisotropic embeddings)?
* **Systems Advantage ("Compute to the Data")**: Can it run in compiled Rust + hardware SIMD (`NEON` / `AVX2`) or native PostgreSQL index access methods (GiST, PRISM, 16-byte register Fastscan, inline `BioVecInt8G128` varlena) inside shared buffers?
* **Theory Proof**: Write a small deterministic unit test inside `#[cfg(any(test, feature = "pg_test"))]` and run `cargo pgrx test pg18 <test_name>` to prove the mathematical invariant (e.g., exact antisymmetry, unit L2 norm, zero TOAST size, or recall gain).

### Stage 4: Isolated Git Worktree Implementation (`pg_bio_feature_implementation`)
Never pollute `main` with experimental code. Create or use an isolated git worktree under `/Users/jonatas/code/pg_bio/.worktrees/<slug>`:
```bash
git -C /Users/jonatas/code/pg_bio worktree add .worktrees/<slug> -b feat/<branch-slug>
```
Follow the 5 phases of [`pg_bio_feature_implementation`](../pg_bio_feature_implementation/SKILL.md):
1. **Rust Implementation (`src/lib.rs`)**: Single-threaded, panic-safe, SIMD-accelerated `#[pg_extern(immutable, parallel_safe)]` functions.
2. **Unit & Regression Testing**: `cargo pgrx test pg18`.
3. **Performance Benchmarking (`scripts/benchmark_<slug>.py`)**: Benchmark both in-engine (`pgrx`) and against real data in `bio_demo` (`localhost:28818`), comparing latency/throughput/accuracy against Python baselines.
4. **Python SDK & Docs (`pgbio-py/` & `docs/<NN>_<slug>.md` & `README.md`)**: Update `pgbio-py/pgbio/client.py`, add a numbered deep-dive doc in `docs/`, and update `README.md`.
5. **Executed Jupyter Notebook (`experiments/<NN>_<slug>.ipynb`)**: Execute via `uv run jupyter nbconvert --to notebook --execute --inplace experiments/<NN>_<slug>.ipynb`.

### Stage 5: Self-Documenting Pull Request (`gh pr create`)
When the theory test passes and benchmarks show clear scientific and systems gains:
1. Stage only specific, verified files (`git add <file>` — NEVER `git add .` or `git add -A`).
2. Push the branch and create a comprehensive, self-documenting Pull Request via `gh pr create --repo jonatas/pg_bio` containing:
   * Mathematical & biological problem statement with paper citations.
   * SQL signatures and architecture details.
   * Full Markdown benchmark tables (latency, throughput, recall/accuracy, TOAST/disk size).
   * Verification checklist (`cargo pgrx test pg18`, `bio_demo` query validation, executed `.ipynb`).
3. Update the attempt status in `pg_bio_research_attempts` to `'PR_Opened'` (or `'Merged'` once merged).

### Stage 6: Companion Blog Post on `../ideia.me` (`pg_bio_blog_writer` & `postgres-rust-series-writer`)
When the feature is merged (or ready for publication alongside the PR):
1. Follow [`pg_bio_blog_writer`](../pg_bio_blog_writer/SKILL.md) to craft an original, technically rigorous post in `/Users/jonatas/code/ideia.me/_posts/YYYY-MM-DD-<slug>.md` tying the new bioinformatics research + PostgreSQL Rust/SIMD internals to real discoveries in `bio_demo`.
2. Download local `.pdb` structures to `/Users/jonatas/code/ideia.me/assets/models/` for interactive side-by-side 3Dmol.js viewers with chemistry-teaching buttons.
3. Register the post in `/Users/jonatas/code/bio-pipelines/base/discovery_routes.py` (`IDEIAME_BLOG_POSTS`) and cross-link it in `/Users/jonatas/code/ideia.me/_includes/pg_bio_promo.md`.

### Stage 7: Enqueue & Scaffold Video Episode (`pg_bio_tutor_and_video`)
1. Insert a backlog task into `public.pg_bio_video_backlog` (and `data/research_ledger.json`) via `manage_research_ledger.py enqueue-video`.
2. Scaffold or render the episode bundle in `/Users/jonatas/code/pg_bio/videos/<episode_slug>/` using `.agents/skills/pg_bio_tutor_and_video/scripts/stage_feature_video.py` and sync it to `/Users/jonatas/code/ideia.me/pg-bio-videos.html` when requested or as part of post-merge automation.
