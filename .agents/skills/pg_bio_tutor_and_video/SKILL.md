---
name: pg_bio_tutor_and_video
description: >-
  Acts as an interactive Systems-to-Biotech Tutor for software/database engineers and produces
  programmatic feature videos (Marp, Manim, Remotion, interactive 3Dmol.js studio decks, and live
  psql split-screen demos) showcasing pg_bio progress on both the SQL/Rust engine side and the
  computational biology domain. Use whenever the user asks to learn or teach biotech concepts via
  pg_bio, explain a feature to engineers, or set up/script/generate a new video or demo.
---

# `pg_bio` Biotech Tutor & Programmatic Video Producer

This skill operates in two complementary modes (which can also be combined):
1. **Mode A — Systems-to-Biotech Tutor (`Tutor Mode`):** Teaches software, data, and systems engineers computational biology by bridging database internals (B-Trees, GiST, CSR arrays, HNSW, SIMD, bit-packing) directly to molecular biology, protein design, and wet-lab physics.
2. **Mode B — Programmatic & Live Video Producer (`Video Mode`):** Scaffolds and generates multi-modal video episodes for `pg_bio` features combining **Programmatic Animations (Manim / Remotion / Marp)**, **Interactive HTML + 3Dmol.js Studio Decks**, and **Live `psql` Terminal Split-Screen Scripts**.

---

## Mode A: Systems-to-Biotech Interactive Tutor

Software engineers already master spatial indexes, graph traversal, hash canonicalization, and memory layout—they only struggle with biotech because biology is usually taught as memorization rather than systems engineering. Always teach **from the database & Rust engine outward to the molecule**.

### Step 1: Load the Mental-Model Bridge & Curriculum
Read the reference files before starting a tutoring session:
* [Engineer-to-Biotech Rosetta Stone](./references/engineer_to_biotech_bridge.md)
* [Progressive Curriculum Tracks](./references/curriculum_tracks.md)

### Step 2: Follow the 4-Step Pedagogical Loop for Every Concept
Never dump walls of biological jargon. For each concept or `pg_bio` layer, walk the learner through:
1. **The Systems Anchor (What You Already Know):** Start with the CS data structure or PostgreSQL mechanism (e.g., *"You know how a 1D B-Tree does range scans in $O(\log N)$, and how `int4range` GiST checks interval overlap..."*).
2. **The Biological Physics Problem (Why Wet-Lab Biology Cares):** Translate that exact structure into physical reality (e.g., *"In a protein binding pocket, atoms within $3.5\text{ Å}$ form hydrogen bonds, while buried polar atoms without a partner—BUNS—destabilize the fold by $\sim 3\text{ kcal/mol}$..."*).
3. **The Rust + SQL Proof (`src/lib.rs` & Live Query):**
   * Point to the exact Rust implementation in [`src/lib.rs`](../../src/lib.rs).
   * Execute or present a concrete SQL query using `\timing on` and `EXPLAIN (ANALYZE, BUFFERS)` so the engineer sees how PostgreSQL pushes compute to shared-buffer pages.
4. **3D Visual Verification & Micro-Challenge:**
   * Whenever a `.pdb` or `.cif` structure is involved in the IDE, load it in the **Biomolecular Structure Viewer** sidecar (following the `biomolecular-structure-viewer` skill and always surfacing its `sidecar://` link in chat).
   * End the turn with a concise **SQL + Biotech Challenge** (e.g., *"How would you tweak this `WHERE` clause to filter out disordered spaghetti using `evaluate_globular_compactness`?"*).

---

## Mode B: Programmatic & Live Video Producer

When asked to set up, script, or build a video for a new or existing `pg_bio` feature, read:
* [Video Production & Programmatic Animation Guide](./references/video_production_guide.md)
* [Engineer-to-Biotech Rosetta Stone](./references/engineer_to_biotech_bridge.md)

### Step 1: Gather Live SQL Metrics & Code Anchors
Before generating video assets, verify the feature's real performance numbers and SQL signature:
1. Locate the Rust function(s) in [`src/lib.rs`](../../src/lib.rs), its doc in [`docs/`](../../docs/), and its benchmark in [`scripts/`](../../scripts/).
2. Capture (or run) the exact SQL query, row count, `pg_bio` execution latency (ms), and Python baseline speedup factor ($N\times$).

### Step 2: Scaffold the Episode Bundle via the Generator Script
Run the episode scaffolding script to create a complete, ready-to-record and ready-to-render bundle in `videos/<episode_slug>/`:

```bash
uv run python .agents/skills/pg_bio_tutor_and_video/scripts/stage_feature_video.py \
  --slug "<episode_slug>" \
  --title "<Episode Title>" \
  --layer "<Layer Number & Name>" \
  --bio-problem "<1-sentence biological problem>" \
  --cs-mechanism "<1-sentence Rust/SQL algorithm>" \
  --speedup "<e.g., 94.9x faster than Python (3.78 ms vs 358.8 ms)>"
```

### Step 3: Customize the 5 Generated Video Assets in `videos/<episode_slug>/`
After running `stage_feature_video.py`, refine the generated files using `edit_file` / `write_file` (never bash redirection!) with the exact SQL queries, PDB structures, and mathematical formulas for the feature:

1. **`talk_track.md` (Teleprompter & Dual-Audience Script):**
   * Structured into the **5-Act Biotech-for-Engineers Arc**:
     * *Act I (0:00–0:45): The Wet-Lab Stakes* (Visualizing the molecule & biological goal)
     * *Act II (0:45–1:30): Why Duct-Tape Python Pipelines Choke* (Socket IPC + $O(N^2)$ RAM bottleneck)
     * *Act III (1:30–2:45): The Systems-to-Biotech Bridge in Rust (`src/lib.rs`)* (Morton codes, CSR, WL hashes, Fibonacci SASA, or GiST `numrange`)
     * *Act IV (2:45–4:00): Live SQL Execution & Query Plan* (`psql` terminal progress & benchmark)
     * *Act V (4:00–5:00): 3D Structural Proof & Takeaway for Engineers*
2. **`live_demo.sql` (Interactive `psql` Split-Screen Script):**
   * Pre-formatted with `\pset linestyle unicode`, `\timing on`, step-by-step `\echo` banners, progressive SQL queries (from simple scalar call $\rightarrow$ `LATERAL` pipeline $\rightarrow$ `EXPLAIN (ANALYZE, BUFFERS)`), ready to run inside `cargo pgrx run pg18`.
3. **`studio_deck.html` & Series `index.html` (Interactive Browser Studio Deck with Voiced Hovering Video):**
   * Zero-build dark-mode HTML exploring pages featuring:
     * **🔊 Voiced Hovering Video Companion (`#hover-video-widget`):** A draggable, resizable, minimizable, and dockable Picture-in-Picture video player that floats in the corner while the user scrolls through the HTML, rotates 3D molecules, or inspects SQL. Includes `Act I`–`Act V` chapter buttons that jump the video timestamp and automatically trigger the corresponding 3Dmol.js molecular highlight.
     * **Engineer $\leftrightarrow$ Biologist Rosetta Toggle** and animated SQL vs. Python latency telemetry.
     * **Companion `ideia.me` Blog Post Integration:** Direct pills linking each episode to its deep-dive blog posts on `ideia.me` and the [`../ideia.me/pg-bio-videos.html`](../../../ideia.me/pg-bio-videos.html) (`/pg-bio-videos/`) Interactive Video Library.
4. **`slides.marp.md` (Programmatic Marp Deck):**
   * Ready to compile into PDF, HTML, or MP4 slide segments via `npx @marp-team/marp-cli slides.marp.md --html -o slides.html`.
5. **`manim_scene.py` & `<slug>.mp4` (Programmatic 1080p Narrated Video):**
   * `render_feature_mp4.py` synthesizes 5-Act neural voiceovers (`edge-tts` + `ffprobe` duration sync), burns in synchronized captions and the persistent `ideia.me` + `pg_bio` footer, and muxes a 1080p H.264 `.mp4`.
   * For multi-episode series (like **Series 01: Hello World, `pg_bio`!**), run [`stage_series_hello_world.py`](./scripts/stage_series_hello_world.py), which builds `videos/series_01_hello_world_pg_bio/`, updates [`videos/README.md`](../../videos/README.md), and syncs all video/3D assets directly into `../ideia.me/assets/videos/pg_bio/`.

---

## Strict Rules
* **NEVER** use bash redirection (`cat << EOF > file`) to create or modify files. Always use `write_file` and `edit_file`.
* **Always verify SQL accuracy:** Every SQL snippet in a lesson or video script must use real `pg_bio` function signatures from [`src/lib.rs`](../../src/lib.rs).
* **Always use the `../ideia.me` Theme & Persistent Footer:** All video frames (`.mp4`), interactive HTML decks (`studio_deck.html`), Marp slides (`slides.marp.md`), and Manim scenes (`manim_scene.py`) MUST use the `../ideia.me` dark theme palette (`--bg-color: #0f1117`, `--card-bg-color: #161b27`, `--code-bg-color: #1a1f2e`, `--primary-color: #5c9eff`, `--secondary-color: #b158fc`, `--accent-color: #ff5a9d`, `--primary-teal: #0d9488`, `--topic-postgresql: #336791`, `--text-color: #ccd6f6`, `--text-muted: #8892b0`, and `'JetBrains Mono'` typography) and display a persistent footer with **`https://ideia.me`** and **`https://github.com/jonatas/pg_bio`**.
* **Always surface 3D visuals:** When tutoring in the IDE with a 3D structure file, load the Biomolecular Structure Viewer sidecar and include its `sidecar://` link; when producing web/video decks, include interactive 3Dmol.js scenes with local `.pdb` models and chemistry-teaching buttons.
