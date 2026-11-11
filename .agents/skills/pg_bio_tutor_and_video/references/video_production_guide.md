# Video Production & Programmatic Animation Guide (`pg_bio`)

Every `pg_bio` video is designed to accomplish two things simultaneously:
1. **Show tangible SQL & Rust engine progress** (latency drop, `EXPLAIN ANALYZE` index usage, zero-IPC shared-buffer architecture).
2. **Expand a software engineer's domain into Biotech** (showing *why* a mathematical invariant in Rust solves a real molecular failure mode in 3D).

---

## 1. The 5-Act Episode Blueprint (4 to 6 Minutes)

| Act | Time | Screen Layout | Narrative Goal |
| :--- | :--- | :--- | :--- |
| **Act I: The Molecular Hook** | `0:00–0:45` | **Interactive 3D Viewer (`studio_deck.html` / Mol* Sidecar)** | Show the 3D protein/ligand immediately. Point out the physical chemistry problem (e.g., *"See this polar oxygen buried inside the hydrophobic core with no H-bond partner? That single BUNS defect unfolds the whole enzyme."*). |
| **Act II: Why Python Pipelines Choke** | `0:45–1:30` | **Manim Animation (`manim_scene.py`) or Architecture Slide** | Visualize the traditional bottleneck: pulling 100k PDB strings or vectors over a socket into Python, $O(N^2)$ loops, and GIL contention vs. keeping compute inside Postgres 8KB shared-buffer pages. |
| **Act III: The Rust & CS Insight** | `1:30–2:45` | **Split Screen: `src/lib.rs` + Manim Math Scene** | Walk through the exact algorithm in `src/lib.rs` using the **Systems-to-Biotech Bridge** (e.g., how `compute_atomic_sasa` bins atoms into $5\text{ Å}$ 3D hash cells + 64-point Fibonacci spheres, or how `z_order_encode` interleaves bits). |
| **Act IV: Live SQL Progress & Benchmark** | `2:45–4:00` | **Live Terminal (`psql` running `live_demo.sql`)** | Execute the progressive SQL queries live with `\timing on`. Show the millisecond latency, the `EXPLAIN (ANALYZE, BUFFERS)` plan, and the exact speedup factor ($16\times$ to $1,513\times$ over Python). |
| **Act V: 3D Biological Proof & Engineer Takeaway** | `4:00–5:00` | **`studio_deck.html` Dual 3Dmol.js Viewer** | Click the interactive teaching buttons in `studio_deck.html` to highlight the conserved triad, cryptic pocket, or verified ligand pose, and summarize the engineering lesson. |

---

## 2. The 3 Programmatic Video Modalities

The scaffolding script [stage_feature_video.py](../scripts/stage_feature_video.py) generates all three modalities inside `videos/<episode_slug>/` so you can record live, render programmatic animations, or combine both in an editor:

### Modality 1: Interactive Browser Studio Deck (`studio_deck.html`)
* **Zero build step:** Open directly in Chrome or serve via `python3 -m http.server 8080`.
* **Built-in features:**
  * Slide-by-slide keyboard navigation (`ArrowRight` / `ArrowLeft` / `Space`).
  * **Engineer $\leftrightarrow$ Biologist Rosetta Toggle:** Flips terminology in real time between CS systems language and structural biology language.
  * **Animated SQL Benchmark Telemetry:** Visualizes `pg_bio` vs. Python wall-clock latency bars and throughput multipliers.
  * **Synchronized 3Dmol.js Teaching Stage:** Loads local `.pdb` structures and provides interactive buttons (`Highlight Catalytic Pocket`, `Show Hydrophobic Core / BUNS`, `Toggle Surface Electrostatics`) that dynamically update a chemistry-teaching callout box below the viewer.
  * Note on 3Dmol.js API: `$3Dmol.viewers` is an Object map, so scripts always iterate via `Object.values(window.$3Dmol.viewers || {})`.

### Modality 2: Programmatic Mathematical Animations (`manim_scene.py`)
* Uses **Manim Community** to render crisp 60fps 1080p/4K MP4 animations of the underlying CS/math transformation:
  * **Scene 1 (`BottleneckVsSharedBuffers`):** Contrasts TCP socket serialization of rows to Python against zero-copy `pgrx` pointers inside PostgreSQL shared buffers.
  * **Scene 2 (`AlgorithmicBridgeScene`):** Animates the core multiomics reduction (3D $(x,y,z) \rightarrow$ 64-bit Morton curve, HNSW $\rightarrow$ Smith-Waterman funnel, or 8-column `numrange` GiST box intersection).
* **Render Command:**
  ```bash
  uv run --with manim manim -pqh videos/<episode_slug>/manim_scene.py FeatureScene
  ```

### Modality 3: Live `psql` Split-Screen Runner (`live_demo.sql` + `slides.marp.md`)
* **`live_demo.sql`** is formatted with clean Unicode box-drawing borders, `\timing on`, and three escalating SQL stages:
  1. **Stage 1 (The Atomic Primitive):** A single-row call showing the input/output types so engineers immediately grasp the SQL signature.
  2. **Stage 2 (The Multi-Constraint Discovery Pipeline):** A realistic `WITH ... CROSS JOIN LATERAL` screening query combining index filtering and physics verification.
  3. **Stage 3 (The Engine Proof):** `EXPLAIN (ANALYZE, BUFFERS)` proving index scan efficiency and sub-10ms execution.
* **`slides.marp.md`** allows instant export to HTML/PDF/MP4 slide decks via Marp CLI:
  ```bash
  npx @marp-team/marp-cli videos/<episode_slug>/slides.marp.md --html -o videos/<episode_slug>/marp_slides.html
  ```

---

## 3. Mandatory `../ideia.me` Theme Palette & Persistent Footer

Every generated video asset (`.mp4`, `studio_deck.html`, `slides.marp.md`, `manim_scene.py`) uses the exact color variables from [`../ideia.me/css/main.css`](../../../../ideia.me/css/main.css), [`../ideia.me/css/custom.css`](../../../../ideia.me/css/custom.css), and [`../ideia.me/assets/css/topic-colors.css`](../../../../ideia.me/assets/css/topic-colors.css):

| Token | Hex | Usage |
| :--- | :--- | :--- |
| `--bg-color` | `#0f1117` | Main canvas & 3Dmol.js viewer background |
| `--card-bg-color` | `#161b27` | Content cards & panels |
| `--code-bg-color` | `#1a1f2e` | SQL code blocks, `psql` terminal, and caption bar |
| `--primary-color` (dark) | `#5c9eff` | Primary blue headings, highlights, and Morton spheres |
| `--secondary-color` (dark) | `#b158fc` | Secondary purple accents & SQL keywords |
| `--accent-color` (dark) | `#ff5a9d` | Accent pink (Python baseline bar & nucleophilic water) |
| `--primary-teal` | `#0d9488` | `pg_bio` Rust speedup bar & catalytic Glu |
| `--topic-postgresql` | `#336791` | Header layer pill badge |
| `--text-color` / `--text-muted` | `#ccd6f6` / `#8892b0` | Main body copy & muted labels (`JetBrains Mono`) |

### Persistent Footer Requirement
Every frame of the `.mp4` video and the bottom of `studio_deck.html` / `slides.marp.md` must render a persistent footer bar containing:
* **Left:** `🌐 Blog & 3D Discoveries: https://ideia.me`
* **Center:** `🧬 pg_bio — Multiomics in PostgreSQL`
* **Right:** `🐙 Repo: https://github.com/jonatas/pg_bio`
