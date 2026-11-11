#!/usr/bin/env python3
"""
Scaffolds a complete pg_bio feature video & interactive tutoring bundle inside `videos/<slug>/`.
Styled with the exact `../ideia.me` dark theme palette and persistent footer linking to:
  - https://ideia.me
  - https://github.com/jonatas/pg_bio

Generates:
  1. talk_track.md     - 5-Act teleprompter script bridging Systems Engineering <-> Biotech
  2. live_demo.sql     - Progressive psql terminal demo with \timing on and EXPLAIN ANALYZE
  3. slides.marp.md    - Programmatic Marp slide deck (exportable to HTML/PDF/MP4)
  4. manim_scene.py    - Programmatic 60fps Manim animation of the Rust/SQL algorithm
  5. studio_deck.html  - Interactive browser studio deck with Engineer<->Biologist Rosetta
                         toggle, neural voiceover bar, SQL benchmark telemetry, 3Dmol.js stage,
                         and persistent ideia.me + pg_bio footer
  6. <slug>.mp4        - 1080p narrated & captioned MP4 video
"""

import argparse
import subprocess
from datetime import date
from pathlib import Path

BLOG_URL = "https://ideia.me"
REPO_URL = "https://github.com/jonatas/pg_bio"


def build_talk_track(args: argparse.Namespace) -> str:
    return f"""# 🎬 Video Talk Track: {args.title}
**Layer:** {args.layer}  
**Date:** {date.today().isoformat()}  
**Benchmark Headline:** {args.speedup}  
**Blog:** [{BLOG_URL}]({BLOG_URL}) | **Repository:** [{REPO_URL}]({REPO_URL})

---

## Act I (`0:00–0:45`) — The Molecular Hook (Why Biotech Cares)
* **Visual:** Open `studio_deck.html` (Slide 1 & 3Dmol.js Stage) or the Mol* Sidecar Viewer.
* **Talking Point (Biotech Context):**
  > "{args.bio_problem}"
* **What to point out on screen:**
  * Rotate the 3D structure so the viewer sees the physical binding pocket / catalytic motif.
  * Explain why a naive sequence match fails in the wet lab and why atomic-level or multi-constraint physics is mandatory.

---

## Act II (`0:45–1:30`) — Why Duct-Tape Python Pipelines Choke
* **Visual:** Play `manim_scene.py` (`FeatureScene` Part 1) or Act II of `{args.slug}.mp4`.
* **Talking Point (Systems Bottleneck):**
  > "In a traditional bioinformatics stack, you pull hundreds of thousands of rows over a TCP socket into Python, deserialize strings into objects, and run $O(N^2)$ loops in user space. 99% of your wall-clock time is spent moving bytes across the wire instead of doing science."

---

## Act III (`1:30–2:45`) — The Systems-to-Biotech Bridge in Rust (`src/lib.rs`)
* **Visual:** Split screen between [`src/lib.rs`](../../src/lib.rs) and `studio_deck.html` (toggle the **Engineer $\\leftrightarrow$ Biologist Rosetta Switch**).
* **Talking Point (The CS Insight):**
  > "{args.cs_mechanism}"
* **Key Translation for Engineers:**
  * **Systems View:** `{args.cs_mechanism}`
  * **Biological View:** `{args.bio_problem}`

---

## Act IV (`2:45–4:00`) — Live SQL Progress & Benchmark (`live_demo.sql`)
* **Visual:** Terminal running `psql` connected to `pg_bio` (`\\i videos/{args.slug}/live_demo.sql`).
* **Talking Point (The SQL Proof):**
  1. Run **Stage 1 (The Atomic Primitive)** to show the clean SQL signature.
  2. Run **Stage 2 (The Multi-Constraint Discovery Pipeline)** using `CROSS JOIN LATERAL` so index filtering prunes the search space before heavy physics runs.
  3. Run **Stage 3 (`EXPLAIN (ANALYZE, BUFFERS)`)** and highlight the execution time: **{args.speedup}**.

---

## Act V (`4:00–5:00`) — 3D Structural Proof & Call to Action
* **Visual:** Return to `studio_deck.html` (Interactive 3Dmol.js Teaching Stage + Footer).
* **Action:** Click the **"1. Highlight Catalytic Triad / Metal"** and **"2. Show 3D Morton Bounding Pocket"** buttons and read the dynamic chemistry callout.
* **Closing Hook:**
  > "Explore all interactive 3D discoveries at {BLOG_URL} and star or contribute to the engine at {REPO_URL}."
"""


def build_live_demo_sql(args: argparse.Namespace) -> str:
    return f"""-- ============================================================================
-- pg_bio Live Video Demo: {args.title}
-- Layer: {args.layer}
-- Target Benchmark: {args.speedup}
-- Blog: {BLOG_URL}  |  Repo: {REPO_URL}
-- Run inside psql: \\i videos/{args.slug}/live_demo.sql
-- ============================================================================

\\pset linestyle unicode
\\pset border 2
\\timing on

\\echo ''
\\echo '======================================================================'
\\echo ' STAGE 1: THE ATOMIC PRIMITIVE (Testing the Rust pgrx Function in SQL)'
\\echo '======================================================================'

SELECT
    molecule_wl_hash('C1=CC=CC=C1') = molecule_wl_hash('c1ccccc1') AS kekule_eq_aromatic,
    molecules_are_identical('[Na+].O=C([O-])C1=CC=CC=C1', 'O=C(O)c1ccccc1') AS salt_stripped_eq,
    round((parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O') % parse_smiles('O=C(O)c1ccccc1'))::numeric, 4) AS tanimoto_sim;

\\echo ''
\\echo '======================================================================'
\\echo ' STAGE 2: MULTI-CONSTRAINT BIOTECH PIPELINE (Index Filter -> Physics)'
\\echo '======================================================================'

WITH candidate_enzymes(design_id, sequence, active_site_coords) AS (
    VALUES
        ('Design_Alpha_ZnO', 'MFEGFERRLVDWPSTDKTKREEVDHCH', '[{{"x":0.0,"y":0.0,"z":0.0,"name":"ZN"}},{{"x":2.05,"y":0.1,"z":0.0,"name":"HOH"}}]'::jsonb)
)
SELECT
    design_id,
    round(isoelectric_point(sequence)::numeric, 2) AS pI,
    round(gravy_score(sequence)::numeric, 3) AS hydropathy,
    round(sequence_entropy(sequence)::numeric, 3) AS shannon_entropy,
    prosite_match(sequence, 'H-x-H') AS has_catalytic_motif
FROM candidate_enzymes;

\\echo ''
\\echo '======================================================================'
\\echo ' STAGE 3: ENGINE EXECUTION PLAN & SHARED-BUFFER TELEMETRY'
\\echo '======================================================================'

EXPLAIN (ANALYZE, BUFFERS)
SELECT
    z_order_encode(12.4, -4.8, 29.1) AS morton_3d_key,
    bio_translate_dna('ATGTTCGAAGGTTTCGAACGTCGTCTGGTTGACTGG') AS translated_peptide;

\\echo ''
\\echo '🌐 Blog: {BLOG_URL}   |   🐙 GitHub: {REPO_URL}'
"""


def build_marp_slides(args: argparse.Namespace) -> str:
    return f"""---
marp: true
theme: gaia
class: invert
paginate: true
backgroundColor: #0f1117
color: #ccd6f6
footer: "🌐 {BLOG_URL}  |  🧬 pg_bio Multiomics Engine  |  🐙 {REPO_URL}"
style: |
  section {{
    font-family: 'JetBrains Mono', 'Fira Code', 'SF Mono', Menlo, Consolas, monospace;
    padding: 48px;
    background-color: #0f1117;
    color: #ccd6f6;
  }}
  h1, h2 {{
    color: #5c9eff;
  }}
  code {{
    background: #1a1f2e;
    color: #ccd6f6;
  }}
  .bridge-box {{
    border: 1px solid #5c9eff;
    border-radius: 12px;
    padding: 20px;
    background: #161b27;
  }}
  footer {{
    color: #8892b0;
    font-size: 0.75rem;
  }}
---

# 🧬 {args.title}
### {args.layer}
**Bridging PostgreSQL Systems Engineering & Computational Biology**

* **Benchmark:** `{args.speedup}`
* **Engine:** Native Rust (`pgrx`) inside PostgreSQL 18

---

## 1. The Wet-Lab Biology Problem
<div class="bridge-box">

**Why Biologists Care:**  
{args.bio_problem}

</div>

* Traditional pipelines export gigabytes of sequences, coordinates, or trajectories over TCP sockets to Python.
* Network serialization + $O(N^2)$ interpreter loops turn millisecond questions into multi-hour batch jobs.

---

## 2. The Systems-to-Biotech Bridge (`src/lib.rs`)

| Systems Engineer Mental Model | Structural Biology Reality |
| :--- | :--- |
| `{args.cs_mechanism}` | `{args.bio_problem}` |
| Zero-IPC Shared-Buffer Pointers (`pgrx`) | Scan millions of molecules/residues in L1/L2 cache |
| `CROSS JOIN LATERAL` Index Short-Circuiting | Prune 99.9% of decoys before running 3D physics |

---

## 3. Declarative SQL Discovery Pipeline

```sql
-- Push multi-constraint biological physics directly to the data:
SELECT design_id,
       round(isoelectric_point(sequence)::numeric, 2) AS pI,
       round(gravy_score(sequence)::numeric, 3) AS gravy
FROM candidate_proteins
WHERE prosite_match(sequence, '[LIVM]-x(2)-H-x-H') = TRUE;
```

**Result:** `{args.speedup}`
"""


def build_manim_scene(args: argparse.Namespace) -> str:
    return f'''"""
Programmatic Manim Animation for: {args.title} ({args.layer})
Styled with the exact `ideia.me` theme palette and persistent footer.
Render with:
    uv run --with manim manim -pqh videos/{args.slug}/manim_scene.py FeatureScene
"""

from manim import (
    Arrow,
    Create,
    FadeIn,
    FadeOut,
    RoundedRectangle,
    Scene,
    Text,
    VGroup,
    Write,
    DOWN,
    LEFT,
    RIGHT,
    UP,
)

# Exact ideia.me Theme Hex Colors
IDEIA_BG = "#0f1117"
IDEIA_CARD = "#161b27"
IDEIA_BLUE = "#5c9eff"
IDEIA_PURPLE = "#b158fc"
IDEIA_PINK = "#ff5a9d"
IDEIA_TEAL = "#0d9488"
IDEIA_TEXT = "#ccd6f6"
IDEIA_MUTED = "#8892b0"


class FeatureScene(Scene):
    def construct(self):
        self.camera.background_color = IDEIA_BG

        title = Text("{args.title}", font_size=32, color=IDEIA_BLUE).to_edge(UP)
        subtitle = Text("{args.layer}  |  {args.speedup}", font_size=19, color=IDEIA_PURPLE).next_to(title, DOWN, buff=0.18)
        footer_bar = Text(
            "{BLOG_URL}   |   pg_bio Multiomics Engine   |   {REPO_URL}",
            font_size=16,
            color=IDEIA_MUTED,
        ).to_edge(DOWN, buff=0.22)

        self.play(Write(title), FadeIn(subtitle), FadeIn(footer_bar))

        py_box = RoundedRectangle(corner_radius=0.2, width=5.5, height=3.2, color=IDEIA_PINK, fill_opacity=0.15)
        py_title = Text("Traditional Python Pipeline", font_size=21, color=IDEIA_PINK).move_to(py_box.get_top() + DOWN * 0.4)
        py_steps = VGroup(
            Text("1. TCP Socket Pull (GBs of rows)", font_size=16, color=IDEIA_TEXT),
            Text("2. PyObject & NumPy Allocation", font_size=16, color=IDEIA_TEXT),
            Text("3. O(N^2) User-Space Loops", font_size=16, color=IDEIA_TEXT),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25).move_to(py_box.get_center() + DOWN * 0.2)
        py_group = VGroup(py_box, py_title, py_steps).to_edge(LEFT, buff=0.6).shift(DOWN * 0.2)

        pg_box = RoundedRectangle(corner_radius=0.2, width=5.5, height=3.2, color=IDEIA_TEAL, fill_opacity=0.15)
        pg_title = Text("pg_bio (Rust inside Postgres)", font_size=21, color=IDEIA_BLUE).move_to(pg_box.get_top() + DOWN * 0.4)
        pg_steps = VGroup(
            Text("1. Zero-IPC Shared Buffer Pointers", font_size=16, color=IDEIA_TEXT),
            Text("2. Index Sieve (B-Tree / GiST / HNSW)", font_size=16, color=IDEIA_TEXT),
            Text("3. LLVM SIMD Physics on Survivors", font_size=16, color=IDEIA_TEXT),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.25).move_to(pg_box.get_center() + DOWN * 0.2)
        pg_group = VGroup(pg_box, pg_title, pg_steps).to_edge(RIGHT, buff=0.6).shift(DOWN * 0.2)

        bridge_arrow = Arrow(py_group.get_right(), pg_group.get_left(), buff=0.2, color=IDEIA_PURPLE)

        self.play(Create(py_box), Write(py_title), FadeIn(py_steps))
        self.wait(0.8)
        self.play(Create(bridge_arrow), Create(pg_box), Write(pg_title), FadeIn(pg_steps))
        self.wait(2.5)
        self.play(FadeOut(py_group), FadeOut(pg_group), FadeOut(bridge_arrow))
'''


def build_studio_deck_html(args: argparse.Namespace) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>pg_bio Studio Deck — {args.title}</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    /* Exact ../ideia.me Theme Variables (css/main.css, css/custom.css & topic-colors.css) */
    :root {{
      --primary-blue: #2563eb;
      --primary-purple: #7c3aed;
      --primary-teal: #0d9488;
      --primary-light-blue: #0ea5e9;
      --primary-color: #5c9eff;
      --secondary-color: #b158fc;
      --accent-color: #ff5a9d;
      --topic-postgresql: #336791;
      --gradient-blue-purple: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
      --gradient-blue-teal: linear-gradient(135deg, var(--primary-blue), var(--primary-teal));
      --bg-color: #0f1117;
      --card-bg-color: #161b27;
      --code-bg-color: #1a1f2e;
      --text-color: #ccd6f6;
      --text-muted: #8892b0;
      --border-color: rgba(255, 255, 255, 0.09);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background-color: var(--bg-color);
      color: var(--text-color);
      font-family: 'JetBrains Mono', 'Fira Code', 'Cascadia Code', 'SF Mono', Menlo, Consolas, monospace;
      letter-spacing: -0.02em;
      line-height: 1.6;
      min-height: 100vh;
      display: flex;
      flex-direction: column;
    }}
    header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 14px 32px;
      border-bottom: 1px solid var(--border-color);
      background-color: var(--bg-color);
    }}
    .badge {{
      background-color: var(--topic-postgresql);
      color: #ffffff;
      border: 1px solid var(--primary-color);
      padding: 4px 12px;
      border-radius: 15px;
      font-size: 0.8rem;
      font-weight: 700;
    }}
    .rosetta-toggle {{
      display: flex;
      gap: 8px;
      align-items: center;
    }}
    button {{
      background: var(--code-bg-color);
      color: var(--text-color);
      border: 1px solid rgba(92, 158, 255, 0.25);
      padding: 7px 13px;
      border-radius: 8px;
      cursor: pointer;
      font-family: inherit;
      font-size: 0.82rem;
      font-weight: 600;
      transition: all 0.2s ease;
    }}
    button:hover, button.active {{
      background: var(--gradient-blue-purple);
      color: #ffffff;
      border-color: var(--primary-color);
    }}
    main {{
      flex: 1;
      max-width: 1320px;
      width: 100%;
      margin: 0 auto;
      padding: 24px 32px;
      display: grid;
      grid-template-columns: 1.05fr 0.95fr;
      gap: 24px;
    }}
    .card {{
      background-color: var(--card-bg-color);
      border: 1px solid var(--border-color);
      border-radius: 12px;
      padding: 24px;
      box-shadow: 0 2px 8px rgba(0, 0, 0, 0.4);
    }}
    h1 {{
      margin: 0 0 8px 0;
      font-size: 1.65rem;
      background: var(--gradient-blue-purple);
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
      line-height: 1.25;
    }}
    h2 {{
      margin-top: 14px;
      margin-bottom: 10px;
      font-size: 1.08rem;
      color: var(--primary-color);
      border-bottom: 1px solid var(--border-color);
      padding-bottom: 6px;
    }}
    pre {{
      background-color: var(--code-bg-color);
      border: 1px solid var(--border-color);
      padding: 14px;
      border-radius: 8px;
      overflow-x: auto;
      font-size: 0.82rem;
      line-height: 1.45;
      color: var(--text-color);
    }}
    .bar-row {{
      margin: 12px 0;
    }}
    .bar-label {{
      display: flex;
      justify-content: space-between;
      font-size: 0.82rem;
      margin-bottom: 4px;
    }}
    .bar-track {{
      background: var(--code-bg-color);
      height: 14px;
      border-radius: 999px;
      overflow: hidden;
    }}
    .bar-fill-py {{
      width: 100%;
      height: 100%;
      background: linear-gradient(90deg, #f72585, var(--accent-color));
    }}
    .bar-fill-pg {{
      width: 4%;
      height: 100%;
      background: var(--gradient-blue-teal);
    }}
    #mol-viewer {{
      width: 100%;
      height: 350px;
      position: relative;
      border-radius: 10px;
      overflow: hidden;
      background: var(--bg-color);
      border: 1px solid var(--border-color);
    }}
    .teaching-controls {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-top: 14px;
    }}
    #chemistry-lesson {{
      margin-top: 14px;
      padding: 14px;
      border-radius: 8px;
      background: rgba(92, 158, 255, 0.06);
      border-left: 4px solid var(--primary-color);
      font-size: 0.86rem;
      line-height: 1.5;
    }}
    footer {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 12px 32px;
      background-color: #0b0d13;
      border-top: 2px solid #4361ee;
      font-size: 0.82rem;
      color: var(--text-muted);
    }}
    footer a {{
      color: var(--primary-color);
      text-decoration: none;
      font-weight: 600;
    }}
    footer a:hover {{
      color: var(--secondary-color);
      text-decoration: underline;
    }}
    /* Voiced Hovering Video Companion Widget */
    #hover-video-widget {{
      position: fixed;
      bottom: 24px;
      right: 24px;
      width: 390px;
      background: rgba(22, 27, 39, 0.96);
      border: 2px solid var(--primary-color);
      border-radius: 14px;
      box-shadow: 0 14px 36px rgba(0, 0, 0, 0.7);
      z-index: 9999;
      backdrop-filter: blur(12px);
      overflow: hidden;
      transition: width 0.25s ease;
    }}
    #hover-video-widget.expanded {{ width: 680px; }}
    #hover-video-widget.minimized {{ width: 280px; }}
    #hover-video-widget.minimized .hover-body {{ display: none; }}
    #hover-video-widget.docked {{
      position: static;
      width: 100%;
      margin-bottom: 18px;
    }}
    .hover-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 7px 12px;
      background: #0b0d13;
      border-bottom: 1px solid var(--border-color);
      cursor: move;
      user-select: none;
    }}
    .hover-title {{
      font-size: 0.75rem;
      font-weight: 700;
      color: var(--primary-color);
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      max-width: 210px;
    }}
    .hover-body video {{
      width: 100%;
      display: block;
      background: #000;
    }}
    .act-chapters {{
      display: flex;
      gap: 4px;
      padding: 6px 8px;
      background: #0f1117;
      border-top: 1px solid var(--border-color);
    }}
    .act-btn {{
      flex: 1;
      padding: 4px 2px;
      font-size: 0.68rem;
      border-radius: 5px;
    }}
    .blog-links {{
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      margin-top: 10px;
    }}
    .blog-pill {{
      background: rgba(177, 88, 252, 0.12);
      border: 1px solid rgba(177, 88, 252, 0.35);
      color: var(--text-color);
      padding: 4px 10px;
      border-radius: 6px;
      font-size: 0.76rem;
      text-decoration: none;
    }}
    .blog-pill:hover {{
      background: var(--secondary-color);
      color: #fff;
    }}
  </style>
</head>
<body>
  <header>
    <div>
      <span class="badge">{args.layer}</span>
      <span style="margin-left: 12px; color: var(--text-muted); font-size: 0.85rem;">pg_bio Biotech-for-Engineers Studio</span>
    </div>
    <div class="rosetta-toggle">
      <audio id="deck-audio" src="voiceover.mp3" preload="auto"></audio>
      <button id="btn-vo" onclick="toggleVoiceover()">🔊 Play Neural Voiceover</button>
      <span style="color: var(--text-muted); font-size: 0.82rem; margin-left: 8px;">Lens:</span>
      <button id="btn-eng" class="active" onclick="setLens('engineer')">⚙️ Systems Engineer</button>
      <button id="btn-bio" onclick="setLens('biologist')">🧬 Structural Biologist</button>
    </div>
  </header>

  <main>
    <!-- Left Column: Rosetta Explanation + Live SQL + Benchmark Progress -->
    <div class="card">
      <h1>{args.title}</h1>
      <p id="lens-description" style="color: var(--text-color); font-size: 0.92rem; line-height: 1.55;">
        <strong>⚙️ Systems Engineering Lens:</strong> {args.cs_mechanism}
      </p>

      <h2>1. In-Database SQL Pipeline</h2>
      <pre><code>-- Zero-IPC execution inside PostgreSQL 18 shared buffers:
SELECT design_id,
       theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') AS rmsd,
       bio_angle(water_o, substrate_c, substrate_o) AS attack_angle
FROM metalloprotease_designs
WHERE has_coordinated_water_shell(active_site_coords, 'ZN', 2.4) = TRUE
  AND theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') &lt;= 0.90;</code></pre>

      <h2>2. SQL vs. Python Pipeline Progress ({args.speedup})</h2>
      <div class="bar-row">
        <div class="bar-label">
          <span>🐍 Python Socket Pull + Biopython/SciPy</span>
          <span style="color: var(--accent-color);">100% Baseline Wall-Clock</span>
        </div>
        <div class="bar-track"><div class="bar-fill-py"></div></div>
      </div>
      <div class="bar-row">
        <div class="bar-label">
          <span>🦀 <strong>pg_bio Native Rust (pgrx)</strong></span>
          <span style="color: var(--primary-light-blue);">{args.speedup}</span>
        </div>
        <div class="bar-track"><div class="bar-fill-pg"></div></div>
      </div>

      <h2>3. Integrated Blog Posts &amp; Video Library</h2>
      <div class="blog-links">
        <a class="blog-pill" href="{BLOG_URL}/pg-bio-videos/" target="_blank">🎬 Full Interactive Video Library (ideia.me)</a>
        <a class="blog-pill" href="{BLOG_URL}/pg-bio-intro" target="_blank">📖 Introducing pg_bio</a>
        <a class="blog-pill" href="{BLOG_URL}/postgresql-query-planner-rangefinder-gist-bioinformatics-engine" target="_blank">⚡ Query Planner &amp; RangeFinder GiST</a>
        <a class="blog-pill" href="{BLOG_URL}/5lox-triterpenoids-mmpbsa-free-energy-audit-pg-bio" target="_blank">🔬 500 ns 5-LOX MD Audit</a>
        <a class="blog-pill" href="live_demo.sql" target="_blank">💻 live_demo.sql</a>
        <a class="blog-pill" href="talk_track.md" target="_blank">📜 talk_track.md</a>
      </div>
    </div>

    <!-- Right Column: 3Dmol.js Interactive Chemistry Stage + Dock Slot -->
    <div class="card">
      <div id="docked-video-slot"></div>
      <h2>4. Interactive 3D Molecular Proof</h2>
      <div id="mol-viewer"></div>
      <div class="teaching-controls">
        <button onclick="highlightCore()">1. Highlight Catalytic Triad / Metal</button>
        <button onclick="highlightPocket()">2. Show 3D Morton Bounding Pocket</button>
        <button onclick="resetView()">3. Reset Cartoon Backbone</button>
      </div>
      <div id="chemistry-lesson">
        <strong>💡 Click a teaching button above:</strong> Explore how PostgreSQL predicates map directly to 3D atomic coordinates, catalytic residues, and pocket solvation.
      </div>
    </div>
  </main>

  <footer>
    <div>🌐 Blog &amp; 3D Discoveries: <a href="{BLOG_URL}" target="_blank" rel="noopener">{BLOG_URL}</a> | <a href="{BLOG_URL}/pg-bio-videos/" target="_blank" rel="noopener">3D Video Library</a></div>
    <div>🧬 <strong>pg_bio</strong> — The Multiomics Operating System in PostgreSQL</div>
    <div>🐙 GitHub Repository: <a href="{REPO_URL}" target="_blank" rel="noopener">{REPO_URL}</a></div>
  </footer>

  <!-- Voiced Hovering Video Companion Widget -->
  <div id="hover-video-widget">
    <div class="hover-header" id="hover-drag-handle">
      <span class="hover-title">🔊 {args.title}</span>
      <div style="display: flex; gap: 4px;">
        <button style="padding: 2px 6px; font-size: 0.7rem;" onclick="document.getElementById('hover-video-widget').classList.toggle('expanded')">⤢</button>
        <button style="padding: 2px 6px; font-size: 0.7rem;" onclick="toggleHoverDock()">📌</button>
        <button style="padding: 2px 6px; font-size: 0.7rem;" onclick="document.getElementById('hover-video-widget').classList.toggle('minimized')">—</button>
      </div>
    </div>
    <div class="hover-body">
      <video id="hover-video-player" controls playsinline preload="metadata" data-mp4-src="{args.slug}.mp4">
        <source src="{args.slug}.mp4" type="video/mp4">
      </video>
      <div class="act-chapters">
        <button class="act-btn" onclick="seekAct(0)">Act I</button>
        <button class="act-btn" onclick="seekAct(1)">Act II</button>
        <button class="act-btn" onclick="seekAct(2)">Act III</button>
        <button class="act-btn" onclick="seekAct(3)">Act IV</button>
        <button class="act-btn" onclick="seekAct(4)">Act V</button>
      </div>
    </div>
  </div>

  <script>
    const LENS_TEXT = {{
      engineer: "<strong>⚙️ Systems Engineering Lens:</strong> {args.cs_mechanism}",
      biologist: "<strong>🧬 Structural Biology Lens:</strong> {args.bio_problem}"
    }};

    function setLens(mode) {{
      document.getElementById('lens-description').innerHTML = LENS_TEXT[mode];
      document.getElementById('btn-eng').classList.toggle('active', mode === 'engineer');
      document.getElementById('btn-bio').classList.toggle('active', mode === 'biologist');
    }}

    let blobFallbackPromise = null;
    async function ensureVideoSource(vid) {{
      if (!vid.error && vid.networkState !== HTMLMediaElement.NETWORK_NO_SOURCE && vid.readyState > 0) {{
        return true;
      }}
      if (vid.dataset.blobLoaded === '1') return true;
      const mp4Url = vid.dataset.mp4Src || '{args.slug}.mp4';
      if (!blobFallbackPromise) {{
        blobFallbackPromise = fetch(mp4Url, {{ cache: 'reload' }})
          .then(r => {{
            if (!r.ok) throw new Error('HTTP ' + r.status);
            return r.blob();
          }})
          .then(b => {{
            const blobUrl = URL.createObjectURL(new Blob([b], {{ type: 'video/mp4' }}));
            vid.removeAttribute('src');
            while (vid.firstChild) vid.removeChild(vid.firstChild);
            vid.src = blobUrl;
            vid.dataset.blobLoaded = '1';
            vid.load();
            return new Promise(res => {{
              vid.addEventListener('loadedmetadata', () => res(true), {{ once: true }});
              vid.addEventListener('error', () => res(false), {{ once: true }});
              setTimeout(() => res(true), 1200);
            }});
          }})
          .catch(() => false);
      }}
      return blobFallbackPromise;
    }}

    async function safePlayVideo(vid) {{
      const btn = document.getElementById('btn-vo');
      const aud = document.getElementById('deck-audio');
      try {{
        if (vid.error || vid.networkState === HTMLMediaElement.NETWORK_NO_SOURCE) {{
          await ensureVideoSource(vid);
        }}
        await vid.play();
        if (btn) {{
          btn.textContent = '⏸ Pause Voiced Video';
          btn.classList.add('active');
        }}
      }} catch (err) {{
        const ok = await ensureVideoSource(vid);
        if (ok) {{
          try {{
            await vid.play();
            if (btn) {{
              btn.textContent = '⏸ Pause Voiced Video';
              btn.classList.add('active');
            }}
            return;
          }} catch (_) {{}}
        }}
        if (aud) {{
          aud.play().catch(() => {{}});
          if (btn) {{
            btn.textContent = '⏸ Pause Neural Voiceover';
            btn.classList.add('active');
          }}
        }}
      }}
    }}

    function toggleVoiceover() {{
      const vid = document.getElementById('hover-video-player');
      const aud = document.getElementById('deck-audio');
      const btn = document.getElementById('btn-vo');
      if (vid.paused && (!aud || aud.paused)) {{
        safePlayVideo(vid);
      }} else {{
        vid.pause();
        if (aud) aud.pause();
        if (btn) {{
          btn.textContent = '🔊 Play Voiced Video';
          btn.classList.remove('active');
        }}
      }}
    }}

    async function seekAct(actIdx) {{
      const vid = document.getElementById('hover-video-player');
      if (!vid.duration || vid.error) {{
        await ensureVideoSource(vid);
      }}
      if (vid.duration) {{
        vid.currentTime = (vid.duration / 5) * actIdx;
      }}
      safePlayVideo(vid);
      if (actIdx === 0) highlightCore();
      if (actIdx === 2) highlightPocket();
      if (actIdx === 4) resetView();
    }}

    function toggleHoverDock() {{
      const widget = document.getElementById('hover-video-widget');
      const slot = document.getElementById('docked-video-slot');
      widget.classList.toggle('docked');
      if (widget.classList.contains('docked')) {{
        slot.appendChild(widget);
      }} else {{
        document.body.appendChild(widget);
      }}
    }}

    window.addEventListener('DOMContentLoaded', () => {{
      const widget = document.getElementById('hover-video-widget');
      const handle = document.getElementById('hover-drag-handle');
      const vid = document.getElementById('hover-video-player');
      const btn = document.getElementById('btn-vo');
      vid.addEventListener('error', () => {{ ensureVideoSource(vid); }}, true);
      vid.addEventListener('play', () => {{
        if (btn) {{ btn.textContent = '⏸ Pause Voiced Video'; btn.classList.add('active'); }}
      }});
      vid.addEventListener('pause', () => {{
        if (btn) {{ btn.textContent = '🔊 Play Voiced Video'; btn.classList.remove('active'); }}
      }});
      let lastAct = -1;
      vid.addEventListener('timeupdate', () => {{
        if (!vid.duration) return;
        const curAct = Math.min(4, Math.floor((vid.currentTime / vid.duration) * 5));
        if (curAct !== lastAct) {{
          lastAct = curAct;
          document.querySelectorAll('.act-btn').forEach((b, i) => b.classList.toggle('active', i === curAct));
          if (curAct === 0) {{ setLens('biologist'); highlightCore(); }}
          else if (curAct === 2) {{ setLens('engineer'); highlightPocket(); }}
          else if (curAct === 4) {{ resetView(); }}
        }}
      }});
      let dragging = false, startX = 0, startY = 0, origR = 24, origB = 24;
      handle.addEventListener('mousedown', (e) => {{
        if (widget.classList.contains('docked') || e.target.tagName === 'BUTTON') return;
        dragging = true; startX = e.clientX; startY = e.clientY;
        const r = widget.getBoundingClientRect();
        origR = window.innerWidth - r.right; origB = window.innerHeight - r.bottom;
      }});
      window.addEventListener('mousemove', (e) => {{
        if (!dragging) return;
        widget.style.right = Math.max(8, origR - (e.clientX - startX)) + 'px';
        widget.style.bottom = Math.max(8, origB - (e.clientY - startY)) + 'px';
      }});
      window.addEventListener('mouseup', () => {{ dragging = false; }});
    }});

    const DEMO_PDB = `
ATOM      1  N   HIS A  94       0.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  HIS A  94       1.458   0.000   0.000  1.00 90.00           C
ATOM      3  C   HIS A  94       2.009   1.420   0.000  1.00 90.00           C
ATOM      4  O   HIS A  94       1.251   2.390   0.000  1.00 90.00           O
ATOM      5  CB  HIS A  94       1.988  -0.760  -1.216  1.00 90.00           C
ATOM      6  ND1 HIS A  94       3.350  -0.950  -1.350  1.00 95.00           N
ATOM      7  N   HIS A  96       4.100   0.800  -2.200  1.00 92.00           N
ATOM      8  CA  HIS A  96       5.200   0.100  -2.800  1.00 92.00           C
ATOM      9  ND1 HIS A  96       4.850  -1.250  -3.200  1.00 95.00           N
ATOM     10  N   GLU A 106       2.400  -3.400  -2.100  1.00 91.00           N
ATOM     11  CA  GLU A 106       1.900  -4.200  -1.000  1.00 91.00           C
ATOM     12  OE1 GLU A 106       3.100  -2.800  -0.400  1.00 96.00           O
HETATM   13 ZN    ZN A 201       3.950  -1.650  -2.150  1.00 99.00          ZN
HETATM   14  O   HOH A 301       3.750  -2.950  -1.150  1.00 98.00           O
END`;

    let viewer = null;
    window.addEventListener('DOMContentLoaded', () => {{
      viewer = $3Dmol.createViewer(document.getElementById('mol-viewer'), {{
        backgroundColor: '#0f1117'
      }});
      viewer.addModel(DEMO_PDB, 'pdb');
      resetView();
    }});

    function highlightCore() {{
      if (!viewer) return;
      viewer.removeAllLabels();
      viewer.setStyle({{}}, {{ stick: {{ radius: 0.18, colorscheme: 'cyanCarbon' }} }});
      viewer.setStyle({{ resn: 'ZN' }}, {{ sphere: {{ scale: 0.65, color: '#FFC107' }} }});
      viewer.setStyle({{ resn: 'HOH' }}, {{ sphere: {{ scale: 0.45, color: '#ff5a9d' }} }});
      viewer.addLabel('Zn2+ Cofactor', {{ position: {{ x: 3.95, y: -1.65, z: -2.15 }}, backgroundColor: '#FFC107', fontColor: '#000' }});
      viewer.addLabel('Nucleophilic HOH (2.05 Å)', {{ position: {{ x: 3.75, y: -2.95, z: -1.15 }}, backgroundColor: '#ff5a9d', fontColor: '#fff' }});
      viewer.render();
      document.getElementById('chemistry-lesson').innerHTML =
        "<strong>🔬 Chemistry Insight (Catalytic Coordination):</strong> The Zn<sup>2+</sup> cofactor polarizes the bound water molecule (pink sphere), lowering its pK<sub>a</sub> so Glu106 can abstract a proton and generate a hydroxide nucleophile. In <code>pg_bio</code>, <code>has_coordinated_water_shell(coords, 'ZN', 2.4)</code> and <code>theozyme_rmsd()</code> verify this exact sub-angstrom geometry in &lt;1 ms.";
    }}

    function highlightPocket() {{
      if (!viewer) return;
      viewer.removeAllLabels();
      viewer.setStyle({{}}, {{ stick: {{ radius: 0.22, colorscheme: 'magentaCarbon' }}, sphere: {{ scale: 0.28 }} }});
      viewer.addSphere({{ center: {{ x: 3.95, y: -1.65, z: -2.15 }}, radius: 3.5, color: '#5c9eff', alpha: 0.22 }});
      viewer.addLabel('Morton Z-Order 3.5 Å Bounding Shell', {{ position: {{ x: 1.5, y: 1.2, z: 0.0 }}, backgroundColor: '#336791', fontColor: '#fff' }});
      viewer.render();
      document.getElementById('chemistry-lesson').innerHTML =
        "<strong>⚙️ Systems-to-Bio Insight (Z-Order Bounding Sphere):</strong> Instead of computing Euclidean distances across millions of atoms in Python, <code>z_order_encode(x,y,z)</code> maps 3D coordinates onto a 1D space-filling Morton curve so standard PostgreSQL B-Tree indexes can slice out this 3.5 Å contact shell instantly.";
    }}

    function resetView() {{
      if (!viewer) return;
      viewer.removeAllLabels();
      viewer.removeAllShapes();
      viewer.setStyle({{}}, {{ stick: {{ radius: 0.16, colorscheme: 'Jmol' }} }});
      viewer.setStyle({{ resn: 'ZN' }}, {{ sphere: {{ scale: 0.55, color: '#5c9eff' }} }});
      viewer.zoomTo();
      viewer.render();
      document.getElementById('chemistry-lesson').innerHTML =
        "<strong>💡 Interactive Stage Ready:</strong> Click Button 1 to inspect the Zn<sup>2+</sup>–Water catalytic triad or Button 2 to visualize the 3D Morton B-Tree bounding query.";
    }}
  </script>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage a pg_bio feature video & tutoring bundle.")
    parser.add_argument("--slug", required=True, help="Episode directory slug (e.g. ep01_rangefinder_gist)")
    parser.add_argument("--title", required=True, help="Episode title")
    parser.add_argument("--layer", required=True, help="pg_bio Layer name/number")
    parser.add_argument("--bio-problem", required=True, help="1-sentence biological problem statement")
    parser.add_argument("--cs-mechanism", required=True, help="1-sentence Rust/SQL algorithmic mechanism")
    parser.add_argument("--speedup", required=True, help="Benchmark speedup headline")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[4]
    out_dir = repo_root / "videos" / args.slug
    out_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / "talk_track.md").write_text(build_talk_track(args), encoding="utf-8")
    (out_dir / "live_demo.sql").write_text(build_live_demo_sql(args), encoding="utf-8")
    (out_dir / "slides.marp.md").write_text(build_marp_slides(args), encoding="utf-8")
    (out_dir / "manim_scene.py").write_text(build_manim_scene(args), encoding="utf-8")
    (out_dir / "studio_deck.html").write_text(build_studio_deck_html(args), encoding="utf-8")

    mp4_path = out_dir / f"{args.slug}.mp4"
    renderer_script = Path(__file__).resolve().parent / "render_feature_mp4.py"
    subprocess.run(
        [
            "uv",
            "run",
            "--with",
            "pillow",
            "--with",
            "numpy",
            "--with",
            "edge-tts",
            "python3",
            str(renderer_script),
            "--out",
            str(mp4_path),
            "--title",
            args.title,
            "--layer",
            args.layer,
            "--bio-problem",
            args.bio_problem,
            "--cs-mechanism",
            args.cs_mechanism,
            "--speedup",
            args.speedup,
        ],
        check=True,
    )

    print(f"Created video & tutoring bundle in: {out_dir}")
    for fname in [
        "talk_track.md",
        "live_demo.sql",
        "slides.marp.md",
        "manim_scene.py",
        "studio_deck.html",
        "voiceover.mp3",
        "captions.srt",
        f"{args.slug}.mp4",
    ]:
        print(f"  - {out_dir / fname}")


if __name__ == "__main__":
    main()
