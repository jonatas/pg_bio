#!/usr/bin/env python3
"""
Builds Series 01: "Hello World, pg_bio! — Intro to Biotech & Basic Multiomics in PostgreSQL"
inside `videos/series_01_hello_world_pg_bio/` and syncs the interactive video library to
`../ideia.me/assets/videos/pg_bio/series_01_hello_world_pg_bio/`.

Includes:
  - Voiced Hovering Video Companion (#hover-video-widget: draggable, dockable, minimizable,
    with 5 Act chapter triggers synchronized to the 3Dmol.js viewer)
  - Direct cross-links to companion `ideia.me` blog posts, `pg_bio` docs, notebooks, and SQL scripts
  - Root `videos/README.md` and `videos/index.html` catalog linking all material
"""

import argparse
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List

BLOG_URL = "https://ideia.me"
REPO_URL = "https://github.com/jonatas/pg_bio"

EPISODES: List[Dict[str, Any]] = [
    {
        "ep_num": 1,
        "slug": "ep01_central_dogma_and_genomics",
        "visual_mode": "central_dogma_gist",
        "title": "Ep 1: Hello World, pg_bio! The Central Dogma in SQL",
        "subtitle": "DNA -> RNA -> Protein Physics & 1D Genomic Safe Harbors (int4range GiST)",
        "layer": "Layers 1 & 6: Genomics (int4range GiST) & Sequence Physics",
        "bio_problem": (
            "Biology's Central Dogma compiles 4 DNA bases (A, T, G, C) into 20 amino acids. "
            "When inserting synthetic genes via CRISPR, you must find a Safe Harbor near a tissue promoter "
            "without intersecting tumor suppressor genes, then verify the protein's charge (pI) and hydropathy."
        ),
        "cs_mechanism": (
            "PostgreSQL 1D int4range GiST overlap (&&) and distance (<->) operators scan billions of base pairs, "
            "while single-pass Rust byte slices compute bio_translate_dna(), isoelectric_point(), and gravy_score() in 200 ns."
        ),
        "speedup": "200x faster than BioPython (~200 ns/protein in Rust)",
        "docs_ref": "../../docs/07_genomics_io.md",
        "notebook_ref": "../../experiments/06_genomics_io_tutorial.ipynb",
        "blogs": [
            {"title": "Introducing pg_bio", "url": f"{BLOG_URL}/pg-bio-intro"},
            {"title": "Genomic Safe-Harbor Architecture", "url": f"{BLOG_URL}/teaflon-synthetic-biology-part-8-genomic-architecture"},
            {"title": "Rust pgrx Memory Contexts & TOAST", "url": f"{BLOG_URL}/postgresql-rust-ffi-pgrx-memory-contexts-varlena-toast-pg-bio"},
        ],
        "bullets": [
            "• Central Dogma in SQL: bio_translate_dna('ATGTTC...') -> 'MF...'",
            "• 1D Safe Harbor GiST: safe.locus <-> promoter.locus < 5000",
            "• Isoelectric Point (pI): 100-step Henderson-Hasselbalch bisection",
            "• Hydropathy (GRAVY) & Shannon Entropy in a single byte pass",
        ],
        "rosetta_mappings": [
            [
                "int4range + GiST Index (&& and <->)",
                "1D spatial interval tree with O(log N) overlap & nearest-neighbor operators",
                "Finds CRISPR genomic Safe Harbors near promoters while avoiding oncogenes",
            ],
            [
                "Compile-Time Match Lookup (bio_translate_dna)",
                "Zero-allocation 3-byte ASCII codon slice -> amino acid byte stream",
                "Translates raw genomic DNA reading frames into functional protein sequences",
            ],
            [
                "Root-Finding Bisection (isoelectric_point & gravy_score)",
                "Single-pass byte histogram + 100-step Henderson-Hasselbalch solver",
                "Predicts whether an enzyme stays soluble or precipitates in acid/salt vats",
            ],
        ],
        "sql_lines": [
            "bio_demo=# SELECT",
            "    bio_translate_dna('ATGTTCGAAGGTTTCGAACGTCGT')       AS peptide_seq,",
            "    round(isoelectric_point('MFEGFERRLVDWPSTDK')::numeric, 2) AS pI_charge,",
            "    round(gravy_score('MFEGFERRLVDWPSTDK')::numeric, 3)       AS hydropathy,",
            "    round(sequence_entropy('MFEGFERRLVDWPSTDK')::numeric, 3)  AS shannon_bits,",
            "    prosite_match('MFEGFERRLVDWPSTDKHCH', 'H-x-H')            AS metal_motif",
            "WHERE int4range(10500, 11800) <-> int4range(8000, 9200) < 5000;",
        ],
        "sql_results": [
            "┌─────────────┬───────────┬────────────┬──────────────┬─────────────┐",
            "│ peptide_seq │ pI_charge │ hydropathy │ shannon_bits │ metal_motif │",
            "├─────────────┼───────────┼────────────┼──────────────┼─────────────┤",
            "│ MFEGFERR    │      4.68 │     -1.124 │        3.585 │ t           │",
            "└─────────────┴───────────┴────────────┴──────────────┴─────────────┘",
            "Time: 0.214 ms  (200x faster than BioPython)",
        ],
        "narrations": [
            (
                "Hello World, and welcome to Episode 1 of pg_bio! "
                "If you are a software or database engineer, biology can look like a wall of jargon. "
                "But the Central Dogma of biology is literally a compiler: DNA is 4-letter source code on disk, "
                "RNA is an ephemeral message queue, and proteins are compiled 3D state machines."
            ),
            (
                "In traditional bioinformatics, engineers export gigabytes of FASTA and VCF files to Python just to translate codons "
                "or check genomic intervals. With pg_bio, that entire compiler runs natively inside PostgreSQL shared buffers in Rust."
            ),
            (
                "Here is your first Systems-to-Biotech bridge: a human chromosome is just a 1D integer coordinate space. "
                "Using PostgreSQL int4range and GiST indexes, finding a CRISPR safe-harbor locus near a promoter without touching "
                "a tumor suppressor gene is a sub-millisecond interval query."
            ),
            (
                "Look at psql: in one query, bio_translate_dna compiles DNA codons into amino acids, "
                "isoelectric_point runs 100-step Henderson-Hasselbalch bisection to find net charge, "
                "and gravy_score measures water-repellency in 200 nanoseconds."
            ),
            (
                "You just executed genomics and protein sequence physics in 0.2 milliseconds inside SQL. "
                "Explore all articles at ideia dot me and check out the pg_bio repository on GitHub."
            ),
        ],
    },
    {
        "ep_num": 2,
        "slug": "ep02_3d_pockets_and_morton_curves",
        "visual_mode": "morton_3d_pocket",
        "title": "Ep 2: From 1D Strings to 3D Atoms with Morton Curves",
        "subtitle": "3D PDB Folding, 64-Bit Z-Order B-Tree Bounding Boxes & Kabsch RMSD",
        "layer": "Layers 2, 6 & 8: 3D Spatial Z-Order Curves & Structural Superposition",
        "bio_problem": (
            "Once translated, a 1D amino acid string folds into a 3D atomic cloud measured in Angstroms. "
            "Catalysis happens inside a tiny 3.5 to 5.0 Angstrom 3D pocket—where residues that are 100 positions "
            "apart in the string snap together in 3D space to coordinate a metal cofactor."
        ),
        "cs_mechanism": (
            "z_order_encode(x, y, z) bit-interleaves 21-bit 3D Cartesian floats into a single 64-bit BIGINT B-Tree key, "
            "while kabsch_rmsd() and theozyme_rmsd() run stack-allocated 3x3 Jacobi SVD to align 3D pockets in microseconds."
        ),
        "speedup": "120x faster 3D pocket extraction (Zero KD-Tree build time)",
        "docs_ref": "../../docs/03_spatial_indexing.md",
        "notebook_ref": "../../experiments/08_pdb_parsing.ipynb",
        "blogs": [
            {"title": "Indexing 3D Space Like a Video Game", "url": f"{BLOG_URL}/indexing-3d-space-like-a-video-game"},
            {"title": "Store & Query 3D Protein Structures in SQL", "url": f"{BLOG_URL}/store-query-3d-protein-structures-postgresql"},
            {"title": "Query Planner & RangeFinder Pocket GiST", "url": f"{BLOG_URL}/postgresql-query-planner-rangefinder-gist-bioinformatics-engine"},
        ],
        "bullets": [
            "• Native PDB Parser: bio_parse_pdb(pdb_text) -> relational atom rows",
            "• 3D -> 1D Morton Key: z_order_encode(x, y, z) -> BIGINT B-Tree index",
            "• Rotational Alignment: kabsch_rmsd() & theozyme_rmsd() via 3x3 SVD",
            "• Bond & Attack Geometry: bio_distance(), bio_angle(), bio_dihedral()",
        ],
        "rosetta_mappings": [
            [
                "64-Bit Morton Space-Filling Curve (z_order_encode)",
                "Interleaves (X, Y, Z) bits via split_by_3 into a 1D B-Tree indexable BIGINT",
                "Slices out 4.0 Å 3D enzyme binding pockets across millions of atoms instantly",
            ],
            [
                "Set-Returning Streaming Parser (bio_parse_pdb)",
                "Zero-copy Rust iterator parsing ATOM/HETATM records into SQL tuples",
                "Turns raw AlphaFold / RCSB PDB files into queryable 3D atomic tables",
            ],
            [
                "Stack-Allocated 3x3 Jacobi SVD (theozyme_rmsd)",
                "Computes optimal 3D rotation matrix without external LAPACK/SciPy calls",
                "Verifies sub-angstrom catalytic triad & Zn2+ water shell superposition",
            ],
        ],
        "sql_lines": [
            "bio_demo=# WITH pocket_bounds AS (",
            "    SELECT z_order_encode(1.0, -4.0, -4.0) AS min_z,",
            "           z_order_encode(6.0,  1.0,  1.0) AS max_z",
            ")",
            "SELECT atom_id, residue_name, x, y, z",
            "FROM teaflon_atoms, pocket_bounds",
            "WHERE z_curve_index BETWEEN min_z AND max_z",
            "  AND theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') <= 0.90;",
        ],
        "sql_results": [
            "┌─────────┬──────────────┬───────┬───────┬───────┬───────────────┐",
            "│ atom_id │ residue_name │   x   │   y   │   z   │ theozyme_rmsd │",
            "├─────────┼──────────────┼───────┼───────┼───────┼───────────────┤",
            "│      13 │ ZN           │  3.95 │ -1.65 │ -2.15 │         0.412 │",
            "│      14 │ HOH          │  3.75 │ -2.95 │ -1.15 │         0.412 │",
            "└─────────┴──────────────┴───────┴───────┴───────┴───────────────┘",
            "Time: 0.490 ms  (B-Tree Z-Order Index Scan + 3x3 Jacobi SVD)",
        ],
        "narrations": [
            (
                "Welcome to Episode 2 of Hello World, pg_bio! "
                "In Episode 1, we treated DNA and proteins as 1D strings. "
                "But in reality, a protein folds into a 3D molecular machine where catalysis happens inside a 4-Angstrom binding pocket."
            ),
            (
                "Normally, finding which atoms sit inside a 3D sphere requires pulling millions of coordinates into Python "
                "and building an in-memory KD-tree or running brute-force Pythagorean distance loops."
            ),
            (
                "pg_bio solves this with Morton coding—also known as a Z-Order space-filling curve. "
                "Our Rust function z_order_encode interleaves the bits of X, Y, and Z coordinates into a single 64-bit integer. "
                "That means a standard 1D PostgreSQL B-Tree index now acts as a 3D spatial bounding box!"
            ),
            (
                "Even better, once the B-Tree narrows down the pocket atoms, theozyme_rmsd runs a zero-allocation "
                "3-by-3 Jacobi SVD right inside SQL to rotationally align the catalytic metal and nucleophilic water."
            ),
            (
                "That is 3D structural biology powered by B-Trees and linear algebra inside PostgreSQL. "
                "Visit ideia dot me and the pg_bio GitHub repo to try the interactive 3D notebooks."
            ),
        ],
    },
    {
        "ep_num": 3,
        "slug": "ep03_scrna_csr_and_attention_maps",
        "visual_mode": "csr_sparse_graph",
        "title": "Ep 3: Transcriptomics & CSR Sparse Attention Graphs",
        "subtitle": "Slashing 20-Billion-Cell scRNA Matrices & Mapping Allosteric Hotspots",
        "layer": "Layers 3 & 9: Single-Cell CSR Arrays, Neural Attention & Allostery",
        "bio_problem": (
            "Every cell in your body shares the exact same DNA, yet a neuron and a T-cell behave completely differently "
            "because of which genes they express (Transcriptomics). Sequencing 1 million cells across 20,000 genes "
            "creates a 20-billion-cell matrix that is 90% zeroes—crushing traditional SQL schemas and RAM."
        ),
        "cs_mechanism": (
            "Compressed Sparse Row (CSR) storage (SparseAttentionMap, create_sparse_map, build_energy_hotspot_map) "
            "stores only non-zero contacts and expression counts, reducing O(L^2) dense scans to O(nnz) nanosecond slice lookups."
        ),
        "speedup": "31.1x faster than Python/NumPy & 90% less RAM footprint",
        "docs_ref": "../../docs/12_mmpbsa_free_energy_auditing.md",
        "notebook_ref": "../../experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb",
        "blogs": [
            {"title": "500 ns 5-LOX Triterpenoid MD Audit in SQL", "url": f"{BLOG_URL}/5lox-triterpenoids-mmpbsa-free-energy-audit-pg-bio"},
            {"title": "Analyzing Hi-C 3D Genomic Contact Maps in SQL", "url": f"{BLOG_URL}/analyzing-hic-3d-genomic-contact-maps-sql"},
        ],
        "bullets": [
            "• Sparse CSR Type: SparseAttentionMap (row_ptr, col_idx, values)",
            "• Transformer Attention in SQL: get_top_interacting_residues(map, res_id, k)",
            "• Cryptic Pocket Discovery: discover_cryptic_hotspots() in 7.97 ms",
            "• Allosteric Mechanism Index: allosteric_mechanism_index() over CSR graphs",
        ],
        "rosetta_mappings": [
            [
                "Compressed Sparse Row (CSR: row_ptr, col_idx, val)",
                "Drops 90%+ zero entries from 20-billion-point matrices and LxL contact maps",
                "Stores single-cell RNA expression & protein attention graphs in L1 cache",
            ],
            [
                "O(nnz) Slice Traversal (get_top_interacting_residues)",
                "Indexes directly into row_ptr[i..i+1] without scanning dense zeroes",
                "Extracts long-range amino acid couplings from ESM-2 transformer heads",
            ],
            [
                "CSR Energy Network Profiler (discover_cryptic_hotspots)",
                "Combines MM/PBSA ledger auditing with sparse graph neighborhood filtering",
                "Uncovers hidden allosteric anchor residues outside the known active site",
            ],
        ],
        "sql_lines": [
            "bio_demo=# SELECT",
            "    system_id,",
            "    get_top_interacting_residues(attention_map, 94, 3) AS top_coupled_residues,",
            "    discover_cryptic_hotspots(",
            "        build_energy_hotspot_map(677, residue_ids, energy_arr, 677, -0.15),",
            "        677, ARRAY[68, 108, 109, 110], 0.35, 3",
            "    ) AS cryptic_allosteric_anchors",
            "FROM protein_attention_and_energy_maps;",
        ],
        "sql_results": [
            "┌────────────┬──────────────────────┬────────────────────────────┐",
            "│ system_id  │ top_coupled_residues │ cryptic_allosteric_anchors │",
            "├────────────┼──────────────────────┼────────────────────────────┤",
            "│ 5LOX_3B4F  │ {96,106,201}         │ {126,29,66}                │",
            "└────────────┴──────────────────────┴────────────────────────────┘",
            "Time: 0.580 ms  (31.14x faster than dense NumPy matrix unpacking)",
        ],
        "narrations": [
            (
                "Welcome to Episode 3 of Hello World, pg_bio! "
                "If every cell in your body has the exact same DNA database, why is a brain cell different from a liver cell? "
                "The answer is Transcriptomics and allosteric regulation—which genes are turned on, and how distant residues talk to each other."
            ),
            (
                "When scientists run Single-Cell RNA sequencing on 1 million cells across 20,000 human genes, "
                "it creates a 20-billion-point matrix. Because a cell only expresses about 2,000 genes at a time, "
                "over 90 percent of that matrix is zeroes. Storing dense columns in SQL or Pandas wastes gigabytes of RAM."
            ),
            (
                "Systems engineers already know the cure: Compressed Sparse Row, or C-S-R. "
                "pg_bio implements a native SparseAttentionMap type using row pointers, column indices, and non-zero values."
            ),
            (
                "Whether you are unpacking single-cell RNA dosage, traversing ESM-2 transformer attention heads with "
                "get_top_interacting_residues, or discovering cryptic allosteric pockets across 500-nanosecond molecular dynamics trajectories, "
                "CSR scans only the non-zero edges in microseconds."
            ),
            (
                "That slashes memory by 90 percent and runs 31 times faster than dense Python matrices. "
                "Read the full deep dive at ideia dot me and explore pg_bio on GitHub."
            ),
        ],
    },
    {
        "ep_num": 4,
        "slug": "ep04_joint_vector_multiomics",
        "visual_mode": "hnsw_twilight_zone",
        "title": "Ep 4: The Multiomics Holy Grail — Vectors & Live UniProt",
        "subtitle": "1280-D HNSW Embeddings, Twilight-Zone Homology & Native UniProt HTTP in SQL",
        "layer": "Layers 4 & 5: Joint Vector Spaces (halfvec HNSW) & Live UniProt FDW",
        "bio_problem": (
            "99% of microbial proteins in extreme environments (deep-sea vents, salt flats, acid lakes) are uncharacterized "
            "'Dark Proteome' orphans. Their sequences have mutated so heavily (<30% identity) that BLAST fails to recognize "
            "that they fold into the exact same 3D enzyme shape."
        ),
        "cs_mechanism": (
            "1280-D halfvec HNSW cosine distance (<=>) narrows 2.34M proteins down to 15 candidates in 0.4 ms, "
            "short-circuiting into O(N) two-row Smith-Waterman DP and enriching winners live via bio_fetch_uniprot()."
        ),
        "speedup": "Scans 2.34M proteins in 0.4 ms + Zero external ETL scripts",
        "docs_ref": "../../docs/01_vector_homology.md",
        "notebook_ref": "../../experiments/09_huggingface_esm2_discovery.ipynb",
        "blogs": [
            {"title": "Why Proteins Are High-Dimensional Vectors", "url": f"{BLOG_URL}/why-proteins-are-just-high-dimensional-vectors"},
            {"title": "Scaling pg_bio Swarm: 8.8M Edges in Dark Proteome", "url": f"{BLOG_URL}/scaling-pg-bio-swarm-8-8-million-edges-all-time-twilight-records-dark-proteome"},
            {"title": "The 27-Hour HNSW Index & halfvec", "url": f"{BLOG_URL}/the-27-hour-index"},
        ],
        "bullets": [
            "• 50% RAM Savings: 1280-D halfvec HNSW graph indexing (<=>)",
            "• Twilight-Zone Bell Curve: 3D fold convergence vs sequence divergence",
            "• O(N) Rolling-Buffer DP: smith_waterman_score() on HNSW survivors",
            "• Live UniProt SQL Client: bio_fetch_uniprot() & bio_search_uniprot()",
        ],
        "rosetta_mappings": [
            [
                "1280-D halfvec HNSW Graph Index (<=>)",
                "16-bit float approximate nearest neighbor graph in shared buffers",
                "Finds proteins with identical 3D folds across 2.34M sequences in 0.4 ms",
            ],
            [
                "Query-Planner Funnel + Two-Row Smith-Waterman DP",
                "HNSW prunes 99.999% first; O(N) rolling-buffer DP scores only 15 survivors",
                "Uncovers extremophile Twilight-Zone orthologs (15%–75% sequence identity)",
            ],
            [
                "Synchronous Rust HTTP SRF (bio_fetch_uniprot)",
                "ureq + serde JSON projection into Postgres UniprotEntry composite types",
                "Enriches orphan discoveries with live organism & taxonomy inside pure SQL",
            ],
        ],
        "sql_lines": [
            "bio_demo=# WITH candidates AS (",
            "    SELECT uniprot_id, sequence, (embedding <=> :bait_vec) AS cos_dist",
            "    FROM dark_proteome ORDER BY embedding <=> :bait_vec LIMIT 10",
            ")",
            "SELECT c.uniprot_id, u.organism, round(c.cos_dist::numeric, 4) AS dist,",
            "       smith_waterman_score(c.sequence, :bait_seq) AS sw_score",
            "FROM candidates c",
            "CROSS JOIN LATERAL bio_fetch_uniprot(c.uniprot_id) u;",
        ],
        "sql_results": [
            "┌────────────┬───────────────────────────┬────────┬──────────┐",
            "│ uniprot_id │         organism          │  dist  │ sw_score │",
            "├────────────┼───────────────────────────┼────────┼──────────┤",
            "│ A0A5B9DCV6 │ Deep-sea Asgard archaeon  │ 0.0682 │      312 │",
            "│ L0JL84     │ Halophilic archaeon       │ 0.0741 │      289 │",
            "└────────────┴───────────────────────────┴────────┴──────────┘",
            "Time: 0.420 ms local HNSW + live UniProt REST enrichment",
        ],
        "narrations": [
            (
                "Welcome to Episode 4, the finale of our Hello World, pg_bio intro series! "
                "Now we unite all layers into the Multiomics Holy Grail: Joint Vector Spaces and live biological database integration."
            ),
            (
                "In extreme environments like deep-sea hydrothermal vents, evolution mutates an enzyme's amino acid sequence "
                "until less than 25 percent of the letters match known proteins. Traditional sequence search tools like BLAST go blind "
                "in this Twilight Zone, even when the 3D catalytic fold is preserved."
            ),
            (
                "By storing 1280-dimensional ESM-2 protein foundation model embeddings as 16-bit halfvecs indexed with HNSW, "
                "PostgreSQL narrows 2.34 million dark proteins down to the top 15 structural twins in 0.4 milliseconds."
            ),
            (
                "Then the PostgreSQL query planner short-circuits: it runs our compiled Rust Smith-Waterman dynamic programming "
                "and Henderson-Hasselbalch pI calculations on only those 15 survivors, and calls bio_fetch_uniprot via LATERAL join "
                "to pull live organism and taxonomy metadata straight from UniProt into your SQL table."
            ),
            (
                "From 1D DNA intervals, to 3D Morton pockets, to sparse CSR graphs, to HNSW vector discovery—you now know "
                "the core multiomics stack in PostgreSQL. Dive deeper at ideia dot me and star pg_bio on GitHub!"
            ),
        ],
    },
]


def build_series_hub_html(episodes: List[Dict[str, Any]]) -> str:
    cards_html = []
    for idx, ep in enumerate(episodes):
        slug = ep["slug"]
        blog_pills = "".join(
            f'<a class="blog-pill" href="{b["url"]}" target="_blank">📖 {b["title"]}</a>'
            for b in ep["blogs"]
        )
        cards_html.append(
            f"""
      <div class="ep-card {'active' if idx == 0 else ''}" id="card-{slug}" onclick="selectEpisode({idx})">
        <div class="ep-header">
          <span class="badge">Episode 0{ep['ep_num']}</span>
          <span class="speedup-pill">{ep['speedup']}</span>
        </div>
        <h3>{ep['title']}</h3>
        <p class="ep-sub">{ep['subtitle']}</p>
        <p class="ep-desc">{ep['bio_problem']}</p>
        <div class="ep-links" onclick="event.stopPropagation();">
          <button onclick="selectEpisode({idx})">🔊 Load in Hovering Video &amp; 3D Stage</button>
          <a class="btn-link" href="{slug}/studio_deck.html" target="_blank">🧪 Standalone Studio Deck</a>
          <a class="btn-link" href="{slug}/talk_track.md" target="_blank">📜 Talk Track</a>
          <a class="btn-link" href="{slug}/live_demo.sql" target="_blank">💻 SQL Script</a>
          <a class="btn-link" href="{ep['docs_ref']}" target="_blank">📘 pg_bio Doc</a>
        </div>
        <div class="ep-links" style="margin-top: 8px;" onclick="event.stopPropagation();">
          {blog_pills}
        </div>
      </div>"""
        )

    first = episodes[0]
    episodes_json = json.dumps(episodes)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>Hello World, pg_bio! — Interactive Series Explorer &amp; Voiced Hovering Video</title>
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <style>
    /* Exact ../ideia.me Theme Variables */
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
      padding-bottom: 80px;
    }}
    header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 14px 36px;
      border-bottom: 1px solid var(--border-color);
      background-color: #0b0d13;
    }}
    .series-badge {{
      background-color: var(--topic-postgresql);
      color: #fff;
      border: 1px solid var(--primary-color);
      padding: 5px 14px;
      border-radius: 999px;
      font-size: 0.82rem;
      font-weight: 700;
    }}
    main {{
      flex: 1;
      max-width: 1400px;
      width: 100%;
      margin: 0 auto;
      padding: 24px 36px;
    }}
    h1 {{
      margin: 0 0 6px 0;
      font-size: 1.9rem;
      background: var(--gradient-blue-purple);
      -webkit-background-clip: text;
      -webkit-text-fill-color: transparent;
    }}
    .interactive-split {{
      display: grid;
      grid-template-columns: 1.05fr 0.95fr;
      gap: 22px;
      margin: 20px 0 28px 0;
    }}
    .stage-card {{
      background: var(--card-bg-color);
      border: 1px solid var(--border-color);
      border-radius: 14px;
      padding: 22px;
    }}
    pre {{
      background: var(--code-bg-color);
      border: 1px solid var(--border-color);
      padding: 12px;
      border-radius: 8px;
      font-size: 0.8rem;
      overflow-x: auto;
      color: var(--text-color);
    }}
    #hub-3d-viewer {{
      width: 100%;
      height: 340px;
      position: relative;
      border-radius: 10px;
      overflow: hidden;
      background: #090c12;
      border: 1px solid var(--border-color);
    }}
    .ep-grid {{
      display: grid;
      grid-template-columns: repeat(2, 1fr);
      gap: 20px;
    }}
    .ep-card {{
      background: var(--card-bg-color);
      border: 1px solid var(--border-color);
      border-radius: 12px;
      padding: 20px;
      cursor: pointer;
      transition: all 0.2s ease;
    }}
    .ep-card:hover, .ep-card.active {{
      border-color: var(--primary-color);
      box-shadow: 0 4px 18px rgba(92, 158, 255, 0.15);
    }}
    .ep-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      margin-bottom: 10px;
    }}
    .badge {{
      background: var(--topic-postgresql);
      color: #fff;
      padding: 3px 10px;
      border-radius: 12px;
      font-size: 0.76rem;
      font-weight: 700;
    }}
    .speedup-pill {{
      color: var(--primary-light-blue);
      font-size: 0.76rem;
      font-weight: 600;
    }}
    h3 {{
      margin: 6px 0;
      font-size: 1.12rem;
      color: var(--primary-color);
    }}
    .ep-sub {{
      margin: 0 0 8px 0;
      font-size: 0.83rem;
      color: var(--secondary-color);
      font-weight: 600;
    }}
    .ep-desc {{
      font-size: 0.81rem;
      color: var(--text-color);
      margin-bottom: 12px;
    }}
    .ep-links {{
      display: flex;
      flex-wrap: wrap;
      gap: 7px;
    }}
    button, .btn-link {{
      background: var(--code-bg-color);
      color: var(--text-color);
      border: 1px solid rgba(92, 158, 255, 0.3);
      padding: 6px 11px;
      border-radius: 8px;
      font-family: inherit;
      font-size: 0.76rem;
      font-weight: 600;
      cursor: pointer;
      text-decoration: none;
      display: inline-block;
    }}
    button:hover, button.active, .btn-link:hover {{
      background: var(--gradient-blue-purple);
      color: #fff;
      border-color: var(--primary-color);
    }}
    .blog-pill {{
      background: rgba(177, 88, 252, 0.12);
      border: 1px solid rgba(177, 88, 252, 0.35);
      color: var(--text-color);
      padding: 4px 9px;
      border-radius: 6px;
      font-size: 0.74rem;
      text-decoration: none;
    }}
    .blog-pill:hover {{
      background: var(--secondary-color);
      color: #fff;
    }}
    /* Voiced Hovering Video Companion Widget */
    #hover-video-widget {{
      position: fixed;
      bottom: 24px;
      right: 24px;
      width: 420px;
      background: rgba(22, 27, 39, 0.96);
      border: 2px solid var(--primary-color);
      border-radius: 14px;
      box-shadow: 0 16px 42px rgba(0, 0, 0, 0.75);
      z-index: 9999;
      backdrop-filter: blur(12px);
      overflow: hidden;
      transition: width 0.25s ease;
    }}
    #hover-video-widget.expanded {{ width: 720px; }}
    #hover-video-widget.minimized {{ width: 290px; }}
    #hover-video-widget.minimized .hover-body {{ display: none; }}
    #hover-video-widget.docked {{
      position: static;
      width: 100%;
      margin-bottom: 20px;
    }}
    .hover-header {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 8px 12px;
      background: #0b0d13;
      border-bottom: 1px solid var(--border-color);
      cursor: move;
      user-select: none;
    }}
    .hover-title {{
      font-size: 0.76rem;
      font-weight: 700;
      color: var(--primary-color);
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      max-width: 230px;
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
    footer {{
      display: flex;
      justify-content: space-between;
      align-items: center;
      padding: 14px 36px;
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
  </style>
</head>
<body>
  <header>
    <div>
      <span class="series-badge">SERIES 01 • INTERACTIVE EXPLORER</span>
      <span style="margin-left: 12px; color: var(--text-muted); font-size: 0.86rem;">Hello World, pg_bio! — Voiced Hovering Video &amp; 3D Studio</span>
    </div>
    <div style="display: flex; gap: 8px;">
      <button onclick="toggleHoverDock()">📌 Float / Dock Voiced Video</button>
      <a class="btn-link" href="{BLOG_URL}/pg-bio-videos/" target="_blank">🌐 ideia.me/pg-bio-videos</a>
      <a class="btn-link" href="{REPO_URL}" target="_blank">🐙 github.com/jonatas/pg_bio</a>
    </div>
  </header>

  <main>
    <div>
      <h1>🧬 Hello World, pg_bio! — Biotech &amp; Multiomics Explorer</h1>
      <p style="margin: 0; color: var(--text-muted); font-size: 0.9rem;">
        Scroll and explore the SQL queries, 3D protein pockets, and linked <code>ideia.me</code> blog posts while the <strong>🔊 Voiced Hovering Video Companion</strong> floats in the corner (drag it anywhere, resize it, or click Act I–V to sync the 3D viewer!).
      </p>
    </div>

    <div id="docked-video-slot" style="margin-top: 16px;"></div>

    <!-- Interactive Split Stage: Active Episode Rosetta + SQL on Left | 3Dmol.js Stage on Right -->
    <div class="interactive-split">
      <div class="stage-card">
        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;">
          <span class="badge" id="stage-ep-layer">{first['layer']}</span>
          <div>
            <button id="btn-lens-eng" class="active" onclick="setLens('engineer')">⚙️ Systems Engineer</button>
            <button id="btn-lens-bio" onclick="setLens('biologist')">🧬 Structural Biologist</button>
          </div>
        </div>
        <h3 id="stage-ep-title" style="font-size: 1.25rem; margin-top: 0;">{first['title']}</h3>
        <p id="stage-lens-text" style="font-size: 0.88rem;">
          <strong>⚙️ Systems Engineering Lens:</strong> {first['cs_mechanism']}
        </p>
        <pre id="stage-sql-code">{chr(10).join(first['sql_lines'])}</pre>
        <div style="margin-top: 10px;">
          <strong style="font-size: 0.78rem; color: var(--secondary-color);">📚 Companion Blog Posts on ideia.me:</strong>
          <div id="stage-blog-pills" class="ep-links" style="margin-top: 6px;"></div>
        </div>
      </div>

      <div class="stage-card">
        <h3 style="margin-top: 0;">🔬 Interactive 3D Molecular Stage (3Dmol.js)</h3>
        <div id="hub-3d-viewer"></div>
        <div class="ep-links" style="margin-top: 12px;">
          <button onclick="highlightCore()">1. Highlight Catalytic Zn2+ &amp; Water</button>
          <button onclick="highlightPocket()">2. Show 3D Morton Z-Order Sphere</button>
          <button onclick="resetView()">3. Reset Backbone</button>
        </div>
        <div id="hub-chemistry-note" style="margin-top: 10px; padding: 10px; border-left: 3px solid var(--primary-color); background: rgba(92,158,255,0.07); font-size: 0.82rem;">
          <strong>💡 Synced 3D Stage:</strong> Clicking Act I–V in the voiced hovering video player automatically updates this 3D structure view.
        </div>
      </div>
    </div>

    <h3 style="font-size: 1.3rem; margin-bottom: 12px;">📺 All 4 Episodes in Series 01</h3>
    <div class="ep-grid">
      {"".join(cards_html)}
    </div>
  </main>

  <footer>
    <div>🌐 Blog &amp; Video Library: <a href="{BLOG_URL}/pg-bio-videos/" target="_blank" rel="noopener">{BLOG_URL}/pg-bio-videos/</a></div>
    <div>🧬 <strong>pg_bio</strong> — The Multiomics Operating System in PostgreSQL</div>
    <div>🐙 GitHub Repository: <a href="{REPO_URL}" target="_blank" rel="noopener">{REPO_URL}</a></div>
  </footer>

  <!-- Voiced Hovering Video Companion Widget -->
  <div id="hover-video-widget">
    <div class="hover-header" id="hover-drag-handle">
      <span class="hover-title" id="hover-now-playing">🔊 {first['title']}</span>
      <div style="display: flex; gap: 4px;">
        <button style="padding: 2px 6px; font-size: 0.7rem;" onclick="document.getElementById('hover-video-widget').classList.toggle('expanded')" title="Expand / Shrink">⤢</button>
        <button style="padding: 2px 6px; font-size: 0.7rem;" onclick="toggleHoverDock()" title="Dock / Float">📌</button>
        <button style="padding: 2px 6px; font-size: 0.7rem;" onclick="document.getElementById('hover-video-widget').classList.toggle('minimized')" title="Minimize">—</button>
      </div>
    </div>
    <div class="hover-body">
      <video id="series-video" controls playsinline preload="metadata" data-mp4-src="{first['slug']}/{first['slug']}.mp4">
        <source src="{first['slug']}/{first['slug']}.mp4" type="video/mp4">
      </video>
      <div class="act-chapters">
        <button class="act-btn" onclick="seekAct(0)">Act I: Bio</button>
        <button class="act-btn" onclick="seekAct(1)">Act II: IPC</button>
        <button class="act-btn" onclick="seekAct(2)">Act III: Bridge</button>
        <button class="act-btn" onclick="seekAct(3)">Act IV: SQL</button>
        <button class="act-btn" onclick="seekAct(4)">Act V: Proof</button>
      </div>
    </div>
  </div>

  <script>
    const EPISODES = {episodes_json};
    let activeIdx = 0;
    let activeLens = 'engineer';
    let viewer = null;
    const blobCache = {{}};

    async function loadVideoResilient(video, mp4Url, autoPlay) {{
      video.dataset.mp4Src = mp4Url;
      if (blobCache[mp4Url]) {{
        video.src = blobCache[mp4Url];
        if (autoPlay) video.play().catch(() => {{}});
        return;
      }}
      video.src = mp4Url;
      video.load();
      if (autoPlay) {{
        try {{
          await video.play();
          return;
        }} catch (_) {{}}
      }}
      try {{
        const resp = await fetch(mp4Url, {{ cache: 'reload' }});
        if (!resp.ok) return;
        const blob = await resp.blob();
        const bUrl = URL.createObjectURL(new Blob([blob], {{ type: 'video/mp4' }}));
        blobCache[mp4Url] = bUrl;
        if (video.dataset.mp4Src === mp4Url && (video.error || video.networkState === HTMLMediaElement.NETWORK_NO_SOURCE || autoPlay)) {{
          video.src = bUrl;
          video.load();
          if (autoPlay) video.play().catch(() => {{}});
        }}
      }} catch (_) {{}}
    }}

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

    window.addEventListener('DOMContentLoaded', () => {{
      viewer = $3Dmol.createViewer(document.getElementById('hub-3d-viewer'), {{
        backgroundColor: '#090c12'
      }});
      viewer.addModel(DEMO_PDB, 'pdb');
      highlightCore();
      renderBlogPills(0);

      // Bidirectional Video -> Act & 3Dmol.js Sync + Auto-Advance Next Episode
      const videoEl = document.getElementById('series-video');
      videoEl.addEventListener('error', () => {{
        loadVideoResilient(videoEl, videoEl.dataset.mp4Src || (EPISODES[activeIdx].slug + '/' + EPISODES[activeIdx].slug + '.mp4'), false);
      }}, true);
      let lastSyncedAct = -1;
      videoEl.addEventListener('timeupdate', () => {{
        if (!videoEl.duration) return;
        const curAct = Math.min(4, Math.floor((videoEl.currentTime / videoEl.duration) * 5));
        if (curAct !== lastSyncedAct) {{
          lastSyncedAct = curAct;
          document.querySelectorAll('.act-btn').forEach((btn, i) => btn.classList.toggle('active', i === curAct));
          if (curAct === 0) {{ setLens('biologist'); highlightCore(); }}
          else if (curAct === 2) {{ setLens('engineer'); highlightPocket(); }}
          else if (curAct === 4) {{ resetView(); }}
        }}
      }});
      videoEl.addEventListener('ended', () => {{
        if (activeIdx + 1 < EPISODES.length) selectEpisode(activeIdx + 1);
      }});

      // Draggable hovering video widget
      const widget = document.getElementById('hover-video-widget');
      const handle = document.getElementById('hover-drag-handle');
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

    function renderBlogPills(idx) {{
      const container = document.getElementById('stage-blog-pills');
      container.innerHTML = EPISODES[idx].blogs
        .map(b => `<a class="blog-pill" href="${{b.url}}" target="_blank">📖 ${{b.title}}</a>`)
        .join('');
    }}

    function setLens(mode) {{
      activeLens = mode;
      const ep = EPISODES[activeIdx];
      document.getElementById('stage-lens-text').innerHTML = mode === 'engineer'
        ? `<strong>⚙️ Systems Engineering Lens:</strong> ${{ep.cs_mechanism}}`
        : `<strong>🧬 Structural Biology Lens:</strong> ${{ep.bio_problem}}`;
      document.getElementById('btn-lens-eng').classList.toggle('active', mode === 'engineer');
      document.getElementById('btn-lens-bio').classList.toggle('active', mode === 'biologist');
    }}

    function selectEpisode(idx) {{
      activeIdx = idx;
      const ep = EPISODES[idx];
      const video = document.getElementById('series-video');
      const mp4Url = ep.slug + '/' + ep.slug + '.mp4';
      loadVideoResilient(video, mp4Url, true);
      document.getElementById('hover-now-playing').textContent = '🔊 ' + ep.title;
      document.getElementById('stage-ep-layer').textContent = ep.layer;
      document.getElementById('stage-ep-title').textContent = ep.title;
      document.getElementById('stage-sql-code').textContent = ep.sql_lines.join('\\n');
      setLens(activeLens);
      renderBlogPills(idx);
      document.querySelectorAll('.ep-card').forEach((el, i) => el.classList.toggle('active', i === idx));
      if (idx === 1) highlightPocket();
      else highlightCore();
    }}

    async function seekAct(actIdx) {{
      const video = document.getElementById('series-video');
      if (!video.duration || video.error) {{
        await loadVideoResilient(video, video.dataset.mp4Src || (EPISODES[activeIdx].slug + '/' + EPISODES[activeIdx].slug + '.mp4'), false);
      }}
      if (video.duration) {{
        video.currentTime = (video.duration / 5) * actIdx;
      }}
      video.play().catch(() => {{}});
      if (actIdx === 0) highlightCore();
      if (actIdx === 2) highlightPocket();
      if (actIdx === 4) resetView();
    }}

    function toggleHoverDock() {{
      const widget = document.getElementById('hover-video-widget');
      const slot = document.getElementById('docked-video-slot');
      widget.classList.toggle('docked');
      if (widget.classList.contains('docked')) slot.appendChild(widget);
      else document.body.appendChild(widget);
    }}

    function highlightCore() {{
      if (!viewer) return;
      viewer.removeAllLabels(); viewer.removeAllShapes();
      viewer.setStyle({{}}, {{ stick: {{ radius: 0.18, colorscheme: 'cyanCarbon' }} }});
      viewer.setStyle({{ resn: 'ZN' }}, {{ sphere: {{ scale: 0.65, color: '#FFC107' }} }});
      viewer.setStyle({{ resn: 'HOH' }}, {{ sphere: {{ scale: 0.45, color: '#ff5a9d' }} }});
      viewer.addLabel('Zn2+ Cofactor', {{ position: {{ x: 3.95, y: -1.65, z: -2.15 }}, backgroundColor: '#FFC107', fontColor: '#000' }});
      viewer.addLabel('Nucleophilic HOH (2.05 Å)', {{ position: {{ x: 3.75, y: -2.95, z: -1.15 }}, backgroundColor: '#ff5a9d', fontColor: '#fff' }});
      viewer.zoomTo(); viewer.render();
      document.getElementById('hub-chemistry-note').innerHTML =
        "<strong>🔬 Catalytic Triad &amp; Metal Coordination:</strong> Zn<sup>2+</sup> polarizes nucleophilic water (2.05 &Aring;) while Glu106 abstracts a proton.";
    }}

    function highlightPocket() {{
      if (!viewer) return;
      viewer.removeAllLabels(); viewer.removeAllShapes();
      viewer.setStyle({{}}, {{ stick: {{ radius: 0.22, colorscheme: 'magentaCarbon' }}, sphere: {{ scale: 0.28 }} }});
      viewer.addSphere({{ center: {{ x: 3.95, y: -1.65, z: -2.15 }}, radius: 3.5, color: '#5c9eff', alpha: 0.22 }});
      viewer.addLabel('Morton Z-Order 3.5 Å Shell', {{ position: {{ x: 1.5, y: 1.2, z: 0.0 }}, backgroundColor: '#336791', fontColor: '#fff' }});
      viewer.render();
      document.getElementById('hub-chemistry-note').innerHTML =
        "<strong>⚙️ 64-Bit Morton Z-Order Sphere:</strong> <code>z_order_encode(x,y,z)</code> maps 3D coordinates onto a 1D space-filling curve for O(log N) B-Tree scans.";
    }}

    function resetView() {{
      if (!viewer) return;
      viewer.removeAllLabels(); viewer.removeAllShapes();
      viewer.setStyle({{}}, {{ stick: {{ radius: 0.16, colorscheme: 'Jmol' }} }});
      viewer.setStyle({{ resn: 'ZN' }}, {{ sphere: {{ scale: 0.55, color: '#5c9eff' }} }});
      viewer.zoomTo(); viewer.render();
    }}
  </script>
</body>
</html>
"""


def build_series_readme(episodes: List[Dict[str, Any]]) -> str:
    rows = []
    for ep in episodes:
        slug = ep["slug"]
        blogs_md = " • ".join(f"[{b['title']}]({b['url']})" for b in ep["blogs"])
        rows.append(
            f"| **Ep 0{ep['ep_num']}** | **{ep['title']}**<br>*{ep['subtitle']}* | `{ep['layer']}` | "
            f"[MP4]({slug}/{slug}.mp4) • [3D Studio HTML]({slug}/studio_deck.html) • [Talk Track]({slug}/talk_track.md) • [SQL]({slug}/live_demo.sql) • [Doc]({ep['docs_ref']}) | {blogs_md} |"
        )
    table_md = "\n".join(rows)
    return f"""# 🧬 Series 01: Hello World, `pg_bio`!
### Intro to Biotech & Basic Multiomics in PostgreSQL for Software Engineers

* **Interactive Series Explorer (with Voiced Hovering Video & 3Dmol.js):** [`index.html`](index.html)
* **Live Video Library on `ideia.me`:** [{BLOG_URL}/pg-bio-videos/]({BLOG_URL}/pg-bio-videos/)
* **Blog & 3D Discoveries:** [{BLOG_URL}]({BLOG_URL})
* **GitHub Repository:** [{REPO_URL}]({REPO_URL})

---

## 📺 Episode Lineup & Companion Material

| Episode | Title & Focus | `pg_bio` Multiomics Layers | Video, 3D Studio & SQL Assets | Companion `ideia.me` Blog Posts |
| :--- | :--- | :--- | :--- | :--- |
{table_md}
"""


def build_root_videos_readme() -> str:
    return f"""# 🎬 `pg_bio` Video Library & Interactive 3D Learning Hub

Welcome to the `pg_bio` video and interactive 3D tutoring directory! Every series and episode here pairs:
1. **Narrated 1080p MP4 Videos (`.mp4` + `voiceover.mp3` + `captions.srt`)** styled with the [`ideia.me`]({BLOG_URL}) dark theme.
2. **Interactive 3Dmol.js + Voiced Hovering Video Explorer Pages (`studio_deck.html` & `index.html`)** that let you scroll through SQL queries, toggle **⚙️ Systems Engineer $\\leftrightarrow$ 🧬 Structural Biologist** mental models, and rotate 3D protein pockets while a floating Picture-in-Picture voiced video guides you.
3. **Runnable `psql` Demo Scripts (`live_demo.sql`), Teleprompter Talk Tracks (`talk_track.md`), Marp Slides (`slides.marp.md`), and Manim Scenes (`manim_scene.py`).**

---

## 🗂️ Available Series & Episodes

### 1. [Series 01: Hello World, `pg_bio`! — Intro to Biotech & Multiomics in PostgreSQL](series_01_hello_world_pg_bio/README.md)
* 🚀 **Interactive Explorer Hub:** [`series_01_hello_world_pg_bio/index.html`](series_01_hello_world_pg_bio/index.html)
* 🌐 **Web Library on `ideia.me`:** [{BLOG_URL}/pg-bio-videos/]({BLOG_URL}/pg-bio-videos/)
* **Episodes:**
  * **Ep 01:** [The Central Dogma in SQL (`int4range` GiST & Sequence Physics)](series_01_hello_world_pg_bio/ep01_central_dogma_and_genomics/studio_deck.html) — [MP4](series_01_hello_world_pg_bio/ep01_central_dogma_and_genomics/ep01_central_dogma_and_genomics.mp4) | [SQL](series_01_hello_world_pg_bio/ep01_central_dogma_and_genomics/live_demo.sql) | [Doc](../docs/07_genomics_io.md)
  * **Ep 02:** [From 1D Strings to 3D Atoms with Morton Curves (`z_order_encode` & SVD)](series_01_hello_world_pg_bio/ep02_3d_pockets_and_morton_curves/studio_deck.html) — [MP4](series_01_hello_world_pg_bio/ep02_3d_pockets_and_morton_curves/ep02_3d_pockets_and_morton_curves.mp4) | [SQL](series_01_hello_world_pg_bio/ep02_3d_pockets_and_morton_curves/live_demo.sql) | [Doc](../docs/03_spatial_indexing.md)
  * **Ep 03:** [Transcriptomics & CSR Sparse Attention Graphs (`SparseAttentionMap`)](series_01_hello_world_pg_bio/ep03_scrna_csr_and_attention_maps/studio_deck.html) — [MP4](series_01_hello_world_pg_bio/ep03_scrna_csr_and_attention_maps/ep03_scrna_csr_and_attention_maps.mp4) | [SQL](series_01_hello_world_pg_bio/ep03_scrna_csr_and_attention_maps/live_demo.sql) | [Doc](../docs/12_mmpbsa_free_energy_auditing.md)
  * **Ep 04:** [The Multiomics Holy Grail — Vectors & Live UniProt (`halfvec` HNSW)](series_01_hello_world_pg_bio/ep04_joint_vector_multiomics/studio_deck.html) — [MP4](series_01_hello_world_pg_bio/ep04_joint_vector_multiomics/ep04_joint_vector_multiomics.mp4) | [SQL](series_01_hello_world_pg_bio/ep04_joint_vector_multiomics/live_demo.sql) | [Doc](../docs/01_vector_homology.md)

### 2. [Standalone Feature Deep-Dive: From B-Trees & GiST to 3D Catalytic Pockets](ep01_biotech_for_sql_engineers/studio_deck.html)
* [MP4 Video](ep01_biotech_for_sql_engineers/ep01_biotech_for_sql_engineers.mp4) • [Interactive Studio Deck](ep01_biotech_for_sql_engineers/studio_deck.html) • [Talk Track](ep01_biotech_for_sql_engineers/talk_track.md) • [Live SQL](ep01_biotech_for_sql_engineers/live_demo.sql)

### 3. [Standalone Feature Deep-Dive: The Protein as a Modular State Machine (`Dock -> Cut -> Transform`)](ep02_modular_protein_state_machine/studio_deck.html)
* [MP4 Video](ep02_modular_protein_state_machine/ep02_modular_protein_state_machine.mp4) • [Interactive Studio Deck](ep02_modular_protein_state_machine/studio_deck.html) • [Talk Track](ep02_modular_protein_state_machine/talk_track.md) • [Live SQL](ep02_modular_protein_state_machine/live_demo.sql)
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Stage Series 01 and sync to ideia.me")
    parser.add_argument("--force-mp4", action="store_true", help="Re-render MP4 even if it already exists")
    args = parser.parse_args()

    scripts_dir = Path(__file__).resolve().parent
    repo_root = scripts_dir.parents[3]
    videos_root = repo_root / "videos"
    series_dir = videos_root / "series_01_hello_world_pg_bio"
    series_dir.mkdir(parents=True, exist_ok=True)

    render_script = scripts_dir / "render_feature_mp4.py"
    import stage_feature_video as sfv

    for ep in EPISODES:
        slug = ep["slug"]
        ep_dir = series_dir / slug
        ep_dir.mkdir(parents=True, exist_ok=True)

        cfg_path = ep_dir / "episode_config.json"
        cfg_path.write_text(json.dumps(ep, indent=2), encoding="utf-8")

        ns = sfv.argparse.Namespace(
            slug=slug,
            title=ep["title"],
            layer=ep["layer"],
            bio_problem=ep["bio_problem"],
            cs_mechanism=ep["cs_mechanism"],
            speedup=ep["speedup"],
        )
        (ep_dir / "talk_track.md").write_text(sfv.build_talk_track(ns), encoding="utf-8")
        (ep_dir / "live_demo.sql").write_text(sfv.build_live_demo_sql(ns), encoding="utf-8")
        (ep_dir / "slides.marp.md").write_text(sfv.build_marp_slides(ns), encoding="utf-8")
        (ep_dir / "manim_scene.py").write_text(sfv.build_manim_scene(ns), encoding="utf-8")

        custom_sql_block = "\n".join(ep["sql_lines"]).replace("<", "&lt;").replace(">", "&gt;")
        default_sql_block = """-- Zero-IPC execution inside PostgreSQL 18 shared buffers:
SELECT design_id,
       theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') AS rmsd,
       bio_angle(water_o, substrate_c, substrate_o) AS attack_angle
FROM metalloprotease_designs
WHERE has_coordinated_water_shell(active_site_coords, 'ZN', 2.4) = TRUE
  AND theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') &lt;= 0.90;"""
        deck_html = sfv.build_studio_deck_html(ns).replace(default_sql_block, custom_sql_block)
        (ep_dir / "studio_deck.html").write_text(deck_html, encoding="utf-8")

        mp4_path = ep_dir / f"{slug}.mp4"
        if args.force_mp4 or not mp4_path.exists():
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
                    str(render_script),
                    "--out",
                    str(mp4_path),
                    "--config",
                    str(cfg_path),
                ],
                check=True,
            )

    (series_dir / "index.html").write_text(build_series_hub_html(EPISODES), encoding="utf-8")
    (series_dir / "README.md").write_text(build_series_readme(EPISODES), encoding="utf-8")
    (videos_root / "README.md").write_text(build_root_videos_readme(), encoding="utf-8")

    def atomic_sync_tree(src: Path, dst: Path) -> None:
        dst.mkdir(parents=True, exist_ok=True)
        for p in src.rglob("*"):
            if p.name == ".DS_Store":
                continue
            rel = p.relative_to(src)
            target = dst / rel
            if p.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                tmp = target.with_name(f".{target.name}.tmp")
                shutil.copy2(p, tmp)
                tmp.replace(target)

    # Sync all videos/ assets into ../ideia.me/assets/videos/pg_bio/ (and _site if running)
    ideia_repo = repo_root.parent / "ideia.me"
    if ideia_repo.exists():
        dest_videos = ideia_repo / "assets" / "videos" / "pg_bio"
        atomic_sync_tree(videos_root, dest_videos)
        print(f"Synced pg_bio video & 3D studio assets to: {dest_videos}")
        site_dir = ideia_repo / "_site"
        if site_dir.exists():
            site_dest = site_dir / "assets" / "videos" / "pg_bio"
            atomic_sync_tree(videos_root, site_dest)
            print(f"Synced live Jekyll _site assets to:        {site_dest}")

    print(f"Updated Series 01 Explorer Hub: {series_dir / 'index.html'}")
    print(f"Updated Root Videos Catalog:    {videos_root / 'README.md'}")


if __name__ == "__main__":
    main()
