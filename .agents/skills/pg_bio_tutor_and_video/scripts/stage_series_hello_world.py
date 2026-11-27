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
        "title": "Ep 1: How Life Compiles DNA & Survives a pH 11 Soda Lake",
        "display_title": "Ep 1: How Life Compiles DNA & Survives a pH 11 Soda Lake",
        "subtitle": "The Central Dogma Compiler, 1D Safe-Harbor Intervals & Extremophile Charge Armor",
        "layer": "Episode 01 • Central Dogma & Surface Charge Heuristics",
        "bio_problem": (
            "If you drop an ordinary human or E. coli enzyme into Lake Magadi—a pH 10.5 volcanic soda lake saturated "
            "with 3.5 Molar salt—the water is ripped off its surface and the protein clumps into dead sludge in seconds. "
            "Yet haloalkaliphilic archaea thrive there using the exact same 20 amino acid building blocks. How?"
        ),
        "cs_mechanism": (
            "By coating the outer protein shell with negatively charged Aspartate and Glutamate carboxylates (lowering "
            "the Isoelectric Point pI to ~4.3 via Henderson-Hasselbalch bisection) while keeping the catalytic core intact—and "
            "moving our compute directly to the database pages where millions of genomes live."
        ),
        "speedup": "0.44 Å 3D core match across 68% sequence divergence (pI 4.38 vs 6.40)",
        "organism_pair": {
            "stage_header": "3D ORGANISM MOTION: NORMAL HALOPHILE vs. LAKE MAGADI SODA-LAKE ARCHAEON",
            "org1_label": "Bait: A0ABD5UZ94 (Standard Hydrolase · 312 aa)",
            "org2_label": "Discovery: D3SWG7 (Natrialba magadii · Lake Magadi pH 10.5)",
            "pdb1": "AF-A0ABD5UZ94-F1-model_v4_aligned.pdb",
            "pdb2": "AF-D3SWG7-F1-model_v4.pdb",
            "pocket_indices": [25, 34, 42, 55, 70, 82],
            "bridge_badge": "0.44 Å Identical Core · Acidic Armor",
            "metric_line1": "3D Core RMSD: 0.44 Å (TM = 0.996)  |  Isoelectric Point pI: 4.38 (Acidic Brine Shield)",
            "metric_line2": "Orange = Negatively Charged Asp/Glu Hydration Cage  |  Emerald = Conserved Catalytic Triad",
        },
        "trap_cards": [
            ("The Naive Assumption", "If two proteins do the exact same job, their amino acid letter sequences should look almost identical."),
            ("The Soda-Lake Paradox", "In 3.5 M salt, hydrophobic surfaces stick together and precipitate. To survive, the organism must mutate 65% of its outer letters!"),
            ("Why Exporting Data Hits a Wall", "In this intro episode, notice one systems truth: pulling millions of genomes over a socket into Python wastes 99% of time moving bytes."),
            ("The Heuristic Solution", "Bring two simple physics heuristics—1D interval trees and Henderson-Hasselbalch charge bisection—straight to the data."),
        ],
        "docs_ref": "../../docs/07_genomics_io.md",
        "notebook_ref": "../../experiments/06_genomics_io_tutorial.ipynb",
        "blogs": [
            {"title": "Introducing pg_bio", "url": f"{BLOG_URL}/pg-bio-intro"},
            {"title": "Genomic Safe-Harbor Architecture", "url": f"{BLOG_URL}/teaflon-synthetic-biology-part-8-genomic-architecture"},
            {"title": "Multi-Pocket PoseBusters & Soda-Lake D3SWG7", "url": f"{BLOG_URL}/pocketvina-posebusters-enzref-pslim-tethers-pg-bio"},
        ],
        "bullets": [
            "• Why does a soda-lake enzyme rewrite 65% of its amino acid letters?",
            "• How does negative surface charge (pI < 4.5) lock onto water in 3.5 M salt?",
            "• How do 1D interval trees find a safe genomic home without hitting cancer genes?",
            "• Why moving the math inside the database lets us test hypotheses in milliseconds",
        ],
        "rosetta_mappings": [
            [
                "Heuristic 1: The Henderson-Hasselbalch Charge Balance (pI)",
                "Every Aspartate (D) and Glutamate (E) carries a negative carboxylate charge at neutral/alkaline pH.",
                "Bisecting the pH where net charge hits zero immediately spots salt-loving brine armor (pI ~ 4.3).",
            ],
            [
                "Heuristic 2: The Kyte-Doolittle Hydropathy Balance (GRAVY)",
                "Oily residues repel water and fold inward; charged residues pull water into a protective hydration shell.",
                "A negative GRAVY score (-0.45) proves the protein stays soluble instead of clumping into sludge.",
            ],
            [
                "Heuristic 3: 1D Spatial Interval Trees (Genomic Safe Harbors)",
                "A chromosome is a 1D number line. To insert a gene safely, it must sit near a promoter but never overlap an oncogene.",
                "Interval overlap (&&) and distance (<->) operators solve safe-harbor placement in one step.",
            ],
        ],
        "sql_lines": [
            "SELECT uniprot_id, organism,",
            "       round(isoelectric_point(sequence)::numeric, 2) AS pI_charge,",
            "       round(gravy_score(sequence)::numeric, 3)       AS hydropathy,",
            "       round(sequence_entropy(sequence)::numeric, 3)  AS complexity",
            "FROM proteins",
            "WHERE uniprot_id IN ('A0ABD5UZ94', 'D3SWG7');",
        ],
        "sql_results": [
            "┌────────────┬──────────────────────────┬───────────┬────────────┬────────────┐",
            "│ uniprot_id │         organism         │ pI_charge │ hydropathy │ complexity │",
            "├────────────┼──────────────────────────┼───────────┼────────────┼────────────┤",
            "│ A0ABD5UZ94 │ Halophilic reference     │      4.52 │     -0.318 │      4.012 │",
            "│ D3SWG7     │ Natrialba magadii (pH11) │      4.38 │     -0.412 │      3.984 │",
            "└────────────┴──────────────────────────┴───────────┴────────────┴────────────┘",
            "Heuristic Verified: Acidic pI 4.38 + Negative GRAVY keeps D3SWG7 soluble at 3.5 M NaCl!",
        ],
        "narrations": [
            (
                "Imagine dropping an ordinary enzyme from your body into Lake Magadi in Kenya, a volcanic soda lake with a pH of 10.5 "
                "that is nearly solid with salt. Within seconds, the salt strips the water molecules off the protein's surface, and it collapses "
                "into useless sludge. Yet microbes called haloalkaliphiles swim happily in that brine, using the exact same 20 amino acids we do. How?"
            ),
            (
                "Look at the two rotating organisms on screen. On the left in cyan is a reference hydrolase; on the right in magenta is D3SWG7 "
                "from Natrialba magadii. Almost two thirds of their outer letters have mutated! In this first introductory episode, notice one foundational "
                "reason we do this inside a database: if you pull millions of sequences over a network socket into external scripts just to check their physics, "
                "you wait hours. When you bring the math straight to the data, you can ask questions at the speed of thought."
            ),
            (
                "And the heuristic that explains how D3SWG7 survives is remarkably intuitive. Oily amino acids hate water and fold inward to form a tight structural core. "
                "Meanwhile, two specific amino acids, Aspartate and Glutamate, carry negatively charged carboxylate groups that grab onto water molecules like tiny magnets. "
                "By running a 100-step bisection on their pKa curves, we find the protein's Isoelectric Point, or pI."
            ),
            (
                "Watch the rotating 3D structures as the orange highlights light up. Those orange spheres are the dense shield of negatively charged Aspartate and Glutamate residues "
                "covering the exterior of D3SWG7, dropping its pI to 4.38. That negative armor holds a pressurized shell of water around the enzyme even in 3.5 Molar salt, "
                "while the glowing emerald catalytic triad inside remains aligned to 0.44 Angstroms!"
            ),
            (
                "Once you understand how surface charge and hydropathy protect a protein's 3D core, expressing that heuristic takes just three lines of SQL. "
                "In Episode 2, we will zoom inside that glowing emerald pocket and see how a 19th-century space-filling curve lets us search 3D atomic space in 1D."
            ),
        ],
    },
    {
        "ep_num": 2,
        "slug": "ep02_3d_pockets_and_morton_curves",
        "visual_mode": "morton_3d_pocket",
        "title": "Ep 2: How a 1D Line Slices 3D Atomic Space (Morton Curves & Kabsch SVD)",
        "display_title": "Ep 2: How a 1D Line Slices 3D Atomic Space (Morton & Kabsch)",
        "subtitle": "Turning 3D Protein Clouds into 1D Space-Filling Keys & Aligning Active Sites",
        "layer": "Episode 02 • 3D Morton Space-Filling Curves & Kabsch Superposition",
        "bio_problem": (
            "In a linear protein sequence, the three amino acids that form a catalytic triad might sit at positions 106, "
            "132, and 224—hundreds of letters apart. Only when the chain folds in 3D space do they snap together within "
            "3.5 Angstroms. How do you find a 3D sphere of atoms when standard indexes only understand 1D sorted lists?"
        ),
        "cs_mechanism": (
            "By interleaving the binary bits of the 3D (X, Y, Z) coordinates into a single 64-bit Morton Z-Order integer "
            "so nearby 3D octants stay consecutive on a 1D line, and then solving a 3x3 Jacobi Singular Value Decomposition "
            "(Kabsch algorithm) to rotate two organisms into exact 3D alignment."
        ),
        "speedup": "74x fewer atoms scanned via Layered 3D Morton Octree + Exact Kabsch RMSD",
        "organism_pair": {
            "stage_header": "3D KABSCH SUPERPOSITION: FLUOROACETATE DEHALOGENASE vs. EXTREMOPHILE TWIN",
            "org1_label": "Bait: Q1JU72 Reference Dehalogenase (C-F Bond Cutter)",
            "org2_label": "Discovery: A0A031LSV2 (Kabsch-Aligned Catalytic Hexad)",
            "pdb1": "AF-A0A031LSV2-F1-model_v4_reference.pdb",
            "pdb2": "AF-A0A031LSV2-aligned-to-Q1JU72.pdb",
            "pocket_indices": [18, 29, 44, 58, 72, 86],
            "bridge_badge": "Sub-Angstrom 3D Kabsch Alignment",
            "metric_line1": "21-Level Morton Octree Pruning: 74.05× fewer atoms scanned (100.0% exact recall)",
            "metric_line2": "3×3 Jacobi SVD rotates both organisms so their catalytic triads lock together in 3D!",
        },
        "trap_cards": [
            ("The 1D vs. 3D Blind Spot", "Looking at a 1D sequence string gives you zero clue which atoms are touching in 3D space."),
            ("The Naive Bounding-Box Trap", "If you sort by X first, then Y, then Z, you have to scan huge slabs of empty space across the whole protein."),
            ("The Boundary Crossing Glitch", "A simple 1D Morton interval [min_z..max_z] jumps across octree cube boundaries and sweeps in thousands of unrelated atoms."),
            ("The Layered Octree Eureka", "Subdivide the 3D cube recursively: skip octants outside the sphere, accept octants fully inside, and only check the thin shell!"),
        ],
        "docs_ref": "../../docs/03_spatial_indexing.md",
        "notebook_ref": "../../experiments/08_pdb_parsing.ipynb",
        "blogs": [
            {"title": "Indexing 3D Space Like a Video Game", "url": f"{BLOG_URL}/indexing-3d-space-like-a-video-game"},
            {"title": "Layered Centroids & 3D Morton Octrees", "url": f"{BLOG_URL}/layered-centroids-register-fastscan-morton-octree-ecfp4-rangefinder-pg-bio"},
            {"title": "Store & Query 3D Protein Structures in SQL", "url": f"{BLOG_URL}/store-query-3d-protein-structures-postgresql"},
        ],
        "bullets": [
            "• How do residues 100 letters apart form a 3.5 Å chemical scissor?",
            "• How does bit-interleaving (X, Y, Z) fold 3D space onto a 1D number line?",
            "• Why does a 21-level Octree scan 74× fewer atoms than a single range?",
            "• How does a 3×3 SVD matrix rotate two organisms into exact alignment?",
        ],
        "rosetta_mappings": [
            [
                "Step 1: Bit-Interleaving 3D Coordinates (Morton Z-Order Curve)",
                "Take the 21 binary bits of X, Y, and Z and shuffle them like a deck of cards: Z20 Y20 X20 ... Z0 Y0 X0.",
                "Any two atoms that sit in the same 3D cube now share the exact same high-order prefix on a 1D number line!",
            ],
            [
                "Step 2: Layered Octree Sphere Decomposition",
                "Walk down the 3D octree: if an octant's closest corner is farther than R, prune it; if its farthest corner is inside R, keep all of it.",
                "Cuts candidate atoms by 74× while guaranteeing 100% mathematical recall without taking square roots on interior atoms.",
            ],
            [
                "Step 3: Kabsch 3×3 Singular Value Decomposition (SVD)",
                "Two proteins in different PDB files face random directions in 3D space. We compute their 3×3 cross-covariance matrix H = P^T Q.",
                "Solving U Σ V^T = H gives the exact optimal 3D rotation matrix R = V U^T that locks the two active sites together!",
            ],
        ],
        "sql_lines": [
            "SELECT a.atom_id, (a.coord).name,",
            "       kabsch_rmsd(bait_pocket_coords, candidate_pocket_coords) AS pocket_rmsd_a",
            "FROM z_order_layered_sphere_ranges(10.0, 0.0, 25.0, 5.0, 64) r",
            "JOIN teaflon_atoms a ON a.z_curve_index BETWEEN r.min_z AND r.max_z",
            "WHERE r.fully_inside",
            "   OR ((a.coord).x - 10.0)^2 + ((a.coord).y)^2 + ((a.coord).z - 25.0)^2 <= 25.0;",
        ],
        "sql_results": [
            "┌─────────┬──────────────┬───────────────┬──────────────────────┐",
            "│ atom_id │ residue_name │ pocket_rmsd_a │ octree_interior_skip │",
            "├─────────┼──────────────┼───────────────┼──────────────────────┤",
            "│     106 │ ASP_OD1      │         0.412 │ t (sqrt skipped!)    │",
            "│     132 │ HIS_NE2      │         0.412 │ t (sqrt skipped!)    │",
            "└─────────┴──────────────┴───────────────┴──────────────────────┘",
            "Heuristic Result: 407 atoms checked instead of 30,181 (74.05x pruning, 100% exact)",
        ],
        "narrations": [
            (
                "Here is a puzzle at the heart of structural biology. In an enzyme's 1D sequence, the three amino acids that act as a chemical "
                "scissor might sit at positions 106, 132, and 224. In a text string, they look completely unrelated. Yet when the protein folds "
                "in 3D space, those three residues snap together within 3.5 Angstroms. How can we search 3D space when our fastest indexes only understand a 1D line?"
            ),
            (
                "In 1966, Guy Macdonald Morton proposed a deceptively simple trick. Take the binary bits of your X, Y, and Z coordinates, "
                "and interleave them—one bit from X, one from Y, one from Z—into a single 64-bit integer. Visually, this traces a 3D fractal zig-zag "
                "called a Z-Order curve. Any atoms that live inside the same 3D sub-cube now share the exact same bit prefix on a 1D number line!"
            ),
            (
                "But there is a catch: whenever a 3D sphere crosses the seam between two large cubes, a single 1D interval sweeps across thousands "
                "of false-positive atoms. We fix this with a Layered Octree heuristic: we recursively split the cube into eight octants, throw away octants "
                "outside the sphere, mark octants completely inside the sphere so we don't even need to compute a square root, and coalesce the remaining intervals."
            ),
            (
                "Once we slice out the 3D pocket from two different organisms—like the cyan reference dehalogenase and the magenta extremophile twin rotating on screen—"
                "they are initially facing random directions in space. By taking their 3-by-3 cross-covariance matrix and running a Jacobi Singular Value Decomposition, "
                "the Kabsch algorithm finds the exact 3D rotation that locks their catalytic triads together to 0.41 Angstroms!"
            ),
            (
                "With bit-interleaving, octree pruning, and 3-by-3 SVD rotations working together, we scan 74 times fewer atoms with 100 percent exact recall. "
                "Next, in Episode 3, we will see what happens when an inhibitor binds far away from the active site and sends a mechanical ripple across the protein."
            ),
        ],
    },
    {
        "ep_num": 3,
        "slug": "ep03_scrna_csr_and_attention_maps",
        "visual_mode": "csr_sparse_graph",
        "title": "Ep 3: Hidden Whispers in Proteins — Allosteric Wires & Sparse Graphs",
        "display_title": "Ep 3: Hidden Whispers in Proteins — Allostery & Sparse Graphs",
        "subtitle": "How Distant Triterpenoids Lock Human 5-Lipoxygenase Over 500 Nanoseconds",
        "layer": "Episode 03 • Allosteric Energy Networks & Thermodynamic Closure",
        "bio_problem": (
            "Most drugs work like a cork shoved directly into an enzyme's active site. But some of nature's most potent "
            "anti-inflammatory molecules—like Boswellic acid (AKBA) and alpha-Amyrin from Amazonian resins—bind at a distant "
            "crevice on human 5-Lipoxygenase (5-LOX) and somehow freeze the catalytic iron center from afar. How does that whisper travel?"
        ),
        "cs_mechanism": (
            "By auditing per-residue free-energy decompositions across 500-nanosecond Molecular Dynamics trajectories for "
            "double-entry thermodynamic closure, and compressing only the strong energetic couplings into a Compressed Sparse Row (CSR) "
            "graph to trace cryptic allosteric anchors outside the known pocket."
        ),
        "speedup": "Uncovers cryptic anchors Ile126, Val29 & Leu66 across 121,860 MD residue frames",
        "organism_pair": {
            "stage_header": "500 ns MOLECULAR DYNAMICS MOTION: HUMAN 5-LOX + AKBA vs. ALPHA-AMYRIN",
            "org1_label": "5-LOX + AKBA Boswellic Acid (Ionic Allosteric Lock)",
            "org2_label": "5-LOX + α-Amyrin Triterpenoid (Hydrophobic Wedge)",
            "pdb1": "5LOX-akba-af7-500ns.pdb",
            "pdb2": "5LOX-alpha-amyrin-500ns.pdb",
            "pocket_indices": [18, 28, 34, 42, 55, 68],
            "bridge_badge": "Cryptic Allosteric Anchors (Ile126 · Val29)",
            "metric_line1": "9 × 500 ns MD Trajectories  |  Thermodynamic Double-Entry Closure Audited",
            "metric_line2": "Comparing AKBA (polar/ionic salt bridge at Arg68/Glu134) vs α-Amyrin (hydrophobic wedge at Ile126)",
        },
        "trap_cards": [
            ("The Rigid-Lock Illusion", "Looking at a static X-ray crystal structure makes a protein look like a solid rock. In reality, it breathes and flexes billions of times a second."),
            ("The 90%-Zero Matrix Trap", "In a 676-residue enzyme like 5-LOX, almost all pairs of distant amino acids have zero direct interaction energy. Storing a dense matrix buries the signal in noise."),
            ("The Missing Cavity Energy Trap", "In Poisson-Boltzmann free-energy calculations (MM/PBSA), the non-polar cavity term (+22 kcal/mol) is often omitted from per-residue tables, skewing rankings!"),
            ("The CSR Graph Diff Eureka", "First audit the thermodynamic ledger for exact closure, then store only significant energy edges in a sparse CSR graph and subtract the known pocket!"),
        ],
        "docs_ref": "../../docs/12_mmpbsa_free_energy_auditing.md",
        "notebook_ref": "../../experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb",
        "blogs": [
            {"title": "Auditing 500 ns 5-LOX Triterpenoid MD in SQL", "url": f"{BLOG_URL}/5lox-triterpenoids-mmpbsa-free-energy-audit-pg-bio"},
            {"title": "Analyzing Hi-C 3D Genomic Contact Maps in SQL", "url": f"{BLOG_URL}/analyzing-hic-3d-genomic-contact-maps-sql"},
        ],
        "bullets": [
            "• How does a molecule bound 20 Å away shut down an enzyme's active site?",
            "• Why do AKBA and α-Amyrin use opposite forces (ionic vs hydrophobic) in the same cleft?",
            "• How does a double-entry thermodynamic audit catch missing solvation penalties?",
            "• How does a sparse graph diff reveal hidden anchor residues Ile126 and Val29?",
        ],
        "rosetta_mappings": [
            [
                "Step 1: Double-Entry Thermodynamic Closure Audit",
                "Just like an accounting ledger, the sum of per-residue gas-phase (MM) and polar solvation energies must equal the global system total.",
                "Immediately exposes the +20.7 to +26.0 kcal/mol Poisson-Boltzmann non-polar cavity discrepancy before you trust a false energy ranking.",
            ],
            [
                "Step 2: Compressed Sparse Row (CSR) Energy Graph",
                "Discard the 95% of residues that have near-zero coupling with the ligand and pack only stabilizing contacts (≤ -0.15 kcal/mol) into CSR arrays.",
                "Turns a noisy 676-residue trajectory into a razor-sharp allosteric wiring diagram.",
            ],
            [
                "Step 3: Cryptic Hotspot Set-Difference & Mechanism Index",
                "Subtract the 9 pre-selected docking-box residues from the CSR neighborhood and compare polar/ionic vs. hydrophobic stabilization.",
                "Discovers unannotated cryptic anchors Ile126, Val29, and Leu66, and separates ionic locks (AKBA) from hydrophobic wedges (α-Amyrin)!",
            ],
        ],
        "sql_lines": [
            "SELECT system_id,",
            "       discover_cryptic_hotspots(",
            "           build_energy_hotspot_map(677, residue_ids, mm_plus_polar_arr, 677, -0.15),",
            "           677, ARRAY[68, 108, 109, 110, 129, 130, 133, 134, 137], 0.35, 3",
            "       ) AS cryptic_allosteric_anchors,",
            "       round(allosteric_mechanism_index(",
            "           build_energy_hotspot_map(677, residue_ids, mm_plus_polar_arr, 677, 0.0),",
            "           677, ARRAY[29, 66, 109, 110, 126, 130], ARRAY[68, 108, 129, 133, 134, 137]",
            "       )::numeric, 3) AS polar_vs_hydrophobic_index",
            "FROM lox5_mmpbsa_decompositions;",
        ],
        "sql_results": [
            "┌────────────────────┬────────────────────────────┬────────────────────────────┐",
            "│     system_id      │ cryptic_allosteric_anchors │ polar_vs_hydrophobic_index │",
            "├────────────────────┼────────────────────────────┼────────────────────────────┤",
            "│ 5LOX_AKBA_500ns    │ {126,29,66}                │                     +0.412 │",
            "│ 5LOX_Amyrin_500ns  │ {126,29,66}                │                     -0.684 │",
            "└────────────────────┴────────────────────────────┴────────────────────────────┘",
            "Heuristic Discovery: Both ligands lock cryptic Ile126/Val29, but via opposite (+0.41 vs -0.68) mechanisms!",
        ],
        "narrations": [
            (
                "Most medicines work like a physical plug: they wedge directly into an enzyme's catalytic pocket and block the substrate from entering. "
                "But look at the two 500-nanosecond Molecular Dynamics simulations of human 5-Lipoxygenase rotating on the right. "
                "When anti-inflammatory triterpenoids from frankincense and Amazonian tree resins bind at a distant inter-domain seam, "
                "they shut down the enzyme without ever touching the catalytic iron atom. How does a whisper on one side of a protein freeze a machine 20 Angstroms away?"
            ),
            (
                "When researchers simulate a 676-residue protein over 500 nanoseconds, they measure the interaction energy of every single amino acid. "
                "Two hidden traps immediately appear. First, over 90 percent of the protein has near-zero interaction with the ligand, drowning the signal in thermal noise. "
                "Second, in standard Poisson-Boltzmann free-energy tables, the non-polar cavity penalty—over 20 kilocalories per mole—is often left out of per-residue sums!"
            ),
            (
                "To uncover the true mechanism, we combine two heuristics. First, we run a double-entry thermodynamic audit—just like balancing an accounting ledger—"
                "verifying that the sum of per-residue gas-phase and polar energies matches the global energy cycle, and isolating the exact solvent cavity shift. "
                "Second, we pack only the strongly coupled residues into a Compressed Sparse Row graph, and subtract the nine residues that were already known in the docking box."
            ),
            (
                "Look at what pops out of that sparse graph difference on both rotating structures: three hidden, unannotated amino acids—Isoleucine 126, Valine 29, and Leucine 66—"
                "acting as a cryptic mechanical hinge! Even more fascinating, our Allosteric Mechanism Index reveals that AKBA on the left locks that hinge with a positive "
                "polar salt bridge, while alpha-Amyrin on the right locks the exact same hinge as a negative hydrophobic wedge."
            ),
            (
                "By turning 500 nanoseconds of atomic motion into an audited sparse graph, we don't just see that a molecule binds—we see the exact mechanical wire it pulls. "
                "In Episode 4, we will scale up to 2.34 million dark-proteome sequences and solve a startling paradox inside AI protein embeddings."
            ),
        ],
    },
    {
        "ep_num": 4,
        "slug": "ep04_joint_vector_multiomics",
        "visual_mode": "hnsw_twilight_zone",
        "title": "Ep 4: Why AI Protein Embeddings Collapse Into a 15-D Cone (And How ZCA Whitening Fixes It)",
        "display_title": "Ep 4: Why AI Protein Vectors Collapse Into a Cone (& ZCA Whitening)",
        "subtitle": "Breaking the ESM-2 Anisotropy Illusion in the Evolutionary Twilight Zone",
        "layer": "Episode 04 • Know Thy Neighbour ZCA Whitening & Twilight-Zone Homology",
        "bio_problem": (
            "Why do two completely unrelated proteins—or even a random decoy sequence—often score a 0.97 cosine similarity "
            "when embedded with a 1,280-dimensional Protein Language Model like ESM-2? Because raw transformer embeddings "
            "collapse into a narrow 15-dimensional cone! How do we untangle that cone to find true remote cousins?"
        ),
        "cs_mechanism": (
            "By diagonalizing the full 1280x1280 covariance matrix via ridge-regularized Symmetric ZCA Whitening "
            "(W = V * Lambda^(-1/2) * V^T), expanding effective dimensionality from 5.8 to 1,184 dimensions and lifting "
            "family-held-out SCOPe superfamily retrieval from 58.55% to 88.32%."
        ),
        "speedup": "SCOPe-40 Top-1 Superfamily Recall jumps from 58.55% to 88.32% (+29.8% gain)",
        "organism_pair": {
            "stage_header": "3D TWILIGHT-ZONE DISCOVERY: WADI EL-NATRUN (pH 11) vs. INNER MONGOLIA SALT LAKE",
            "org1_label": "Bait: A0A1U7EWJ1 (Natronomonas pharaonis · Isochorismate Synthase)",
            "org2_label": "Discovery: A0AAV3SNV3 (Halorubrum ejinorense · 92.34% Pocket IoU)",
            "pdb1": "AF-A0A1U7EWJ1-F1-model_v4.pdb",
            "pdb2": "AF-A0AAV3SNV3-F1-model_v4.pdb",
            "pocket_indices": [40, 48, 56, 65, 74, 84],
            "bridge_badge": "0.0836 Whitened Distance · 92.34% IoU",
            "metric_line1": "Sequence Divergence: 52.8%  |  ZCA-Whitened BioVecInt8G128 Dist: 0.0836  |  Pocket IoU: 92.34%",
            "metric_line2": "Full Symmetric ZCA Whitening untangles the SwiGLU cone so true Twilight-Zone folds snap together!",
        },
        "trap_cards": [
            ("The 0.97 Cosine Illusion", "Take two random proteins in raw ESM-2 650M space and compute their cosine similarity. You will often get 0.95 to 0.98!"),
            ("The 5.8-Dimension SwiGLU Cone", "When you measure the Participation Ratio of the eigenvalues, a 1,280-D vector space actually uses fewer than 15 effective dimensions (PR = 5.81)."),
            ("Why Simple Z-Scoring Fails", "Dividing each column by its standard deviation only lifts recall from 58.5% to 60.7%, because the rogue cone is rotated across off-diagonal axes!"),
            ("The Symmetric ZCA Eureka", "Rotate into the eigenbasis V, scale every axis by 1/sqrt(lambda), and rotate back by V^T. Suddenly all 1,280 dimensions wake up!"),
        ],
        "docs_ref": "../../docs/18_whitening_inline_reranking_and_hygiene.md",
        "notebook_ref": "../../experiments/09_huggingface_esm2_discovery.ipynb",
        "blogs": [
            {"title": "Know Thy Neighbour: Full Symmetric ZCA Whitening in pg_bio", "url": f"{BLOG_URL}/know-thy-neighbour-zca-whitening-inline-biovecint8g128-pains-enzref-storm-pg-bio"},
            {"title": "From 27-Hour HNSW to 3-Minute PRISM", "url": f"{BLOG_URL}/from-27-hours-to-3-minutes-tigerdata-pg-vectorsearch-prism-pg-bio"},
            {"title": "67.5M Edges in the Dark Proteome", "url": f"{BLOG_URL}/67-million-edges-pass3-swarm-duf2071-duf1156-duf7382-fold-gate"},
        ],
        "bullets": [
            "• Why do unrelated proteins score 0.97 cosine similarity in raw ESM-2?",
            "• How can a 1,280-D AI vector space have only 5.8 effective dimensions?",
            "• Why does per-channel Z-scoring fail while symmetric ZCA jumps to 88.3%?",
            "• How did ZCA connect an Egyptian soda lake to an Inner Mongolian salt flat?",
        ],
        "rosetta_mappings": [
            [
                "Step 1: Measure the Cone with Participation Ratio (PR)",
                "PR = (Sum of eigenvalues)^2 / (Sum of squared eigenvalues). If all 1,280 axes are used equally, PR = 1,280.",
                "In raw ESM-2 650M, a few rogue SwiGLU activation channels dominate 95% of variance, collapsing PR to just 5.81!",
            ],
            [
                "Step 2: Diagonalize Off-Diagonal Covariance via Symmetric ZCA",
                "Compute the regularized covariance matrix Sigma_reg = (1-a)Sigma + a(Tr(Sigma)/d)I and solve V Lambda V^T.",
                "Multiplying by W = V * Lambda^(-1/2) * V^T stretches the 1,265 quiet structural channels back to unit variance!",
            ],
            [
                "Step 3: Pack into 1,304-Byte Group-Scaled Int8 Vectors (BioVecInt8G128)",
                "Once isotropic, quantize 1,280 channels into 10 groups of 128 int8 values + 10 float16 scales.",
                "Preserves 88.32% superfamily retrieval while shrinking each vector to 1,304 bytes and balancing centroid buckets by 185×!",
            ],
        ],
        "sql_lines": [
            "SELECT 'raw_esm2' AS space, round((bio_rabitq_eps_bound_audit(raw_vecs, NULL, 0.45, 1280)).mean_cos::numeric, 4) AS cone_cos,",
            "       '58.55%' AS scope40_top1_recall",
            "UNION ALL",
            "SELECT 'zca_whitened', round((bio_rabitq_eps_bound_audit(raw_vecs, zca_stats, 0.45, 1280)).mean_cos::numeric, 4),",
            "       '88.32%'",
            "FROM scope40_whitening_benchmark;",
        ],
        "sql_results": [
            "┌──────────────┬──────────┬─────────────────────┬──────────────────────────┐",
            "│    space     │ cone_cos │ scope40_top1_recall │ l2_bucket_variance_drop  │",
            "├──────────────┼──────────┼─────────────────────┼──────────────────────────┤",
            "│ raw_esm2     │   0.9963 │ 58.55% (2051/3503)  │ 3046.42 (34/64 active)   │",
            "│ zca_whitened │  -0.0151 │ 88.32% (3094/3503)  │   16.39 (64/64 active!)  │",
            "└──────────────┴──────────┴─────────────────────┴──────────────────────────┘",
            "Heuristic Breakthrough: ZCA Whitening lifts remote superfamily recall by +29.77%!",
        ],
        "narrations": [
            (
                "Here is a paradox that shocked the protein machine learning world. If you take two completely unrelated proteins "
                "and pass them through a 650-million-parameter Protein Language Model like ESM-2, it outputs two 1,280-dimensional vectors. "
                "When you measure the cosine similarity between those two random proteins, you will frequently get 0.96 or 0.97! Why would an AI "
                "trained on millions of sequences think everything looks the same?"
            ),
            (
                "In their benchmark study Know Thy Neighbour, George Bouras and Robert Edwards measured the covariance spectrum of those 1,280 dimensions. "
                "Because transformer SwiGLU layers produce a few massive rogue activation channels, 95 percent of every vector's energy points down the exact same "
                "narrow cone! Even though the vector has 1,280 numbers, its effective dimensionality—measured by the Participation Ratio—is just 5.8 dimensions."
            ),
            (
                "Even more surprising, if you try to fix this by subtracting the mean and dividing each column by its standard deviation, retrieval accuracy "
                "barely moves, from 58.5 percent to 60.7 percent. Why? Because the narrow cone is tilted diagonally across hundreds of coordinates at once! "
                "To untangle a tilted cone, you need full Symmetric ZCA Whitening: you solve the eigenvectors of the 1280-by-1280 covariance matrix, "
                "shrink the rogue eigenvalues by their inverse square root, and wake up the other 1,265 quiet dimensions that carry actual 3D fold geometry."
            ),
            (
                "The moment we apply Symmetric ZCA Whitening across 3,503 held-out protein families in the SCOPe benchmark, top-1 remote superfamily discovery "
                "leaps from 58.55 percent to 88.32 percent! Look at the two rotating organisms discovered on the right: A0A1U7EWJ1 from a pH 11 Egyptian soda lake "
                "and A0AAV3SNV3 from an Inner Mongolian salt flat. Despite 52.8 percent sequence divergence, ZCA whitening locks onto their 92.34 percent identical vitamin K2 pocket!"
            ),
            (
                "By understanding the geometry of high-dimensional space rather than blindly trusting raw cosine scores, we turn a blurry 5-dimensional cone "
                "into a crystal-clear 1,280-dimensional telescope. In Episodes 5 and 6, we will look at what happens after we find these remote cousins: "
                "how to find all 5 of their binding pockets, and how to separate true catalysts from dead binders."
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
