# Biotech-for-Engineers Curriculum Tracks (`pg_bio`)

When a software or database engineer asks to learn biotech through `pg_bio`, guide them through these 4 progressive tracks. Each module pairs a **Database/Rust Concept** with a **Biotech Domain Milestone** and a **Hands-On Lab**.

---

## Track 1: From Strings & B-Trees to the Central Dogma (Genomics $\rightarrow$ Sequence Physics)
**Target Audience:** Backend / SQL engineers taking their first step into bioinformatics.

### Module 1.1: The Genome as a 1D Spatial Interval Index
* **CS Concept:** `int4range`, GiST indexes, overlap (`&&`) and distance (`<->`) operators, streaming VCF parsing (`TableIterator`).
* **Biotech Concept:** Chromosomes, promoters, enhancers, tumor suppressors, and "Safe Harbor" CRISPR knock-in sites.
* **Code & Docs:** [`docs/07_genomics_io.md`](../../docs/07_genomics_io.md) | [`experiments/06_genomics_io_tutorial.ipynb`](../../experiments/06_genomics_io_tutorial.ipynb)
* **Hands-On SQL Lab:** Parse a `.vcf` variant block via `parse_vcf()`, translate coding DNA via `bio_translate_dna()`, and query a safe-harbor locus near a tissue promoter that avoids tumor suppressor intervals.

### Module 1.2: Sequence Physics in 200 Nanoseconds (Why Extremophiles Survive)
* **CS Concept:** Single-pass ASCII byte histograms, Shannon entropy (`sequence_entropy`), Henderson-Hasselbalch root-finding bisection (`isoelectric_point`), and finite-state token matching (`prosite_match`).
* **Biotech Concept:** Amino acid pKa charge curves, hydropathy (`gravy_score`), low-complexity disordered regions, and catalytic PROSITE motifs.
* **Code & Docs:** [`docs/04_native_chemistry.md`](../../docs/04_native_chemistry.md) | [`docs/11_sequence_analysis.md`](../../docs/11_sequence_analysis.md)
* **Hands-On SQL Lab:** Compare a mesophilic enzyme against a halophilic (salt-loving) archaeon ortholog using `isoelectric_point(seq)` and `gravy_score(seq)` to observe the acidic surface shift ($\Delta\text{pI}$).

---

## Track 2: From HNSW Vectors & Morton Curves to 3D Structural Biology
**Target Audience:** AI/ML and spatial-database engineers learning protein folding and structural homology.

### Module 2.1: Mining the Dark Proteome in the "Twilight Zone"
* **CS Concept:** `halfvec` 1280-D HNSW cosine distance (`<=>`) short-circuiting into $O(N)$ rolling-buffer Smith-Waterman dynamic programming (`smith_waterman_score`).
* **Biotech Concept:** Why 3D structure is more conserved than 1D sequence; discovering uncharacterized extremophile enzymes with <30% sequence identity.
* **Code & Docs:** [`docs/01_vector_homology.md`](../../docs/01_vector_homology.md) | [`docs/08_native_alignment.md`](../../docs/08_native_alignment.md) | [`examples/swarm_discovery/swarm_discovery.py`](../../examples/swarm_discovery/swarm_discovery.py)
* **Hands-On SQL Lab:** Run the Twilight-Zone bell-curve query combining `<=>` and `smith_waterman_score` to find a divergent structural twin, then fetch its live metadata via `bio_fetch_uniprot()`.

### Module 2.2: 3D Atomic Space via 64-Bit Morton Curves & Jacobi SVD
* **CS Concept:** Bit-interleaving 3D $(x,y,z)$ floats into 1D `BIGINT` keys (`z_order_encode`) for B-Tree bounding boxes + stack-allocated $3\times 3$ Jacobi SVD (`kabsch_rmsd`).
* **Biotech Concept:** Angstrom ($\text{Å}$) atomic coordinates, binding pockets, $\text{C}\alpha$ backbones, and rotational superposition.
* **Code & Docs:** [`docs/03_spatial_indexing.md`](../../docs/03_spatial_indexing.md) | [`docs/09_pdb_parsing.md`](../../docs/09_pdb_parsing.md) | [`experiments/08_pdb_parsing.ipynb`](../../experiments/08_pdb_parsing.ipynb)
* **Hands-On SQL Lab:** Parse a PDB structure with `bio_parse_pdb()`, extract a $5\text{ Å}$ pocket around a catalytic metal ion using `z_order_encode`, and visualize the pocket in the **Biomolecular Structure Viewer** sidecar.

---

## Track 3: Cheminformatics & High-Throughput Virtual Screening in SQL
**Target Audience:** Graph algorithm and search-engine engineers exploring small-molecule drug discovery.

### Module 3.1: Molecular Graphs, Hückel $4n+2$ Aromaticity & VF2 Isomorphism
* **CS Concept:** Cycle perception via DFS, $\pi$-electron bitmask checks, 64-bit Weisfeiler-Lehman canonical hashing (`molecule_wl_hash`), VF2 backtracking subgraph isomorphism (`@>`), and hardware `POPCNT` Tanimoto (`%`).
* **Biotech Concept:** SMILES strings, Kekulé resonance, counter-ion salts, keto-enol tautomers, pharmacophore warheads, and SAR similarity.
* **Code & Docs:** [`docs/14_aromaticity_wl_vf2_cheminformatics.md`](../../docs/14_aromaticity_wl_vf2_cheminformatics.md) | [`experiments/05_cheminformatics_tutorial.ipynb`](../../experiments/05_cheminformatics_tutorial.ipynb)
* **Hands-On SQL Lab:** Prove canonical equivalence across Kekulé, aromatic, and salted SMILES using `molecule_wl_hash`, then screen substrates with `@>` and `%`.

### Module 3.2: 100M+ Ligand Screening with RangeFinder Pocket-Interval GiST
* **CS Concept:** 8-column `numrange` GiST index (`@>`) and composite B-Tree range pruning driven by pocket-point-cloud inference (`predict_pocket_ranges`).
* **Biotech Concept:** Physicochemical complementarity (`mw`, `logp`, `tpsa`, `hbd`, `hba`, `rotb`, `q`, `fsp3`) and reverse target fishing across the human proteome.
* **Code & Docs:** [`docs/RANGEFINDER_GIST_INDEXING.md`](../../docs/RANGEFINDER_GIST_INDEXING.md) | [`experiments/rangefinder_pocket_gist_screening.ipynb`](../../experiments/rangefinder_pocket_gist_screening.ipynb)
* **Hands-On SQL Lab:** Predict compatible ligand intervals from a $10\text{ Å}$ pocket and run both forward B-Tree ligand filtering and reverse GiST target fishing in sub-millisecond SQL.

---

## Track 4: De Novo Enzyme Design, Chaperone Pipelines & Wet-Lab Compilation
**Target Audience:** Systems architects building end-to-end generative AI-to-wet-lab protein engineering pipelines.

### Module 4.1: Theozyme Geometry, SASA/BUNS Auditing & RosettaCommons Chaperone
* **CS Concept:** $5\text{ Å}$ 3D voxel grid hashing + 64-point Fibonacci sphere SASA (`bio_rasa`), CSR H-bond networks (`bio_hbond_network`), BUNS polar burial detection (`audit_hydrogen_bonds`), and lineage hashing (`chaperone_lineage_hash`).
* **Biotech Concept:** Why generative backbones (`RFdiffusion3 -> Foundry-MPNN -> Boltz-2/AF3`) fail without coordinated $\text{Zn}^{2+}$ water shells, Bürgi–Dunitz attack angles ($100^\circ\text{–}115^\circ$), and zero Buried Unsatisfied Polar Atoms (BUNS).
* **Code & Docs:** [`docs/15_chaperone_protein_design_workflows.md`](../../docs/15_chaperone_protein_design_workflows.md) | [`experiments/10_metalloprotease_theozyme_design.ipynb`](../../experiments/10_metalloprotease_theozyme_design.ipynb) | [`experiments/13_chaperone_protein_design_workflows.ipynb`](../../experiments/13_chaperone_protein_design_workflows.ipynb)
* **Hands-On SQL Lab:** Audit all 8 Official RosettaCommons Chaperone Core Examples in a single `CROSS JOIN LATERAL evaluate_chaperone_core_example_block(...)` query.

### Module 4.2: Thermodynamic Ledgers, RK4 Cascades & Wet-Lab DNA Cassettes
* **CS Concept:** Double-entry MM/PBSA reconciliation (`audit_mmpbsa_decomposition`), CSR cryptic pocket discovery (`discover_cryptic_hotspots`), 4th-order Runge-Kutta ODE integration (`simulate_cascade_kinetics`), and Golden Gate DNA cassette compilation (`bio_design_expression_cassette`).
* **Biotech Concept:** Cavity/dispersion solvation penalties, allosteric networks, multi-enzyme product inhibition sinks, and host-specific codon optimization (*E. coli*, *B. subtilis*, *H. salinarum*).
* **Code & Docs:** [`docs/12_mmpbsa_free_energy_auditing.md`](../../docs/12_mmpbsa_free_energy_auditing.md) | [`docs/13_wetlab_synbio_and_teaflon_validation.md`](../../docs/13_wetlab_synbio_and_teaflon_validation.md)
* **Hands-On SQL Lab:** Compile an expression-ready, `BsaI`-free Golden Gate DNA construct for a validated enzyme design and simulate its 48-hour reaction kinetics inside PostgreSQL.
