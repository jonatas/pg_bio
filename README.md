# 🧬 pg_bio: The Multiomics Operating System for Synthetic Biology

**`pg_bio`** is a native PostgreSQL extension written in highly optimized Rust (`pgrx`). 

It is engineered for a single purpose: to completely eradicate the data bottlenecks in modern biology by pushing complex multiomics mathematics *directly into the database engine*. By unifying Genomics, Proteomics, and Transcriptomics into a shared computational space, `pg_bio` allows scientists to execute entire drug discovery pipelines and synthetic biology designs in milliseconds natively in SQL.

> 🔗 **The `pg_bio` Synthetic Biology & Multiomics Ecosystem:**
> - **Core Database Engine (Rust `pgrx` + PostgreSQL 18):** [**`jonatas/pg_bio`**](https://github.com/jonatas/pg_bio) *(this repository — 3D Morton Octrees, PRISM + Inline `BioVecInt8G128` Vector Homology, RangeFinder 10 Å Pocket GiST, Cheminformatics, and RK4 Solvers)*
> - **Declarative Bioprocess Compiler, Autonomous Protein Correlation Swarm & Unified Web UI:** [**`jonatas/bio-pipelines`**](https://github.com/jonatas/bio-pipelines) — Multi-feedstock 4-Module ($M_1 \rightarrow M_2 \rightarrow M_3 \rightarrow M_4$) bioprocess compiler, autonomous dark-proteome discovery swarm ([`base/swarm_discovery.py`](https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_discovery.py)), and live 3D Kabsch / RangeFinder / CETSA-TPP Explorer UI ([`base/swarm_web.py`](https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_web.py))
> - **Original PTFE $\rightarrow$ Fluorapatite Case Study:** [**`jonatas/teaflon`**](https://github.com/jonatas/teaflon) ([`ideia.me/teaflon-intro`](https://ideia.me/teaflon-intro))

> **🎬 Interactive Video Series & 3D Studio: "Hello World, `pg_bio`!"**
> New to computational biology or multiomics in PostgreSQL? Explore our narrated 1080p video series paired with interactive **3Dmol.js** structures and **Systems Engineer $\leftrightarrow$ Structural Biologist** Rosetta Stone lenses:
> * 🌐 **Live 3D Video Library on `ideia.me`:** [https://ideia.me/pg-bio-videos/](https://ideia.me/pg-bio-videos/)
> * 🗂️ **Local Video & 3D Studio Catalog:** [`videos/README.md`](videos/README.md) · [`videos/series_01_hello_world_pg_bio/index.html`](videos/series_01_hello_world_pg_bio/index.html)
> * **Episodes:**
>   1. [**Ep 01: The Central Dogma in SQL (`int4range` GiST & Sequence Physics)**](videos/series_01_hello_world_pg_bio/ep01_central_dogma_and_genomics/studio_deck.html) — [MP4](videos/series_01_hello_world_pg_bio/ep01_central_dogma_and_genomics/ep01_central_dogma_and_genomics.mp4) · [SQL](videos/series_01_hello_world_pg_bio/ep01_central_dogma_and_genomics/live_demo.sql)
>   2. [**Ep 02: From 1D Strings to 3D Atoms with Morton Curves (`z_order_encode` & SVD)**](videos/series_01_hello_world_pg_bio/ep02_3d_pockets_and_morton_curves/studio_deck.html) — [MP4](videos/series_01_hello_world_pg_bio/ep02_3d_pockets_and_morton_curves/ep02_3d_pockets_and_morton_curves.mp4) · [SQL](videos/series_01_hello_world_pg_bio/ep02_3d_pockets_and_morton_curves/live_demo.sql)
>   3. [**Ep 03: Transcriptomics & CSR Sparse Attention Graphs (`SparseAttentionMap`)**](videos/series_01_hello_world_pg_bio/ep03_scrna_csr_and_attention_maps/studio_deck.html) — [MP4](videos/series_01_hello_world_pg_bio/ep03_scrna_csr_and_attention_maps/ep03_scrna_csr_and_attention_maps.mp4) · [SQL](videos/series_01_hello_world_pg_bio/ep03_scrna_csr_and_attention_maps/live_demo.sql)
>   4. [**Ep 04: The Multiomics Holy Grail — Vectors & Live UniProt (`halfvec` HNSW)**](videos/series_01_hello_world_pg_bio/ep04_joint_vector_multiomics/studio_deck.html) — [MP4](videos/series_01_hello_world_pg_bio/ep04_joint_vector_multiomics/ep04_joint_vector_multiomics.mp4) · [SQL](videos/series_01_hello_world_pg_bio/ep04_joint_vector_multiomics/live_demo.sql)

---

## 🛑 The Data Bottleneck in Modern Biology

The ultimate frontier of biology is **Multiomics**—the fusion of DNA, RNA, and Protein data. Today, this is done using duct-tape pipelines:
1. Downloading terabytes of sequence data from UniProt.
2. Building an in-memory KD-Tree in Python to measure 3D distances.
3. Using custom tools like `bedtools` for genomic overlaps.
4. Loading 20-billion-row RNA matrices into massive Pandas dataframes.

**The result:** The network I/O and RAM overhead crushes most laptops, and forces cloud deployments into hours of agonizing processing. 

## 🚀 The Solution: The Multiomics Stack inside Postgres

`pg_bio` eliminates this entire nightmare. It provides four distinct computational layers natively inside PostgreSQL. We move the compute directly to the data, resulting in order-of-magnitude speedups.

---

## 🔬 Deep Dive: The 4 Multiomics Layers

### Layer 1: Genomics (1D Spatial Targeting)
When engineering synthetic genes (e.g., bio-ceramics or CRISPR therapies), you must find a "Safe Harbor" in the genome. It must be near the correct tissue-specific promoter (to turn it on) but cannot overlap with any tumor suppressor genes (which could cause cancer).

* **The Engine:** `pg_bio` treats the human genome as a 1D spatial environment utilizing native PostgreSQL **Range Types (`int4range`)** and **`GiST` (Generalized Search Tree) indexes**.
* **The Advantage:** You can use native overlap operators (`&&`) and adjacency operators (`<->`) to scan billions of base pairs in milliseconds.
* **SQL Example:**
  ```sql
  -- Instantly find a Safe Harbor insertion site next to the Nail Promoter, avoiding Cancer genes.
  SELECT safe.genomic_range FROM human_genome safe
  JOIN human_genome promoter ON safe.chromosome = promoter.chromosome
  WHERE safe.genomic_range <-> promoter.genomic_range < 5000
    AND NOT EXISTS (
      SELECT 1 FROM human_genome cancer_genes
      WHERE cancer_genes.feature_type = 'tumor_suppressor' 
        AND safe.genomic_range && cancer_genes.genomic_range
    );
  ```

### Layer 2: Proteomics & Epigenomics (3D Spatial Curves & Layered Morton Octrees)
Once transcribed, linear sequences fold into 3D atomic structures. Finding which atoms are within a 4.0–10.0 Ångstrom radius of a binding pocket usually requires extracting millions of coordinates into Python for brute-force Pythagorean math—and a naive single-interval 1D Z-Order range (`z_index BETWEEN min_z AND max_z`) sweeps across thousands of false-positive atoms whenever the sphere crosses an octree boundary.

* **The Engine:** `pg_bio` pairs `z_order_encode(x, y, z)` (64-bit **Morton Coding / Z-Order curves** interleaving 21-bit 3D coordinates) with `z_order_layered_sphere_ranges(cx, cy, cz, r, max_ranges)`—a **21-level 3D Octree Centroid Router** that prunes disjoint octants ($d_{\min}^2 > R^2$), flags interior cubes ($d_{\max}^2 \le R^2$) to skip `sqrt()` rechecks, and coalesces the remaining leaf intervals into $\le 64$ tight B-tree Z-ranges.
* **The Advantage:** Across 200,000 atoms, `z_order_layered_sphere_ranges` scans **74.05× fewer B-tree candidate rows** (`407.6` vs. `30,181` atoms) and runs **10.25× faster (`0.96 ms` vs. `9.89 ms`)** with **100.0% exact mathematical recall**.
* **Reference Documentation & Benchmark:** [`docs/17_layered_centroids_multi_index_experiments.md`](docs/17_layered_centroids_multi_index_experiments.md) | [`scripts/benchmark_layered_centroids_multi_index.py`](scripts/benchmark_layered_centroids_multi_index.py)
* **SQL Example:**
  ```sql
  -- Extracting a 10.0 Å 3D Binding Pocket in 0.96 ms (74.05x fewer atoms scanned, 100% exact recall)
  SELECT a.atom_id, (a.coord).name
  FROM z_order_layered_sphere_ranges(10.0, 0.0, 25.0, 5.0, 64) r
  JOIN teaflon_atoms a ON a.z_curve_index BETWEEN r.min_z AND r.max_z
  WHERE r.fully_inside
     OR ((a.coord).x - 10.0)^2 + ((a.coord).y - 0.0)^2 + ((a.coord).z - 25.0)^2 <= 25.0;
  ```

### Layer 3: Transcriptomics (Single-Cell CSR Arrays)
To understand drug dosage, scientists use Single-Cell RNA-Seq. Humans have 20,000 genes. Sequencing 1 million cells creates a 20-Billion point matrix. Traditional SQL schemas grind to a halt.

* **The Engine:** Because a single cell only expresses ~2,000 genes at a time, the matrix is 85%+ zeroes. `pg_bio` implements **Compressed Sparse Row (CSR)** architecture natively using parallel PostgreSQL Arrays (`expressed_gene_ids INT[]` and `expression_counts REAL[]`).
* **The Advantage:** Storage overhead is slashed by 90%. Using parallel `LATERAL UNNEST() WITH ORDINALITY`, PostgreSQL unpacks the sparse data dynamically in execution memory to calculate highly precise RNA dosage levels in milliseconds.

### Layer 4: Joint Vector Spaces (`pg_vectorsearch` PRISM + `pgvector` `halfvec(1280) STORAGE PLAIN`)
Traditionally, DNA sequences and 3D Protein folds live in different worlds. `pg_bio` breaks down that wall.

* **The Engine:** By combining Tiger Data's [`pg_vectorsearch`](https://github.com/timescale/pg_vectorsearch) **PRISM** index access method (1-bit RaBitQ quantization with theoretical error bounds, SOAR spill-orthogonal assignment, and hierarchical centroid routing) with `pgvector`'s `halfvec(1280) STORAGE PLAIN` columns and `pg_bio`'s `<~>` hybrid Smith-Waterman operator, we map Genomic Foundation Models (like Enformer) and Proteomic Foundation Models (like ESM-2) into a **shared 1280-dimensional mathematical space**.
* **The Advantage:** Across **2,340,766** UniProt proteins, PRISM shrinks the vector index by **9.72×** (`8,772 MB` HNSW $\rightarrow$ `902 MB` PRISM), saves **7.65 GiB** of total disk/RAM, builds in **3.08 minutes** (down from 27 hours), and accelerates `<=>` cosine and `<~>` hybrid homology searches by **7.5× to 11.3×** (`26–37 ms` vs `197–417 ms`).
* **Reference Documentation & Benchmark:** [`docs/16_pg_vectorsearch_vs_pgvector_benchmark.md`](docs/16_pg_vectorsearch_vs_pgvector_benchmark.md) | [`scripts/setup_and_benchmark_parallel_2m_proteome.py`](scripts/setup_and_benchmark_parallel_2m_proteome.py)
* **SQL Example:**
  ```sql
  -- Build a 902 MB PRISM index over 2.34M 1280-D protein embeddings in ~3 minutes:
  CREATE INDEX idx_proteins_prism_embedding ON proteins USING prism (embedding vec16_cosine_ops);

  -- Find which DNA enhancer matches the structural embedding of the P0DPB7 Protein:
  SELECT dna.locus, (dna.embedding <=> protein.embedding) as cross_omic_distance
  FROM dna_embeddings dna, proteins protein
  WHERE protein.uniprot_id = 'P0DPB7'
  ORDER BY cross_omic_distance ASC LIMIT 1;
  ```

### Layer 5: Direct Biological Database Integration (UniProt API inside SQL)
Instead of relying on external Python scripts to fetch biological metadata, `pg_bio` contains a native Rust-powered synchronous HTTP client (`ureq`). It can fetch and parse JSON from the UniProt API directly into specialized PostgreSQL Composite Types dynamically.

* **The Engine:** The `bio_fetch_uniprot()` function queries `rest.uniprot.org`, parses the JSON via `serde`, and projects it instantly into a Postgres `UniprotEntry` type (`id`, `organism`, `taxonomy`, `sequence`).
* **The Advantage:** Zero external scripts. You can enrich millions of orphan protein IDs dynamically inside pure SQL pipelines. 
* **SQL Example (Multi-Row Batch Fetch & Insert):**

  Using a `CROSS JOIN LATERAL`, you can map the HTTP fetch over a batch of IDs, evaluate it exactly once per row (saving network calls), and `INSERT` the results directly into your metadata table in a single atomic SQL transaction!

  ```sql
  -- 1. Create a table for the newly discovered metadata
  CREATE TABLE dark_proteome_metadata (
      uniprot_id VARCHAR(20) PRIMARY KEY,
      organism TEXT,
      taxonomy TEXT[],
      sequence TEXT
  );

  -- 2. The Batch Fetch & Insert Pipeline
  WITH new_discoveries(uniprot_id) AS (
      VALUES 
          ('A0A5B9DCV6'), -- Deep-sea Asgard archaeon (RuBisCO)
          ('L0JL84')      -- Halophilic archaeon
  ),
  enriched_data AS (
      -- The LATERAL join ensures the HTTP request evaluates exactly once per ID
      SELECT u.id, u.organism, u.taxonomy, u.sequence
      FROM new_discoveries d
      CROSS JOIN LATERAL bio_fetch_uniprot(d.uniprot_id) u
  )
  INSERT INTO dark_proteome_metadata (uniprot_id, organism, taxonomy, sequence)
  SELECT id, organism, taxonomy, sequence 
  FROM enriched_data
  RETURNING uniprot_id, organism;
  ```

* **SQL Example (Native API Search via Virtual Tables):**

  We also provide a **Set Returning Function (SRF)**, which acts mechanically identical to a read-only Foreign Data Wrapper (FDW). You can pass standard UniProt Lucene query parameters directly into the function and query it exactly like a local table.
  
  ```sql
  -- Find the top 5 longest curated proteins for the BRCA1 gene in Humans
  SELECT id, organism, length(sequence) AS seq_len 
  FROM bio_search_uniprot('gene:BRCA1 AND taxonomy_id:9606 AND reviewed:true')
  ORDER BY seq_len DESC 
  LIMIT 5;

  -- Filter archaeal proteins by exact length and sequence
  SELECT id, sequence
  FROM bio_search_uniprot('taxonomy_id:2594042 AND length:[1 TO 50]')
  WHERE sequence LIKE 'M%'; -- standard SQL filters apply to the returning set
  ```
  
  **Supported UniProt Query Parameters:**
  The `query` argument supports all standard UniProtKB REST API search fields, including:
  * `gene:BRCA1` (Filter by gene name)
  * `taxonomy_id:9606` (Filter by taxonomic ID, e.g., 9606 = Homo sapiens)
  * `organism_name:"Homo sapiens"` (Filter by exact organism string)
  * `length:[1 TO 200]` (Filter by amino acid sequence length)
  * `reviewed:true` (Only return manually curated Swiss-Prot entries)
  * `keyword:KW-0002` (Filter by specific biological keywords, e.g., 3D-structure)

---

### Layer 6: Dry Lab Tooling (Structural & Sequence Analysis)
We have fully integrated classic dry-lab functions directly into the SQL engine, removing the need for Biopython or PyMOL in many scenarios.
* **Native Sequence Translation:** `bio_translate_dna('ATGCGT')` instantly translates DNA to Amino Acids using standard codon tables.
* **Native PDB Parsing:** `bio_parse_pdb(pdb_string)` parses raw PDB files into relational SQL rows.
* **Structural Alignment:** `bio_rmsd(coords_a, coords_b)` calculates the RMSD between two protein conformations directly in PostgreSQL.

### Layer 7: Cheminformatics (Hückel $4n+2$ Aromaticity Perception, WL/ECFP4 Hashing, VF2 `@>` & Salt/Tautomer Normalization)
To discover novel drugs and enzyme substrates, you must traverse massive libraries of chemical compounds—where the exact same molecule may be ingested in Kekulé form (`C1=CC=CC=C1`), queried in aromatic form (`c1ccccc1`), or stored as a ChEMBL/PubChem salt (`[Na+].O=C([O-])C1=CC=CC=C1`).

* **The Engine:** `pg_bio` implements a zero-dependency Rust molecular graph normalizer and indexer inside PostgreSQL:
  1. **Salt Stripping & Prototropic Tautomer Normalization:** Splits disconnected `.` fragments into connected components, retains the largest heavy-atom parent molecule, neutralizes counter-ion salt charges (`CC(=O)[O-].[Na+]` $\rightarrow$ `CC(=O)O`), and canonicalizes aliphatic 1,3-prototropic keto-enol / imidic-amide tautomers (`C=CO` $\rightarrow$ `CC=O`, `CC(=N)O` $\rightarrow$ `CC(=O)N`).
  2. **Hückel’s $4n+2$ Aromaticity Perception:** Traces fundamental and fused cycles ($3..=10$), sums $\pi$-electron contributions across $sp^2$ conjugated ring atoms (`(pi_sum >= 6) && ((pi_sum - 2) & 3) == 0`), and rewrites aromatic rings to `BOND_AROMATIC` so `C1=CC=CC=C1` and `c1ccccc1` compile to the **exact same `Molecule` struct**.
  3. **Weisfeiler-Lehman (WL) 64-Bit Canonical Hash & 1024-Bit ECFP4 Bitset:** Computes vertex-permutation-invariant 64-bit `wl_hash` (`molecule_wl_hash`) and `16 × u64` bitsets using `fmix64` and sorted neighbor hashes, enabling $O(1)$ B-tree lookups (**`0.673 ms` across 10,000 molecules — 1,513.7x faster than Python pull**) and hardware-popcount (`count_ones`) Tanimoto similarity (`Molecule % Molecule`).
  4. **VF2 Subgraph Isomorphism Recheck (`@>`):** Exact backtracking subgraph matcher (`target_mol @> query_mol` and `smiles_contains`) verifying atom numbers, aromaticity, and bond types across cyclic and fused scaffolds.
* **Reference Documentation & Benchmark:** [`docs/14_aromaticity_wl_vf2_cheminformatics.md`](docs/14_aromaticity_wl_vf2_cheminformatics.md) | [`scripts/benchmark_aromaticity_wl_vf2.py`](scripts/benchmark_aromaticity_wl_vf2.py)
* **SQL Example:**
  ```sql
  -- 1. Canonical equivalence across Kekulé, Aromatic, and Salted representations (100% recall)
  SELECT molecule_wl_hash('C1=CC=CC=C1') = molecule_wl_hash('c1ccccc1') AS benzene_eq,
         molecules_are_identical('[Na+].O=C([O-])C1=CC=CC=C1', 'O=C(O)c1ccccc1') AS salt_normalized,
         parse_smiles('CC(=O)OC1=CC=CC=C1C(=O)O') @> parse_smiles('c1ccccc1') AS vf2_substructure,
         parse_smiles('C1=CC=CC=C1') % parse_smiles('c1ccccc1') AS tanimoto;

  -- 2. Screen Swarm Discovery enzyme substrates using B-tree WL hash, VF2 (@>), and ECFP4 Tanimoto (%)
  SELECT family_name, substrate_name, wl_hash
  FROM swarm_enzyme_substrates
  WHERE mol @> parse_smiles('c1ccc(C(=O)O)cc1')
    AND (mol % parse_smiles('O=C(OCCO)c1ccc(C(=O)O)cc1')) >= 0.75;
  ```

### Layer 8: De Novo Enzyme Design, Theozyme Geometry & Conformational Preorganization
Scaling *de novo* enzyme design success rates from **14.7%** (`Zn` series) to **87.5%** (`ZnO` series) (Chen, Hilvert, Baker et al., *Nature* 2026) requires moving beyond coarse $\text{C}\alpha$ proximity to explicitly model metal cofactors ($\text{Zn}^{2+}$), coordinated nucleophilic water ($\text{HOH}$), strict Bürgi–Dunitz attack angles ($100^\circ - 115^\circ$), two-sided $\beta$-strand pairing, and conformational preorganization across Molecular Dynamics (MD) / MLFF ensembles.

* **The Engine:** `pg_bio` natively indexes `HETATM` cofactors and catalytic water molecules in Morton-coded `ResidueCoord` arrays, executes rotationally invariant 3D **Kabsch / Horn Quaternion Superposition** (`theozyme_rmsd`, `kabsch_rmsd`) against reference active-site motifs (such as `'4qhp_zn_pentad'`), computes bond/attack angles (`bio_angle`) and torsions (`bio_dihedral`), and evaluates trajectory histograms via the `PreorganizationProfile` composite type (`ensemble_fraction_below`).
* **The Scientist's Advantage:** Screen thousands of RFdiffusion / AlphaFold3 enzyme-substrate complexes and MD ensembles in **< 7 milliseconds** directly inside PostgreSQL—without exporting gigabytes of PDB trajectories to Python/PyRosetta.
* **SQL Example:**
  ```sql
  -- Filter de novo metalloprotease constructs for coordinated Zn(II)-water shells,
  -- sub-angstrom 4qhp theozyme RMSD, Burgi-Dunitz attack angles, and >80% MD preorganization:
  SELECT design_id,
         theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') AS rmsd,
         bio_angle(water_o, substrate_c, substrate_o) AS burgi_dunitz_angle,
         ensemble_fraction_below(proton_transfer_dist_hist, 3.0) AS preorg_fraction
  FROM metalloprotease_designs
  WHERE has_coordinated_water_shell(active_site_coords, 'ZN', 2.4) = TRUE
    AND theozyme_rmsd(active_site_coords, '4qhp_zn_pentad') <= 0.90
    AND bio_angle(water_o, substrate_c, substrate_o) BETWEEN 100.0 AND 115.0
    AND bio_distance(glu_oe, water_h) < 3.0
    AND ensemble_fraction_below(proton_transfer_dist_hist, 3.0) >= 0.80
  ORDER BY rmsd ASC;
  ```

### Layer 9: MM/GBSA & MM/PBSA Thermodynamic Ledger Auditing & CSR Allosteric Hotspot Discovery
When analyzing explicit-solvent Molecular Dynamics (MD) trajectories and per-residue end-state free-energy decompositions (such as the UFAM / COARACI **$9 \times 500\text{ ns}$ 5-lipoxygenase triterpenoid benchmark**, [Zenodo 23565306](https://doi.org/10.5281/zenodo.23565306)), medicinal chemists need to verify whether per-residue decompositions close mathematically against global thermodynamic cycles (`GGAS`, `GSOLV`, `DELTA TOTAL`) and discover stabilizing cryptic sub-pocket residues outside pre-selected binding pockets.

* **The Engine:** `pg_bio` provides native Rust double-entry thermodynamic auditors (`audit_mmpbsa_decomposition`, `solvent_shift_penalty`) and Compressed Sparse Row (CSR) energetic network profilers (`build_energy_hotspot_map`, `discover_cryptic_hotspots`, `allosteric_mechanism_index`).
* **The Scientist's Advantage:** Screens **180 protein-ligand MD systems (`121,860` residue rows) in `7.97 ms` (31.14x faster than Python/NumPy)**, automatically flagging the $+20.72\text{ to }+26.05\text{ kcal/mol}$ non-decomposed Poisson-Boltzmann cavity/dispersion penalty and uncovering unannotated cryptic hydrophobic anchors (`Ile126`, `Val29`, `Leu66` in human 5-LOX).
* **Interactive Tutorial & Blog Post:**
  * 📓 **Executed Jupyter Notebook:** [`experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb`](experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb)
  * 📖 **Reference Documentation:** [`docs/12_mmpbsa_free_energy_auditing.md`](docs/12_mmpbsa_free_energy_auditing.md)
  * 🌐 **Technical Deep-Dive & 3Dmol.js Blog Post:** [Auditing 500 ns Molecular Dynamics in SQL: Triterpenoids, 5-Lipoxygenase, and Cryptic Allosteric Pockets with pg_bio](https://ideia.me/5lox-triterpenoids-mmpbsa-free-energy-audit-pg-bio)
* **SQL Example:**
  ```sql
  -- Build a 676-residue CSR energetic hotspot map and discover cryptic anchors outside the 9-residue pocket
  SELECT
      system_id,
      discover_cryptic_hotspots(
          build_energy_hotspot_map(677, residue_ids, mm_plus_polar_arr, 677, -0.15),
          677,
          ARRAY[68, 108, 109, 110, 129, 130, 133, 134, 137],
          0.35,
          5
      ) AS cryptic_hotspots,
      allosteric_mechanism_index(
          build_energy_hotspot_map(677, residue_ids, mm_plus_polar_arr, 677, 0.0),
          677,
          ARRAY[29, 66, 109, 110, 126, 130],
          ARRAY[68, 108, 129, 133, 134, 137]
      ) AS mechanism_index
  FROM lox5_mmpbsa_decompositions;
  ```

### Layer 12: Wet-Lab Synthetic Biology, Multi-Host Codon Cassettes, 3D Globular Gate & Defluorination Digital Twin
Transitioning an *in silico* multiomics design (such as the **TeaFlon** PTFE-degrading & fluorapatite-mineralizing TriFusion enzyme across **7,338 proteins and 930 organisms**) into physical wet-lab validation requires host-specific codon optimization, plasmid cassette assembly, 3D globular folding quality control, and enzymatic assay digital twins.

* **The Engine:** `pg_bio` provides native Rust wet-lab compilers (`bio_reverse_translate`, `bio_codon_adaptation_index`, `bio_design_expression_cassette`), 3D gyration tensor globular gates (`evaluate_globular_compactness`), built-in Fluoroacetate Dehalogenase theozyme templates (`'q1ju72_dehalogenase_hexad'`), and Michaelis-Menten + Fluorapatite $\text{Ca}_5(\text{PO}_4)_3\text{F}$ kinetic simulators (`simulate_defluorination_kinetics`).
* **Reference Documentation & Python Friction Roadmap:** [`docs/13_wetlab_synbio_and_teaflon_validation.md`](docs/13_wetlab_synbio_and_teaflon_validation.md)
* **SQL Example:**
  ```sql
  -- Design wet-lab expression cassettes across 9 bacterial, fungal, archaeal & mammalian hosts,
  -- audit 3D globular compactness, and simulate 48h Fluoride ISE + Fluorapatite yield:
  SELECT host,
         (bio_design_expression_cassette('MFEGFERRLVDWPSTDKTKREEVD', host, 'T7')).gc_content AS gc,
         (bio_design_expression_cassette('MFEGFERRLVDWPSTDKTKREEVD', host, 'T7')).cai_score AS cai,
         (bio_design_expression_cassette('MFEGFERRLVDWPSTDKTKREEVD', host, 'T7')).restriction_sites_free AS bsa1_free
  FROM unnest(ARRAY['e_coli', 'b_subtilis', 'p_putida', 's_cerevisiae', 'h_salinarum', 'h_sapiens']) AS host;
  ```

### Layer 13: RosettaCommons Chaperone Protein Design Workflows & Official Core Examples (`RFD3 -> Foundry-MPNN -> Boltz-2 / AF3`)
Supporting [**RosettaCommons Chaperone**](https://rosettacommons.org/2026/10/19/introducing-chaperone-a-unified-framework-for-protein-design-workflows/) (`RosettaCommons/chaperone`) and its [8 Official Core Examples](https://rosettacommons.github.io/chaperone/core-examples/) (`nanoparticle-design`, `peptide-binder`, `cortisol-binder`, `dna-binder`, `enzyme-glycosidase`, `enzyme-protease`, `ribozyme-hammerhead`, `multi-mpnn`) requires pushing post-oracle evaluation blocks (**RMSD**, **RASA**, **Hydrogen Bonds**, **BUNS**, and **`PipelineDataHandler` Lineage**) directly into PostgreSQL (**`209.04 ms` for 400 complexes — `36.9x` faster than Python pull**).

* **The Engine:** `pg_bio` implements a fused Rust Chaperone evaluation engine (`bio_rasa` using 64-point Shrake–Rupley Fibonacci spheres and Tien 2013 RASA normalization, `bio_hbond_network` projecting 3D H-bonds into CSR `SparseAttentionMap` matrices, `audit_hydrogen_bonds` counting cross-interface H-bonds and **Buried Unsatisfied Polar Atoms / BUNS**, `evaluate_chaperone_block` / `evaluate_chaperone_design`, `chaperone_lineage_hash`, and `chaperone_core_examples` / `evaluate_chaperone_core_example_block`).
* **Reference Documentation & Benchmark:** [`docs/15_chaperone_protein_design_workflows.md`](docs/15_chaperone_protein_design_workflows.md) | [`scripts/benchmark_chaperone_workflows.py`](scripts/benchmark_chaperone_workflows.py) | [`experiments/13_chaperone_protein_design_workflows.ipynb`](experiments/13_chaperone_protein_design_workflows.ipynb)
* **SQL Example:**
  ```sql
  -- 1. Audit all 8 Official RosettaCommons Chaperone Core Examples in milliseconds:
  SELECT ex.example_name, ex.category,
         audit.sc_rmsd, audit.mean_plddt, audit.interface_delta_sasa,
         audit.cross_interface_hbonds, audit.buried_unsatisfied_polars AS buns,
         audit.active_site_or_motif_rmsd, audit.ligand_or_cofactor_verified,
         audit.passes_example_gate
  FROM chaperone_core_examples() AS ex
  CROSS JOIN LATERAL evaluate_chaperone_core_example_block(
      ex.example_name, ex.candidate_pdb, ex.reference_pdb, ex.ligand_smiles
  ) AS audit;

  -- 2. Filter Chaperone RFdiffusion3 -> Foundry-MPNN -> Boltz-2/AF3 candidates by scRMSD, pLDDT,
  -- interface delta-SASA, cross-interface H-bonds, zero BUNS traps, and compile lineage hashes:
  SELECT d.design_id,
         chaperone_lineage_hash('boltz2', d.parent_mpnn_id, d.sequence, d.cycle_num) AS lineage_id,
         eval.sc_rmsd,
         eval.mean_plddt,
         eval.interface_delta_sasa,
         eval.cross_interface_hbonds,
         eval.buried_unsatisfied_polars AS buns,
         eval.design_composite_score
  FROM chaperone_benchmark_designs d
  CROSS JOIN LATERAL evaluate_chaperone_block(d.oracle_pdb, d.rfd3_backbone_pdb, 'A', 'B') AS eval
  WHERE eval.passes_chaperone_gate = TRUE
  ORDER BY eval.design_composite_score DESC;
  ```

### Layer 19: PocketVina Multi-Pocket SAS Conditioning & 21-Check PoseBusters Physical Validity Gate
Inspired by [**PocketVina Atlas**](https://pocketvina-atlas.mdc-berlin.de/about/) ([Sarigun et al., *arXiv:2506.20043*](https://arxiv.org/abs/2506.20043)) and [**PoseBusters**](https://doi.org/10.1039/D3SC04185A) (Buttenschoen et al., *Chem. Sci.* 2024), proteome-scale virtual screening requires moving beyond single guessed docking boxes to **Multi-Pocket Conditioning** (`pk1..pk5`) while verifying that top poses pass **21 physical, geometric, steric, and salt-hygiene checks**—without silently hiding salted vendor entries (`[Na+].O=C([O-])...`) or teleporting sterically clashing poses into false positives.

* **The Engine:** `pg_bio` implements three native Rust functions directly inside PostgreSQL (**`197.92 ms` for 120 complexes (`1.65 ms/complex`) — `20.96×` faster than Python pull**):
  1. **`bio_detect_top_pockets(pdb_text, max_pockets DEFAULT 5)`:** Probes the protein surface with 24-ray Fibonacci Connolly Solvent Accessible Surface (SAS) spheres ($r_{\text{probe}} = 1.40\,\text{\AA}$) accelerated by a $5.0\,\text{\AA}$ 3D spatial voxel hash grid, filters for buried concave clefts, clusters SAS points via $4.0\,\text{\AA}$ BFS flood-fill ($\ge 7.5\,\text{\AA}$ inter-pocket separation), and emits `pk1..pk5` along with their 64-bit Morton `z_order_center`, 3D SAS coordinates (`sas_coords_flat`), 9D `RangeFinderPocket` intervals, and 16-byte CPU register Fastscan sketch (`fastscan_16b`).
  2. **`bio_posebusters_audit(protein_pdb, ligand_sdf_or_smiles, cx, cy, cz)` & `bio_posebusters_audit_row(...)`:** Audits all **21 PoseBusters & chemical hygiene checks** (valence, connectivity, charge, salt counterion exclusion, covalent bond lengths/angles, internal clash, $3\times 3$ Jacobi SVD aromatic ring flatness $\le 0.25\,\text{\AA}$, $sp^2$ flatness, $sp^3$ tetrahedral volume, pairwise protein–ligand VdW ratio $d_{ij}/(r_{\text{vdW},i}+r_{\text{vdW},j}) \ge 0.75$, volume overlap, pocket burial, cofactor clearance, H-bond geometry, BUNS, RangeFinder compatibility, and PAINS/aggregator hygiene) with a **Strict Non-Hiding Dual-Ledger Provenance** contract (`raw_passes_all_21`, `was_desalted`, `rescue_possible` $\le 1.25\,\text{\AA}$, `rescue_rmsd_shift_a`, `provenance_tag`).
  3. **`bio_multipocket_posebusters_screen(protein_pdb, ligand_smiles_or_sdf, max_pockets DEFAULT 5)`:** Fuses top-$K$ pocket detection, 9D RangeFinder selectivity, and 21-check PoseBusters validation across all candidate pockets of a protein in a single SQL query, flagging the winning pocket with `is_best_pocket = TRUE` (`★`).
* **Reference Documentation & Benchmark:** [`docs/19_pocketvina_multipocket_and_posebusters.md`](docs/19_pocketvina_multipocket_and_posebusters.md) | [`scripts/benchmark_pocketvina_posebusters.py`](scripts/benchmark_pocketvina_posebusters.py) | [`experiments/16_pocketvina_multipocket_posebusters.ipynb`](experiments/16_pocketvina_multipocket_posebusters.ipynb)
* **SQL Example:**
  ```sql
  -- 1. Detect the Top 5 P2Rank-style SAS pockets of an enzyme and emit 9D RangeFinder + 16-byte Fastscan sketches:
  SELECT pocket_rank, pocket_id, p2rank_score, probability,
         center_x, center_y, center_z, sas_points_count, residue_count,
         mean_plddt, residue_ids[1:4] AS top_residues,
         octet_length(fastscan_16b) AS fastscan_bytes
  FROM bio_detect_top_pockets(protein_pdb, 5);

  -- 2. Audit a salted vendor drug SMILES with Strict Non-Hiding Dual-Ledger Provenance:
  --    (passes_all_21 = TRUE on the desalted drug, while raw_passes_all_21 = FALSE and was_desalted = TRUE)
  SELECT passes_all_21, raw_passes_all_21, was_desalted,
         rescue_possible, rescue_rmsd_shift_a, provenance_tag,
         passed_count, total_checks, failed_checks,
         min_vdw_ratio, max_ring_rmsd_a
  FROM bio_posebusters_audit_row(
      protein_pdb,
      '[Na+].O=C([O-])c1ccccc1O', -- Sodium salicylate salt -> desalted to Salicylic acid
      10.0, 10.0, 10.0
  );

  -- 3. Multi-Pocket Conditioning + 21-Check PoseBusters Screen for Diclofenac (Pipeline #6 feedstock):
  SELECT pocket_rank,
         CASE WHEN is_best_pocket THEN '★ ' || pocket_id ELSE pocket_id END AS pocket_badge,
         p2rank_score, pocket_probability, estimated_vina_kcal_mol,
         rangefinder_compatibility, posebusters_valid, posebusters_passed_checks
  FROM bio_multipocket_posebusters_screen(
      protein_pdb,
      'O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl',
      5
  );
  ```

### Layer 20: SIMD `EnzRef` Catalytic Role-Contrast Kinetics, Reaction SMILES Splitter & Disordered Microprotein `pSLiM` Window Pooler
Inspired by [**`EnzRef`**](https://doi.org/10.64898/2026.11.15.766718) (Ding et al., *bioRxiv* 2026 — Reference-Centered Multimodal Adaptation of Protein Language Models for Enzyme Turnover Prediction) and [**`pSLiM`** / **`SLiMNet`**](https://doi.org/10.64898/2026.05.04.722713) (Poretsky et al., *BMC Bioinformatics* 2026; McFee & Kim 2026), Stage-2.5 enzymatic ranking and synthetic metabolon design require distinguishing true forward catalytic turnover from dead-end substrate binding ($\mathbf{r}_{\text{prod}} = \mathbf{r}_{\text{react}}$) and indexing disordered microprotein Short Linear Motifs (`sORFs <= 110 aa`) without diluting them via full-sequence mean pooling or cleaving hydrophobic domain cores across synthetic fusion linkers.

* **The Engine:** `pg_bio` implements native ARM64 NEON / x86_64 AVX2 role-contrast and motif-pooling primitives directly inside PostgreSQL 18 (**`2,092,090` 1280-D `EnzRef` triplets/sec/core** and **`617,845` ZCA-whitened SLiM windows/sec/core**):
  1. **Reaction SMILES Agent-Stripping Splitter & Role Vectorizers (`bio_split_reaction_roles`, `bio_reaction_role_vectors`, `bio_batch_reaction_role_vectors`):** Decomposes `A.B>>C.D` or `A>agent>B` (stripping intermediate agents/solvents) and vectorizes dot-separated multi-component reactants and products into 1024-bit Morgan ECFP4 + 9D `RangeFinder` descriptors projected to any target dimension $D$ (e.g. `1280`).
  2. **Reference-Centered SIMD `EnzRef` Kinetics (`bio_enzref_kinetics`, `bio_enzref_flux_rank`, `bio_enzref_kinetics_g128`, `bio_enzref_flux_rank_g128`, `bio_enzref_kinetics_smiles`):** Computes centered deviations $\mathbf{d}_{\text{react}} = \mathbf{r}_{\text{react}} - \mathbf{v}_{\text{bg}}$ and $\mathbf{d}_{\text{prod}} = \mathbf{r}_{\text{prod}} - \mathbf{v}_{\text{bg}}$, symmetric mean response $\mathbf{m} = \frac{1}{2}(\mathbf{d}_{\text{react}} + \mathbf{d}_{\text{prod}})$, directed forward role contrast $\mathbf{c}_{\text{fwd}} = \mathbf{d}_{\text{prod}} - \mathbf{d}_{\text{react}}$, and product-retention `backpressure_risk` ($K_p$ stall) in a single fused SIMD pass over `float4[]` or inline 1,304-byte zero-TOAST `BioVecInt8G128` varlenas—guaranteeing exact antisymmetry under reaction reversal ($\text{FluxRank}(e, p, r) = -\text{FluxRank}(e, r, p)$) and exact zeroing on dead binders.
  3. **Disordered Microprotein `pSLiM` Window Pooler & Bivalent Tether Matcher (`bio_pslim_window_pool`, `bio_match_dual_slim_tether`, `bio_audit_fusion_junction_sanity`):** Extracts strided 5–12 aa motif windows gated by disorder/loop masks, applies in-memory corpus ZCA whitening + exact L2 unit normalization ($\|\mathbf{w}\|_2 = 1.0$, emitting 68-byte `vec16` vectors compatible with 2-level layered-centroid `fastscan` tables), and matches non-overlapping dual-SLiM metabolon couplers (`end_a + min_linker <= start_b`) while enforcing a hydrophobic-core linker junction sanity gate (`Single_Polypeptide_Fusion` vs. `Bicistronic_Operon_Required`).
* **Reference Documentation:** [`docs/20_enzref_role_contrast_and_pslim_tethers.md`](docs/20_enzref_role_contrast_and_pslim_tethers.md)
* **SQL Example:**
  ```sql
  -- 1. Evaluate EnzRef catalytic flux rank and backpressure risk directly from reaction SMILES:
  SELECT reactants_smiles, products_smiles,
         forward_flux_score, backpressure_risk,
         role_contrast_norm, mean_catalytic_response,
         cosine_alignment, norm_ratio_react_over_prod
  FROM bio_enzref_kinetics_smiles(
      ARRAY(SELECT sin(g::float8)::float4 FROM generate_series(1, 1280) g),
      'CC(=O)Oc1ccccc1C(=O)O.O>[H+]>CC(=O)O.O=C(O)c1ccccc1O'
  );

  -- 2. Match a disordered microprotein as a Dual-SLiM Metabolon Tether with Linker Junction Sanity:
  SELECT start_a, end_a, subseq_a, sim_a, l2_centroid_a,
         start_b, end_b, subseq_b, sim_b, l2_centroid_b,
         linker_len, linker_gravy, junction_sanity_pass,
         recommended_operon_mode, dual_tether_score
  FROM bio_match_dual_slim_tether(
      'MKKLVLVAVILGLAVAGSAQAPVGRDTVPDPLRERAQELGVKVTVEGPGGGSGGGGSLLPAIVHIEQDKTKREEVD',
      ARRAY[0.45, -0.20, 0.30, 0.15, -0.40, 0.25, 0.10, -0.35, 0.20, 0.10, -0.15, 0.30, 0.05, -0.25, 0.40, 0.12]::real[],
      ARRAY[-0.30, 0.50, -0.15, 0.40, 0.20, -0.35, 0.25, 0.10, -0.45, 0.30, 0.15, -0.20, 0.35, 0.10, -0.25, 0.18]::real[],
      8, 2, 6, NULL, NULL
  );
  ```

### Layer 21: Directional `ipSAE` & $5 \times 5$ `MiniPAE` Co-Folding Confidence Gate
Inspired by [**`ipSAE`**](https://doi.org/10.1101/2025.02.10.637595) (Dunbrack, *bioRxiv* 2025/2026 — *"Rēs ipSAE loquuntur: What's wrong with AlphaFold's ipTM score and how to fix it"*) and [**`MiniPAE`**](https://doi.org/10.1093/bib/bbaf501) (Veinstein et al., *Briefings in Bioinformatics* 2025; Omidi et al., 2024), scoring multi-chain complexes or synthetic metabolon tethers ($M_1 \rightarrow M_2 \rightarrow M_3$) with global `ipTM` causes severe false negatives whenever one partner contains an Intrinsically Disordered Region (IDR) or flexible linker (`GGGSGGGGS`). Because global `ipTM` sets $d_0(L_A + L_B)$ using total sequence length and averages over all partner residues, an 80-aa microprotein with a 10-aa docked SLiM ($\text{PAE} = 0.45\,\text{\AA}$) and a 70-aa free disordered tail ($\text{PAE} = 28.5\,\text{\AA}$) drops to $\text{ipTM} \approx 0.125$.

* **The Engine:** `pg_bio` implements Dunbrack's contact-restricted directional `ipSAE` ($d_0(N_i) = 1.0$ for $N_i < 27$ else $1.24(N_i - 15)^{1/3} - 1.8$, evaluated strictly over qualifying contacts $\text{PAE}_{i,j} < 10.0\,\text{\AA}$) paired with Veinstein's 2D-prefix-sum $5 \times 5$ sliding-window `MiniPAE` and a `< 2,000 B` zero-TOAST `"PGPA"` quantized CSR `bytea` container (`bio_ipsae_minipae_audit`, `bio_ipsae_minipae_from_pae_flat`, `bio_pack_sparse_pae_csr`, `bio_peek_pae_csr_header`).
* **Reference Documentation:** [`docs/21_ipsae_minipae_cofolding_gate.md`](docs/21_ipsae_minipae_cofolding_gate.md)
* **SQL Example:**
  ```sql
  -- Audit a 2-chain receptor + disordered SLiM tether complex and rescue it from global ipTM IDR dilution:
  SELECT chain_a_id, chain_b_id,
         ipsae_a_to_b, ipsae_b_to_a, ipsae_min, ipsae_max,
         iptm_global_diluted, idr_dilution_rescue_delta,
         minipae_5x5_a, best_window_a_start, best_window_b_start,
         passes_cofold_gate, gate_verdict, packed_csr_bytes
  FROM bio_ipsae_minipae_audit(complex_pdb_text, NULL, 10.0, 3.5);
  ```

---
## 🛠️ Scale & Optimization Infrastructure

`pg_bio` isn't just mathematically clever; it is heavily optimized for multi-million scale deployment on commercial hardware.

### 1. Inline Float16 Vector Compression (`halfvec(1280) STORAGE PLAIN`)
High-dimensional 1280-D ESM-2 vectors (`2,568 bytes`) exceed PostgreSQL's 2 KB TOAST threshold by default, spawning 6.11 GB of out-of-line TOAST chunks across 2.34 million proteins. `pg_bio` pairs `halfvec(1280)` (16-bit floating-point) with `STORAGE PLAIN` to eliminate external TOAST tables and random chunk lookups.

### 2. Tiered Vector & Multi-Index Routing (`pg_vectorsearch` PRISM + Layered Centroids + Register Fastscan)
* **`pg_vectorsearch` (`USING prism` — Primary 1M+ Proteome Engine):** Uses 1-bit RaBitQ quantization with theoretical error bounds, SOAR spill-orthogonal secondary cluster assignment, and hierarchical centroid routing. Across **2,340,766** proteins, PRISM shrinks the vector index by **9.72× (`8,772 MB` $\rightarrow$ `902 MB`)**, builds in **3.08 minutes** (vs. 27 hours for HNSW), and accelerates homology queries by **7.5×–11.3×** (`1.22 ms` to `26 ms`).
* **100% Lossless Layered Centroids & Register-Direct Fastscan (`PR #14`):** Ports `pg_vectorsearch`'s hierarchical centroid routing, SOAR boundary assignment, and register-direct Fastscan upper bounds to `pg_bio`'s other three biological index families:
  * **3D Morton Octree (`z_order_layered_sphere_ranges`):** `74.05×` fewer atoms scanned, **`10.25×` faster (`0.96 ms`)**, `100.0%` exact recall.
  * **1024-Bit ECFP4 Tanimoto (`ecfp4_fastscan_ub_raw`):** Prunes **`99.89%`** of molecules in 64-bit CPU registers, **`31.19×` faster (`16.98 ms`)** at `100.0%` exact recall (or **`158.15×` faster at `3.35 ms`** with SOAR dual-centroid routing).
  * **9D RangeFinder Active-Site Pockets (`pocket_quantize_fastscan` & `pocket_fastscan_iou_ub_raw`):** Compresses 350-byte pockets by **22× into two 64-bit CPU registers (`16 bytes`)**, running **`13.24×` faster (`7.79 ms`)** at `100.0%` exact Recall@10 (or **`20.42×` faster at `5.05 ms`** with 37 balanced Level-2 centroid envelopes + SOAR).

---

## 💻 Python SDK & Tooling Integration

You don't need to write raw SQL to leverage this power.

* **`pgbio-py`:** A seamless Python SDK that connects directly to the Rust engine.
* **Jupyter Integration:** Full compatibility with Pandas and Matplotlib. View our complete tutorial notebook demonstrating all four layers in action: `https://github.com/jonatas/teaflon/blob/main/pg_bio_multiomics_pipeline/teaflon_multiomics_tutorial.ipynb`.

* **Dry Lab Tutorial:** Learn how to use these new structural and sequence analysis tools natively in SQL: [examples/bio_dry_lab_tutorial.ipynb](examples/bio_dry_lab_tutorial.ipynb)


```python
from pgbio import PgBioClient

client = PgBioClient("postgresql://localhost:28818/bio_demo")

# Query the Joint Vector Space natively from Python
dna_matches = client.find_joint_enhancers(protein_id="P0DPB7")
print(f"Strongest regulatory match: {dna_matches[0].locus}")
```

## 🐳 Quickstart

```bash
# 1. Clone the repository
git clone https://github.com/jonatas/pg_bio.git
cd pg_bio

# 2. Spin up the Rust-optimized Postgres Engine
docker-compose up -d

# 3. Compile the Database and Interactive Notebooks
uv run jupyter nbconvert --execute https://github.com/jonatas/teaflon/blob/main/pg_bio_multiomics_pipeline/teaflon_multiomics_tutorial.ipynb

# 4. Run the Declarative Multi-Feedstock Bioprocess Compiler & Autonomous Swarm UI (jonatas/bio-pipelines)
git clone https://github.com/jonatas/bio-pipelines.git ../bio-pipelines
cd ../bio-pipelines && uv sync
uv run python base/compiler.py
uv run python base/swarm_web.py
```

*(For an extended deep-dive into the architectural theory, read [`docs/MULTIOMICS_ARCHITECTURE.md`](docs/MULTIOMICS_ARCHITECTURE.md))*

---

## 📊 Master Catalog of `pg_bio` Correlation Types & Performance Architecture

Below is the complete reference of correlation engines built into [`src/lib.rs`](src/lib.rs), how each correlation works intuitively, how `pg_bio` executes it natively inside PostgreSQL, why it outperforms Python pipelines, and where to find the official documentation and runnable examples.

| # | Correlation Type | Simple Explanation (How It Works) | Technical `pg_bio` Approach (`src/lib.rs` & SQL) | Why `pg_bio` (Rust in Postgres) Is Faster Than Python | Official Docs & Examples |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | **Twilight-Zone Homology Index**<br>*(3D Convergence vs. Sequence Divergence)* | Finds pairs of proteins that fold into the **exact same 3D shape** even though their amino acid sequences have mutated heavily (`15%–75%` similarity), filtering out both 100% clones and random k-mer noise. | Combines HNSW cosine distance (`<=>`), [`smith_waterman_score`](src/lib.rs#L1964-L1995), and [`sequence_entropy`](src/lib.rs#L1780-L1805) in a bell-curve filter:<br>`(1 - d_cos)^2 * [4 * s * (1 - s)] * (L_min / L_max)^2 * min(H / 4.0, 1.0)`. | **Zero-IPC & $O(N)$ Two-Row DP:** [`smith_waterman_score`](src/lib.rs#L1964-L1995) runs in compiled Rust using two flat `Vec<i32>` rolling buffers directly on Postgres `TEXT` pointers—no $M \times N$ Python matrix allocation and no socket serialization of millions of sequences. | 📖 [`01_vector_homology.md`](docs/01_vector_homology.md)<br>📖 [`08_native_alignment.md`](docs/08_native_alignment.md)<br>🐝 [`jonatas/bio-pipelines` (`swarm_discovery.py`)](https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_discovery.py)<br>📓 [`07_native_alignment.ipynb`](experiments/07_native_alignment.ipynb) |
| **2** | **Extremophile Niche Shift Index**<br>*(Electrostatic & Hydropathy Remodeling: $\Delta\text{pI}$ / $\Delta\text{GRAVY}$)* | Measures how a protein rewires its surface charge (acidic vs. basic) and water-repellency to survive extreme habitats (deep-sea brine, acid vents, soda lakes) while keeping the same 3D core. | Computes [`isoelectric_point`](src/lib.rs#L1828-L1905) via 100-step bisection on Henderson-Hasselbalch $\text{p}K_a$ curves and [`gravy_score`](src/lib.rs#L1806-L1827) (Kyte-Doolittle hydropathy) inline with `<=>`. | **Single-Pass Byte-Slice Histograms:** Counts amino acid frequencies in a single CPU register pass (`seq.bytes()`) and solves $\text{pI}$ in stack-allocated `f64` math (~200 ns/protein vs. ~40 µs in BioPython). | 📖 [`11_sequence_analysis.md`](docs/11_sequence_analysis.md)<br>📖 [`04_native_chemistry.md`](docs/04_native_chemistry.md)<br>🐝 [`jonatas/bio-pipelines` (`swarm_web.py`)](https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_web.py) |
| **3** | **Layer 8 Interface Complementarity**<br>*(Protein–Protein / Complex Docking Compatibility)* | Checks whether two proteins can physically dock together—like a lock and key—where a positively charged or convex patch on one matches a negatively charged or concave pocket on the other. | [`interface_complementarity`](src/lib.rs#L3041-L3070) splits embeddings into **Steric/Shape** (first half: anti-correlated $p_i \cdot (-q_i)$) and **Electrostatic/Hydrophobic** (second half: complementary product), mapped via `tanh`. | **Fused SIMD Vector Loop:** Operates directly on Postgres `real[]` / `vector` slices in L1 cache without copying arrays across Python/NumPy FFI boundaries. | 📖 [`MULTIOMICS_ARCHITECTURE.md`](docs/MULTIOMICS_ARCHITECTURE.md)<br>🐝 [`swarm_discovery/README.md`](examples/swarm_discovery/README.md)<br>📓 [`09_huggingface_esm2_discovery.ipynb`](experiments/09_huggingface_esm2_discovery.ipynb) |
| **4** | **Active-Site Motif & Metal-Cage Gate**<br>*(Catalytic Residue Conservation)* | Checks if a newly discovered dark protein still has the exact chemical "teeth" (e.g., a `C-x(2)-C-x(12)-H` zinc finger or `[CH]-x-[CH]` iron-sulfur cage) needed to run a chemical reaction. | [`prosite_match`](src/lib.rs#L276-L365) compiles Prosite syntax (`[LIVM]-x(2,4)-{P}-C`) into a native Rust finite-state token matcher; [`has_coordinated_water_shell`](src/lib.rs#L2710-L2786) checks metal-water geometry. | **Zero-Regex-Engine Overhead:** Avoids Python `re` backtracking by stepping through fixed `PrositeToken` enums (`Exact`, `AnyOf`, `NotOf`, `Wildcard(min,max)`) directly on ASCII bytes. | 📖 [`04_native_chemistry.md`](docs/04_native_chemistry.md)<br>📓 [`10_metalloprotease_theozyme_design.ipynb`](experiments/10_metalloprotease_theozyme_design.ipynb)<br>📓 [`bio_dry_lab_tutorial.ipynb`](examples/bio_dry_lab_tutorial.ipynb) |
| **5** | **3D Z-Order (Morton) Pocket & Theozyme Correlation**<br>*(Spatial Active-Site Superposition)* | Compares the 3D arrangement of catalytic atoms inside an enzyme's pocket against known catalytic templates ("theozymes"), regardless of sequence order. | [`z_order_encode`](src/lib.rs#L34-L46) interleaves 21-bit $(x,y,z)$ coordinates via bit-magic (`split_by_3`) into a 64-bit `BIGINT` B-tree key; [`kabsch_rmsd`](src/lib.rs#L2823-L2935) & [`theozyme_rmsd`](src/lib.rs#L3002-L3040) run Jacobi SVD. | **3D-to-1D B-Tree Indexing + Stack $3\times 3$ SVD:** Turns a 3D bounding-box scan into a 1D B-tree range scan (`i64`), and runs a zero-allocation $3\times 3$ Jacobi SVD in Rust (no LAPACK/SciPy call overhead). | 📖 [`03_spatial_indexing.md`](docs/03_spatial_indexing.md)<br>📖 [`09_pdb_parsing.md`](docs/09_pdb_parsing.md)<br>📓 [`08_pdb_parsing.ipynb`](experiments/08_pdb_parsing.ipynb)<br>📓 [`10_metalloprotease_theozyme_design.ipynb`](experiments/10_metalloprotease_theozyme_design.ipynb) |
| **6** | **Substrate & Catalytic Warhead Correlation**<br>*(Enzyme $\leftrightarrow$ Small-Molecule Pairing)* | Verifies that an enzyme family's chemical substrate contains the reactive "warhead" bond that the enzyme cleaves or oxidizes, even if drawn in different tautomer, salt, or Kekulé forms. | [`parse_smiles`](src/lib.rs#L1596-L1623) runs DFS cycle detection, Hückel $4n+2$ $\pi$-electron aromaticity ([`pgbio_perceive_aromaticity`](src/lib.rs#L1272-L1338)), 64-bit Weisfeiler-Lehman hashing ([`molecule_wl_hash`](src/lib.rs#L1624-L1630)), VF2 subgraph isomorphism (`@>`, [`pgbio_vf2_is_subgraph`](src/lib.rs#L1437-L1595)), and ECFP4 Tanimoto (`%`). | **Bit-Packed Graph & Native Operators (`@>`, `%`):** Stores molecules as compact `PgBioAtom`/`PgBioBond` structs inside Postgres pages and executes recursive VF2 subgraph matching and 64-bit `POPCNT` Tanimoto (**1,513.7x faster than Python pull**). | 📖 [`14_aromaticity_wl_vf2_cheminformatics.md`](docs/14_aromaticity_wl_vf2_cheminformatics.md)<br>📖 [`05_cheminformatics.md`](docs/05_cheminformatics.md)<br>📖 [`06_smiles_rules.md`](docs/06_smiles_rules.md)<br>📓 [`05_cheminformatics_tutorial.ipynb`](experiments/05_cheminformatics_tutorial.ipynb) |
| **7** | **Sparse CSR Allosteric & Cryptic Pocket Correlation**<br>*(Long-Range Neural Attention & Energy Networks)* | Tracks how binding a molecule at a distant surface site sends a ripple through the protein's internal contact network to open or close the active site. | [`create_sparse_map`](src/lib.rs#L142-L157), [`discover_cryptic_hotspots`](src/lib.rs#L3461-L3504), [`allosteric_mechanism_index`](src/lib.rs#L3505-L3526), and [`sparse_cross_reactivity_score`](src/lib.rs#L3106-L3148) use Compressed Sparse Row (CSR) graphs. | **$O(\text{nnz})$ CSR Traversal vs. $O(L^2)$ Dense Matrices:** Stores only non-zero contacts (`row_ptr`, `col_idx`, `values`), cutting memory by >95% and scanning row slices in nanoseconds (**31.14x faster than Python/NumPy**). | 📖 [`02_attention_traversal.md`](docs/02_attention_traversal.md)<br>📖 [`12_mmpbsa_free_energy_auditing.md`](docs/12_mmpbsa_free_energy_auditing.md)<br>📓 [`11_5lox_triterpenoid_mmpbsa_validation.ipynb`](experiments/11_5lox_triterpenoid_mmpbsa_validation.ipynb) |
| **8** | **Globular Compactness & Hydrophobic Core Correlation**<br>*(Foldability vs. Disordered Spaghetti)* | Distinguishes proteins that fold into a tight, stable 3D ball with a water-free core from floppy, unfolded "spaghetti" predictions. | [`evaluate_globular_compactness`](src/lib.rs#L3773-L3908) and [`bio_rasa`](src/lib.rs#L4222-L4304) compute Radius of Gyration ($R_g$ vs. Flory $R_0 = 2.2 N^{0.38}$) and 64/92-point Fibonacci-sphere Shrake-Rupley SASA/RASA. | **Spatial Voxel Grid Hashing for SASA:** Bins atoms into $5\text{ Å}$ 3D hash cells (`HashMap<(i32,i32,i32), Vec<usize>>`) in [`compute_atomic_sasa`](src/lib.rs#L4142-L4221), turning $O(N^2)$ sphere overlap checks into $O(N)$ local neighborhood lookups. | 📖 [`13_wetlab_synbio_and_teaflon_validation.md`](docs/13_wetlab_synbio_and_teaflon_validation.md)<br>📖 [`15_chaperone_protein_design_workflows.md`](docs/15_chaperone_protein_design_workflows.md)<br>📖 [`07_on_the_fly_folding.md`](docs/07_on_the_fly_folding.md) |
| **9** | **Hydrogen-Bond Network & Buried Unsatisfied Polar (BUNS) Audit**<br>*(Thermodynamic Pocket Stability)* | Finds polar atoms (N, O) trapped deep inside a water-free protein core without a hydrogen-bond partner—a major cause of misfolding or unstable designs. | [`bio_hbond_network`](src/lib.rs#L4366-L4480) and [`audit_hydrogen_bonds`](src/lib.rs#L4481-L4682) combine Shrake-Rupley SASA burial ($\text{SASA} < 0.5\text{ Å}^2$) with donor–acceptor distance ($\le 3.5\text{ Å}$) and antecedent angle ($\ge 100^\circ$). | **Fused Single-Pass SASA + Angle Geometry:** Shares the precomputed atomic SASA buffer with the H-bond angular validator in a single Rust call ([`evaluate_chaperone_block`](src/lib.rs#L4807-L4853)), eliminating redundant PDB parsing. | 📖 [`15_chaperone_protein_design_workflows.md`](docs/15_chaperone_protein_design_workflows.md)<br>📓 [`13_chaperone_protein_design_workflows.ipynb`](experiments/13_chaperone_protein_design_workflows.ipynb) |
| **10** | **Conformational Preorganization Correlation**<br>*(MD Trajectory & Ensemble Gating)* | Measures what fraction of the time a flexible enzyme naturally spends in its "ready-to-react" state before the substrate even arrives. | [`create_preorganization_profile`](src/lib.rs#L3167-L3229), [`ensemble_fraction_below`](src/lib.rs#L3230-L3258), and [`ensemble_fraction_between`](src/lib.rs#L3259-L3275) compress millions of MD trajectory frames into quantized histogram bins inside Postgres. | **$O(B)$ Quantized Histogram Interpolation:** Replaces loading gigabytes of `.xtc`/`.dcd` trajectory files in Python with a 20-bin linear interpolation over `PreorganizationProfile` in <1 µs per row. | 📓 [`10_metalloprotease_theozyme_design.ipynb`](experiments/10_metalloprotease_theozyme_design.ipynb)<br>📖 [`12_mmpbsa_free_energy_auditing.md`](docs/12_mmpbsa_free_energy_auditing.md) |
| **11** | **MS/MS Proteomics Isotope & Fragment Correlation**<br>*(Experimental Validation of Dark Proteins)* | Matches predicted dark-proteome sequences against real-world mass spectrometry peaks ($b/y$ ions and $^{13}\text{C}$ isotope envelopes) to prove the protein is actually expressed in cells. | [`generate_theoretical_fragments`](src/lib.rs#L2176-L2252), [`is_mass_within_ppm`](src/lib.rs#L2161-L2175), and [`isotope_fit_score`](src/lib.rs#L2253-L2323) compute Poisson $^{13}\text{C}$ envelopes ($e^{-\lambda}\lambda^k/k!$) and cosine similarity against `Spectrum`. | **Set-Returning Rust Iterators (`TableIterator`):** Streams $b/y/c/z$ fragments and `.mgf` spectra ([`parse_mgf`](src/lib.rs#L2386-L2466)) directly into the PostgreSQL executor pipeline with zero temporary Python list materialization. | 📖 [`10_mass_spectrometry.md`](docs/10_mass_spectrometry.md)<br>📓 [`bio_dry_lab_tutorial.ipynb`](examples/bio_dry_lab_tutorial.ipynb) |
| **12** | **Genomic Safe-Harbor & SynBio Host Expression Correlation**<br>*(Discovery $\rightarrow$ Wet-Lab DNA Construct)* | Links a discovered protein to genomic variants (`.vcf`), safe-harbor insertion loci, and optimal DNA codons/GC%/RBS for expression in *E. coli*, *B. subtilis*, or *H. salinarum*. | [`parse_vcf`](src/lib.rs#L2467-L2538), [`bio_reverse_translate`](src/lib.rs#L3634-L3649), [`bio_codon_adaptation_index`](src/lib.rs#L3650-L3704), and [`bio_design_expression_cassette`](src/lib.rs#L3705-L3772) emit complete Golden-Gate/Gibson expression cassettes in SQL. | **Compile-Time Lookup Tables (`match (host, aa)`):** Translates and scores CAI in a single zero-allocation loop over static Rust `(&'static str, f64)` codon tables. | 📖 [`07_genomics_io.md`](docs/07_genomics_io.md)<br>📖 [`13_wetlab_synbio_and_teaflon_validation.md`](docs/13_wetlab_synbio_and_teaflon_validation.md)<br>📓 [`06_genomics_io_tutorial.ipynb`](experiments/06_genomics_io_tutorial.ipynb) |
| **13** | **Layer 14 Industrial Bioprocess Cascade Kinetics (RK4), Pareto Domain Scoring & Biomass Stoichiometry**<br>*(Multi-Feedstock Upcycling & Equilibrium Sinks)* | Ranks `(Anchor × Scissor)` domain pairs across hydropathy, pI, and sequence entropy, integrates 48-hour two-stage enzymatic cascades with competitive product inhibition ($K_i$) and surface adsorption boost via 4th-order Runge-Kutta (RK4), and solves 3-enzyme biomass stoichiometry ($w_1+w_2+w_3=1.0$). | [`pareto_score_domain_pair`](src/lib.rs#L5482-L5501), [`simulate_cascade_kinetics`](src/lib.rs#L5343-L5410), [`optimize_consortium_stoichiometry`](src/lib.rs#L5414-L5478), and [`translate_and_merge_pdb_chains`](src/lib.rs#L5506-L5558) execute natively inside PostgreSQL 18. | **In-Database RK4 & Simplex Grid in < 1 ms:** Eliminates shipping thousands of candidate pairs and ODE time-series arrays to Python/SciPy; evaluates 12 full 97-step RK4 trajectories in ~24 ms directly inside SQL. | 🚀 [**`jonatas/bio-pipelines`**](https://github.com/jonatas/bio-pipelines)<br>🧪 [**`jonatas/teaflon`**](https://github.com/jonatas/teaflon) |
| **14** | **Layer 15 RangeFinder Pocket-Interval GiST & B-Tree Virtual Screening**<br>*(100M+ Ligand Filtering & Reverse Target Fishing, `arXiv:2611.03379v1`)* | Predicts 9 compatible ligand property intervals ($\hat{I}_{pj} = [\hat{\mu}_{pj} - 2|\hat{\sigma}_{pj}|, \hat{\mu}_{pj} + 2|\hat{\sigma}_{pj}|]$) from a $10\text{ \AA}$ protein pocket point cloud to filter 100M+ molecules via B-Tree or reverse-screen 20,000+ pockets via Multi-Column `numrange` GiST before touching a GPU. | [`predict_pocket_ranges`](src/lib.rs#L6057-L6269), [`compute_ligand_properties_row`](src/lib.rs#L6497-L6525), [`pocket_numranges`](src/lib.rs#L6590-L6619), [`pocket_matches_ligand_8prop`](src/lib.rs#L6410-L6427), and [`evaluate_rangefinder_screening`](src/lib.rs#L6433-L6492). | **0.167 ms/pocket Inference + 94.9x B-Tree & 16.6x GiST Speedup:** Screens 250,000 ligands in **3.78 ms** (94.9x faster than Python) and reverse-fishes 20,000 pockets via 8-column `numrange` GiST (`@>`) in **0.68 ms**. | 📖 [`RANGEFINDER_GIST_INDEXING.md`](docs/RANGEFINDER_GIST_INDEXING.md)<br>🐝 [`rangefinder_pocket_screening/`](examples/rangefinder_pocket_screening/)<br>⚡ [`benchmark_rangefinder_gist.py`](scripts/benchmark_rangefinder_gist.py)<br>📓 [`rangefinder_pocket_gist_screening.ipynb`](experiments/rangefinder_pocket_gist_screening.ipynb) |
| **15** | **Anisotropy Whitening (`Know Thy Neighbour`), Inline `BioVecInt8G128` Varlena, Structure/PAINS Hygiene & Stage-2.5 (`EnzRef` / `STORM`) Re-Rankers** | Eliminates `ESM-2`/`ProstT5` SwiGLU cone anisotropy (`0.9963 -> -0.0151`) to tighten RaBitQ error bounds (`eps_bound`) by **2.78×**, keeps 1280-D Stage-2 reranking vectors 100% inline on the Postgres heap (`1,304 B < 2,000 B` TOAST threshold), repairs missing PDB/mmCIF atoms via Kabsch SVD + `< 2.2 Å` `CYX` disulfides + `HID`/`HIE`/`HIP` tautomers, filters PAINS/aggregators at **11,995 mols/s**, and re-ranks enzyme/ligand candidates via `EnzRef` SIMD role-contrast and `STORM` 64-hot `ID128` + Tanimoto cliff-preserving smoothing. | [`bio_fit_whitening_stats`](src/lib.rs#L8372-L8380), [`bio_whiten_embedding`](src/lib.rs#L8392-L8465), [`BioVecInt8G128`](src/lib.rs#L8657-L8963), [`bio_sanitize_structure`](src/lib.rs#L9615-L9638), [`bio_is_pains_or_aggregator`](src/lib.rs#L9640-L9935), [`bio_enzref_role_contrast`](src/lib.rs#L9941-L10075), [`bio_id128_hash`](src/lib.rs#L10077-L10122), [`bio_smooth_assay_series`](src/lib.rs#L10124-L10265). | **2.78× Tighter RaBitQ Bound, Zero TOAST Reads (`1,304 B` Inline Varlena), 2.00× Faster SIMD `i32` Reranking & 52,970 rows/s `EnzRef` SIMD:** Fuses ARM64 NEON / x86_64 AVX2 whitening, 128-channel group `i8+f16` `i32` accumulation, and 3-vector role-contrast inside Postgres buffer pages. | 📖 [`18_whitening_inline_reranking_and_hygiene.md`](docs/18_whitening_inline_reranking_and_hygiene.md)<br>📖 [`17_layered_centroids_multi_index_experiments.md`](docs/17_layered_centroids_multi_index_experiments.md) |
| **16** | **Layer 19 PocketVina Multi-Pocket SAS Conditioning & 21-Check PoseBusters Physical Validity Gate**<br>*(Multi-Pocket `pk1..pk5` Concavity + 21-Check Steric/Chemical Audit)* | Detects and ranks the **Top 5 binding pockets** per protein using 24-ray Fibonacci Connolly SAS probes + $5.0\text{ \AA}$ voxel hashing (`P2Rank`-inspired), emits 9D `RangeFinderPocket` intervals + 16-byte register Fastscan sketches per pocket, and audits any ligand pose (SDF, PDBQT, or desalted SMILES) across all **21 PoseBusters physical, geometric, salt-hygiene, and steric clash checks** inside SQL. | [`bio_detect_top_pockets`](src/lib.rs#L12047-L12545), [`bio_posebusters_audit`](src/lib.rs#L13159-L13718), [`bio_posebusters_audit_row`](src/lib.rs#L13719-L13783), and [`bio_multipocket_posebusters_screen`](src/lib.rs#L13784-L13915). | **Zero Python/RDKit IPC + Fused Salt Rescue & $3\times 3$ Jacobi SVD Ring Planarity:** Automatically rescues salted vendor SMILES (`[Na+].O=C([O-])...`), rejects bare salt counterions (`[Na+]`) and PAINS motifs, and audits VdW clashes ($d/r_{\text{vdW}} \ge 0.75$) and ring flatness ($\le 0.25\text{ \AA}$) in single-digit milliseconds in SQL. | 📖 [`19_pocketvina_multipocket_and_posebusters.md`](docs/19_pocketvina_multipocket_and_posebusters.md)<br>⚡ [`benchmark_pocketvina_posebusters.py`](scripts/benchmark_pocketvina_posebusters.py)<br>📓 [`16_pocketvina_multipocket_posebusters.ipynb`](experiments/16_pocketvina_multipocket_posebusters.ipynb) |
| **17** | **Layer 20 SIMD `EnzRef` Catalytic Role-Contrast Kinetics & Disordered Microprotein `pSLiM` Window Pooler**<br>*(Exact Antisymmetric $k_{\text{cat}}$ Flux, Backpressure & Dual-SLiM Metabolon Tethers)* | Decomposes reaction SMILES (`A.B>>C.D`, stripping `>agent>`), computes reference-centered substrate/product deviations ($\mathbf{d}_{\text{react}} = \mathbf{r} - \mathbf{e}$, $\mathbf{c}_{\text{fwd}} = \mathbf{p} - \mathbf{r}$) with exact antisymmetry and product-retention backpressure ($K_p$), and extracts 16-D ZCA-whitened unit-L2 SLiM windows (`vec16`) to match bivalent microprotein tethers with hydrophobic-core linker junction sanity gating. | [`bio_split_reaction_roles`](src/lib.rs#L1749-L1761), [`bio_reaction_role_vectors`](src/lib.rs#L1834-L1868), [`bio_enzref_kinetics`](src/lib.rs#L10933-L10948), [`bio_enzref_flux_rank`](src/lib.rs#L10954-L10960), [`bio_enzref_kinetics_g128`](src/lib.rs#L10964-L10971), [`bio_enzref_kinetics_smiles`](src/lib.rs#L10987-L11022), [`bio_pslim_window_pool`](src/lib.rs#L11188-L11298), [`bio_match_dual_slim_tether`](src/lib.rs#L11310-L11515), [`bio_audit_fusion_junction_sanity`](src/lib.rs#L11522-L11575). | **2,092,090 `EnzRef` 1280-D Triplets/sec/core (`0.48 µs/row`) & 617,845 ZCA-Whitened SLiM Windows/sec/core:** Computes all 4 directional moments in a single streaming NEON/AVX2 pass over `float4[]` or inline 1,304-byte zero-TOAST `BioVecInt8G128` varlenas without heap allocations. | 📖 [`20_enzref_role_contrast_and_pslim_tethers.md`](docs/20_enzref_role_contrast_and_pslim_tethers.md)<br>📖 [`18_whitening_inline_reranking_and_hygiene.md`](docs/18_whitening_inline_reranking_and_hygiene.md) |
| **18** | **Layer 21 Directional `ipSAE` & $5 \times 5$ `MiniPAE` Co-Folding Confidence Gate**<br>*(IDR Dilution Rescue for Co-Folded Complexes & Disordered SLiM Tethers)* | Rescues disordered microprotein SLiM tethers (`sORFs <= 110 aa`) and multi-domain chimeras from AlphaFold/Boltz-2 global `ipTM` IDR dilution by evaluating Dunbrack's contact-restricted $d_0(N_i)$ directional `ipSAE` ($\text{PAE}_{i,j} < 10.0\,\text{\AA}$) and Veinstein's $5 \times 5$ contiguous `MiniPAE` window over zero-TOAST `< 2,000 B` `"PGPA"` quantized CSR `bytea` containers. | [`bio_ipsae_minipae_audit`](src/lib.rs#L14585-L14710), [`bio_ipsae_minipae_from_pae_flat`](src/lib.rs#L14515-L14582), [`bio_pack_sparse_pae_csr`](src/lib.rs#L14460-L14483), [`bio_peek_pae_csr_header`](src/lib.rs#L14488-L14512). | **Zero-TOAST `< 2,000 B` `"PGPA"` CSR Storage + $O(L_A L_B)$ 2D-Prefix `MiniPAE` + $O(1)$ 16-Byte Header Peek:** Compresses 1.44 MB dense PAE matrices by **>750×** into inline heap tuples and evaluates directional `ipSAE` + $5 \times 5$ `MiniPAE` in microseconds without Python JSON parsing. | 📖 [`21_ipsae_minipae_cofolding_gate.md`](docs/21_ipsae_minipae_cofolding_gate.md)<br>⚡ [`benchmark_ipsae_minipae_gate.py`](scripts/benchmark_ipsae_minipae_gate.py)<br>📓 [`18_ipsae_minipae_cofolding_gate.ipynb`](experiments/18_ipsae_minipae_cofolding_gate.ipynb) |

---

### ⚡ Why `pg_bio` Is 10x–1,000x Faster Than a Python Pipeline

When executing a multi-modal discovery query (such as [`jonatas/bio-pipelines` (`base/swarm_discovery.py`)](https://github.com/jonatas/bio-pipelines/blob/main/base/swarm_discovery.py)), `pg_bio` fuses all filtering layers inside a single PostgreSQL backend process:

1. **Zero Data Movement (Compute Pushed to Shared-Buffer Pages):**
   * **Python Stack (`psycopg` + `BioPython` + `SciPy` + `RDKit`):** Must serialize sequences, 1280-D vectors, or 3,000-line PDB strings over a socket from PostgreSQL, allocate Python `str`/`list`/`numpy.ndarray` objects, and garbage-collect them.
   * **`pg_bio` Engine:** Reads PostgreSQL shared-buffer memory pages directly in Rust via `pgrx` pointers. Only the final winning rows ever cross the wire to the client.
2. **Query-Planner Short-Circuiting (HNSW Index $\rightarrow$ B-Tree Anti-Join $\rightarrow$ DP Verification):**
   * PostgreSQL uses the **HNSW vector index** (`<=>`) to narrow **2,340,573 proteins down to 15 candidates in ~0.4 ms**, applies the `NOT EXISTS` B-tree anti-join, and executes [`smith_waterman_score`](src/lib.rs#L1964-L1995) and [`interface_complementarity`](src/lib.rs#L3041-L3070) on **only those 5–15 survivors**.
3. **LLVM-Vectorized Native Rust Loops vs. Python Interpreter Overhead:**
   * Inner loops for [`smith_waterman_score`](src/lib.rs#L1964-L1995) ($O(M \cdot N)$ dynamic programming), [`compute_atomic_sasa`](src/lib.rs#L4142-L4221) (Fibonacci sphere integration), and [`pgbio_vf2_is_subgraph`](src/lib.rs#L1437-L1595) (recursive subgraph isomorphism) compile with `rustc -C opt-level=3`, auto-vectorizing into ARM64 NEON / x86 AVX2 SIMD instructions with zero GIL contention.
4. **Domain-Specific Spatial & Graph Reductions:**
   * **3D Atomic Coords $\rightarrow$ 64-bit Morton Integer ([`z_order_encode`](src/lib.rs#L34-L46)):** Replaces 3D Euclidean scans with 1D B-tree index lookups ([`docs/03_spatial_indexing.md`](docs/03_spatial_indexing.md)).
   * **Dense $L \times L$ Contact Matrix $\rightarrow$ CSR Graph ([`create_sparse_map`](src/lib.rs#L142-L157)):** Reduces $O(L^2)$ attention and H-bond matrices to $O(\text{nnz})$ linear arrays ([`docs/02_attention_traversal.md`](docs/02_attention_traversal.md)).
   * **All-vs-All SASA $\rightarrow$ $5\text{ Å}$ Voxel Hash Grid ([`compute_atomic_sasa`](src/lib.rs#L4142-L4221)):** Drops surface-area and H-bond burial complexity from $O(N_{\text{atoms}}^2)$ to $O(N_{\text{atoms}})$ ([`docs/15_chaperone_protein_design_workflows.md`](docs/15_chaperone_protein_design_workflows.md)).

