# 🧬 pg_bio: The Multiomics Operating System for Synthetic Biology

**`pg_bio`** is a native PostgreSQL extension written in highly optimized Rust (`pgrx`). 

It is engineered for a single purpose: to completely eradicate the data bottlenecks in modern biology by pushing complex multiomics mathematics *directly into the database engine*. By unifying Genomics, Proteomics, and Transcriptomics into a shared computational space, `pg_bio` allows scientists to execute entire drug discovery pipelines and synthetic biology designs in milliseconds natively in SQL.

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

### Layer 2: Proteomics & Epigenomics (3D Spatial Curves)
Once transcribed, linear sequences fold into 3D atomic structures. Finding which atoms are within a 4.0 Ångstrom radius of a binding pocket usually requires extracting millions of coordinates into Python for brute-force Pythagorean math.

* **The Engine:** `pg_bio` contains a native Rust function `z_order_encode(x, y, z)`. This uses **Morton Coding (Z-Order curves)** to mathematically interleave 3D Cartesian coordinates into a single 1-Dimensional integer.
* **The Advantage:** Standard, ultra-fast 1D PostgreSQL B-Tree indexes now act as 3D bounding boxes. Extract exact 3D binding pockets or detect cross-protein molecular clashes across millions of atoms instantly—with zero setup time.
* **SQL Example:**
  ```sql
  -- Extracting a 3D Binding Pocket in milliseconds
  SELECT atom_id, residue_name FROM teaflon_atoms
  WHERE z_curve_index BETWEEN z_order_encode(5.0, -5.0, 20.0) AND z_order_encode(15.0, 5.0, 30.0)
    AND x BETWEEN 5.0 AND 15.0 AND y BETWEEN -5.0 AND 5.0 AND z BETWEEN 20.0 AND 30.0;
  ```

### Layer 3: Transcriptomics (Single-Cell CSR Arrays)
To understand drug dosage, scientists use Single-Cell RNA-Seq. Humans have 20,000 genes. Sequencing 1 million cells creates a 20-Billion point matrix. Traditional SQL schemas grind to a halt.

* **The Engine:** Because a single cell only expresses ~2,000 genes at a time, the matrix is 85%+ zeroes. `pg_bio` implements **Compressed Sparse Row (CSR)** architecture natively using parallel PostgreSQL Arrays (`expressed_gene_ids INT[]` and `expression_counts REAL[]`).
* **The Advantage:** Storage overhead is slashed by 90%. Using parallel `LATERAL UNNEST() WITH ORDINALITY`, PostgreSQL unpacks the sparse data dynamically in execution memory to calculate highly precise RNA dosage levels in milliseconds.

### Layer 4: Joint Vector Spaces (The Cross-Omic Holy Grail)
Traditionally, DNA sequences and 3D Protein folds live in different worlds. `pg_bio` breaks down that wall.

* **The Engine:** By utilizing `pgvector` alongside `pg_bio`, we map the outputs of Genomic Foundation Models (like Enformer) and Proteomic Foundation Models (like ESM-2) into a **shared 1280-dimensional mathematical space**.
* **The Advantage:** For the first time, you can execute cross-domain queries. You can ask the database to mathematically link a physical DNA regulatory enhancer directly to a 3D Protein structure using Cosine Distance (`<=>`).
* **SQL Example:**
  ```sql
  -- Find which DNA enhancer matches the structural embedding of the P0DPB7 Protein
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

---
## 🛠️ Scale & Optimization Infrastructure

`pg_bio` isn't just mathematically clever; it is heavily optimized for multi-million scale deployment on commercial hardware.

### 1. Vector Compression (`halfvec`)
High-dimensional vectors require massive amounts of RAM. `pg_bio` natively supports casting embedding columns to `halfvec` (16-bit floating-point). This cuts RAM and SSD footprint in half with zero loss in structural accuracy, allowing graphs of 2.5 million proteins to build natively in-memory on standard laptops.

### 2. Intelligent Graph Indexing (HNSW vs IVFFlat)
Vector indexing is fully tunable to hardware profiles. 
* Use **`HNSW` (Hierarchical Navigable Small World)** for the highest recall accuracy, utilizing optimized `maintenance_work_mem` ceilings to build massive interconnected graphs. 
* Use **`IVFFlat` (K-Means Clustering)** to bypass severe RAM limitations, dropping index build times from days to minutes on localized machines by sorting tuples into hyper-clusters (lists=1500).

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
git clone https://github.com/your-org/pg_bio.git
cd pg_bio

# 2. Spin up the Rust-optimized Postgres Engine
docker-compose up -d

# 3. Compile the Database and Interactive Notebooks
uv run jupyter nbconvert --execute https://github.com/jonatas/teaflon/blob/main/pg_bio_multiomics_pipeline/teaflon_multiomics_tutorial.ipynb
```

*(For an extended deep-dive into the architectural theory, read [`docs/MULTIOMICS_ARCHITECTURE.md`](docs/MULTIOMICS_ARCHITECTURE.md))*

---

## 📊 Master Catalog of `pg_bio` Correlation Types & Performance Architecture

Below is the complete reference of correlation engines built into [`src/lib.rs`](src/lib.rs), how each correlation works intuitively, how `pg_bio` executes it natively inside PostgreSQL, why it outperforms Python pipelines, and where to find the official documentation and runnable examples.

| # | Correlation Type | Simple Explanation (How It Works) | Technical `pg_bio` Approach (`src/lib.rs` & SQL) | Why `pg_bio` (Rust in Postgres) Is Faster Than Python | Official Docs & Examples |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **1** | **Twilight-Zone Homology Index**<br>*(3D Convergence vs. Sequence Divergence)* | Finds pairs of proteins that fold into the **exact same 3D shape** even though their amino acid sequences have mutated heavily (`15%–75%` similarity), filtering out both 100% clones and random k-mer noise. | Combines HNSW cosine distance (`<=>`), [`smith_waterman_score`](src/lib.rs#L1964-L1995), and [`sequence_entropy`](src/lib.rs#L1780-L1805) in a bell-curve filter:<br>`(1 - d_cos)^2 * [4 * s * (1 - s)] * (L_min / L_max)^2 * min(H / 4.0, 1.0)`. | **Zero-IPC & $O(N)$ Two-Row DP:** [`smith_waterman_score`](src/lib.rs#L1964-L1995) runs in compiled Rust using two flat `Vec<i32>` rolling buffers directly on Postgres `TEXT` pointers—no $M \times N$ Python matrix allocation and no socket serialization of millions of sequences. | 📖 [`01_vector_homology.md`](docs/01_vector_homology.md)<br>📖 [`08_native_alignment.md`](docs/08_native_alignment.md)<br>🐝 [`swarm_discovery.py`](examples/swarm_discovery/swarm_discovery.py)<br>📓 [`07_native_alignment.ipynb`](experiments/07_native_alignment.ipynb) |
| **2** | **Extremophile Niche Shift Index**<br>*(Electrostatic & Hydropathy Remodeling: $\Delta\text{pI}$ / $\Delta\text{GRAVY}$)* | Measures how a protein rewires its surface charge (acidic vs. basic) and water-repellency to survive extreme habitats (deep-sea brine, acid vents, soda lakes) while keeping the same 3D core. | Computes [`isoelectric_point`](src/lib.rs#L1828-L1905) via 100-step bisection on Henderson-Hasselbalch $\text{p}K_a$ curves and [`gravy_score`](src/lib.rs#L1806-L1827) (Kyte-Doolittle hydropathy) inline with `<=>`. | **Single-Pass Byte-Slice Histograms:** Counts amino acid frequencies in a single CPU register pass (`seq.bytes()`) and solves $\text{pI}$ in stack-allocated `f64` math (~200 ns/protein vs. ~40 µs in BioPython). | 📖 [`11_sequence_analysis.md`](docs/11_sequence_analysis.md)<br>📖 [`04_native_chemistry.md`](docs/04_native_chemistry.md)<br>🐝 [`swarm_web.py`](examples/swarm_discovery/swarm_web.py) |
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

---

### ⚡ Why `pg_bio` Is 10x–1,000x Faster Than a Python Pipeline

When executing a multi-modal discovery query (such as [`examples/swarm_discovery/swarm_discovery.py`](examples/swarm_discovery/swarm_discovery.py) or [`examples/continuous_discovery/continuous_night_discovery.py`](examples/continuous_discovery/continuous_night_discovery.py)), `pg_bio` fuses all filtering layers inside a single PostgreSQL backend process:

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

