# Layer 12: Wet-Lab Synthetic Biology, Multi-Host Codon Cassettes, 3D Globular Gate & Defluorination Digital Twin

While Layers 1–11 push genomics, 3D Z-order spatial indexing, single-cell CSR matrices, cross-omic vectors, UniProt APIs, PDB/VCF parsing, SMILES cheminformatics, de novo theozyme geometry, MS/MS proteomics, and MM/PBSA thermodynamic auditing into PostgreSQL, transitioning a computational design into a **physical wet-lab experiment** historically required dropping back into Python (`Biopython`, `PyMOL`, `NumPy`, `SciPy`) for four recurring bottlenecks:

1. **Multi-Host Reverse Translation & Codon Adaptation Index (CAI)** across diverse bacterial, fungal, archaeal, and mammalian expression hosts.
2. **Wet-Lab Expression Cassette Assembly & Cloning Safety Audits** (Promoter + RBS/Kozak + CDS + Stop + Terminator, GC% window, and Golden Gate `BsaI` restriction site screening).
3. **3D Globular Compactness & Folding Quality Gating** (`pLDDT`, Radius of Gyration $R_g$, Flory globular ratio $R_g / [2.25 \cdot N^{0.38}]$, and gyration tensor principal-axis elongation) to filter out disordered spaghetti and coiled-coil rods without loading PDBs into Python/NumPy.
4. **Enzymatic Defluorination & Biomineralization Digital Twin Kinetics** (simulating Fluoride Ion-Selective Electrode `[F⁻]` ppm release and stoichiometric $\text{Ca}_5(\text{PO}_4)_3\text{F}$ fluorapatite precipitation directly in SQL).

Layer 12 eliminates these Python bottlenecks natively in Rust (`pgrx`).

---

## 1. Multi-Host Reverse Translation & Codon Adaptation Index (`bio_reverse_translate`, `bio_codon_adaptation_index`)

Different expression hosts exhibit drastically different tRNA pools and GC content biases—from AT-rich *Saccharomyces cerevisiae* (~38.6% GC) to hyper-halophilic Archaea like *Halobacterium salinarum* (~71.6% GC):

```sql
SELECT host,
       bio_reverse_translate('MFEGFERRLVDWPSTDKTKREEVD', host) AS cds,
       round((bio_gc_content(bio_reverse_translate('MFEGFERRLVDWPSTDKTKREEVD', host)) * 100)::numeric, 1) AS gc_pct,
       round(bio_codon_adaptation_index(bio_reverse_translate('MFEGFERRLVDWPSTDKTKREEVD', host), host)::numeric, 3) AS cai,
       bio_translate_dna(bio_reverse_translate('MFEGFERRLVDWPSTDKTKREEVD', host)) AS roundtrip_aa
FROM unnest(ARRAY[
    'e_coli', 'b_subtilis', 'p_putida', 'c_necator',
    's_cerevisiae', 'p_pastoris', 'h_salinarum', 's_acidocaldarius', 'h_sapiens'
]) AS host;
```

---

## 2. Wet-Lab Plasmid Expression Cassette Compiler (`bio_design_expression_cassette`)

Constructs the complete synthetic plasmid cassette (`ExpressionCassette` composite type) with host-matched promoter (`T7`, `GAL1`, `AOX1`, `PFDA`, `KRT14`), ribosome binding site / Kozak sequence, optimal stop codon, terminator, GC content, CAI score, and Golden Gate `BsaI` (`GGTCTC` / `GAGACC`) restriction site verification:

```sql
SELECT (bio_design_expression_cassette(
    'MFEGFERRLVDWPSTDKTKREEVD',
    'e_coli',
    'T7'
))::text;
```

---

## 3. 3D Globular Compactness Gate (`evaluate_globular_compactness`)

Computes 3D structural folding metrics (`GlobularFoldingProfile`) in a single pass over PDB `ATOM ... CA` records:
- Auto-normalizes `pLDDT` B-factors whether scaled in `[0, 1]` (local ESMFold) or `[0, 100]` (AlphaFold DB).
- Computes Radius of Gyration $R_g = \sqrt{\frac{1}{N}\sum \|\mathbf{r}_i - \bar{\mathbf{r}}\|^2}$ and Flory Compactness Ratio $R_g / (2.25 \cdot N^{0.38})$.
- Solves the $3 \times 3$ gyration covariance tensor via Jacobi eigendecomposition to measure principal-axis elongation $\sqrt{\lambda_1 / \lambda_2}$.

```sql
SELECT evaluate_globular_compactness(pg_read_file('/path/to/teaflon_trifusion_local.pdb'))::text;
-- Returns: {"ca_count":418,"mean_plddt":82.89,"radius_of_gyration":22.96,"compactness_ratio":1.030,"elongation":1.62,"is_globular":true}
```

---

## 4. Fluoroacetate Dehalogenase Theozyme (`q1ju72_dehalogenase_hexad`) & Defluorination ISE Digital Twin (`simulate_defluorination_kinetics`)

`theozyme_rmsd` now includes built-in Kabsch superposition templates for **Fluoroacetate Dehalogenase (`Q1JU72` / PDB `3BJX`)**:
- `'q1ju72_dehalogenase_triad'`: Catalytic triad (`Asp104:OD2` nucleophile, `His271:NE2` general base, `Asp128:OD2` charge-relay acid).
- `'q1ju72_dehalogenase_hexad'`: Full 6-residue catalytic triad + halide-binding hole (`Asp104:OD2`, `His271:NE2`, `Asp128:OD2`, `Trp150:NE1`, `Tyr212:OH`, `Arg108:NH2`).

And `simulate_defluorination_kinetics` models the coupled Michaelis-Menten + Hydrophobin surface adsorption + $\text{Ca}_5(\text{PO}_4)_3\text{F}$ Fluorapatite biomineralization time course natively in SQL:

```sql
SELECT g AS time_hr,
       wt.fluoride_ppm AS wt_f_ppm,
       tf.fluoride_ppm AS teaflon_f_ppm,
       tf.fluorapatite_mg_per_l AS teaflon_fluorapatite_mg_l
FROM generate_series(0, 48, 12) AS g
CROSS JOIN LATERAL simulate_defluorination_kinetics(1.5, 450.0, 200.0, 0.25, 10.0, g::float8) wt
CROSS JOIN LATERAL simulate_defluorination_kinetics(12.0, 45.0, 1000.0, 5.0, 45.0, g::float8) tf;
```

---

## 5. Python Friction Log & Next Steps Roadmap for `pg_bio`

During the end-to-end validation of the **TeaFlon** multiomics compiler across **7,338 proteins and 930 organisms**, we recorded every computation that previously forced scientists out of PostgreSQL into Python, what we migrated into Rust in Layer 12, and the highest-leverage next steps for `pg_bio`:

| Bottleneck Previously Done in Python | How Layer 12 Solved It Natively in `pg_bio` | Next Step / Future Evolution for `pg_bio` |
| :--- | :--- | :--- |
| **1. `halfvec` vs `real[]` Cast Friction in SDK** (`find_templates.py` failed on `embedding_cosine_distance(halfvec, real[])`) | Updated `PgBioClient.find_homologues` to use native HNSW `<=>` (`embedding <=> get_esm_embedding(%s)::halfvec`). | Add overloaded `embedding_cosine_distance(halfvec, halfvec)` and native `get_esm_embedding_halfvec(text)` in Rust so users never have to cast manually. |
| **2. Single-Host Python Codon Tables** (`orca_dna_prep.py` used Biopython + a hardcoded 20-codon *E. coli* dict) | Added `bio_reverse_translate`, `bio_codon_adaptation_index`, and `bio_design_expression_cassette` supporting 9 bacterial, fungal, archaeal, and mammalian hosts. | Add mRNA secondary-structure hairpins ($\Delta G$ folding energy near the 5' RBS) and full Golden Gate / Gibson overhang assembly planners in SQL. |
| **3. PyMOL Headless Geometry & NumPy Rg/pLDDT Loops** (`measure.py` & `amazonzyme_swarm.py` parsed PDB text in Python) | Added `evaluate_globular_compactness(pdb_text)` with 3x3 Jacobi gyration eigendecomposition and `q1ju72_dehalogenase_hexad` in `theozyme_rmsd`. | Add native mmCIF (`.cif`) string parsing alongside PDB string parsing in `evaluate_globular_compactness`, plus solvent-accessible surface area (SASA / Shrake-Rupley) in Rust. |
| **4. Python ODE / Kinetic Curves** (`digital_twin_defluorination.py` computed hardcoded exponential curves in NumPy) | Added `simulate_defluorination_kinetics` SRF modeling Michaelis-Menten turnover, Hydrophobin surface adsorption boost, and stoichiometric $\text{Ca}_5(\text{PO}_4)_3\text{F}$ precipitation. | Generalize to a multi-step metabolic pathway ODE solver (`bio_simulate_pathway_ode`) accepting stoichiometric matrices and $k_{\text{cat}}/K_m$ vectors directly in SQL. |
| **5. Rigid-Body Multi-Domain Fusion Coordinate Stitching** (`export_trifusion_pdb.py` used PyMOL `cmd.translate` + `cmd.fab`) | Verified folded 418-aa TriFusion coordinates (`teaflon_trifusion_local.pdb`) natively via `bio_parse_pdb` and `evaluate_globular_compactness`. | Add `bio_assemble_fusion_pdb(pdb_domains text[], linkers text[])` in Rust using internal coordinate (NERF) backbone placement to stitch multi-domain chimeras without PyMOL. |
