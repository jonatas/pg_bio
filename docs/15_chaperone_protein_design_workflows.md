# Tutorial 15: RosettaCommons Chaperone Workflows & Official Core Examples in SQL — Shrake–Rupley RASA, 3D Hydrogen-Bond CSR Networks, BUNS & Lineage Auditing

## 1. The Scientist's Problem: The Post-Oracle Bottleneck in Agentic Protein Design

Modern computational protein design pipelines—standardized by [**RosettaCommons Chaperone**](https://rosettacommons.org/2026/10/19/introducing-chaperone-a-unified-framework-for-protein-design-workflows/) (`RosettaCommons/chaperone`)—chain together modular blocks via a Polars-backed `PipelineDataHandler`:

$$\text{RFdiffusion3 / RFD3NA (Backbone)} \longrightarrow \text{Foundry-MPNN / LigandMPNN (Sequence)} \longrightarrow \text{AlphaFold3 / Boltz-2 (Oracle)} \longrightarrow \text{Evaluation Blocks (RMSD, RASA, H-Bond)}$$

While deep learning inference runs on GPUs, **evaluating, filtering, and tracking lineage across 10,000–100,000 predicted structures** creates a massive CPU/I/O bottleneck when each evaluation block (`rmsd`, `rasa`, `hbond`) spawns separate Python processes to parse loose PDB/mmCIF files from disk:
1. **Surface Exposure & Aggregation Traps (`rasa`):** De novo designs often fold with high pLDDT in silico yet aggregate in wet-lab expression because bulky hydrophobic sidechains (`Val, Ile, Leu, Met, Phe, Trp, Tyr`) remain solvent-exposed, or the binding interface buries insufficient $\Delta\text{SASA}$.
2. **Buried Unsatisfied Polar Atoms (`hbond` / BUNS):** Burying a polar $\text{N}$ or $\text{O}$ atom at a binder–target interface ($\text{SASA} < 1.0\text{ \AA}^2$) without a compensating hydrogen bond or metal/water coordination incurs a $+3\text{ to }+5\text{ kcal/mol}$ desolvation penalty—destroying experimental binding affinity.
3. **Multi-Cycle Lineage Fragmentation (`PipelineDataHandler`):** Iterative design loops (`MPNN <-> AF3/Boltz-2`) require deterministic parent–child lineage tracking so surviving hits can be traced back to their exact `RFdiffusion3` seed backbone and compiled directly into wet-lab expression cassettes.

---

## 2. The 5-Stage `pg_bio` Chaperone Engine

### Stage 1: Shrake–Rupley Golden-Section Fibonacci Sphere SASA & Tien (2013) RASA (`bio_rasa`)
- Evaluates 64 Golden-Ratio spiral test points per heavy atom ($R_i = r_{\text{vdW},i} + r_{\text{probe}}$ with $r_{\text{probe}} = 1.4\text{ \AA}$) using 3D bounding-box pruning and nearest-neighbor distance sorting so occluded interior points short-circuit on the first 1–2 neighbor checks.
- Normalizes per-residue Solvent Accessible Surface Area ($\text{SASA}_i$, in $\text{\AA}^2$) by the **Tien et al. (2013) empirical Gly-X-Gly maximum ASA scale** (plus nucleotide ASA for DNA/RNA `DA, DG, DC, DT, A, U, G, C`) to return `RasaProfile`:
  - `total_sasa`, `polar_sasa`, `apolar_sasa` ($\text{\AA}^2$)
  - `mean_rasa` ($\in [0, 1]$)
  - `core_buried_fraction`: Fraction of residues with $\text{RASA}_i < 0.20$
  - `surface_hydrophobic_fraction`: Fraction of solvent-exposed residues ($\text{RASA}_i \ge 0.25$) that are bulky hydrophobic residues (`VAL, ILE, LEU, MET, PHE, TRP, TYR`)
  - `per_residue_rasa`: Full `FLOAT8[]` array for residue-level filtering

### Stage 2: 3D Hydrogen-Bond CSR Network (`bio_hbond_network`) & BUNS Auditor (`audit_hydrogen_bonds`)
- Identifies all backbone, sidechain, nucleobase, and `HETATM` ligand donor/acceptor heavy atoms ($\text{N}, \text{O}, \text{S}, \text{F}$), enforces donor–acceptor distance ($2.20\text{ \AA} \le d_{\text{DA}} \le 3.50\text{ \AA}$) and zero-allocation dot-product antecedent angular geometry ($\angle(\text{Ant}-\text{D}\cdots\text{A}) \ge 90^\circ \iff (\mathbf{x}_{\text{Ant}} - \mathbf{x}_{\text{D}})\cdot(\mathbf{x}_{\text{A}} - \mathbf{x}_{\text{D}}) \le 0$), and weights interactions via $w = \exp(-((d_{\text{DA}} - 2.90)/0.45)^2)$.
- `bio_hbond_network(pdb_text, 3.5)` projects the 3D hydrogen-bond graph directly into `pg_bio`'s Compressed Sparse Row (CSR) `SparseAttentionMap` type.
- `audit_hydrogen_bonds(pdb_text, binder_chain, target_chain)` computes `total_hbonds`, `backbone_hbonds`, `sidechain_hbonds`, `cross_interface_hbonds`, `buried_unsatisfied_polars` (**BUNS**), and `hbond_satisfaction_rate`.

### Stage 3: Single-Pass Fused Chaperone Evaluator (`evaluate_chaperone_block` / `evaluate_chaperone_design`)
- Fuses Kabsch self-consistency $\text{C}\alpha / \text{C4}'$ RMSD (`sc_rmsd`), `mean_plddt`, `mean_rasa`, `core_buried_fraction`, `surface_hydrophobic_fraction`, buried interface $\Delta\text{SASA} = \text{SASA}_{\text{binder}} + \text{SASA}_{\text{target}} - \text{SASA}_{\text{complex}}$, `cross_interface_hbonds`, `buried_unsatisfied_polars` (BUNS), and 3D globular compactness (`is_globular`) into a single pass.

### Stage 4: Content-Addressable Lineage Tracking (`chaperone_lineage_hash`)
- Generates deterministic 64-bit FNV-1a / MurmurHash3 (`fmix64`) lineage node IDs (`chp_<stage>_c<cycle>_<hash16>`) across `rfd3 -> mpnn -> af3 / boltz2` blocks.

### Stage 5: Official RosettaCommons Chaperone Core Examples (`chaperone_core_examples`, `evaluate_chaperone_core_example`)
- Built-in support for all **8 official Chaperone Core Examples** ([`rosettacommons.github.io/chaperone/core-examples/`](https://rosettacommons.github.io/chaperone/core-examples/)):

| Official Core Example | Modality | Chaperone Pipeline Blocks | `pg_bio` In-Database Verification Gate |
| :--- | :--- | :--- | :--- |
| **`nanoparticle-design`** | `symmetric_assembly` | `rfd3(C3) -> dock(I3) -> tied_mpnn -> af3` | Symmetric oligomer $\Delta\text{SASA} \ge 100\text{ \AA}^2$, $\text{scRMSD} \le 2.0\text{ \AA}$, $\text{BUNS} \le 3$ |
| **`peptide-binder`** | `protein_peptide` | `rfd3(peptide) -> mpnn -> af3/boltz2 -> cycle` | Cross-interface $\beta$-strand H-bonds, $\Delta\text{SASA}$, $\text{scRMSD} \le 2.0\text{ \AA}$ |
| **`cortisol-binder`** | `protein_ligand` | `rfd3(cortisol) -> ligand_mpnn -> boltz2/af3` | `cortisol_binder` pocket theozyme RMSD + steroid VF2 subgraph `@>` + `HETATM` H-bonds |
| **`dna-binder`** | `protein_dna` | `rfd3na(dsDNA) -> mpnn -> af3` | B-DNA (`DA, DG, DC, DT`) major-groove H-bonds + $\Delta\text{SASA}$ |
| **`enzyme-glycosidase`** | `enzyme_catalysis` | `rfd3(dyad) -> ligand_mpnn -> branch(apo/cov/holo)` | Retaining `GLU/GLU/TRP/TYR` catalytic dyad `theozyme_rmsd <= 0.90` + glycan H-bonds |
| **`enzyme-protease`** | `enzyme_catalysis` | `rfd3(triad) -> ligand_mpnn -> branch(monomer/complex)` | Catalytic triad `HIS/ASP/SER/GLY` `theozyme_rmsd <= 0.90` + substrate H-bonds |
| **`ribozyme-hammerhead`** | `rna_catalysis` | `rfd3na(hammerhead) -> rna_mpnn -> af3` | RNA `C4'` $\text{scRMSD}$, `G12/C17/A9/MG` core theozyme RMSD + `has_coordinated_water_shell(..., 'MG', 2.5)` |
| **`multi-mpnn`** | `inverse_folding_benchmark` | `rfd3 -> submit(Protein/Ligand/Soluble/FoundryMPNN) -> gather` | Multi-backend solubility (`surface_hydrophobic_fraction`), $\text{scRMSD}$, and `pLDDT` ranking |

---

## 3. SQL Examples

```sql
-- 1. Audit all 8 Official RosettaCommons Chaperone Core Examples in a single SQL query
SELECT
    ex.example_name,
    ex.category,
    ROUND(audit.sc_rmsd::numeric, 3) AS sc_rmsd,
    ROUND(audit.mean_plddt::numeric, 1) AS plddt,
    ROUND(audit.interface_delta_sasa::numeric, 1) AS delta_sasa,
    audit.cross_interface_hbonds AS iface_hb,
    audit.buried_unsatisfied_polars AS buns,
    ROUND(audit.active_site_or_motif_rmsd::numeric, 3) AS motif_rmsd,
    audit.ligand_or_cofactor_verified AS verified,
    ROUND(audit.composite_score::numeric, 3) AS score,
    audit.passes_example_gate AS pass_gate
FROM chaperone_core_examples() AS ex
CROSS JOIN LATERAL evaluate_chaperone_core_example_block(
    ex.example_name, ex.candidate_pdb, ex.reference_pdb, ex.ligand_smiles
) AS audit;

-- 2. Filter a 400-complex Chaperone campaign & compile lineage hashes
SELECT
    d.design_id,
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

## 4. Benchmark Summary (`scripts/benchmark_chaperone_workflows.py`)

| Benchmark Metric (PostgreSQL 18 on Apple Silicon) | Python Pull + Evaluation | `pg_bio` Native Rust Engine | Speedup / Result |
| :--- | :--- | :--- | :--- |
| **400 Multi-Chain Complexes (`144` Heavy Atoms Each)** | `7,709.11 ms` | **`209.04 ms`** (`0.52 ms`/complex) | **`36.9x` faster** |
| **4-Archetype Failure Mode Discrimination** (Pass vs Drift/Hydro/BUNS) | Partial (no BUNS/RASA) | **`100 / 400` exact recall** | **100% precision** |
| **8 Official Chaperone Core Examples Suite** (`chaperone_core_examples`) | Manual multi-script setup | **`~4.5 ms` total (`8 / 8` verified)** | **Zero external scripts** |
