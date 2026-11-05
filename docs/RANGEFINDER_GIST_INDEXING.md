# RangeFinder Pocket-Interval GiST & B-Tree Indexing in `pg_bio`

**Reference:** Ruvalcaba, A. N., Alvin, C., Mukhopadhyay, S., & Bess, A. (2026). *RangeFinder: Predicting Ligand Property Ranges from Protein Pockets.* NeurIPS 2026 AI for Drug Discovery (AI4DD) Workshop (`arXiv:2611.03379v1`).

---

## 1. The Scientist's Advantage: Filtering 100M+ Molecules Before Touching a GPU

Structure-based drug discovery pipelines typically rely on either:
1. **3D Generative Models** (`Pocket2Mol`, `TargetDiff`, `3D-SBDD`) that sample individual molecules atom-by-atom inside a binding site, or
2. **Docking & Scoring Oracles** (`GNINA`, `AutoDock Vina`, `AlphaFold 3`, `Boltz-2`) that evaluate one protein–ligand complex at a time.

Both paradigms operate on **individual candidate molecules** rather than characterizing the **compatible physicochemical hyper-volume** of the binding pocket up front.

**RangeFinder (`arXiv:2611.03379v1`)** reformulates pocket-to-ligand modeling as a **9-dimensional interval prediction problem**. From a $10\text{ \AA}$ protein pocket point cloud ($5\text{ \AA}$ radius graph + 46 atomic chemistry descriptors + 1280-d ESM-2 embeddings), RangeFinder predicts the empirical mean $\hat{\mu}_{pj}$ and scale $\hat{\sigma}_{pj}$ across **nine ligand physicochemical properties**, parameterizing each interval with guaranteed positive width:

$$\hat{I}_{pj} = \left[\hat{\mu}_{pj} - 2|\hat{\sigma}_{pj}|,\; \hat{\mu}_{pj} + 2|\hat{\sigma}_{pj}|\right]$$

| # | Property | Column / Field | Definition | Degeneracy Note (Appendix A.1) |
|---|---|---|---|---|
| 1 | **MW** | `mw` | Molecular weight (Da) | Non-degenerate |
| 2 | **logP** | `logp` | Wildman-Crippen octanol–water partition coefficient | Non-degenerate |
| 3 | **HBD** | `hbd` | Hydrogen-bond donors ($\text{N-H} + \text{O-H} + \text{S-H}$) | Non-degenerate |
| 4 | **HBA** | `hba` | Hydrogen-bond acceptors ($\text{N} + \text{O}$) | Non-degenerate |
| 5 | **TPSA** | `tpsa` | Ertl (2000) topological polar surface area ($\text{\AA}^2$) | Non-degenerate |
| 6 | **Rotatable Bonds** | `rotatable_bonds` | Non-ring single bonds between non-terminal heavy atoms (excl. amides) | Non-degenerate |
| 7 | **Formal Charge** | `formal_charge` | Net formal charge | **Degenerate ($\sigma < 0.1$) in 81% of BioLiP2 pockets** |
| 8 | **Aromatic Rings** | `aromatic_rings` | Hückel $4n+2$ aromatic ring count (+ auxiliary BCE logit) | Non-degenerate |
| 9 | **$\text{Fsp}^3$** | `fsp3` | Fraction of $sp^3$-hybridized carbons ($N_{C,sp^3} / N_{C,\text{total}}$) | Non-degenerate |

---

## 2. Database Architecture: When to Use Composite B-Tree vs. Multi-Column GiST

In PostgreSQL 18 + `pg_bio`, virtual screening operates in two complementary directions:

### Direction A: Forward Virtual Screening (1 Target Pocket $\to$ 100M+ Molecules)
Here, the `molecules` table stores scalar `float8` columns (`mw`, `logp`, `tpsa`, `hba`, `hbd`, `rotatable_bonds`, `aromatic_rings`, `fsp3`), and the query supplies a single pocket's interval bounds $[\hat{\mu}_{pj} - 2|\hat{\sigma}_{pj}|, \hat{\mu}_{pj} + 2|\hat{\sigma}_{pj}|]$.
- **Optimal Index:** **Covering Composite B-Tree Index** on the highest-selectivity continuous descriptors `(mw, logp, tpsa, hba) INCLUDE (mol_id, chembl_id, hbd, rotatable_bonds, aromatic_rings, fsp3)` paired with IEEE-754 `BETWEEN` or `pg_bio`'s native `pocket_matches_ligand_8prop(pocket, props)` filter.
- **Why B-Tree Beats `numrange` Here:** Scalar `float8 BETWEEN` uses hardware IEEE-754 comparisons and a single B-Tree range scan, achieving **3.78 ms** across 250,000 compounds (**94.9× faster** than Python client-side filtering and **21.8× faster** than `numrange @> numeric` BCD arithmetic).

### Direction B: Reverse Target Fishing / Polypharmacology (1 Query Molecule $\to$ 20,000+ Protein Pockets)
Here, the `protein_pockets` table stores **intervals** (`[min, max]`) for every pocket in the proteome (e.g., all human AlphaFold/BioLiP2 pockets), and a scientist asks: *"Which protein pockets in the database have predicted property ranges compatible with my lead molecule?"*
- **Why B-Tree Fails on Intervals:** A B-Tree on `(mw_min, mw_max)` executing `WHERE mw_min <= 206.28 AND mw_max >= 206.28` must scan every pocket with `mw_min <= 206.28` (potentially half the table) before filtering `mw_max`.
- **Optimal Index:** **Multi-Column GiST Index over `numrange`** (`USING gist (mw_range, logp_range, tpsa_range, hba_range, hbd_range, rot_range, arom_range, fsp3_range)`). PostgreSQL's GiST R-tree splits 1D/multi-D intervals into bounding boxes, executing a 20,000-pocket reverse screen in **0.68 ms** (**16.6× faster** than sequential scan).

---

## 3. Live PostgreSQL 18 Benchmark Summary (`scripts/benchmark_rangefinder_gist.py`)

| Benchmark Stage | Workload | Method | Latency / Metric | Speedup / Notes |
| :--- | :--- | :--- | :--- | :--- |
| **1. Pocket Inference** | $10\text{ \AA}$ PDB Pocket (`6N9Q`, $138$ atoms, $2,570$ $5\text{ \AA}$ edges) | **`pg_bio` Native Rust (`predict_pocket_ranges`)**<br>RangeFinder $SE(3)$-Transformer (A10G GPU, Table 1)<br>Foldseek Baseline (Table 1)<br>APoc Baseline (Table 1) | **0.167 ms / pocket**<br>`12.200 ms / pocket`<br>`127.000 ms / pocket`<br>`26,300.000 ms / pocket` | **73.1×** faster than A10G GPU<br>**760.5×** faster than Foldseek<br>**157,485×** faster than APoc |
| **2. Forward Screen** | `1 Pocket -> 250,000 Molecules` ($198$ hits, $1,262.6\times$ reduction) | Python Client-Side Pull + Filter<br>PG Seq Scan (`numrange @> numeric` BCD)<br>PG Seq Scan (`float8 BETWEEN` IEEE-754)<br>**PG Composite B-Tree Index Scan** | `359.11 ms`<br>`82.35 ms`<br>`11.98 ms`<br>**`3.78 ms`** | `1.0×` baseline<br>`4.4×` speedup<br>`30.0×` speedup<br>**`94.9×` speedup** (`11.84 MB` index) |
| **3. Reverse Screen** | `1 Query Ligand -> 20,000 Pockets` (Ibuprofen, $10$ compatible pockets) | PG Sequential Scan<br>**PG 8-Column `numrange` GiST Scan (`@>`)** | `11.35 ms`<br>**`0.68 ms`** | `1.0×` baseline<br>**`16.6×` speedup** (`7.95 MB` GiST index) |
| **4. Charge Degeneracy** | Appendix A.1 9-Property vs. 8-Property Audit (`6N9Q_active`) | 9-Property (with `formal_charge`)<br>**8-Property (no `formal_charge`)** | Range-IoU `0.546` \| RRF `0.000`<br>**Range-IoU `0.615` \| RRF `0.615`** | Eliminates zero-IoU collapse on neutral pockets |

---

## 4. SQL Quickstart & Runnable Examples

See [`examples/rangefinder_pocket_screening/`](../examples/rangefinder_pocket_screening/) and [`experiments/rangefinder_pocket_gist_screening.ipynb`](../experiments/rangefinder_pocket_gist_screening.ipynb) for complete runnable examples.

```sql
-- 1. Compute the 9 RangeFinder properties from any SMILES string in pure Rust
SELECT * FROM compute_ligand_properties_row('CC(=O)Oc1ccccc1C(=O)O');

-- 2. Predict the 9 [mu - 2|sigma|, mu + 2|sigma|] intervals directly from a 10 Å PDB pocket
SELECT * FROM pocket_numranges(
    predict_pocket_ranges('6N9Q_active', pdb_text, 0.0, 0.0, 0.0, 10.0)
);

-- 3. Reverse Target Fishing via Multi-Column GiST Index on numrange
CREATE INDEX idx_pockets_rangefinder_gist ON protein_pockets
USING gist (mw_range, logp_range, tpsa_range, hba_range, hbd_range, rot_range, arom_range, fsp3_range);

WITH lig AS (
    SELECT * FROM compute_ligand_properties_row('CC(C)Cc1ccc(C(C)C(=O)O)cc1')
)
SELECT p.pocket_id, p.uniprot_id
FROM protein_pockets p, lig
WHERE p.mw_range   @> lig.mw::numeric
  AND p.logp_range @> lig.logp::numeric
  AND p.tpsa_range @> lig.tpsa::numeric
  AND p.hba_range  @> lig.hba::numeric
  AND p.hbd_range  @> lig.hbd::numeric
  AND p.rot_range  @> lig.rotatable_bonds::numeric
  AND p.arom_range @> lig.aromatic_rings::numeric
  AND p.fsp3_range @> lig.fsp3::numeric;
```
