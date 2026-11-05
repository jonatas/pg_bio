# RangeFinder (`arXiv:2611.03379v1`) Pocket-Interval GiST & B-Tree Screening Example

This example demonstrates **`pg_bio` Layer 15** — an in-database Rust implementation of **RangeFinder** (*Predicting Ligand Property Ranges from Protein Pockets*, `arXiv:2611.03379v1`) paired with PostgreSQL 18 **Covering Composite B-Tree** and **8-Column `numrange` GiST** indexes.

## Files in This Example

- [`schema_and_queries.sql`](./schema_and_queries.sql): Self-contained PostgreSQL 18 DDL and queries for:
  1. Creating `rf_example_pockets` with an 8-column `numrange` GiST index (`idx_rf_example_pockets_gist`).
  2. Creating `rf_example_ligands` with a covering composite B-Tree index (`idx_rf_example_ligands_btree`).
  3. Computing 9 physicochemical descriptors (`mw`, `logp`, `hbd`, `hba`, `tpsa`, `rotatable_bonds`, `formal_charge`, `aromatic_rings`, `fsp3`) from SMILES in Rust via `compute_ligand_properties_row()`.
  4. Forward virtual screening (`1 Pocket -> Ligand Library`) and Reverse target fishing (`1 Ligand -> Proteome Pockets`).
- [`run_rangefinder_example.py`](./run_rangefinder_example.py): End-to-end Python + `PgBioClient` + `psycopg` runner.

## Quickstart

```bash
# 1. Run the self-contained example
python examples/rangefinder_pocket_screening/run_rangefinder_example.py

# 2. Run the 250,000-compound B-Tree & 20,000-pocket GiST benchmark suite
python scripts/benchmark_rangefinder_gist.py
```

## Live PostgreSQL 18 Benchmark Results (`bio_demo`)

| Benchmark Stage | Workload | Method | Latency | Speedup |
| :--- | :--- | :--- | :--- | :--- |
| **1. Pocket Inference** | $10\text{ \AA}$ PDB Pocket ($138$ atoms, $2,570$ $5\text{ \AA}$ edges) | **`pg_bio` Native Rust (`predict_pocket_ranges`)** | **0.167 ms / pocket** | **73.1×** vs. A10G GPU (`12.2 ms`)<br>**760.5×** vs. Foldseek (`127 ms`) |
| **2. Forward Virtual Screen** | `1 Pocket -> 250,000 Ligands` ($1,262.6\times$ library reduction, $198$ hits) | Python Client Pull + Filter<br>PG Seq Scan (`numrange @> numeric` BCD)<br>PG Seq Scan (`float8 BETWEEN` IEEE-754)<br>**PG Covering Composite B-Tree Scan** | `359.11 ms`<br>`82.35 ms`<br>`11.98 ms`<br>**`3.78 ms`** | `1.0×` baseline<br>`4.4×` speedup<br>`30.0×` speedup<br>**`94.9×` speedup** |
| **3. Reverse Target Fishing** | `1 Query Ligand -> 20,000 Pockets` ($10$ compatible pockets) | PG Sequential Scan<br>**PG 8-Column `numrange` GiST Scan (`@>`)** | `11.35 ms`<br>**`0.68 ms`** | `1.0×` baseline<br>**`16.6×` speedup** |
| **4. Appendix A.1 Audit** | 9-Property vs. 8-Property (dropping degenerate `formal_charge`) | 9-Property (`include_formal_charge = true`)<br>**8-Property (`include_formal_charge = false`)** | Range-IoU `0.546`, RRF `0.000`<br>**Range-IoU `0.615`, RRF `0.615`** | Restores recall on neutral pockets ($\sigma_{\text{charge}} < 0.1$) |
