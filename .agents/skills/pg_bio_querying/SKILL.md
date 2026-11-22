---
name: pg_bio_querying
description: Teaches the agent how to construct SQL queries using the custom pg_bio PostgreSQL extension for Z-Order spatial searches and Vector embeddings.
---

# Querying the pg_bio Extension

When asked to perform bioinformatics searches, structural queries, or spatial alignment inside the database, use the `pg_bio` extension capabilities.

## 1. Z-Order Spatial Indexing (Layered 3D Morton Octree & Bounding Boxes)
To find atoms within a radius of a specific target, do NOT use brute-force `sqrt()` math or a naive single-interval `z_index BETWEEN min_z AND max_z` scan (which suffers Morton Z-jumps across octree boundaries). Use `z_order_layered_sphere_ranges(cx, cy, cz, radius, max_ranges)` to decompose the 3D sphere into $\le 64$ tight octree Z-intervals (74× fewer atoms scanned, 10.25× faster, 100% exact recall):

**SQL Pattern (Layered 3D Morton Octree — Recommended):**
```sql
SELECT a.*
FROM z_order_layered_sphere_ranges(TGT_X, TGT_Y, TGT_Z, RADIUS, 64) r
JOIN your_table a ON a.z_index BETWEEN r.min_z AND r.max_z
WHERE r.fully_inside
   OR ((a.coord).x - TGT_X)^2 + ((a.coord).y - TGT_Y)^2 + ((a.coord).z - TGT_Z)^2 <= (RADIUS * RADIUS);
```

## 2. ESM-2 Vector Embeddings (`pg_vectorsearch` PRISM + `pgvector` `halfvec(1280) STORAGE PLAIN` + `<~>`)
Store 1280-D protein embeddings as `halfvec(1280) STORAGE PLAIN` indexed with Tiger Data's `pg_vectorsearch` (`USING prism (embedding vec16_cosine_ops)`), and re-rank top candidates with `pg_bio`'s `<~>` hybrid Smith-Waterman operator (or `pg_bio_search_homologs_prism`).

**SQL Pattern:**
```sql
SET LOCAL prism.nprobe = 12;
SELECT p.uniprot_id, p.name, (p.embedding <=> q.embedding) AS cosine_dist
FROM proteins p, (SELECT embedding FROM proteins WHERE uniprot_id = 'P0DPB7') q
WHERE p.uniprot_id <> 'P0DPB7'
ORDER BY p.embedding <=> q.embedding ASC
LIMIT 10;
```

## 3. 1024-Bit ECFP4 & 9D RangeFinder Register-Direct Fastscan
- **1024-Bit ECFP4 Tanimoto Fastscan (100% Exact Recall, 31× Faster):** Pre-filter using `popcount` and `ecfp4_fastscan_ub_raw(m.popcount, m.block_counts, m.sketch_lo, m.sketch_hi, q.popcount, q.block_counts, q.sketch_lo, q.sketch_hi) >= THRESHOLD` before calling `ecfp4_tanimoto_words_raw(m.ecfp_words, q.ecfp_words, m.popcount, q.popcount)`.
- **9D RangeFinder Active-Site Pocket Fastscan (100% Exact Recall, 13.2×–20.4× Faster):** Pre-filter using 16-byte register sketches (`pocket_quantize_fastscan` $\rightarrow$ `lo_code`, `hi_code`) via `pocket_fastscan_iou_ub_raw(p.lo_code, p.hi_code, q.lo_code, q.hi_code) >= MIN_IOU` before calling `pocket_iou_8prop_raw(...)`.

## 3. Sparse Attention Maps
To find which amino acids have the highest Attention weights against a target residue, use the `get_top_interacting_residues` function against a `SparseAttentionMap` type.

**SQL Pattern:**
```sql
SELECT get_top_interacting_residues(attention_data, TARGET_RESIDUE_ID, LIMIT_INT) 
FROM protein_attention_maps;
```

## Testing
Always run queries using standard `psql` connected to the `pgrx` instance (e.g., via `cargo pgrx run pg18` or `psycopg` scripts).
