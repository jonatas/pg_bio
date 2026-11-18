# Research & Benchmark: Should `pg_bio` Replace `pgvector` with `pg_vectorsearch` or Combine Them?

**Branch / Worktree:** `/Users/jonatas/code/pg_bio_vectorsearch` (`feat/pg-vectorsearch-benchmark`)  
**PostgreSQL Engine:** PostgreSQL 18.3 (`aarch64-apple-darwin`, Apple Accelerate BLAS + ARM64 NEON SIMD)  
**Extensions Tested:** `pg_bio` `0.0.0` + `vector` (`pgvector`) `0.8.6` + `pg_vectorsearch` `0.2.0-dev` (`5bb4e169`)

---

## 1. Executive Verdict: **Combine (`pgvector` + `pg_vectorsearch` + `pg_bio`), Do NOT Replace**

It makes **overwhelming architectural and performance sense to COMBINE `pgvector` and `pg_vectorsearch` inside `pg_bio`**, rather than replacing `pgvector` outright:

1. **Zero-Copy Binary Interoperability by Design:**  
   `pg_vectorsearch` (`vec32` / `vec16`) uses the **exact same varlena binary layout** (`varlena header + int16 dim + int16 unused + float[]`) as `pgvector` (`vector` / `halfvec`). When both extensions are installed, `pg_vectorsearch_setup_pgvector_compat()` automatically installs `WITHOUT FUNCTION` (`RelabelType`) zero-overhead binary casts and registers `pgvector`'s `<->`, `<#>`, and `<=>` operators directly into `prism`'s operator families (`vec32_*_ops`, `vec16_*_ops`).
   - **Result:** Existing `pg_bio` tables (`proteins.embedding halfvec(1280)`), Python SDK queries (`pgbio-py`), and `pg_bio`'s `<~>` hybrid distance operator (`bio_feature`) work **without rewriting tables or changing SQL queries**.
2. **Massive Proteome Scale (PRISM) vs. Small In-Memory Assays (HNSW) & Bit/Sparse Vectors (`pgvector`):**
   - **Use `pg_vectorsearch` (`USING prism`) for large embedding tables ($\ge 10\text{k}$ to $100\text{M}+$ rows):** On `bio_demo`'s `proteins` table (2.34M rows of `halfvec(1280)`), `pgvector`'s HNSW index consumes **8,772 MB (8.77 GB)**—8.3× larger than the main table heap (1,054 MB). In contrast, `pg_vectorsearch` PRISM uses **1-bit RaBitQ quantization with theoretical error bounds + hierarchical centroid routing**, shrinking the index by **9.8× (default SOAR + boundary replication)** to **17.2× (compact)**, building **20× faster**, and querying **7.1×–9.8× faster** at 98.4%–99.2% Recall@10.
   - **Keep `pgvector` for column types, ecosystem compatibility, small tables, and bit/sparse vectors:** `pg_vectorsearch` is currently an **alpha** release targeting **PostgreSQL 18+** only, and only implements dense float32/float16 ANN (`vec32`/`vec16`). `pgvector` remains essential for `pg13–pg17` compatibility, `vector`/`halfvec` type definitions (`requires = 'vector'` in `pg_bio.control`), `bit` vectors (Hamming/Jaccard distances), `sparsevec`, and sub-10k-row tables where in-memory HNSW needs zero k-means training.

---

## 2. License Audit & Compatibility Matrix

All three components in the combined stack use **permissive, non-copyleft open-source licenses** with zero SSPL/BSL/AGPL restrictions:

| Component | Repository / Author | License | Key Terms & Commercial / Open-Source Compatibility |
| :--- | :--- | :--- | :--- |
| **`pg_vectorsearch`** | `timescale/pg_vectorsearch` (Timescale, Inc. d/b/a Tiger Data) | **The PostgreSQL License** | Fully permissive (MIT/BSD-equivalent). Grants permission to use, copy, modify, and distribute for any purpose without fee. All C/SIMD/RaBitQ source files in `src/` carry `Licensed under the PostgreSQL License`. *(Note: Distinct from old Timescale License / TSL—this is pure OSI PostgreSQL License).* |
| **`pgvector`** | `pgvector/pgvector` (Andrew Kane) | **The PostgreSQL License** | Fully permissive (MIT/BSD-equivalent). 100% compatible with `pg_vectorsearch` and `pg_bio`. |
| **`pg_bio` (Rust / `pgrx`)** | `jonatas/pg_bio` (`pgrx`, `bio`, `noodles`, `pdbtbx`, `purr`) | **MIT / Apache-2.0 / PostgreSQL Compatible** | All Rust dependencies (`pgrx = 0.17.0`, `bio = 4.0.1`, `noodles = 0.116.0`, `pdbtbx = 0.12.0`, `purr = 0.9.0`) use permissive MIT/Apache-2.0 licenses. |

---

## 3. Empirical Benchmarks on Real 1280-D ESM-2 Proteome Embeddings

All benchmarks were executed on the live `bio_demo` database (`PostgreSQL 18.3` on Apple Silicon M-series with ARM64 NEON + Apple Accelerate BLAS) using real 1280-dimensional `halfvec(1280)` ESM-2 protein embeddings from UniProt (`scripts/benchmark_pg_vectorsearch_vs_pgvector.py`).

### 3.1 Critical Discovery #1: `STORAGE PLAIN` Eliminates the 1280-D TOAST Bottleneck

A 1280-dimensional `halfvec(1280)` occupies $4 + 1280 \times 2 = 2,564\text{ bytes}$, which exceeds PostgreSQL's `TOAST_TUPLE_THRESHOLD` (~2,000 bytes).
- In `bio_demo.proteins` (2.34M rows), default `STORAGE EXTERNAL` forces PostgreSQL to push every embedding into a **6,112 MB TOAST table** (`main_heap = 1,054 MB`, `toast_heap = 6,112 MB`, `hnsw_index = 8,772 MB`).
- Because PRISM's Stage-2 RaBitQ error-bounded reranker fetches full-precision vectors from the heap for the top candidates, TOASTed embeddings incur an extra B-tree + TOAST page lookup per candidate:

| Storage Mode (`50,000` proteins, `halfvec(1280)`) | Main Heap | TOAST Heap | PRISM Build Time | PRISM Query Latency (`nprobe=20`, p50) | Speedup |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`STORAGE EXTERNAL` (PostgreSQL Default)** | 21.73 MB | 130.36 MB | 6.87 s | 4.83 ms (mean 5.33 ms) | 1.00× |
| **`STORAGE PLAIN` (Inline Heap Tuple)** | 130.21 MB | 22.71 MB *(seq only)* | **3.18 s** | **1.67 ms** (mean 1.67 ms) | **2.89× faster query / 2.16× faster build** |

---

### 3.2 Index Build Speed & Storage Footprint (`20,000` and `50,000` Proteins, `1280-D halfvec`)

| Index Access Method | `20k` Build Time | `20k` Index Size | `50k` Build Time | `50k` Index Size | Projected `2.34M` (`proteins`) Size | Memory / Disk Reduction vs. HNSW |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`pgvector` HNSW** (`m=16, ef_construction=64`) | 22.94 s | 78.03 MB | 63.23 s | 195.16 MB | **8,772 MB** *(measured live)* | 1.00× (Baseline) |
| **`pgvector` IVFFlat** (`lists=224 / 316`) | 22.24 s | 53.28 MB | 44.89 s | 131.88 MB | ~6,170 MB | 1.48× smaller |
| **`pg_vectorsearch` PRISM (Default: SOAR + Boundary + RaBitQ)** | **1.10 s** | **8.21 MB** | **3.18 s** | **19.87 MB** | **~930 MB** | **9.82× smaller, 19.9× faster build** |
| **`pg_vectorsearch` PRISM (Compact: `soar_lambda=0, boundary_epsilon=0`)** | **0.94 s** | **5.52 MB** | **2.70 s** | **11.38 MB** | **~533 MB** | **17.15× smaller, 23.4× faster build** |

---

### 3.3 Query Latency, Throughput (QPS) & Recall (`50,000` Proteins, `1280-D halfvec`)

Brute-force sequential scan baseline (`enable_indexscan = off`): **p50 = 118.44 ms** (mean = 117.51 ms).

| Engine & Index Configuration | Parameter (`nprobe` / `ef_search` / `probes`) | p50 Latency (ms) | Mean Latency (ms) | Throughput (QPS) | Recall@10 (%) | Recall@50 (%) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`pg_vectorsearch` PRISM (Default)** | `nprobe = 5` | **1.51 ms** | 1.57 ms | **636.3** | 56.40% | 42.08% |
| **`pg_vectorsearch` PRISM (Default)** | `nprobe = 20` | **1.67 ms** | 1.67 ms | **599.9** | 88.40% | 78.08% |
| **`pg_vectorsearch` PRISM (Default)** | `nprobe = 40` | **1.90 ms** | 2.16 ms | **462.6** | **98.40%** | **93.52%** |
| **`pg_vectorsearch` PRISM (Default)** | `nprobe = 60` | **2.75 ms** | 2.95 ms | **339.2** | **99.20%** | **97.92%** |
| **`pg_vectorsearch` PRISM (Compact)** | `nprobe = 40` | 3.85 ms | 4.17 ms | 240.0 | 85.20% | 76.64% |
| **`pgvector` HNSW** | `ef_search = 10` | 3.78 ms | 3.65 ms | 274.3 | 86.80% | 18.08% |
| **`pgvector` HNSW** | `ef_search = 40` *(default)* | 8.75 ms | 8.15 ms | 122.7 | 96.00% | 77.84% |
| **`pgvector` HNSW** | `ef_search = 100` | 13.56 ms | 13.22 ms | 75.6 | **98.40%** | **94.80%** |
| **`pgvector` HNSW** | `ef_search = 200` | 26.83 ms | 26.76 ms | 37.4 | **99.60%** | **98.72%** |
| **`pgvector` IVFFlat** | `probes = 10` | 12.91 ms | 12.84 ms | 77.9 | 50.40% | 37.20% |
| **`pgvector` IVFFlat** | `probes = 40` | 22.43 ms | 23.09 ms | 43.3 | 80.80% | 70.80% |

**Key Takeaway at Iso-Recall (98.40% Recall@10):**
- **`pg_vectorsearch` PRISM (`nprobe = 40`)** achieves **98.40% Recall@10** at **1.90 ms p50 (462.6 QPS)**.
- **`pgvector` HNSW (`ef_search = 100`)** achieves the same **98.40% Recall@10** at **13.56 ms p50 (75.6 QPS)**.
- **PRISM is 7.14× faster in p50 latency and delivers 6.12× higher QPS while using 9.82× less disk/RAM.**

---

### 3.4 Combined Pipeline Benchmark: Stage-1 PRISM + Stage-2 `pg_bio` `<~>` Hybrid Smith-Waterman Re-Ranking

In `pg_bio`, biological homology discovery (`pg_bio_search_homologs` / `pg_bio_search_homologs_prism`) fuses:
1. **Stage 1 (Evolutionary Vector Oversampling):** Retrieve top-50 candidates via `<=>` cosine distance on `halfvec(1280)`.
2. **Stage 2 (Sequence Physics & Smith-Waterman Hybrid Operator `<~>`):** Filter out low-complexity regions (`sequence_entropy(seq) >= 2.0`), enforce length tolerance ($\pm 30\%$), and re-rank by:
   $$d_{\text{hybrid}}(A, B) = 0.7 \cdot (A_{\text{emb}} \iff B_{\text{emb}}) + 0.3 \cdot \frac{100}{100 + \text{SmithWaterman}(A_{\text{seq}}, B_{\text{seq}})}$$

| Combined Pipeline Configuration (`50,000` Proteins) | End-to-End p50 Latency (ANN + Entropy + Smith-Waterman `<~>`) | Hybrid Top-5 Recall (%) | Speedup vs. HNSW (`ef=100`) |
| :--- | :---: | :---: | :---: |
| **`PRISM (nprobe=auto)` + `pg_bio` `<~>`** | **11.97 ms** | 64.00% | **2.39× faster** |
| **`PRISM (nprobe=30)` + `pg_bio` `<~>`** | **17.59 ms** | **78.67%** | **1.63× faster** |
| **`HNSW (ef_search=40)` + `pg_bio` `<~>`** | 28.02 ms | 81.33% | 1.02× |
| **`HNSW (ef_search=100)` + `pg_bio` `<~>`** | 28.63 ms | 85.33% | 1.00× (Baseline) |
| **`HNSW (ef_search=200)` + `pg_bio` `<~>`** | 42.35 ms | 86.67% | 0.68× |
| **`IVFFlat (probes=10)` + `pg_bio` `<~>`** | 24.87 ms | 40.00% | 1.15× (poor recall) |

*(Note: Of the 17.59 ms total latency in `PRISM (nprobe=30) + pg_bio <~>`, only ~1.8 ms is spent in PRISM vector retrieval; the remaining ~15.8 ms is spent computing exact $O(L_1 L_2)$ Smith-Waterman dynamic programming alignments across the 50 oversampled protein sequences in Rust!)*

---

## 4. Full 2.34-Million Proteome Parallel Benchmark (`proteins` vs. `proteins_prism`)

Both tables are now live and running side-by-side in `bio_demo` (`2,340,766` UniProt proteins, `1280-D halfvec`):
- **`proteins` (Original, untouched):** `embedding halfvec(1280)` (`STORAGE EXTERNAL`) + `idx_protein_embedding` (`USING hnsw (embedding halfvec_cosine_ops)`)
- **`proteins_prism` (Parallel copy):** `embedding halfvec(1280) STORAGE PLAIN` + `idx_proteins_prism_embedding` (`USING prism (embedding vec16_cosine_ops)`, built in **184.68 s / 3.08 min** with `nlist=9216, fan_out=96, nlevels=2, auto_nprobe=48`)

### 4.1 Full-Proteome Storage Comparison (`2,340,766` Rows)

| Table & Index Stack | Main Heap | TOAST Heap | Vector Index Size | Total Footprint (Table + TOAST + Indexes) | Savings vs. HNSW Stack |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **`proteins` (`pgvector` HNSW + `EXTERNAL`)** | 1,053.66 MB | 6,112.26 MB | **8,772.13 MB** | **16,101.81 MB (15.72 GiB)** | Baseline (1.00×) |
| **`proteins_prism` (`pg_vectorsearch` PRISM + `PLAIN`)** | 6,095.75 MB | 1,105.83 MB | **902.26 MB** | **8,266.57 MB (8.07 GiB)** | **9.72× smaller vector index; 7.65 GiB total disk saved** |

### 4.2 Full-Proteome ANN Query Latency, Buffer Efficiency & Cosine Distance Quality

| Parallel Table & Index | Tuning Parameter | p50 Latency | Mean Latency | p95 Latency | QPS | Avg Top-10 Cosine Distance *(lower is better)* | Shared Buffer Blocks / Query | Speedup vs. Iso-Quality HNSW |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`proteins_prism` (PRISM)** | `nprobe = 10` | **9.27 ms** | 8.56 ms | 12.16 ms | **116.8** | 0.07118 | **525.2** | **21.28× faster** than HNSW `ef=40` |
| **`proteins_prism` (PRISM)** | `nprobe = 40` | **20.26 ms** | 18.57 ms | 27.00 ms | **53.8** | 0.06815 | 1,764.0 | **9.73× faster** than HNSW `ef=40` |
| **`proteins_prism` (PRISM)** | `nprobe = auto (48)` | **26.33 ms** | 25.21 ms | 39.61 ms | **39.7** | **0.06794** | 1,994.4 | **7.49× faster** than HNSW `ef=40` *(better dist)* |
| **`proteins_prism` (PRISM)** | `nprobe = 80` | **37.05 ms** | 33.31 ms | 47.70 ms | **30.0** | **0.06732** | 2,932.4 | **11.26× faster** than HNSW `ef=100` *(equal dist)* |
| **`proteins` (HNSW)** | `ef_search = 40` | 197.22 ms | 186.97 ms | 232.60 ms | 5.3 | 0.06798 | 3,304.2 | 1.00× (`ef=40` baseline) |
| **`proteins` (HNSW)** | `ef_search = 100` | 417.31 ms | 397.82 ms | 522.27 ms | 2.5 | 0.06733 | 6,706.2 | 1.00× (`ef=100` baseline) |
| **`proteins` (HNSW)** | `ef_search = 200` | 784.42 ms | 773.63 ms | 926.89 ms | 1.3 | 0.06682 | 12,123.4 | 1.00× (`ef=200` baseline) |

### 4.3 Full-Proteome Combined `pg_bio` `<~>` Hybrid Smith-Waterman Re-Ranking (`2.34M` Rows)

| Parallel Pipeline (`2,340,766` Proteins, Oversample $K=50$) | End-to-End p50 Latency | End-to-End Mean Latency | Avg Top-1 Hybrid `<~>` Score *(lower is better)* | Speedup & Quality vs. HNSW |
| :--- | :---: | :---: | :---: | :---: |
| **`proteins_prism` (`PRISM nprobe=40` + `pg_bio` `<~>`)** | **141.56 ms** | 151.67 ms | **0.27931** | **2.22× faster** than HNSW `ef=40` & **better hybrid match** |
| **`proteins_prism` (`PRISM auto=48` + `pg_bio` `<~>`)** | **141.71 ms** | 153.65 ms | **0.27906** | **4.16× faster** than HNSW `ef=100` & **better hybrid match** |
| **`proteins` (`HNSW ef_search=40` + `pg_bio` `<~>`)** | 313.60 ms | 295.14 ms | 0.28165 | 1.00× (`ef=40` baseline) |
| **`proteins` (`HNSW ef_search=100` + `pg_bio` `<~>`)** | 589.52 ms | 594.75 ms | 0.28068 | 1.00× (`ef=100` baseline) |

---

## 5. New `pg_bio` Features Implemented in This Worktree (`src/lib.rs`)

1. **`evaluate_vector_index_tier(row_count, dim, use_halfvec, ram_budget_mb) -> VectorIndexRecommendation`**:
   - Calculates exact TOAST threshold requirements (`requires_storage_plain`), compares estimated index sizes across HNSW, IVFFlat, PRISM Default, and PRISM Compact, and recommends the optimal index access method (`hnsw` for $<10\text{k}$ rows; `prism` with `vec16_cosine_ops` / `vec32_cosine_ops` for $\ge 10\text{k}$ rows).
2. **`pg_bio_search_homologs_prism(...)`**:
   - Fixes a subtle left-hand cast issue in `pg_bio_search_homologs` (keeping `v_target_emb halfvec` during Stage-1 index scan so `halfvec` HNSW and `vec16` PRISM indexes are never bypassed by an implicit `::vector` cast, and casting `c.embedding::vector` only on the 50 oversampled rows in Stage-2 `<~>` re-ranking).
   - Dynamically sets transaction-local `prism.nprobe` via `set_config('prism.nprobe', p_prism_nprobe::text, true)` when `pg_vectorsearch` is loaded.
