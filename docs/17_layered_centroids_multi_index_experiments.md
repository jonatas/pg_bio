# Applying `pg_vectorsearch` Engineering to `pg_bio`’s Other Indices: Layered Centroids, SOAR, and Register-Direct Error-Bounded Fastscan

**Branch / Worktree:** `/Users/jonatas/code/pg_bio_vectorsearch` (`feat/layered-centroids-multi-index`)  
**Benchmark Script:** [`scripts/benchmark_layered_centroids_multi_index.py`](../scripts/benchmark_layered_centroids_multi_index.py)  
**Benchmark Data:** [`scripts/benchmark_layered_centroids_multi_index_results.json`](../scripts/benchmark_layered_centroids_multi_index_results.json)

---

## 1. Core Systems Engineering Learnings Extracted from `pg_vectorsearch`

Inspecting the C23 source of `timescale/pg_vectorsearch` (`src/index/centroid_search.c`, `src/index/centroid_build.c`, `src/quant/rabitq.c`, `src/quant/fastscan.c`, and `src/pg/amcache.c`) reveals four universal index-engine principles that generalize beyond dense float vectors to **every spatial, chemical, and structural index in `pg_bio`**:

1. **Hierarchical Centroid / Tree Beam Routing (`O(log N)` vs. Flat Scans):**  
   Instead of scanning a flat 1D interval or full table, organize the domain space into a shallow 2-level or $L$-level tree of centroids/cells and route queries top-down with a configurable beam (`nprobe` / `max_ranges`).
2. **SOAR Spill-Orthogonal Secondary Cluster Assignment (`soar_lambda`):**  
   Single-cluster partitioning always suffers from the **boundary problem** (a query near a hyperplane boundary misses neighbors just across the boundary). By projecting items onto orthogonal Walsh-Hadamard axes and flipping the coordinate with the smallest margin to the decision hyperplane (`0.0`), every boundary item gets a deterministic secondary `soar_centroid` that eliminates boundary false negatives at low `nprobe`.
3. **Two-Stage Quantized Fastscan with Provable Theoretical Error Bounds (RaBitQ Pattern):**  
   Compress complex high-dimensional structures into tiny fixed-width bit/byte codes (`16–24 bytes`) that carry a **strict mathematical bound** ($d_{\text{exact}} \ge d_{\text{lower\_bound}}$ or $\text{Sim}_{\text{exact}} \le \text{Sim}_{\text{upper\_bound}}$). Prune `90%–99.9%` of candidates in Stage 1 without false negatives before running expensive Stage-2 exact math.
4. **Register-Direct Zero-Deserialization Scoring (`STORAGE PLAIN` / Native `bigint` Registers):**  
   Avoid varlena TOASTing and CBOR/JSON composite struct deserialization inside hot scan loops by packing quantized codes into native 64-bit `bigint` CPU registers passed by value across the C/Rust FFI boundary.

---

## 2. Experiment 1: 3D Morton Z-Order Atomic Spatial Index (`z_order_layered_sphere_ranges`)

### The Problem with Naive 1-Range Morton Scans
In `pg_bio`, `z_order_encode(x, y, z)` maps `[-500, +500] Å` to `[0, 1000000]` (21 bits per axis) and interleaves bits into a 63-bit Morton code. Whenever a 3D binding-pocket sphere `(cx, cy, cz, R)` crosses a major power-of-2 octree plane (such as `24.288 Å` where `scaled = 524,288 = 2^19`, or `(0, 0, 0)`), `min_z` and `max_z` sit in different high-order octants. A single `WHERE z_index BETWEEN min_z AND max_z` scan traverses massive regions of dead 3D space (**Morton Z-jumps**), scanning `30,181` candidate atoms on average out of `200,000` just to find `~110` true pocket atoms.

### The `pg_vectorsearch`-Inspired Solution
`z_order_layered_sphere_ranges(center_x, center_y, center_z, radius, max_ranges)` descends the 21-level 3D Morton octree top-down:
- **Lower Error Bound Pruning ($d_{\min}^2 > R^2$):** Immediately prunes any 3D octant cube whose closest corner/face is farther than `radius` (`0%` false negatives, `100%` recall).
- **Upper Error Bound Interior Guarantee ($d_{\max}^2 \le R^2$):** Marks any octant cube completely inside the sphere with `fully_contained = true`, allowing SQL to skip the Pythagorean `sqrt()` recheck on interior atoms!
- **Contiguous Run Merging & Gap Coalescing:** Merges adjacent Z-intervals down to `max_ranges` (`24` or `64`) for a compact multi-range B-tree scan.

### Live Benchmark Results (`200,000` 3D Protein Atoms, `20` Pocket Queries)

| 3D Spatial Index Strategy | p50 Latency | Mean Latency | Avg Candidate Atoms Scanned | Candidate Reduction | Mean Speedup | Exact Recall |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Flat 1-Range Morton Scan (`BETWEEN min_z AND max_z`)** | 6.25 ms | 9.89 ms | 30,181.0 atoms | 1.00× | 1.00× | 100.0% |
| **Layered Octree Centroid Router (`max_ranges = 24`)** | **1.64 ms** | **2.33 ms** | **545.8 atoms** | **55.30× fewer** | **4.25× faster** | **100.0%** |
| **Layered Octree Centroid Router (`max_ranges = 64`)** | **0.98 ms** | **0.96 ms** | **407.6 atoms** | **74.05× fewer** | **10.25× faster** | **100.0%** |

---

## 3. Experiment 2: 1024-Bit ECFP4 Cheminformatics Tanimoto Index (`ecfp4_register_descriptor` & `ecfp4_fastscan_ub_raw`)

### The Problem with Full `Molecule % Molecule` Scans
Evaluating `mol % query_mol >= 0.70` across `10,000` molecules requires a sequential scan (`529.65 ms`) because every row detoasts and deserializes the full `Molecule` struct (`Vec<MolAtom>`, `Vec<MolBond>`, `Vec<u64>`).

### The `pg_vectorsearch`-Inspired Solution
`ecfp4_register_descriptor(smiles)` extracts a **24-byte register-direct Fastscan sketch** (`popcount: i16`, `block_counts: i64` packing eight 8-bit sub-block popcounts $[c_0..c_7]$, and `sketch_lo: i64, sketch_hi: i64` packing a 128-bit 8:1 OR-folded bit sketch) plus 2-level Walsh-Hadamard bit-density centroids (`l1_centroid`, `l2_centroid`, `soar_centroid`).

`ecfp4_fastscan_ub_raw(...)` computes a **100% mathematically guaranteed Upper Bound on 1024-bit Tanimoto similarity** purely in CPU registers:
$$U_{\cap}(A, B) = \min\!\left(\sum_{k=0}^7 \min(c_{A,k}, c_{B,k}),\; |A| - \text{popcnt}(S_A \ \& \ \neg S_B),\; |B| - \text{popcnt}(S_B \ \& \ \neg S_A)\right)$$
$$T_{\text{exact}}(A, B) \le T_{\text{ub}}(A, B) = \frac{U_{\cap}(A, B)}{|A| + |B| - U_{\cap}(A, B)}$$

### Live Benchmark Results (`10,000` Diverse Molecules, Tanimoto $\tau = 0.70$)

| Cheminformatics Tanimoto Strategy | p50 Latency | Mean Latency | Avg Survivors Rechecked | Candidates Pruned | p50 Speedup | Recall |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Baseline Full `Molecule % Molecule` Scan** | 529.65 ms | 528.65 ms | 10,000.0 mols | 0.00% | 1.00× | 100.0% |
| **Register Fastscan Error Bound (`ecfp4_fastscan_ub_raw`, Exact)** | **16.98 ms** | **16.55 ms** | **11.2 mols** | **99.89% pruned!** | **31.19× faster** | **100.0%** |
| **Layered Centroid + SOAR + Register Fastscan Bound** | **3.35 ms** | **3.35 ms** | **6.7 mols** | **99.93% pruned!** | **158.15× faster** | 75.17% |

---

## 4. Experiment 3: 9D RangeFinder Active-Site Pocket Index (`pocket_quantize_fastscan` & `pocket_envelope_iou_upper_bound`)

### The Problem with Top-$K$ 9D `RangeFinderPocket` Scans
While GiST (`idx_swarm_pockets_8prop_gist`) supports boolean `@>` containment, finding the top-10 pockets ordered by continuous 8-property `pocket_range_iou` across `4,856` 3D AlphaFold pockets in `swarm_pocket_ranges` requires an $O(N)$ scan (`103.13 ms`) that deserializes `RangeFinderPocket` on every row.

### The `pg_vectorsearch`-Inspired Solution
1. **16-Byte (`2 × i64`) Register Interval Quantization (`pocket_quantize_fastscan`):** Quantizes all 8 lower bounds into `lo_code: i64` (floor byte bins `0..255`) and all 8 upper bounds into `hi_code: i64` (ceil byte bins `0..255`), enabling a 100% guaranteed register-direct IoU upper bound (`pocket_fastscan_iou_ub_raw`).
2. **2-Level Balanced Pocket Centroid Envelopes (`pocket_envelope_iou_upper_bound`) + SOAR:** Partitions pockets into 37 balanced Level-2 clusters via Walsh-Hadamard rotation + SOAR secondary assignment, storing each cluster's bounding envelope `(env_mins, env_maxs, env_min_widths)` in a tiny 37-row centroid table and routing top-$K$ queries to the top `nprobe` clusters.

### Live Benchmark Results (`4,856` Real 3D AlphaFold Pockets from `swarm_pocket_ranges`)

| 9D RangeFinder Pocket Top-10 Strategy | p50 Latency | Mean Latency | p50 Speedup | Recall@10 |
| :--- | :---: | :---: | :---: | :---: |
| **Baseline Composite `RangeFinderPocket` Full Scan** | 103.13 ms | 103.67 ms | 1.00× | 100.0% |
| **16B Register Fastscan Bound + Raw Float8 IoU (Full Table)** | **7.79 ms** | **7.52 ms** | **13.24× faster** | **100.0%** |
| **Layered Centroid Envelope + SOAR (`nprobe = 2`)** | **5.05 ms** | **4.78 ms** | **20.42× faster** | **96.5%** |
| **Layered Centroid Envelope + SOAR (`nprobe = 4`)** | **5.79 ms** | **5.64 ms** | **17.82× faster** | **97.5%** |
| **Layered Centroid Envelope + SOAR (`nprobe = 8`)** | **6.80 ms** | **6.70 ms** | **15.16× faster** | **99.5%** |
