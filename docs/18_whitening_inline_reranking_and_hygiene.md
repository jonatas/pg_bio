# Phase 18: Anisotropy Whitening (`Know Thy Neighbour`), Inline `BioVecInt8G128` Varlena, Structure/PAINS Ingestion Hygiene & Stage-2.5 (`EnzRef` / `STORM`) Re-Rankers

This document details the four high-leverage pre-conditioning, inline storage, ingestion hygiene, and Stage-2.5 post-filter re-ranking primitives added to [`src/lib.rs`](../src/lib.rs) on top of `pg_bio`'s 2-level layered centroid, 32-wide SIMD register `fastscan` (`VPSHUFB` / ARM64 NEON `vtbl`), 3D Morton octree, ECFP4, and `RangeFinder` `numrange` GiST stack.

---

## 1. Anisotropy Centering & Whitening Pre-Transform (`Know Thy Neighbour`) + Heap Inline `BioVecInt8G128`

### The Biological & Mathematical Problem
Raw protein language model embeddings (`ESM-2`, `ProstT5`) occupy a narrow anisotropic cone in $\mathbb{R}^{1280}$ where a handful of SwiGLU channels dominate variance (mean pairwise cosine $> 0.95$). In 1-bit RaBitQ quantization, the alignment factor $\langle \hat{r}, \hat{q}_b \rangle = \frac{\|r\|_1}{\sqrt{D}\|r\|_2}$ drops sharply when energy concentrates in a few rogue dimensions, loosening the theoretical error bound:

$$\varepsilon_{\text{bound}}(r) = \frac{\sqrt{1 - \langle \hat{r}, \hat{q}_b \rangle^2}}{\langle \hat{r}, \hat{q}_b \rangle} \cdot \|r\|_2$$

Additionally, a raw 1280-D `f32` vector (`5,124 B`) or `f16` `halfvec` (`2,564 B`) exceeds PostgreSQL's `~2,000 B` `TOAST_TUPLE_THRESHOLD`, causing out-of-line TOAST I/O during Stage-2 exact re-ranking.

### The Native Rust / SIMD Solution
1. **Whitening Pre-Transform (`bio_fit_whitening_stats`, `bio_fit_zscore_whiten`, `bio_whiten_embedding`, `bio_prepare_fastscan_vector`):**
   - Computes exact per-channel mean $\mu \in \mathbb{R}^D$, top-$k$ principal cone directions $\{u_j\}_{j=1}^k$ via deterministic power iteration + Gram-Schmidt deflation, and inverse standard deviation scales $s_i = (\sigma_i^2 + \varepsilon)^{-1/2}$, packed into a deterministic `PGBW` binary `bytea`.
   - Applies $x_{\text{out}} = \text{normalize}\!\left(s \odot \left((x - \mu) - \sum_{j=1}^k \alpha_j \langle x - \mu, u_j \rangle u_j\right)\right)$ using hardware SIMD (`std::arch::aarch64` NEON `vfmaq_f32`/`vfmsq_f32` and `std::arch::x86_64` AVX2 `_mm256_sub_ps`/`_mm256_mul_ps`).
2. **Zero-CBOR Fixed-Header C-Varlena `BioVecInt8G128` (`1,300 B` payload / `1,304 B` varlena):**
   - Divides 1280-D vectors into 10 groups of 128 `i8` channels with one IEEE 754 `f16` scale per group (`10 * 2B + 1280 * 1B = 1,300 B` payload + `4 B` `SET_VARSIZE` header = `1,304 B` on disk).
   - Stays strictly below the `~2,000 B` `TOAST_TUPLE_THRESHOLD` (`pg_column_size = 1304`), keeping heap tuples 100% inline with zero TOAST reads and evaluating dot/cosine distance in hardware SIMD `i32` accumulators (`vmull_s8` + `vpadalq_s16` on ARM64 NEON; `_mm256_cvtepi8_epi16` + `_mm256_madd_epi16` on AVX2).

### SQL Example
```sql
-- 1. Fit whitening statistics on a corpus of embeddings and audit RaBitQ error bounds
WITH corpus AS (
    SELECT ARRAY(
        SELECT CASE WHEN d = 1 THEN 45.0 + sin(i * d)::float4
                    WHEN d <= 4 THEN 18.0 + cos(i * d)::float4
                    ELSE 0.35 * sin(i * 17 + d)::float4 END
        FROM generate_series(1, 40) i, generate_series(1, 128) d
    ) AS flat_vecs
),
stats AS (
    SELECT bio_fit_whitening_stats(flat_vecs, 128) AS w_stats, flat_vecs FROM corpus
)
SELECT 'raw' AS mode, (bio_rabitq_eps_bound_audit(flat_vecs, NULL, 0.45, 128)).* FROM stats
UNION ALL
SELECT 'whitened' AS mode, (bio_rabitq_eps_bound_audit(flat_vecs, w_stats, 0.45, 128)).* FROM stats;

-- 2. Encode 1280-D embeddings into 1,304-byte inline BioVecInt8G128 and rerank in SIMD i32
SELECT
    pg_column_size(bio_vec_int8_g128_encode(ARRAY(SELECT sin(g::float8)::float4 FROM generate_series(1, 1280) g))) AS inline_varlena_bytes,
    bio_vec_int8_g128_cosine_distance(
        bio_vec_int8_g128_encode(ARRAY(SELECT sin(g::float8)::float4 FROM generate_series(1, 1280) g)),
        bio_vec_int8_g128_encode(ARRAY(SELECT sin(g::float8 + 0.15)::float4 FROM generate_series(1, 1280) g))
    ) AS simd_i32_cosine_dist;
```

---

## 2. Ingestion Hygiene: Structure Sanitization & OpenSMILES PAINS / Aggregator Screen

### The Biological & Chemical Problem
1. **Defective Experimental & Predicted PDB/mmCIF Structures:** Missing backbone carbonyl oxygens (`O`), missing $\beta$-carbons (`CB`), or truncated catalytic sidechains (`SG`, `HIS` imidazole ring), unannotated `< 2.2 Å` disulfide bonds (`CYS` instead of `CYX`), and ambiguous histidine protonation states (`HIS` instead of `HID`, `HIE`, or `HIP`) corrupt downstream 3D Morton octree active-site and hydrogen-bond audits.
2. **Pan-Assay Interference Compounds (PAINS) & Colloidal Aggregators:** Redox-cycling quinones/catechols, covalent rhodanines/maleimides/isothiazolones, azo chromophores, acylhydrazones, and high-$\text{cLogP}$ planar polyaromatic aggregators produce false-positive assay hits if not rejected at `INSERT` or `WHERE` filter time.

### The Native Rust Solution
- **`bio_sanitize_structure(pdb_or_mmcif_text)` & `bio_sanitize_structure_audit(pdb_or_mmcif_text)`:**
  - Single-threaded, panic-safe, zero-Rayon structure sanitizer supporting both PDB and mmCIF `_atom_site` records.
  - Reconstructs missing heavy atoms via Kabsch SVD / Horn 4×4 quaternion superposition (`det(R) = +1`) against canonical Engh & Huber local residue frames.
  - Pairs `< 2.2 Å` `CYS:SG`–`CYS:SG` bridges via Morton Z-order sorting into `CYX` (`SSBOND`).
  - Assigns `HID`, `HIE`, or `HIP` based on explicit hydrogens, $< 2.65\text{ Å}$ metal coordination (`ZN`, `FE`, `CU`, `MN`, `CO`, `NI`, `MG`), and $< 3.5\text{ Å}$ carboxylate salt-bridge / H-bond acceptor geometry.
- **`bio_is_pains_or_aggregator(smiles)`:**
  - Pre-screens a 16-bit chemical feature bitmask (`N`, `O`, `S`, `=`, 5-ring, 6-ring, aromaticity) in $O(|V| + |E|)$ before evaluating deterministic ring/bond structural alerts and Shoichet colloidal aggregation thresholds.

### SQL Example
```sql
-- 1. Audit and repair a raw PDB structure with missing atoms, a 2.04 A disulfide, and Zn/salt-bridged Histidines
SELECT repaired_atoms, cyx_disulfide_bonds, hid_count, hie_count, hip_count
FROM bio_sanitize_structure_audit($$
ATOM      1  N   CYS A   1      -1.458   0.000   0.000  1.00  0.00           N
ATOM      2  CA  CYS A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   CYS A   1       0.530   1.425   0.000  1.00  0.00           C
ATOM      4  SG  CYS A   1      -0.625  -2.548  -0.924  1.00  0.00           S
ATOM      5  N   CYS A   2       2.500   0.000   0.000  1.00  0.00           N
ATOM      6  CA  CYS A   2       3.958   0.000   0.000  1.00  0.00           C
ATOM      7  C   CYS A   2       4.488   1.425   0.000  1.00  0.00           C
ATOM      8  O   CYS A   2       4.587   2.388   0.771  1.00  0.00           O
ATOM      9  CB  CYS A   2       3.440  -0.774  -1.206  1.00  0.00           C
ATOM     10  SG  CYS A   2       1.415  -2.548  -0.924  1.00  0.00           S
$$);

-- 2. Filter out PAINS compounds before RangeFinder GiST or ECFP4 Tanimoto joins
SELECT name, smiles, bio_is_pains_or_aggregator(smiles) AS is_pains
FROM (VALUES
    ('p-Benzoquinone', 'O=C1C=CC(=O)C=C1'),
    ('Rhodanine',      'O=C1CSC(=S)N1'),
    ('Azobenzene',     'c1ccc(N=Nc2ccccc2)cc1'),
    ('Aspirin',        'CC(=O)Oc1ccccc1C(=O)O'),
    ('Caffeine',       'Cn1cnc2c1c(=O)n(C)c(=O)n2C')
) AS t(name, smiles);
```

---

## 3. Stage-2.5 CPU Functional & Affinity Re-Rankers (`EnzRef` + `STORM` L0)

### The Native Rust / SIMD Solution
1. **`EnzRef` Role-Contrast Scorer (`bio_enzref_role_contrast(enzyme_bg_vec, reactant_vec, product_vec) -> float4`):**
   - Computes reference-centered substrate deviation $\delta_{\text{ref}} = r - e$ and directed catalytic reaction contrast $\delta_{\text{cat}} = p - r$ in a single fused SIMD pass (`aarch64` NEON / `x86_64` AVX2):
     $$\text{EnzRef}(e, r, p) = \frac{\langle r - e, p - r \rangle}{\|r - e\|_2 \, \|p - r\|_2} \cdot \ln\!\left(1 + \frac{\|r - e\|_2 \, \|p - r\|_2}{\sqrt{D}}\right)$$
   - Vanishes (`0.0`) for non-catalytic dead binders ($p = r$) and reverses sign under reaction reversal ($r \leftrightarrow p$).
2. **`STORM`-Style Fast Affinity Descriptor Kernel (`bio_id128_hash`, `bio_smooth_assay_series`):**
   - `bio_id128_hash(accession_or_neighbor_id)` generates a deterministic 16-byte (128-bit) target code with **exact 64-hot weight** (`popcount == 64`) via seeded SplitMix64 Fisher-Yates bit selection.
   - `bio_smooth_assay_series(fp_array, raw_scores, sim_threshold)` performs hardware-`POPCNT` bitwise Tanimoto neighborhood smoothing across congeneric series while keeping **both `raw_score` and `smoothed_score`** (plus `neighbor_count` and `activity_cliff_delta`) exposed so SAR activity cliffs are never erased.

### SQL Example
```sql
-- 1. EnzRef Role-Contrast pairwise log10(kcat) ranking
SELECT
    bio_enzref_role_contrast(ARRAY[1.0, 0.0, 0.0, 0.0]::real[], ARRAY[2.0, 1.0, 0.0, 0.0]::real[], ARRAY[3.5, 2.5, 0.0, 0.0]::real[]) AS forward_catalysis,
    bio_enzref_role_contrast(ARRAY[1.0, 0.0, 0.0, 0.0]::real[], ARRAY[2.0, 1.0, 0.0, 0.0]::real[], ARRAY[2.0, 1.0, 0.0, 0.0]::real[]) AS dead_binder,
    bio_enzref_role_contrast(ARRAY[1.0, 0.0, 0.0, 0.0]::real[], ARRAY[2.0, 1.0, 0.0, 0.0]::real[], ARRAY[0.5, -0.5, 0.0, 0.0]::real[]) AS reverse_reaction;

-- 2. STORM bitwise Tanimoto smoothing with activity cliff preservation
SELECT * FROM bio_smooth_assay_series(
    ARRAY['Cc1ccccc1', 'CCc1ccccc1', 'CCCc1ccccc1', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'],
    ARRAY[5.1, 9.2, 5.3, 4.0]::float8[],
    0.25
);
```

---

## 4. `pgrx-bench` Verification Results (PostgreSQL 18.3, Apple Silicon ARM64 NEON)

Reproduce anytime inside PostgreSQL via `SELECT * FROM bio_run_phase4_benchmarks();`:

| Benchmark | Raw / Baseline | Optimized (`pg_bio`) | Gain / Ratio | Unit |
| :--- | :---: | :---: | :---: | :--- |
| **1a. `ESM-2` Cone Anisotropy (Mean Pairwise Cosine)** | `0.9963` | `-0.0151` | **65.98×** isotropy gain | mean pairwise cosine |
| **1b. RaBitQ 1-Bit Error Bound (`eps_bound`)** | `2.8108` | `1.0118` | **2.78× tighter** | `eps_bound` (lower is tighter) |
| **1c. 32-Wide SIMD `fastscan` Pruning Ratio** | `0.0000` (0.0%) | `0.5617` (56.17%) | **+56.17% pruned** | fraction of candidates pruned |
| **2a. 1280-D On-Disk Varlena Size (`f16` vs `BioVecInt8G128`)** | `2,564 B` (TOASTed) | `1,304 B` (100% inline) | **1.97× smaller** (`0` TOAST reads) | bytes (`pg_column_size`, threshold ~2,000 B) |
| **2b. Stage-2 100-Candidate Rerank Latency (`f32` vs `BioVecInt8G128` SIMD `i32`)** | `1,629.70 µs` | `813.97 µs` | **2.00× faster** | µs / 100 candidates (1280-D) |
| **3a. Stage-2.5 `bio_enzref_role_contrast` (1280-D SIMD)** | — | **52,970 rows/sec** | **18.8 µs / row** | rows / sec (1280-D triplets) |
| **3b. OpenSMILES `bio_is_pains_or_aggregator` Screen** | — | **11,995 molecules/sec** | **83.4 µs / molecule** | full SMILES parse + aromaticity + PAINS screen / sec |
