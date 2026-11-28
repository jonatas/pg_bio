# Layer 23: Sequence-Encoded `3Di` Conformational Entropy, `FlexProfile` Cryptic Pocket Breathing Gate, 128-Bit 64-Hot `3Di` `POPCNT` Sketch & `RangeFinder` 8-Property GiST Carbohydrate/PFAS Calibration

In structural enzymology, allosteric drug discovery, and extremophile bioprocess engineering, a static crystal or AlphaFold coordinate snapshot hides two critical physical truths:

1. **Rigid Orthosteric Active Sites vs. Breathing Cryptic Pockets vs. Disordered Loops (*ProtProfileMD* & *AlloDyn* 2026):**
   As demonstrated by Lüth et al. (*"ProtProfileMD: Predicting Conformational Dynamics and 3Di-FlexProfiles Directly from Sequence and Structure"*, *bioRxiv* Jan 2026, [`doi:10.64898/2026.01.21.700698`](https://doi.org/10.64898/2026.01.21.700698)) and *AlloDyn* (*bioRxiv* May 2026, [`doi:10.64898/2026.05.22.727284`](https://doi.org/10.64898/2026.05.22.727284)), a catalytic active site must pre-organize its nucleophilic triad in a **rigid conformational state** ($H_i^{3\text{Di}} < 0.28$, $\Delta G_{\text{open}} \approx 0.0\text{ kcal/mol}$) to avoid paying an entropic penalty during transition-state stabilization. Conversely, a **cryptic allosteric pocket** (such as the 5-LOX triterpenoid hinge or a kinase DFG-out pocket) breathes between closed and open Foldseek `3Di` states inside an intermediate **cryptic breathing band** ($0.30 \le H_i^{3\text{Di}} \le 0.70$, $\Delta G_{\text{open}} \in [0.8, 2.4]\text{ kcal/mol}$), while unstructured surface coils exceed $H_i^{3\text{Di}} > 0.70$.
2. **Twilight-Zone Structural Homology Pre-Filtering via 128-Bit 64-Hot `3Di` `POPCNT` (*ProstT5* / *Foldseek* / *SaProt*):**
   Below 20% sequence identity, 1D amino-acid $k$-mer overlap collapses to noise. Projecting entropy-weighted 20-letter `3Di` structural 3-mers ($w_i = 1 - 0.55 H_i^{3\text{Di}}$, so rigid fold cores dominate over noisy disordered tails) into a **strictly 64-hot 128-bit (`16-byte`) sketch** (`popcount == 64`) turns structural homology pre-filtering into **2 hardware `POPCNT` instructions** (`>50,000,000 comparisons/sec/core`).
3. **RangeFinder 8-Property GiST Carbohydrate (`Fsp3 = 1.00`) & Perfluoroalkyl (`Arom = 0, Fsp3 >= 0.80`) Calibration:**
   Decoupling aliphatic hydrophobicity (`f_hydro`) from mandatory aromatic ring floors (`arom_min <= 0.0`) and incorporating polar-polyol H-bonding density (`polar_polyol_drive`) in `predict_pocket_ranges` ensures 100% 8-property GiST containment recall on fully saturated carbohydrates (`Fsp3 = 1.00`, e.g., glucose, glucosamine) and non-aromatic perfluoroalkyls (`Arom = 0.0, Fsp3 = 0.80`, e.g., PFPeA/PFOA).

---

## 🔬 The Five New SQL Primitives in Layer 23

### 1. `bio_encode_3di_flex_profile(pdb_or_seq, md_b_factor_scale DEFAULT 1.0) -> bytea`
Encodes any 3D PDB or 1D amino-acid sequence into an inline zero-TOAST `"PGFP"` (`pg_bio` FlexProfile) `bytea` container (`32 + 4 * L` bytes; e.g., `1,232 B < 2,000 B` for a 300-aa protein, **0 B TOAST**), storing a 32-byte header (including the 16-byte 64-hot `3Di` `POPCNT` sketch at bytes `16..32`) followed by 4-byte/residue quantized `(dom_3di, h3di_q8, top2_3di, top2_prob_q8)` tuples.

### 2. `bio_peek_flex3di_header(pgfp_bytea)`
Inspects the 32-byte `"PGFP"` header in **$O(1)$ time** without unpacking the per-residue array:

```sql
SELECT n_residues, mean_h3di, cryptic_fraction, rigid_fraction,
       idr_fraction, source_mode, byte_size, octet_length(sketch_128b) AS sketch_bytes
FROM bio_peek_flex3di_header(
    bio_encode_3di_flex_profile('MIVLFWYCMLIVFWYCAAIVLFWYCMLIVFWYC', 1.0)
);
```

### 3. `bio_3di_structural_tokens(pdb_or_seq)`
Emits the 20-letter Foldseek/ProstT5 `seq_3di` string, the interleaved 400-state **SaProt** token string (`UPPERCASE AA` + `lowercase 3di`, e.g., `MaIdVk...`), `mean_h3di`, and the 16-byte 64-hot `sketch_128b`.

### 4. `bio_3di_sketch_similarity(sketch_or_pgfp_a, sketch_or_pgfp_b)`
Computes `shared_bits`, `cosine_64hot`, `tanimoto_64hot`, `hamming_dist`, and `is_twilight_homolog` in 2 hardware `POPCNT` instructions (`u64::count_ones`), accepting either raw 16-byte sketches or full `"PGFP"` blobs.

### 5. `bio_detect_cryptic_flex_pockets(pdb_text, flex_profile_bytea DEFAULT NULL, max_pockets DEFAULT 5)`
Couples Layer 19 multi-pocket SAS detection (`bio_detect_top_pockets`) with per-residue `3Di` conformational entropy ($H_i^{3\text{Di}}$) to compute each pocket's `mean_lining_h3di`, `max_hinge_h3di`, `cryptic_breathing_score`, `delta_g_open_kcal` ($\Delta G_{\text{open}}$), and `pocket_thermodynamic_regime` (`RIGID_ORTHOSTERIC_CLEFT`, `CRYPTIC_BREATHING_POCKET`, or `DISORDERED_NON_POCKET_LOOP`).
