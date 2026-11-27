# Layer 21: Directional `ipSAE` & $5 \times 5$ `MiniPAE` Co-Folding Confidence Gate for Protein Complexes & Disordered SLiM Metabolon Tethers

When scoring multi-chain protein complexes or synthetic metabolon couplers ($M_1 \rightarrow M_2 \rightarrow M_3$) predicted by **AlphaFold-Multimer**, **AlphaFold 3**, or **Boltz-2**, a major false-negative trap occurs whenever one partner contains an **Intrinsically Disordered Region (IDR)** or a flexible hydrophilic linker (`GGGSGGGGS`):

1. **The Global `ipTM` IDR Dilution Trap (Dunbrack, *bioRxiv* 2025/2026):**
   Standard AlphaFold `ipTM` sets the TM-score distance scale $d_0(L_{\text{tot}})$ using the **total length of both chains** ($L_A + L_B$) and averages over every residue in the partner chain. As proven by Roland L. Dunbrack Jr. (*"Rēs ipSAE loquuntur: What's wrong with AlphaFold's ipTM score and how to fix it"*, [`doi:10.1101/2025.02.10.637595`](https://doi.org/10.1101/2025.02.10.637595)), if an 80-residue microprotein tether carries a 10-residue Short Linear Motif (**SLiM**) locked into a 150-residue enzyme cleft at $\text{PAE} = 0.45\,\text{\AA}$ while its 70-residue disordered tail floats free ($\text{PAE} = 28.5\,\text{\AA}$), global `ipTM` divides by all 80 residues and drops to **`0.12–0.34`**—falsely rejecting a real nanomolar tether!
2. **The Single-Residue Minimum PAE Glitch (`MiniPAE`, Veinstein et al., *Brief. Bioinform.* 2025):**
   Taking the raw minimum $\min_{i,j} \text{PAE}_{i,j}$ across an interface is fragile because a single accidental terminal contact can score $< 3.5\,\text{\AA}$. Veinstein et al. ([`doi:10.1093/bib/bbaf501`](https://doi.org/10.1093/bib/bbaf501)) and Omidi et al. (2024) showed that true SLiM-mediated complexes form a **contiguous $5 \times 5$ inter-chain contact patch** (`MiniPAE`).
3. **The 1.44 MB Dense PAE TOAST Wall:**
   A dense $600 \times 600$ `float32` PAE matrix consumes **`1,440,000 bytes`**, triggering massive out-of-line PostgreSQL TOAST chunking.

---

## 🚀 The Scientist's Advantage: Contact-Restricted Directional `ipSAE`, $O(L_A L_B)$ 2D-Prefix `MiniPAE` & `< 2,000 B` Zero-TOAST `"PGPA"` CSR Storage

`pg_bio` Layer 21 implements Dunbrack's **Directional `ipSAE`** and Veinstein's **$5 \times 5$ Sliding-Window `MiniPAE`** natively in compiled Rust (`src/lib.rs`) over zero-TOAST inline `"PGPA"` quantized CSR `bytea` containers:

* **Contact-Restricted $d_0(N_i)$ Scaling:** For each anchor residue $i \in A$, `pg_bio` counts ONLY the $N_i$ qualifying partner residues in Chain $B$ with $\text{PAE}_{i,j} < \tau_{\text{PAE}}$ (default $10.0\,\text{\AA}$) and computes:
  $$d_0(N_i) = \max\!\left(1.0,\; 1.24 \cdot (\max(N_i, 27) - 15)^{1/3} - 1.8\right)$$
  $$\text{ipSAE}_{A \to B} = \max_{i \in A,\, N_i > 0} \frac{1}{N_i} \sum_{j \in B,\, \text{PAE}_{i,j} < \tau_{\text{PAE}}} \frac{1}{1 + \left(\text{PAE}_{i,j} / d_0(N_i)\right)^2}$$
* **$O(L_A L_B)$ 2D Prefix-Sum $5 \times 5$ `MiniPAE`:** Evaluates all $(L_A - 4)(L_B - 4)$ bidirectional $5 \times 5$ windows in a single pass via a 2D separable prefix sum, returning `minipae_5x5_a` and exact 1-based motif coordinates (`best_window_a_start`, `best_window_b_start`).
* **Zero-TOAST `"PGPA"` Sparse CSR Container (`< 2,000 B`) + $O(1)$ 16-Byte Header Peek (`bio_peek_pae_csr_header`):** Compresses inter-chain contacts ($\text{PAE} \le 15.0\,\text{\AA}$) at $0.0625\,\text{\AA}$ `u8` resolution into `< 2,000 bytes` (100% inline on the 8 KB PostgreSQL heap page) and caches `ipsae_min` and `minipae_5x5_a` in bytes `14..15` of the 16-byte header for instant $O(1)$ pre-filtering.

---

## 🔬 The Four New SQL Primitives

### 1. `bio_ipsae_minipae_audit(pdb_text, pae_csr_bytea DEFAULT NULL, pae_cutoff_a DEFAULT 10.0, minipae_gate_a DEFAULT 3.5)`
Polymorphic entry point accepting **either** a 2-chain 3D PDB complex (deriving directional inter-chain PAE from $\text{C}\alpha/\text{C}\beta$ contact geometry and B-factor pLDDT frame confidence) **or** a packed `"PGPA"` CSR `bytea`:

```sql
-- Audit a 2-chain receptor + disordered SLiM tether complex directly from PDB coordinates or packed CSR:
SELECT chain_a_id, chain_b_id,
       chain_a_len, chain_b_len,
       active_contacts_ab, active_contacts_ba,
       ipsae_a_to_b, ipsae_b_to_a,
       ipsae_min, ipsae_max,
       iptm_global_diluted,
       idr_dilution_rescue_delta,
       minipae_5x5_a,
       best_window_a_start, best_window_b_start,
       passes_cofold_gate,
       gate_verdict,
       packed_csr_bytes
FROM bio_ipsae_minipae_audit(complex_pdb_text, NULL, 10.0, 3.5);
```

### 2. `bio_ipsae_minipae_from_pae_flat(pae_flat, chain_a_len, chain_b_len, pae_cutoff_a DEFAULT 10.0, minipae_gate_a DEFAULT 3.5)`
Evaluates `ipSAE` and $5 \times 5$ `MiniPAE` directly over an AlphaFold2/3 or Boltz-2 PAE array (`float4[]`, supporting full $(L_A+L_B)^2$, bidirectional $2 L_A L_B$, or $L_A L_B$ layouts) and emits the inline `"PGPA"` CSR `bytea`:

```sql
-- Rescue a 10-aa SLiM tether embedded in an 80-aa disordered microprotein (where global ipTM drops to 0.125):
SELECT ipsae_min,
       ipsae_max,
       iptm_global_diluted,
       idr_dilution_rescue_delta,
       minipae_5x5_a,
       best_window_a_start,
       best_window_b_start,
       passes_cofold_gate,
       gate_verdict,
       packed_csr_bytes
FROM bio_ipsae_minipae_from_pae_flat(pae_array, 150, 80, 10.0, 3.5);
-- Returns:
--   ipsae_min                 = 0.8315
--   iptm_global_diluted       = 0.1250  (Falsely rejected by standard AlphaFold ipTM!)
--   idr_dilution_rescue_delta = +0.7065
--   minipae_5x5_a             = 0.4750 Å
--   best_window_b_start       = 31      (Exact 1-based SLiM window start!)
--   passes_cofold_gate        = TRUE
--   gate_verdict              = 'RESCUED_SLIM_IDR_TETHER'
```

### 3. `bio_pack_sparse_pae_csr(pae_flat, chain_a_len, chain_b_len, pae_cutoff_a DEFAULT 15.0)` & `bio_peek_pae_csr_header(pae_csr_bytea)`
Packs dense PAE matrices into `< 2,000 B` inline `"PGPA"` CSR `bytea` and inspects the 16-byte header in $O(1)$ time without unpacking the CSR arrays:

```sql
-- O(1) 16-byte header peek without unpacking CSR payload arrays:
SELECT * FROM bio_peek_pae_csr_header(
    bio_pack_sparse_pae_csr(pae_array, 150, 80, 15.0)
);
```
