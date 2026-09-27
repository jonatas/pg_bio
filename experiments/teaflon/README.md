# The Teaflon Multiomics Project

## Overview
The **Teaflon** project is a synthetic biology initiative aimed at reinforcing human biological matrices (like the keratin in nails and hydroxyapatite in bones) using a synthetic bio-ceramic protein: **Syn-TFLN** (Bio-Teflon). 

By leveraging the multiomics capabilities natively built into `pg_bio` (PostgreSQL), we simulate the complete pipeline of engineering this human augmentation.

---

## Step 1: Genomics (The Insertion Problem)
**Goal:** Find a safe location in the human genome to insert the synthetic `Syn-TFLN` gene so it only activates in specific tissues (nails/bones) without disrupting critical biological functions (like tumor suppressors).

**pg_bio Implementation:** 
* Utilizes PostgreSQL native `int4range` types to map the 1-Dimensional coordinates of chromosomes.
* Uses a `GiST` index to execute lightning-fast spatial overlap queries (`&&`) and adjacency queries (`<->`).
* **Validation Criteria:** The gene must be inserted adjacent to the `KRT14` (Keratin 14) promoter and must have zero overlap with `TP53` or `BRCA1`.

---

## Step 2: Proteomics / Spatial Epigenomics (The Folding Problem)
**Goal:** Once the DNA is transcribed, the resulting amino acid sequence must fold into a dense beta-barrel structure with outward-facing Arginine residues that can chemically latch onto the human nail matrix.

**pg_bio Implementation:**
* Utilizes **Morton Coding (Z-Order Curves)** via the `z_order_encode(x, y, z)` Rust function in `pg_bio`.
* 3D atomic point clouds are mathematically unrolled into a 1D B-Tree index.
* **Validation Criteria:** A 3D bounding box spatial query must successfully isolate the `ARG_BIND` pocket residues at coordinates `(X:5-15, Y:-5-5, Z:20-30)`.

---

## Step 3: Transcriptomics (The Dosage Problem)
**Goal:** If the body produces too much Teaflon, the nail/bone becomes brittle. If too little, it provides no reinforcement. We must simulate and analyze the single-cell RNA-Seq expression levels across thousands of cells to ensure dosage hits the "Goldilocks Zone".

**pg_bio Implementation:**
* Utilizes **Compressed Sparse Row (CSR)** architecture using parallel PostgreSQL arrays (`gene_indices int[]`, `expression_values real[]`) to store massive, sparse single-cell matrices without blowing up storage.
* **Validation Criteria:** The average expression dosage of the `Syn-TFLN` gene (Gene ID: 9999) across all keratinocyte cells must fall within the Goldilocks threshold of **15 to 50** RNA copies per cell.
