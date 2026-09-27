# pg_bio: Multiomics Architecture

`pg_bio` extends PostgreSQL to natively support multiomics data processing, enabling the integration of Genomics, Proteomics, and Transcriptomics within a single database engine.

## 1. Genomics (DNA) via Range Types
Genomic annotations (genes, promoters, variants) are mapped to 1-Dimensional chromosomes.
*   **Infrastructure:** PostgreSQL `int4range` (Integer Ranges).
*   **Indexing:** `GiST` (Generalized Search Tree).
*   **Capabilities:** Execute millisecond overlap queries (`&&`) to map mutations to genes, or find adjacent genomic features (`<->`) across billions of base pairs.

## 2. Proteomics & Epigenomics (3D Structural Space)
Biological molecules fold into 3D space. Querying spatial atomic relationships is critical for drug discovery and binding pocket analysis.
*   **Infrastructure:** Custom Rust functions (`z_order_encode`).
*   **Indexing:** B-Tree indexing on 1D Morton Codes (Z-Order curves).
*   **Capabilities:** Compresses 3D point clouds (X, Y, Z) into 1D integers, allowing standard B-Trees to execute sub-millisecond 3D bounding-box spatial queries.

## 3. Transcriptomics (Single-Cell Expression)
Single-cell RNA sequencing generates massive, sparse matrices (e.g., 20,000 genes x 1,000,000 cells), where 85%+ of data points are zeroes.
*   **Infrastructure:** Compressed Sparse Row (CSR) via PostgreSQL parallel arrays (`int[]` for gene indices, `real[]` for expression counts).
*   **Capabilities:** Bypasses billions of zero-values to store sparse expression profiles. Utilizes `UNNEST() WITH ORDINALITY` to dynamically unpack and aggregate RNA dosage on the fly.
