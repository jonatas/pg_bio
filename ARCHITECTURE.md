# pg_bio Architecture

`pg_bio` is a specialized PostgreSQL extension designed to natively understand protein physics, 3D spatial coordinates, and evolutionary embeddings. By bringing bioinformatics directly into the database engine, we eliminate the need for fragile pipelines of flat files (FASTA, PDB) and fragmented scripts, enabling scientists to run massive protein Computational models and spatial queries on local hardware.

## Core Concepts & Mechanisms

### 1. 3D Spatial Queries via Z-Order Indexing (`spiral` integration)
Currently, finding atoms within a certain Ångstrom distance requires loading the entire protein structure into memory. 
- **Implementation:** We will introduce a custom `3d_residue` type. By applying Z-order curves (Morton coding) similar to the `spiral` extension, we can interleave the 3D coordinates into a 1D index.
- **Result:** PostgreSQL can natively perform ultra-fast bounding-box queries, instantly finding spatial clashes, binding pockets, or calculating molecular docking surfaces using standard SQL `WHERE` clauses.

### 2. High-Dimensional Vector Storage for ESM Embeddings (`pgvector` + `pg_vectorsearch` PRISM)
Language models like ESM-2 output dense vectors (e.g., 1280 dimensions) representing the biological meaning of a sequence.
- **Implementation:** We combine `pgvector`'s `halfvec(1280) STORAGE PLAIN` and `vector` types with Tiger Data's [`pg_vectorsearch`](https://github.com/timescale/pg_vectorsearch) **PRISM** index access method (`USING prism (embedding vec16_cosine_ops)`), leveraging 1-bit RaBitQ quantization with theoretical error bounds, SOAR spill-orthogonal cluster replication, and hierarchical centroid routing.
- **Result:** Across 2.34 million 1280-D UniProt proteins, PRISM shrinks the vector index by **9.72x** (`8,772 MB` HNSW $\rightarrow$ `902 MB` PRISM), builds in **3.08 minutes** (down from 27 hours), and accelerates `<=>` cosine and `<~>` hybrid Smith-Waterman homology searches by **7.5x to 11.3x** (`26–37 ms` vs `197–417 ms`).

### 3. Recursive Query Plans for Molecular Graphing
A folded protein is fundamentally a graph, with amino acids as nodes and chemical bonds (covalent, hydrogen) as edges.
- **Implementation:** We will optimize `WITH RECURSIVE` queries in Postgres through custom extensions and operators that traverse these biological edges efficiently.
- **Result:** Tracing allosteric pathways—how a change on one side of a protein physically propagates to the other—becomes a graph traversal problem solved in milliseconds by the database planner.

### 4. Sparse Attention Matrix Storage (Custom TAM)
The massive memory requirements of models like ESMFold (e.g., 11GB) stem largely from the $N \times N$ attention maps.
- **Implementation:** We will build a custom Table Access Method specifically designed for Sparse Attention Matrices. Instead of loading the entire matrix into RAM, Postgres will intelligently fetch only the highly correlated blocks.
- **Result:** This drastically reduces the active memory footprint, allowing large protein models and complex structural alignments to run locally on machines with limited RAM (like 16GB laptops).

## Vision
To build a database that doesn't just store biological data—it **understands** it. `pg_bio` will turn PostgreSQL into a local, high-performance computational biology workstation.
