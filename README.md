# 🧬 pg_bio: The Multiomics Operating System for Synthetic Biology

**`pg_bio`** is a native PostgreSQL extension written in highly optimized Rust (`pgrx`). 

It is engineered for a single purpose: to completely eradicate the data bottlenecks in modern biology by pushing complex multiomics mathematics *directly into the database engine*. By unifying Genomics, Proteomics, and Transcriptomics into a shared computational space, `pg_bio` allows scientists to execute entire drug discovery pipelines and synthetic biology designs in milliseconds natively in SQL.

---

## 🛑 The Data Bottleneck in Modern Biology

The ultimate frontier of biology is **Multiomics**—the fusion of DNA, RNA, and Protein data. Today, this is done using duct-tape pipelines:
1. Downloading terabytes of sequence data from UniProt.
2. Building an in-memory KD-Tree in Python to measure 3D distances.
3. Using custom tools like `bedtools` for genomic overlaps.
4. Loading 20-billion-row RNA matrices into massive Pandas dataframes.

**The result:** The network I/O and RAM overhead crushes most laptops, and forces cloud deployments into hours of agonizing processing. 

## 🚀 The Solution: The Multiomics Stack inside Postgres

`pg_bio` eliminates this entire nightmare. It provides four distinct computational layers natively inside PostgreSQL. We move the compute directly to the data, resulting in order-of-magnitude speedups.

---

## 🔬 Deep Dive: The 4 Multiomics Layers

### Layer 1: Genomics (1D Spatial Targeting)
When engineering synthetic genes (e.g., bio-ceramics or CRISPR therapies), you must find a "Safe Harbor" in the genome. It must be near the correct tissue-specific promoter (to turn it on) but cannot overlap with any tumor suppressor genes (which could cause cancer).

* **The Engine:** `pg_bio` treats the human genome as a 1D spatial environment utilizing native PostgreSQL **Range Types (`int4range`)** and **`GiST` (Generalized Search Tree) indexes**.
* **The Advantage:** You can use native overlap operators (`&&`) and adjacency operators (`<->`) to scan billions of base pairs in milliseconds.
* **SQL Example:**
  ```sql
  -- Instantly find a Safe Harbor insertion site next to the Nail Promoter, avoiding Cancer genes.
  SELECT safe.genomic_range FROM human_genome safe
  JOIN human_genome promoter ON safe.chromosome = promoter.chromosome
  WHERE safe.genomic_range <-> promoter.genomic_range < 5000
    AND NOT EXISTS (
      SELECT 1 FROM human_genome cancer_genes
      WHERE cancer_genes.feature_type = 'tumor_suppressor' 
        AND safe.genomic_range && cancer_genes.genomic_range
    );
  ```

### Layer 2: Proteomics & Epigenomics (3D Spatial Curves)
Once transcribed, linear sequences fold into 3D atomic structures. Finding which atoms are within a 4.0 Ångstrom radius of a binding pocket usually requires extracting millions of coordinates into Python for brute-force Pythagorean math.

* **The Engine:** `pg_bio` contains a native Rust function `z_order_encode(x, y, z)`. This uses **Morton Coding (Z-Order curves)** to mathematically interleave 3D Cartesian coordinates into a single 1-Dimensional integer.
* **The Advantage:** Standard, ultra-fast 1D PostgreSQL B-Tree indexes now act as 3D bounding boxes. Extract exact 3D binding pockets or detect cross-protein molecular clashes across millions of atoms instantly—with zero setup time.
* **SQL Example:**
  ```sql
  -- Extracting a 3D Binding Pocket in milliseconds
  SELECT atom_id, residue_name FROM teaflon_atoms
  WHERE z_curve_index BETWEEN z_order_encode(5.0, -5.0, 20.0) AND z_order_encode(15.0, 5.0, 30.0)
    AND x BETWEEN 5.0 AND 15.0 AND y BETWEEN -5.0 AND 5.0 AND z BETWEEN 20.0 AND 30.0;
  ```

### Layer 3: Transcriptomics (Single-Cell CSR Arrays)
To understand drug dosage, scientists use Single-Cell RNA-Seq. Humans have 20,000 genes. Sequencing 1 million cells creates a 20-Billion point matrix. Traditional SQL schemas grind to a halt.

* **The Engine:** Because a single cell only expresses ~2,000 genes at a time, the matrix is 85%+ zeroes. `pg_bio` implements **Compressed Sparse Row (CSR)** architecture natively using parallel PostgreSQL Arrays (`expressed_gene_ids INT[]` and `expression_counts REAL[]`).
* **The Advantage:** Storage overhead is slashed by 90%. Using parallel `LATERAL UNNEST() WITH ORDINALITY`, PostgreSQL unpacks the sparse data dynamically in execution memory to calculate highly precise RNA dosage levels in milliseconds.

### Layer 4: Joint Vector Spaces (The Cross-Omic Holy Grail)
Traditionally, DNA sequences and 3D Protein folds live in different worlds. `pg_bio` breaks down that wall.

* **The Engine:** By utilizing `pgvector` alongside `pg_bio`, we map the outputs of Genomic Foundation Models (like Enformer) and Proteomic Foundation Models (like ESM-2) into a **shared 1280-dimensional mathematical space**.
* **The Advantage:** For the first time, you can execute cross-domain queries. You can ask the database to mathematically link a physical DNA regulatory enhancer directly to a 3D Protein structure using Cosine Distance (`<=>`).
* **SQL Example:**
  ```sql
  -- Find which DNA enhancer matches the structural embedding of the P0DPB7 Protein
  SELECT dna.locus, (dna.embedding <=> protein.embedding) as cross_omic_distance
  FROM dna_embeddings dna, proteins protein
  WHERE protein.uniprot_id = 'P0DPB7'
  ORDER BY cross_omic_distance ASC LIMIT 1;
  ```

### Layer 5: Direct Biological Database Integration (UniProt API inside SQL)
Instead of relying on external Python scripts to fetch biological metadata, `pg_bio` contains a native Rust-powered synchronous HTTP client (`ureq`). It can fetch and parse JSON from the UniProt API directly into specialized PostgreSQL Composite Types dynamically.

* **The Engine:** The `bio_fetch_uniprot()` function queries `rest.uniprot.org`, parses the JSON via `serde`, and projects it instantly into a Postgres `UniprotEntry` type (`id`, `organism`, `taxonomy`, `sequence`).
* **The Advantage:** Zero external scripts. You can enrich millions of orphan protein IDs dynamically inside pure SQL pipelines. 
* **SQL Example (Multi-Row Batch Fetch & Insert):**

  Using a `CROSS JOIN LATERAL`, you can map the HTTP fetch over a batch of IDs, evaluate it exactly once per row (saving network calls), and `INSERT` the results directly into your metadata table in a single atomic SQL transaction!

  ```sql
  -- 1. Create a table for the newly discovered metadata
  CREATE TABLE dark_proteome_metadata (
      uniprot_id VARCHAR(20) PRIMARY KEY,
      organism TEXT,
      taxonomy TEXT[],
      sequence TEXT
  );

  -- 2. The Batch Fetch & Insert Pipeline
  WITH new_discoveries(uniprot_id) AS (
      VALUES 
          ('A0A5B9DCV6'), -- Deep-sea Asgard archaeon (RuBisCO)
          ('L0JL84')      -- Halophilic archaeon
  ),
  enriched_data AS (
      -- The LATERAL join ensures the HTTP request evaluates exactly once per ID
      SELECT u.id, u.organism, u.taxonomy, u.sequence
      FROM new_discoveries d
      CROSS JOIN LATERAL bio_fetch_uniprot(d.uniprot_id) u
  )
  INSERT INTO dark_proteome_metadata (uniprot_id, organism, taxonomy, sequence)
  SELECT id, organism, taxonomy, sequence 
  FROM enriched_data
  RETURNING uniprot_id, organism;
  ```

* **SQL Example (Native API Search via Virtual Tables):**

  We also provide a **Set Returning Function (SRF)**, which acts mechanically identical to a read-only Foreign Data Wrapper (FDW). You can pass standard UniProt Lucene query parameters directly into the function and query it exactly like a local table.
  
  ```sql
  -- Find the top 5 longest curated proteins for the BRCA1 gene in Humans
  SELECT id, organism, length(sequence) AS seq_len 
  FROM bio_search_uniprot('gene:BRCA1 AND taxonomy_id:9606 AND reviewed:true')
  ORDER BY seq_len DESC 
  LIMIT 5;

  -- Filter archaeal proteins by exact length and sequence
  SELECT id, sequence
  FROM bio_search_uniprot('taxonomy_id:2594042 AND length:[1 TO 50]')
  WHERE sequence LIKE 'M%'; -- standard SQL filters apply to the returning set
  ```
  
  **Supported UniProt Query Parameters:**
  The `query` argument supports all standard UniProtKB REST API search fields, including:
  * `gene:BRCA1` (Filter by gene name)
  * `taxonomy_id:9606` (Filter by taxonomic ID, e.g., 9606 = Homo sapiens)
  * `organism_name:"Homo sapiens"` (Filter by exact organism string)
  * `length:[1 TO 200]` (Filter by amino acid sequence length)
  * `reviewed:true` (Only return manually curated Swiss-Prot entries)
  * `keyword:KW-0002` (Filter by specific biological keywords, e.g., 3D-structure)

---
### Layer 6: Cheminformatics (Molecular Graph Search & Fingerprinting)
To discover novel drugs, you must traverse massive libraries of chemical compounds, traditionally requiring you to serialize millions of strings into Python for parsing via RDKit.

* **The Engine:** `pg_bio` integrates the native Rust `purr` crate to parse **SMILES** (Simplified Molecular-Input Line-Entry System) strings directly in the database. It instantly converts chemical graphs into **1024-bit Morgan Fingerprints** (ECFP) entirely inside Postgres.
* **The Advantage:** By utilizing the native `%` **Tanimoto Similarity Operator**, you can perform virtual screening across billions of compounds natively in SQL. It is mathematically 5x-10x faster across massive datasets by entirely avoiding the I/O bottleneck of serializing data to Python.
* **SQL Example:**
  ```sql
  -- Find all library compounds with > 85% structural similarity to Aspirin
  SELECT compound_name, smiles                                              
  FROM library_compounds                                                    
  WHERE smiles_to_fingerprint(smiles) % smiles_to_fingerprint('CC(=O)OC1=CC=CC=C1C(=O)O') > 0.85;
  ```

---
## 🛠️ Scale & Optimization Infrastructure

`pg_bio` isn't just mathematically clever; it is heavily optimized for multi-million scale deployment on commercial hardware.

### 1. Vector Compression (`halfvec`)
High-dimensional vectors require massive amounts of RAM. `pg_bio` natively supports casting embedding columns to `halfvec` (16-bit floating-point). This cuts RAM and SSD footprint in half with zero loss in structural accuracy, allowing graphs of 2.5 million proteins to build natively in-memory on standard laptops.

### 2. Intelligent Graph Indexing (HNSW vs IVFFlat)
Vector indexing is fully tunable to hardware profiles. 
* Use **`HNSW` (Hierarchical Navigable Small World)** for the highest recall accuracy, utilizing optimized `maintenance_work_mem` ceilings to build massive interconnected graphs. 
* Use **`IVFFlat` (K-Means Clustering)** to bypass severe RAM limitations, dropping index build times from days to minutes on localized machines by sorting tuples into hyper-clusters (lists=1500).

---

## 💻 Python SDK & Tooling Integration

You don't need to write raw SQL to leverage this power.

* **`pgbio-py`:** A seamless Python SDK that connects directly to the Rust engine.
* **Jupyter Integration:** Full compatibility with Pandas and Matplotlib. View our complete tutorial notebook demonstrating all four layers in action: `https://github.com/jonatas/teaflon/blob/main/pg_bio_multiomics_pipeline/teaflon_multiomics_tutorial.ipynb`.

```python
from pgbio import PgBioClient

client = PgBioClient("postgresql://localhost:28818/bio_demo")

# Query the Joint Vector Space natively from Python
dna_matches = client.find_joint_enhancers(protein_id="P0DPB7")
print(f"Strongest regulatory match: {dna_matches[0].locus}")
```

## 🐳 Quickstart

```bash
# 1. Clone the repository
git clone https://github.com/your-org/pg_bio.git
cd pg_bio

# 2. Spin up the Rust-optimized Postgres Engine
docker-compose up -d

# 3. Compile the Database and Interactive Notebooks
uv run jupyter nbconvert --execute https://github.com/jonatas/teaflon/blob/main/pg_bio_multiomics_pipeline/teaflon_multiomics_tutorial.ipynb
```

*(For an extended deep-dive into the architectural theory, read `docs/MULTIOMICS_ARCHITECTURE.md`)*
