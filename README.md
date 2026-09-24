# pg_bio: The High-Performance Bioinformatics Database Engine

**`pg_bio`** is a native PostgreSQL extension written in highly optimized Rust. It is engineered specifically for computational biologists, structural bioinformaticians, and drug discovery scientists who need to accelerate their pipelines.

## 🧬 The Data Bottleneck
Modern biological pipelines (like AlphaFold, ESM, and Hi-C genomic folding) generate billions of 3D coordinates and interaction matrices. Traditional workflows require extracting terabytes of this data out of a database and into Python (NumPy/Pandas) just to run basic spatial searches or homology queries. This creates a massive I/O bottleneck that makes scanning the biosphere extremely slow and expensive.

## 🚀 The Solution: Bring Compute to the Data
`pg_bio` eliminates the network bottleneck by pushing advanced bioinformatics mathematics directly into the database engine. You can now execute entire drug discovery pipelines in milliseconds, natively in SQL or via our Python SDK.

### Core Scientific Capabilities

#### 1. Vector Homology (High-Dimensional Sequence Embeddings)
Replace slow string alignments (BLAST) with sequence semantics. `pg_bio` natively stores and indexes high-dimensional vectors (like K-mer feature hashes or ESM-2 embeddings).
* **The Math:** Natively computes Cosine Distances (`embedding_cosine_distance`) in Rust.
* **The Result:** Instantly identify structural and functional protein twins across the entire database, even if their textual sequences look completely different.

#### 2. Sparse Attention Traversal (CSR Interaction Networks)
Protein language models output massive $N \times N$ attention matrices representing evolutionary coupling. Storing these densely consumes gigabytes per protein.
* **The Math:** `pg_bio` uses a custom `SparseAttentionMap` type based on Compressed Sparse Row (CSR) logic. It natively filters and stores only significant interaction weights.
* **The Result:** Traverse allosteric communication pathways and pinpoint active sites in microseconds.

#### 3. Z-Order Spatial Indexing (Absolute 3D Pockets)
Finding atoms within a 4.0 Ångstrom radius of a drug target typically requires brute-force Pythagorean math ($d = \sqrt{x^2 + y^2 + z^2}$) against millions of records.
* **The Math:** `pg_bio` maps 3D Cartesian coordinates (X, Y, Z) into 1D integers using Morton Coding (Z-Order curves).
* **The Result:** Standard PostgreSQL B-Tree indexes can now be used for 3D bounding boxes. Extract precise 3D binding pockets or detect cross-protein molecular clashes across 5.6+ million atoms in under 200 milliseconds.

---

## 💻 Zero-Friction Python SDK
We know scientists work in Python. You do not need to write raw SQL to use `pg_bio`. The `pgbio-py` SDK elegantly bridges Python to the Rust database extension.

```python
from pgbio import PgBioClient

client = PgBioClient("postgresql://localhost:28818/bio_demo")

# 1. Find Structural Homologues instantly
homologues = client.find_homologues(sequence="MFEGFERRLVD", limit=1)
best_match = homologues[0].uniprot_id

# 2. Traverse the Attention Network to find the active site
# (Which residues communicate heavily with residue #50?)
active_site = client.find_interacting_residues(best_match, target_residue_index=50)

# 3. Fetch the 3D Binding Pocket instantly via Z-Order Indexing
atoms = client.find_atoms_in_radius(x=13.69, y=44.14, z=12.82, radius=4.0)
print(f"Extracted {len(atoms)} atoms in the binding pocket!")
```

## 🐳 Quickstart (Zero-Click Deployment)
We have containerized the entire ecosystem for immediate research use.

```bash
# 1. Clone the repository
git clone https://github.com/your-org/pg_bio.git
cd pg_bio

# 2. Spin up the Rust-optimized Postgres Engine & FastAPI Backend
docker-compose up -d

# 3. (Optional) Run the local Demo Seeder to pull structures directly from RCSB PDB
uv run scripts/seed_bio_demo.py
```

## 📚 Resources & Next Steps

Whether you are a researcher looking to run benchmarks or an engineer deploying the engine, everything you need is linked below.

### 📖 Tutorials & Manuals
* **[Tutorial 1: The Magic of Vector Homology](docs/01_vector_homology.md)** - Understanding biological sequence embeddings.
* **[Tutorial 2: Attention Traversal & Allostery](docs/02_attention_traversal.md)** - Traversing massive interaction networks instantly.
* **[Tutorial 3: Z-Order Spatial Indexing](docs/03_spatial_indexing.md)** - Finding 3D binding pockets and collisions.
* **[Python SDK Manual](pgbio-py/README.md)** - Full documentation for the `pgbio-py` Python client.
* **[FastAPI Backend Manual](backend/README.md)** - Setup instructions for the web API engine.
* *(For an animated visual introduction to the project, open `docs/intro_blog_post.html` in your web browser!)*

### 🔬 Research & Experiments
We have provided raw SQL notes and Python benchmarks proving `pg_bio`'s capabilities across different biological domains:
* **[Experiment 1: Z-Order vs B-Tree](experiments/experiment_1.sql)** - Proof of spatial indexing performance.
* **[Experiment 2: Sparse Attention](experiments/experiment_2.sql)** - Creating and querying compressed neural network weights.
* **[Experiment 3: Hi-C ORCA Genomics](experiments/experiment_3_orca_genomics.sql)** - Using `pg_bio` to query 3D genome architecture and chromatin contact maps.

### 🧪 Scripts & Pipelines
We have heavily organized our Python utilities to help you run local workflows:
* **`scripts/`**: Contains core database utilities, seeding scripts (like `seed_bio_demo.py`), embeddings migration tools, and rigorous benchmark scripts (`benchmark_pgbio.py` which proves in-database Cosine Similarity is **6.4x faster** than fetching vectors to Python).
* **`examples/`**: Contains ready-to-run research pipelines demonstrating Off-Target Prediction, Cross-Species Analysis, and Discovery Workflows natively on your machine!

### 🛠️ Setting Up the Database Locally
Choose the method that best fits your workflow:

1. **Zero-Click Docker (Recommended)**
   Spin up the pre-configured environment in seconds:
   ```bash
   docker-compose up -d
   ```

2. **Native Rust Compilation (For Developers)**
   If you want to compile the extension from source against your local PostgreSQL 18 instance:
   ```bash
   cargo pgrx install --release
   # Then in psql: CREATE EXTENSION pg_bio;
   ```

3. **Seeding Live Data (RCSB PDB)**
   Once your database is running, you can seed it with real 3D coordinates directly from the Protein Data Bank:
   ```bash
   uv run scripts/seed_bio_demo.py
   ```
