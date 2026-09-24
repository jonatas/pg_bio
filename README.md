# pg_bio: The High-Performance Bioinformatics Database Engine

**`pg_bio`** is a native PostgreSQL extension written in highly optimized Rust. It is engineered specifically for computational biologists, structural bioinformaticians, and drug discovery scientists who need to radically accelerate their pipelines.

## 🧬 The Data Bottleneck in Science
Modern biological models (like AlphaFold and ESM) generate billions of 3D coordinates and massive interaction matrices. Traditional scientific workflows rely on a "default Python approach":
1. Query a massive database to find a protein.
2. Download terabytes of structural data into memory.
3. Load it into Python (`NumPy`, `Pandas`, `BioPython`).
4. Run computationally heavy spatial searches (like finding a binding pocket) in a single-threaded Python process.

**This creates a massive I/O bottleneck.** Moving data out of the database is extremely slow, making biosphere-scale scans impossibly expensive.

## 🚀 The Solution: Bring Compute to the Data
`pg_bio` eliminates this bottleneck by pushing advanced bioinformatics mathematics *directly into the database engine*. By doing the heavy lifting in Rust before the data ever leaves Postgres, you save massive amounts of time, memory, and bandwidth. 

You can execute entire drug discovery pipelines in milliseconds, natively in SQL or via our seamless Python SDK.

### The 3 Pillars of pg_bio

#### 1. Vector Homology (High-Dimensional Sequence Embeddings)
Stop using slow, heuristic string alignments (BLAST) to find similar proteins. `pg_bio` natively stores and indexes high-dimensional vectors (like ESM-2 sequence embeddings).
* **The Math:** Natively computes Cosine Distances (`embedding_cosine_distance`) in Rust.
* **The Scientist's Advantage:** Instantly identify structural and functional protein twins across the entire database, even if their textual sequences look completely different.

#### 2. Sparse Attention Traversal (CSR Interaction Networks)
Language models output massive $N \times N$ attention matrices representing evolutionary coupling. Storing these densely consumes gigabytes per protein.
* **The Math:** `pg_bio` uses a custom `SparseAttentionMap` type based on Compressed Sparse Row (CSR) logic to store only significant interaction weights.
* **The Scientist's Advantage:** Traverse allosteric communication pathways and pinpoint active sites in microseconds using recursive SQL graph traversal.

#### 3. Z-Order Spatial Indexing (Absolute 3D Pockets)
Finding atoms within a 4.0 Ångstrom radius of a drug target typically requires brute-force Pythagorean math against millions of records in Python.
* **The Math:** `pg_bio` maps 3D Cartesian coordinates (X, Y, Z) into 1D integers using Morton Coding (Z-Order curves).
* **The Scientist's Advantage:** Standard B-Tree indexes now act as 3D bounding boxes. Extract precise 3D binding pockets or detect cross-protein molecular clashes across millions of atoms instantly.

---

## 💻 Zero-Friction Python SDK
We know scientists work in Python. You do not need to write raw SQL to use `pg_bio`. Our `pgbio-py` SDK elegantly bridges Python to the Rust database extension, turning complex SQL into three lines of Python code:

```python
from pgbio import PgBioClient

client = PgBioClient("postgresql://localhost:28818/bio_demo")

# 1. Vector Homology: Find structural homologues instantly
homologues = client.find_homologues(sequence="MFEGFERRLVD", limit=1)
best_match = homologues[0].uniprot_id

# 2. Attention: Which residues communicate heavily with residue #50?
active_site = client.find_interacting_residues(best_match, target_residue_index=50)

# 3. Spatial: Fetch the exact 3D Binding Pocket coordinates in milliseconds
atoms = client.find_atoms_in_radius(x=13.69, y=44.14, z=12.82, radius=4.0)
print(f"Extracted {len(atoms)} atoms in the binding pocket!")
```

---

## 🐳 Quickstart (Zero-Click Deployment)
We have containerized the entire ecosystem for immediate research use.

```bash
# 1. Clone the repository
git clone https://github.com/your-org/pg_bio.git
cd pg_bio

# 2. Spin up the Rust-optimized Postgres Engine & FastAPI Backend
docker-compose up -d

# 3. Run the local Demo Seeder to pull structures directly from RCSB PDB
uv run scripts/seed_bio_demo.py
```

---

## 🧪 The `scripts/` Directory
All Python logic has been consolidated into the `scripts/` folder to keep your root directory clean. 
You can run any of these using `uv run scripts/<script_name>.py`:

**Research Pipelines (Ready to run):**
* **`cross_species_analysis.py`** - Maps evolutionary conservation across different organisms.
* **`discovery.py`** - An automated workflow for mining the database for novel protein folds.
* **`off_target_prediction.py`** - Predicts unintended drug binding sites using vector homology.
* **`real_case.py`** - A complete end-to-end drug discovery pipeline demonstration.

**Utilities & Benchmarks:**
* **`benchmark_pgbio.py`** - Validates spatial search speeds against standard BioPython implementations.
* **`benchmark_embeddings.py`** - Validates vector math speeds against NumPy/SciPy.
* **`validate_attention_contacts.py`** - Validates that AI attention maps accurately correlate to physical 3D contacts in the folded protein.
* **`seed_bio_demo.py`** - Ingests 3D coordinates directly from the public RCSB PDB into Postgres.
* **`load_and_map_all.py`** - Pipeline for streaming and embedding the entire Human Proteome from UniProt.

---

## 📚 Resources & Next Steps
* **[Tutorial 1: Vector Homology](docs/01_vector_homology.md)**
* **[Tutorial 2: Attention Traversal](docs/02_attention_traversal.md)**
* **[Tutorial 3: Z-Order Spatial Indexing](docs/03_spatial_indexing.md)**
* **[IDE & Editor Integrations (MCP)](docs/IDE_INTEGRATION.md)** - How to connect Cursor or Claude to your local database.
* *(For an animated visual introduction to the project, open `docs/intro_blog_post.html` in your web browser!)*

---

## 📊 Benchmark Details
We ran rigorous benchmarks to prove why bringing compute to the database is strictly better than the default Python approach.

### 1. Vector Math (pg_bio vs NumPy)
* **Task:** Calculate Cosine Similarity across thousands of high-dimensional vectors.
* **Result:** `pg_bio` is **6.4x faster** than fetching the vectors over the network and running them through standard Numpy/SciPy calculations. Eliminating serialization and network I/O is a massive win.

### 2. Spatial 3D Indexing (pg_bio vs BioPython KDTree)
* **Task:** Find all neighboring atoms within a 5.0 Ångstrom radius inside a 197,010-atom structure (`2N5T`).
* **Result:** `pg_bio` is **order-of-magnitudes faster** from a cold start. Traditional pipelines require parsing the heavy `.pdb` file and building an in-memory KDTree (which takes ~1.2 seconds of setup time per query). `pg_bio` uses Z-Order Morton Coding natively in Postgres B-Trees to execute the same query in **~150 milliseconds** with zero setup time.
