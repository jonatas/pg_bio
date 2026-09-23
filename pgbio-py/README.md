# pgbio-py

The official Python SDK for the `pg_bio` PostgreSQL extension.

`pgbio-py` allows computational biologists and engineers to seamlessly interact with `pg_bio` directly from Python. It abstracts away the complex Z-Order spatial SQL and vector distance queries, letting you focus on the science.

## Installation

```bash
pip install -e .
```

## Quickstart

```python
from pgbio import PgBioClient

# Connect to the local pg_bio instance
client = PgBioClient("postgresql://localhost:28818/bio_demo")

# 1. Structural/Functional Homology via Vector Embeddings
print("Searching for homologues...")
homologues = client.find_homologues("MFEGFERRLVD", limit=3)
for h in homologues:
    print(f"Match: {h.name} (Distance: {h.embedding_distance})")

# 2. 3D Spatial Radius Search
print("\nFinding atoms within 5.0 Ångstroms of coordinate (0, 0, 0)...")
atoms = client.find_atoms_in_radius(target_x=0.0, target_y=0.0, target_z=0.0, radius=5.0)
print(f"Found {len(atoms)} atoms instantly using Z-Order Indexing.")
```
