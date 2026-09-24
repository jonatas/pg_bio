# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "pgbio @ file://./pgbio-py",
# ]
# ///
from pgbio import PgBioClient

# Connect to the local pg_bio instance
client = PgBioClient("postgresql://localhost:28818/bio_demo")

print("--- DISCOVERY 1: Vector Homology Mock ---")
# 1. Structural/Functional Homology via Vector Embeddings
homologues = client.find_homologues("MFEGFERRLVD", limit=3)
for h in homologues:
    print(f"Match: {h.name} (Distance: {h.embedding_distance})")
print("(Wait, why are all distances 0? Let's check the Rust code...)")

print("\n--- DISCOVERY 2: Instant Spatial Search ---")
# We know 12AD exists and has tons of atoms. Let's find one of its atoms.
# Since pgbio client only has find_atoms_in_radius, let's just pick a coordinate
# Let's get an atom from 12AD to use as the center.
with client._get_conn() as conn:
    with conn.cursor() as cur:
        cur.execute("SELECT coord::text FROM protein_atoms WHERE uniprot_id = '12AD' LIMIT 1")
        coord_raw = cur.fetchone()[0]
        print(f"Target Atom coord raw: {coord_raw}")
        import json
        try:
            coord = json.loads(coord_raw)
        except:
            # Let's assume it's JSON-like if it parses, otherwise we might have to clean it
            print("Failed to parse coord JSON. Let's fallback to coordinate (0,0,0)")
            coord = {"x": 0.0, "y": 0.0, "z": 0.0}

print(f"Searching for atoms within 5.0 A of {coord['x']}, {coord['y']}, {coord['z']}...")

import time
start_time = time.time()
atoms = client.find_atoms_in_radius(target_x=coord['x'], target_y=coord['y'], target_z=coord['z'], radius=5.0)
duration = time.time() - start_time

print(f"Found {len(atoms)} atoms in {duration*1000:.2f} ms using Z-Order Indexing.")
for a in atoms[:5]:
    print(f"  Atom ID: {a.atom_id} ({a.name}) from {a.uniprot_id} at ({a.x}, {a.y}, {a.z})")

