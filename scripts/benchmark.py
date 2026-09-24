# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "psycopg>=3.1",
#     "numpy",
# ]
# ///

import time
import random
import math
import psycopg
import numpy as np
import sys

# Connect to the pgrx-managed Postgres 18 instance (default port 28818)
DB_URI = "postgresql://localhost:28818/postgres"

NUM_ATOMS = 100_000
TARGET_X, TARGET_Y, TARGET_Z = 0.0, 0.0, 0.0
RADIUS = 5.0 # Ångstroms

def generate_atoms(count):
    print(f"Generating {count} synthetic atoms in a 200x200x200 Å box...")
    atoms = []
    for _ in range(count):
        # Random coords between -100 and 100
        x = random.uniform(-100, 100)
        y = random.uniform(-100, 100)
        z = random.uniform(-100, 100)
        atoms.append((x, y, z, "ALA"))
    return atoms

def python_benchmark(atoms):
    print("\n--- Python Brute-Force (Flat List) ---")
    start_time = time.time()
    
    matches = []
    for (x, y, z, name) in atoms:
        dx = x - TARGET_X
        dy = y - TARGET_Y
        dz = z - TARGET_Z
        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
        if dist <= RADIUS:
            matches.append((x, y, z, name, dist))
            
    end_time = time.time()
    elapsed = (end_time - start_time) * 1000 # ms
    print(f"Found {len(matches)} atoms within {RADIUS}Å")
    print(f"Time taken: {elapsed:.2f} ms")
    return elapsed

def numpy_benchmark(atoms):
    print("\n--- Python Numpy (Vectorized) ---")
    start_time = time.time()
    
    # Convert to numpy array
    coords = np.array([[x, y, z] for x, y, z, _ in atoms])
    target = np.array([TARGET_X, TARGET_Y, TARGET_Z])
    
    # Calculate distances
    dists = np.linalg.norm(coords - target, axis=1)
    matches = np.where(dists <= RADIUS)[0]
    
    end_time = time.time()
    elapsed = (end_time - start_time) * 1000 # ms
    print(f"Found {len(matches)} atoms within {RADIUS}Å")
    print(f"Time taken (including array setup): {elapsed:.2f} ms")
    return elapsed

def pg_bio_benchmark(atoms):
    print("\n--- PostgreSQL pg_bio (Z-Order B-Tree) ---")
    
    try:
        with psycopg.connect(DB_URI, autocommit=True) as conn:
            with conn.cursor() as cur:
                print("Setting up pg_bio database...")
                cur.execute("CREATE EXTENSION IF NOT EXISTS pg_bio;")
                cur.execute("DROP TABLE IF EXISTS atoms;")
                cur.execute("""
                    CREATE TABLE atoms (
                        id serial PRIMARY KEY,
                        coord ResidueCoord,
                        z_index bigint
                    );
                """)
                
                # Insert data
                print(f"Loading {len(atoms)} atoms into PostgreSQL...")
                with cur.copy("COPY atoms (coord, z_index) FROM STDIN") as copy:
                    for x, y, z, name in atoms:
                        # Serialize to our custom Postgres type format (JSON-ish string)
                        coord_str = f'{{"x":{x},"y":{y},"z":{z},"name":"{name}"}}'
                        
                        # We calculate z-index in postgres, but for COPY speed, we can compute it 
                        # or just insert the coord and update. Let's do an INSERT pipeline.
                        pass # Copy is tricky with custom types and function calls.
                
                # Better to use execute_batch
                print("Inserting data (this might take a few seconds)...")
                # We use execute_values or standard inserts. To be fast, we'll insert just the coords
                # and let Postgres calculate the z_index!
                cur.execute("ALTER TABLE atoms DROP COLUMN z_index")
                cur.execute("ALTER TABLE atoms ADD COLUMN z_index bigint GENERATED ALWAYS AS (residue_z_index(coord)) STORED;")
                
                # Fast batch insert
                records = [(f'{{"x":{x},"y":{y},"z":{z},"name":"{name}"}}',) for x, y, z, name in atoms]
                with cur.copy("COPY atoms (coord) FROM STDIN") as copy:
                    for record in records:
                        copy.write_row(record)
                
                print("Building B-Tree Index on Z-Order curve...")
                cur.execute("CREATE INDEX idx_atoms_z ON atoms (z_index);")
                
                # --- BENCHMARK ---
                print("Running spatial query...")
                start_time = time.time()
                
                # To search within RADIUS, we calculate the bounding box Z-Order bounds
                # Box: [TARGET_X - RADIUS, TARGET_Y - RADIUS, TARGET_Z - RADIUS] to [TARGET_X + RADIUS, ...]
                
                min_coord = f'{{"x": {TARGET_X - RADIUS}, "y": {TARGET_Y - RADIUS}, "z": {TARGET_Z - RADIUS}, "name": ""}}'
                max_coord = f'{{"x": {TARGET_X + RADIUS}, "y": {TARGET_Y + RADIUS}, "z": {TARGET_Z + RADIUS}, "name": ""}}'
                tgt_coord = f'{{"x": {TARGET_X}, "y": {TARGET_Y}, "z": {TARGET_Z}, "name": ""}}'
                
                cur.execute(f"""
                    WITH bounds AS (
                        SELECT 
                            residue_z_index('{min_coord}') as min_z,
                            residue_z_index('{max_coord}') as max_z
                    )
                    SELECT count(*) FROM atoms, bounds 
                    WHERE z_index BETWEEN min_z AND max_z
                    AND distance_angstroms(coord, '{tgt_coord}') <= {RADIUS};
                """)
                
                matches = cur.fetchone()[0]
                end_time = time.time()
                elapsed = (end_time - start_time) * 1000
                
                print(f"Found {matches} atoms within {RADIUS}Å")
                print(f"Time taken (End-to-End Query): {elapsed:.2f} ms")
                return elapsed
                
    except Exception as e:
        print(f"Database error: {e}")
        return None

if __name__ == "__main__":
    atoms = generate_atoms(NUM_ATOMS)
    
    t_py = python_benchmark(atoms)
    t_np = numpy_benchmark(atoms)
    t_pg = pg_bio_benchmark(atoms)
    
    if t_pg:
        print("\n=== SUMMARY ===")
        print(f"Python (Brute Force) : {t_py:.2f} ms")
        print(f"Python (Numpy Vector): {t_np:.2f} ms")
        print(f"Postgres (pg_bio)    : {t_pg:.2f} ms")
        print(f"\npg_bio speedup vs Python: {t_py/t_pg:.1f}x")
        print(f"pg_bio speedup vs Numpy: {t_np/t_pg:.1f}x")
