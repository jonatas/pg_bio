# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "psycopg>=3.1",
#     "requests",
# ]
# ///

import psycopg
import requests
import random
import time

# Using the pgrx default local port for Postgres 18
DB_URI_DEFAULT = "postgresql://localhost:28818/postgres"
DB_URI_DEMO = "postgresql://localhost:28818/bio_demo"

def create_database():
    print("Connecting to local pgrx cluster to create 'bio_demo' database...")
    try:
        with psycopg.connect(DB_URI_DEFAULT, autocommit=True) as conn:
            conn.execute("DROP DATABASE IF EXISTS bio_demo;")
            conn.execute("CREATE DATABASE bio_demo;")
        print("Database 'bio_demo' created successfully.")
    except Exception as e:
        print(f"Error creating database: {e}")

def setup_schema_and_seed():
    print("Connecting to 'bio_demo'...")
    with psycopg.connect(DB_URI_DEMO, autocommit=True) as conn:
        with conn.cursor() as cur:
            print("Installing our custom pg_bio extension...")
            cur.execute("CREATE EXTENSION IF NOT EXISTS pg_bio;")
            
            # 1. Create Tables
            cur.execute("""
                CREATE TABLE proteins (
                    uniprot_id VARCHAR(20) PRIMARY KEY,
                    name TEXT,
                    sequence TEXT,
                    embedding REAL[]
                );
            """)
            
            cur.execute("""
                CREATE TABLE protein_atoms (
                    atom_id SERIAL PRIMARY KEY,
                    uniprot_id VARCHAR(20) REFERENCES proteins(uniprot_id),
                    coord ResidueCoord,
                    z_index BIGINT GENERATED ALWAYS AS (residue_z_index(coord)) STORED
                );
            """)

            # 2. Download Real Data from UniProt (Swiss-Prot Human Proteome)
            print("Streaming real Human Proteome data directly from UniProt API...")
            url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:9606)"
            response = requests.get(url, stream=True)
            
            proteins = []
            current_id, current_name, current_seq = None, None, []
            
            # Parse FASTA on the fly, limiting to 2,000 proteins for a quick demo
            for line in response.iter_lines(decode_unicode=True):
                if not line: continue
                if line.startswith(">"):
                    if current_id and len(proteins) < 2000:
                        proteins.append((current_id, current_name, "".join(current_seq)))
                    if len(proteins) >= 2000:
                        break
                    
                    parts = line.split("|")
                    current_id = parts[1] if len(parts) > 1 else "UNKNOWN"
                    current_name = parts[2].split(" OS=")[0] if len(parts) > 2 else "Unknown Protein"
                    current_seq = []
                else:
                    current_seq.append(line.strip())

            # 3. Insert Proteins and Generate Embeddings natively in Postgres
            print(f"Inserting {len(proteins)} human proteins into the database...")
            for pid, name, seq in proteins:
                # We use our custom get_esm_embedding function built in Rust
                cur.execute(
                    "INSERT INTO proteins (uniprot_id, name, sequence, embedding) VALUES (%s, %s, %s, get_esm_embedding(%s))",
                    (pid, name, seq, seq)
                )

            # 4. Generate Spatial 3D atomic coordinates
            print("Generating 200,000 synthetic 3D atomic coordinates attached to the real proteins...")
            records = []
            for pid, _, _ in proteins:
                # 100 atoms per protein
                for _ in range(100):
                    x, y, z = random.uniform(-100, 100), random.uniform(-100, 100), random.uniform(-100, 100)
                    coord_str = f'{{"x":{x},"y":{y},"z":{z},"name":"ALA"}}'
                    records.append((pid, coord_str))
            
            print("Bulk loading spatial data via COPY STDIN...")
            with cur.copy("COPY protein_atoms (uniprot_id, coord) FROM STDIN") as copy:
                for r in records:
                    copy.write_row(r)
            
            # 5. Optimize Indexes
            print("Building massive Z-Order B-Tree spatial index...")
            cur.execute("CREATE INDEX idx_spatial_z_order ON protein_atoms (z_index);")
            
            print("\n🎉 Database 'bio_demo' successfully seeded and heavily optimized!")
            print("You can now connect to it via `psql -p 28818 -d bio_demo` to run spatial queries.")

if __name__ == "__main__":
    create_database()
    setup_schema_and_seed()
