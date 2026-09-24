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
import argparse
import sys
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Configure requests session with retries to handle connection resets
session = requests.Session()
retry = Retry(connect=5, read=5, backoff_factor=0.5, status_forcelist=[ 500, 502, 503, 504 ])
adapter = HTTPAdapter(max_retries=retry)
session.mount('http://', adapter)
session.mount('https://', adapter)

# Using the pgrx default local port for Postgres 18
DB_URI_DEFAULT = "postgresql://localhost:28818/postgres"
DB_URI_DEMO = "postgresql://localhost:28818/bio_demo"

ORGANISM_MAP = {
    "human": "9606",
    "mouse": "10090",
    "ecoli": "83333",
    "yeast": "4932",
    "arabidopsis": "3702"
}

def interactive_prompt():
    print("🌟 Welcome to the pg_bio Database Seeder 🌟")
    print("This script will help you set up your local bio_demo database with data from public databases.\n")
    
    # Source
    print("Select a data source:")
    print("1) UniProt (Proteins with sequence only - synthetic 3D coords)")
    print("2) PDB (Proteins with real 3D atomic coordinates)")
    source_choice = input("Enter 1 or 2 [1]: ").strip()
    source = "pdb" if source_choice == "2" else "uniprot"
    
    # Organism
    print("\nSelect an organism:")
    print("1) human (9606)")
    print("2) mouse (10090)")
    print("3) ecoli (83333)")
    print("4) custom taxonomy ID")
    org_choice = input("Enter 1-4 [1]: ").strip()
    
    organism = "human"
    if org_choice == "2": organism = "mouse"
    elif org_choice == "3": organism = "ecoli"
    elif org_choice == "4": 
        organism = input("Enter custom taxonomy ID: ").strip()

    # Size
    print("\nHow many proteins do you want to download?")
    print("Enter a number (e.g., 10, 100, 1000) or -1 for ALL available.")
    size_choice = input("Limit [100]: ").strip()
    try:
        limit = int(size_choice) if size_choice else 100
    except ValueError:
        limit = 100

    # Reset
    print("\nDo you want to reset (drop) the existing database first? (y/N)")
    reset_choice = input("Reset [N]: ").strip().lower()
    reset = reset_choice == "y" or reset_choice == "yes"

    return source, organism, limit, reset

def create_database(reset=False):
    print("Connecting to local pgrx cluster...")
    try:
        with psycopg.connect(DB_URI_DEFAULT, autocommit=True) as conn:
            if reset:
                print("Resetting 'bio_demo' database...")
                conn.execute("DROP DATABASE IF EXISTS bio_demo;")
            
            # Create if not exists
            try:
                conn.execute("CREATE DATABASE bio_demo;")
                print("Database 'bio_demo' created successfully.")
            except psycopg.errors.DuplicateDatabase:
                print("Database 'bio_demo' already exists.")
    except Exception as e:
        print(f"Error creating database: {e}")
        sys.exit(1)

def get_pdb_ids_for_organism(organism_name, limit):
    print(f"Querying RCSB PDB for {organism_name} structures...")
    # Basic graph query to get PDB IDs for the given organism
    query = {
      "query": {
        "type": "terminal",
        "service": "text",
        "parameters": {
          "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
          "operator": "exact_match",
          "value": organism_name.capitalize() if not organism_name.isdigit() else "Homo sapiens" # Simplification for demo
        }
      },
      "return_type": "entry",
      "request_options": {
        "paginate": {
          "start": 0,
          "rows": limit if limit > 0 else 500
        }
      }
    }
    
    # If organism is a tax ID, use tax_id search
    if organism_name.isdigit():
         query["query"]["parameters"]["attribute"] = "rcsb_entity_source_organism.taxonomy_lineage.id"
         query["query"]["parameters"]["value"] = str(organism_name)
         
    resp = session.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query)
    if resp.status_code != 200:
        print("Failed to query RCSB PDB.")
        return []
        
    data = resp.json()
    if "result_set" not in data:
        return []
        
    return [item["identifier"] for item in data["result_set"]]

def fetch_pdb_data(pdb_id):
    # 1. Fetch FASTA
    fasta_resp = session.get(f"https://www.rcsb.org/fasta/entry/{pdb_id}")
    if fasta_resp.status_code != 200:
        return None, None
        
    lines = fasta_resp.text.splitlines()
    if not lines: return None, None
    
    sequence = ""
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip()
            
    # 2. Fetch PDB coordinates
    pdb_resp = session.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
    if pdb_resp.status_code != 200:
        return sequence, []
        
    atoms = []
    for line in pdb_resp.text.splitlines():
        if line.startswith("ATOM  "):
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append({"name": atom_name, "x": x, "y": y, "z": z})
            except ValueError:
                continue
                
    return sequence, atoms

def setup_schema_and_seed(source, organism_id, organism_raw, limit):
    print(f"\nConnecting to 'bio_demo' to seed {limit if limit > 0 else 'all'} proteins from {source.upper()}...")
    with psycopg.connect(DB_URI_DEMO, autocommit=True) as conn:
        with conn.cursor() as cur:
            print("Installing our custom pg_bio extension...")
            cur.execute("CREATE EXTENSION IF NOT EXISTS pg_bio;")
            
            # 1. Create Tables
            cur.execute("""
                CREATE TABLE IF NOT EXISTS proteins (
                    uniprot_id VARCHAR(20) PRIMARY KEY,
                    name TEXT,
                    sequence TEXT,
                    embedding REAL[]
                );
            """)
            
            cur.execute("""
                CREATE TABLE IF NOT EXISTS protein_atoms (
                    atom_id SERIAL PRIMARY KEY,
                    uniprot_id VARCHAR(20) REFERENCES proteins(uniprot_id),
                    coord ResidueCoord,
                    z_index BIGINT GENERATED ALWAYS AS (residue_z_index(coord)) STORED
                );
            """)

            # 2. Download Data
            inserted_count = 0
            
            if source == "uniprot":
                url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(reviewed:true)+AND+(model_organism:{organism_id})"
                print(f"Streaming data from UniProt API: {url}")
                response = session.get(url, stream=True)
                
                proteins = []
                current_id, current_name, current_seq = None, None, []
                
                for line in response.iter_lines(decode_unicode=True):
                    if not line: continue
                    if line.startswith(">"):
                        if current_id and (limit <= 0 or len(proteins) < limit):
                            proteins.append((current_id, current_name, "".join(current_seq)))
                        
                        if limit > 0 and len(proteins) >= limit:
                            break
                        
                        parts = line.split("|")
                        current_id = parts[1] if len(parts) > 1 else "UNKNOWN"
                        current_name = parts[2].split(" OS=")[0] if len(parts) > 2 else "Unknown Protein"
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                
                if current_id and (limit <= 0 or len(proteins) < limit):
                    proteins.append((current_id, current_name, "".join(current_seq)))
                
                print(f"Inserting {len(proteins)} proteins into the database...")
                new_proteins = []
                for pid, name, seq in proteins:
                    try:
                        cur.execute(
                            """
                            INSERT INTO proteins (uniprot_id, name, sequence, embedding) 
                            VALUES (%s, %s, %s, get_esm_embedding(%s))
                            ON CONFLICT (uniprot_id) DO NOTHING
                            RETURNING uniprot_id;
                            """,
                            (pid, name, seq, seq)
                        )
                        if cur.fetchone():
                            inserted_count += 1
                            new_proteins.append((pid, name, seq))
                    except Exception as e:
                        print(f"Error inserting {pid}: {e}")
                
                print(f"Successfully inserted {inserted_count} new proteins.")

                if new_proteins:
                    print(f"Generating synthetic 3D atomic coordinates for {len(new_proteins)} new proteins...")
                    records = []
                    for pid, _, _ in new_proteins:
                        for _ in range(100):
                            x, y, z = random.uniform(-100, 100), random.uniform(-100, 100), random.uniform(-100, 100)
                            coord_str = f'{{"x":{x},"y":{y},"z":{z},"name":"ALA"}}'
                            records.append((pid, coord_str))
                    
                    print("Bulk loading spatial data via COPY STDIN...")
                    with cur.copy("COPY protein_atoms (uniprot_id, coord) FROM STDIN") as copy:
                        for r in records:
                            copy.write_row(r)
                            
            elif source == "pdb":
                pdb_ids = get_pdb_ids_for_organism(organism_id, limit)
                print(f"Found {len(pdb_ids)} PDB entries. Processing...")
                
                for pdb_id in pdb_ids:
                    # Check if it already exists in the database
                    cur.execute("SELECT 1 FROM proteins WHERE uniprot_id = %s", (pdb_id,))
                    if cur.fetchone():
                        print(f"Skipping {pdb_id} (Already exists in DB).")
                        continue

                    print(f"Fetching {pdb_id}...", end=" ", flush=True)
                    seq, atoms = fetch_pdb_data(pdb_id)
                    
                    if not seq:
                        print("Failed (No sequence found).")
                        continue
                        
                    # Insert protein
                    cur.execute(
                        """
                        INSERT INTO proteins (uniprot_id, name, sequence, embedding) 
                        VALUES (%s, %s, %s, get_esm_embedding(%s))
                        ON CONFLICT (uniprot_id) DO NOTHING
                        RETURNING uniprot_id;
                        """,
                        (pdb_id, f"PDB Structure {pdb_id}", seq, seq)
                    )
                    
                    if not cur.fetchone():
                        print("Skipped (Already exists).")
                        continue
                        
                    inserted_count += 1
                    
                    # Insert atoms
                    if atoms:
                        print(f"Success ({len(atoms)} real atoms). Bulk inserting...", end="")
                        records = []
                        for a in atoms:
                            coord_str = f'{{"x":{a["x"]},"y":{a["y"]},"z":{a["z"]},"name":"{a["name"]}"}}'
                            records.append((pdb_id, coord_str))
                            
                        with cur.copy("COPY protein_atoms (uniprot_id, coord) FROM STDIN") as copy:
                            for r in records:
                                copy.write_row(r)
                        print(" Done.")
                    else:
                        print(" Done (No atoms found).")
                    
            print(f"\nCompleted! Inserted {inserted_count} new entries.")
            
            print("Ensuring spatial Z-Order index exists...")
            cur.execute("CREATE INDEX IF NOT EXISTS idx_spatial_z_order ON protein_atoms (z_index);")
                
            print("\n🎉 Database 'bio_demo' successfully seeded and highly optimized!")
            print("You can now connect to it via `psql -p 28818 -d bio_demo` to run spatial queries.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Running in CLI mode
        parser = argparse.ArgumentParser(description="Seed the bio_demo PostgreSQL database with protein data.")
        parser.add_argument("--source", type=str, default="uniprot", choices=["uniprot", "pdb"], help="Public database source to download from")
        parser.add_argument("--organism", type=str, default="human", help="Organism to filter by (e.g. human, mouse, ecoli). Can also be a taxonomy ID.")
        parser.add_argument("--limit", type=int, default=100, help="Number of proteins to download. Set to -1 for all.")
        parser.add_argument("--reset", action="store_true", help="Drop and recreate the database from scratch.")
        
        args = parser.parse_args()
        
        source = args.source
        organism_raw = args.organism
        limit = args.limit
        reset = args.reset
    else:
        # Interactive mode
        source, organism_raw, limit, reset = interactive_prompt()
        
    organism_id = ORGANISM_MAP.get(organism_raw.lower(), organism_raw)
    
    create_database(reset=reset)
    setup_schema_and_seed(source=source, organism_id=organism_id, organism_raw=organism_raw, limit=limit)
