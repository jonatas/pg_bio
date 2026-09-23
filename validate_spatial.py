import psycopg
import requests
import math
import json

DB_URI = "postgresql://localhost:28818/bio_demo"
PDB_ID = "10EC"

TGT_X = 170.835
TGT_Y = 152.248
TGT_Z = 207.561
RADIUS = 5.0

MIN_X = TGT_X - RADIUS
MAX_X = TGT_X + RADIUS
MIN_Y = TGT_Y - RADIUS
MAX_Y = TGT_Y + RADIUS
MIN_Z = TGT_Z - RADIUS
MAX_Z = TGT_Z + RADIUS

print(f"Validation for PDB {PDB_ID} within {RADIUS}A of ({TGT_X}, {TGT_Y}, {TGT_Z})")

# 1. Database Query using pg_bio spatial index
print("\n--- 1. Querying PostgreSQL with pg_bio Z-Order Index ---")
try:
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            query = f"""
            WITH bounds AS (
                SELECT 
                    residue_z_index('{{"x": {MIN_X}, "y": {MIN_Y}, "z": {MIN_Z}, "name": ""}}') as min_z,
                    residue_z_index('{{"x": {MAX_X}, "y": {MAX_Y}, "z": {MAX_Z}, "name": ""}}') as max_z
            )
            SELECT coord FROM protein_atoms, bounds 
            WHERE uniprot_id = '{PDB_ID}'
            AND z_index BETWEEN min_z AND max_z
            AND distance_angstroms(coord, '{{"x": {TGT_X}, "y": {TGT_Y}, "z": {TGT_Z}, "name": ""}}') <= {RADIUS};
            """
            cur.execute(query)
            db_results = cur.fetchall()
            db_atoms = [json.loads(row[0]) for row in db_results]
            print(f"Database found {len(db_atoms)} atoms.")
except Exception as e:
    print(f"DB Error: {e}")

# 2. Fetch from RCSB API and calculate manually
print("\n--- 2. Fetching from RCSB PDB API and calculating manually ---")
resp = requests.get(f"https://files.rcsb.org/download/{PDB_ID}.pdb")
api_atoms = []

if resp.status_code == 200:
    for line in resp.text.splitlines():
        if line.startswith("ATOM  "):
            try:
                atom_name = line[12:16].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                
                # Euclidean distance
                dist = math.sqrt((x - TGT_X)**2 + (y - TGT_Y)**2 + (z - TGT_Z)**2)
                if dist <= RADIUS:
                    api_atoms.append({"name": atom_name, "x": x, "y": y, "z": z, "dist": dist})
            except ValueError:
                continue
    print(f"API + Manual Math found {len(api_atoms)} atoms.")
else:
    print("Failed to fetch PDB file")

# 3. Compare Results
print("\n--- 3. Comparison ---")
db_coords = set((a['x'], a['y'], a['z']) for a in db_atoms)
api_coords = set((a['x'], a['y'], a['z']) for a in api_atoms)

if db_coords == api_coords:
    print("✅ PERFECT MATCH! The pg_bio spatial index returns the exact same atoms as manual API calculation.")
else:
    print("❌ MISMATCH!")
    print(f"In DB but not in API: {db_coords - api_coords}")
    print(f"In API but not in DB: {api_coords - db_coords}")

print("\nDB Atoms:")
for a in db_atoms:
    print(a)
    
print("\nAPI Atoms:")
for a in api_atoms:
    print(a)
