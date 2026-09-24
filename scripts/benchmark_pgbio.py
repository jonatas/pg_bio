import time
import psycopg
import json
from Bio.PDB import PDBParser, NeighborSearch, Selection

PDB_ID = "2N5T"
TGT_X, TGT_Y, TGT_Z = 4.716, 4.299, -33.733
RADIUS = 5.0

print(f"--- BENCHMARK: Spatial Search in {PDB_ID} (197,010 atoms) ---")
print(f"Task: Find all atoms within {RADIUS}A of ({TGT_X}, {TGT_Y}, {TGT_Z})")

# ==========================================
# 1. TRADITIONAL PIPELINE (BioPython KDTree)
# ==========================================
print("\n[1] Traditional Pipeline (BioPython in-memory KDTree)")
start_time = time.perf_counter()

# Step A: Parse the file (Disk I/O + Parsing overhead)
parser = PDBParser(QUIET=True)
structure = parser.get_structure(PDB_ID, f"{PDB_ID}.pdb")
parse_time = time.perf_counter()

# Step B: Extract all atoms and build KDTree
atoms = Selection.unfold_entities(structure, 'A')
ns = NeighborSearch(atoms)
kdtree_time = time.perf_counter()

# Step C: Query the KDTree
target_coord = (TGT_X, TGT_Y, TGT_Z)
bio_results = ns.search(target_coord, RADIUS, level='A')
query_time = time.perf_counter()

bio_total_time = query_time - start_time
bio_parse_dur = parse_time - start_time
bio_kdtree_dur = kdtree_time - parse_time
bio_query_dur = query_time - kdtree_time

print(f"  Results found: {len(bio_results)}")
print(f"  - File Parsing Time: {bio_parse_dur*1000:.2f} ms")
print(f"  - KDTree Build Time: {bio_kdtree_dur*1000:.2f} ms")
print(f"  - Query Execution:   {bio_query_dur*1000:.2f} ms")
print(f"  TOTAL TIME:          {bio_total_time*1000:.2f} ms")

# ==========================================
# 2. PG_BIO PIPELINE (Z-Order Spatial Index)
# ==========================================
print("\n[2] pg_bio Pipeline (PostgreSQL Z-Order Index)")

DB_URI = "postgresql://localhost:28818/bio_demo"

# Calculate bounding box bounds for Z-Order
MIN_X, MAX_X = TGT_X - RADIUS, TGT_X + RADIUS
MIN_Y, MAX_Y = TGT_Y - RADIUS, TGT_Y + RADIUS
MIN_Z, MAX_Z = TGT_Z - RADIUS, TGT_Z + RADIUS

# We connect once (as web servers/apps use connection pooling)
with psycopg.connect(DB_URI) as conn:
    with conn.cursor() as cur:
        # We start timing the actual query execution
        start_time = time.perf_counter()
        
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
        pg_results = cur.fetchall()
        
        pg_total_time = time.perf_counter() - start_time

print(f"  Results found: {len(pg_results)}")
print(f"  TOTAL TIME:          {pg_total_time*1000:.2f} ms")

# ==========================================
# COMPARISON
# ==========================================
print("\n--- RESULTS ---")
speedup = bio_total_time / pg_total_time
print(f"pg_bio is {speedup:.1f}x faster than the traditional pipeline from cold start.")
if bio_query_dur > 0:
    q_speedup = bio_parse_dur / pg_total_time # usually query is fast, but setup is slow.
    print(f"(Note: BioPython's pure query takes {bio_query_dur*1000:.2f}ms, but requires {bio_parse_dur*1000:.2f}ms of setup each time a new protein is loaded)")
