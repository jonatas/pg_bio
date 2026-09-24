# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "pgbio @ file://./pgbio-py",
# ]
# ///
import psycopg
import json
from pgbio import PgBioClient

DB_URI = "postgresql://localhost:28818/bio_demo"
client = PgBioClient(DB_URI)

with psycopg.connect(DB_URI) as conn:
    with conn.cursor() as cur:
        # Step 1: Get 1BAY's sequence and find homologues
        print("--- STEP 1: Vector Homology ---")
        cur.execute("SELECT sequence FROM proteins WHERE uniprot_id = '1BAY';")
        seq = cur.fetchone()[0]
        
        cur.execute("""
            SELECT uniprot_id, name, embedding_cosine_distance(embedding, get_esm_embedding(%s)) as dist
            FROM proteins 
            WHERE uniprot_id != '1BAY' 
            ORDER BY dist ASC 
            LIMIT 3;
        """, (seq,))
        homologues = cur.fetchall()
        for h in homologues:
            print(f"Homologue: {h[0]} ({h[1]}) - Distance: {h[2]:.4f}")
            
        # Step 2: Attention Traversal
        print("\n--- STEP 2: Attention Traversal ---")
        target_residue = 50
        cur.execute("""
            SELECT get_top_interacting_residues(attention_data, %s, 5)
            FROM protein_attention_maps
            WHERE uniprot_id = '1BAY';
        """, (target_residue,))
        interacting = cur.fetchone()[0]
        print(f"Top 5 interacting residues for residue {target_residue}: {interacting}")
        
        # Step 3: Spatial Indexing
        print("\n--- STEP 3: Spatial Indexing ---")
        # We need to find the coordinates of an atom near residue 50 (let's say the 50th atom)
        cur.execute("""
            SELECT atom_id, coord::text
            FROM protein_atoms 
            WHERE uniprot_id = '1BAY'
            ORDER BY atom_id ASC
            OFFSET %s LIMIT 1;
        """, (target_residue,))
        target_atom = cur.fetchone()
        atom_id = target_atom[0]
        coord_raw = target_atom[1]
        coord = json.loads(coord_raw)
        
        print(f"Target Atom {atom_id} ({coord['name']}) coords: {coord['x']}, {coord['y']}, {coord['z']}")
        
        import time
        start = time.time()
        atoms_in_radius = client.find_atoms_in_radius(coord['x'], coord['y'], coord['z'], radius=4.0)
        dur = time.time() - start
        print(f"Found {len(atoms_in_radius)} atoms in a 4.0A radius in {dur*1000:.2f} ms")
        for a in atoms_in_radius[:3]:
            print(f"  -> Atom ID: {a.atom_id} ({a.name}) from {a.uniprot_id} at {a.x:.2f}, {a.y:.2f}, {a.z:.2f}")

