# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "psycopg",
#     "pgbio @ file://./pgbio-py",
# ]
# ///

import psycopg
import json
import math

DB_URI = "postgresql://localhost:28818/bio_demo"

def setup_attention():
    print("Setting up Ground-Truth Geometric Attention Maps...")
    
    with psycopg.connect(DB_URI) as conn:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS protein_attention_maps (
                    id SERIAL PRIMARY KEY,
                    uniprot_id VARCHAR(20) REFERENCES proteins(uniprot_id),
                    attention_data SparseAttentionMap
                );
            """)
            
            cur.execute("TRUNCATE protein_attention_maps RESTART IDENTITY;")
            
            # Fetch proteins that have actual atoms
            cur.execute("""
                SELECT DISTINCT uniprot_id FROM protein_atoms;
            """)
            proteins = cur.fetchall()
            
            if not proteins:
                print("No atoms found in database. Run 'uv run scripts/seed_bio_demo.py' and select 'pdb' source first.")
                return

            print(f"Generating geometric contact maps for {len(proteins)} structures...")
            
            for (pid,) in proteins:
                print(f"  -> Processing {pid}")
                
                # Fetch atoms for this protein
                cur.execute("SELECT atom_id, coord::text FROM protein_atoms WHERE uniprot_id = %s ORDER BY atom_id ASC;", (pid,))
                atoms = cur.fetchall()
                
                if not atoms:
                    continue
                    
                # To keep it fast for demo, we take the first 100 atoms
                atoms = atoms[:100]
                
                targets_list = []
                weights = []
                # Compute distance matrix and create sparse weights for < 8.0 Angstroms
                for i in range(len(atoms)):
                    c1 = json.loads(atoms[i][1])
                    for j in range(len(atoms)):
                        if i == j:
                            continue
                        c2 = json.loads(atoms[j][1])
                        dx = c1['x'] - c2['x']
                        dy = c1['y'] - c2['y']
                        dz = c1['z'] - c2['z']
                        dist = math.sqrt(dx*dx + dy*dy + dz*dz)
                        
                        if dist < 8.0:
                            # Convert distance into a normalized attention weight (closer = higher weight)
                            weight = max(0.01, 1.0 - (dist / 8.0))
                            targets_list.append(j + 1)
                            weights.append(weight)
                
                if not targets_list:
                    continue
                    
                # Store sparse interactions
                targets_pg = "{" + ",".join(map(str, targets_list)) + "}"
                weights_pg = "{" + ",".join(f"{w:.4f}" for w in weights) + "}"
                
                query = f"""
                    INSERT INTO protein_attention_maps (uniprot_id, attention_data)
                    VALUES (%s, row({len(atoms)}, '{targets_pg}', '{weights_pg}')::SparseAttentionMap)
                """
                cur.execute(query, (pid,))
                
            conn.commit()
            print("\nSuccessfully built ground-truth biological attention maps from real 3D coordinates!")

if __name__ == '__main__':
    setup_attention()
